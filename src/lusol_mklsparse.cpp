#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <sycl/sycl.hpp>
#include <oneapi/mkl/spblas.hpp>
#include <omp.h>

/*
########################################################################
c Copyright 2022 Predictive Science Inc.
c
c Licensed under the Apache License, Version 2.0 (the "License");
c you may not use this file except in compliance with the License.
c You may obtain a copy of the License at
c
c    http://www.apache.org/licenses/LICENSE-2.0
c
c Unless required by applicable law or agreed to in writing, software
c distributed under the License is distributed on an "AS IS" BASIS,
c WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
c implied.
c See the License for the specific language governing permissions and
c limitations under the License.
c#######################################################################
*/

namespace mkl        = oneapi::mkl;
namespace mkl_sparse = oneapi::mkl::sparse;

/* ================================================================ */
/* Global state                                                      */
/* ================================================================ */

static sycl::queue*                sycl_q     = nullptr;
static mkl_sparse::matrix_handle_t L_handle   = nullptr;
static mkl_sparse::matrix_handle_t U_handle   = nullptr;
static float*                      x_32       = nullptr;
static float*                      y_32       = nullptr;
static int                         N_global   = 0;
static const float                 alpha      = 1.0f;

/* ================================================================ */
/* Host ILU0 (ported from pot3d.F90 ilu0 — Saad algorithm)          */
/* CSR arrays use 1-based Fortran-style indexing in IA/JA values.     */
/* ================================================================ */

static void compute_diag_ptr(int N, const int* IA, const int* JA, int* A_da)
{
    for (int i = 1; i <= N; i++) {
        A_da[i - 1] = 0;
        for (int ij = IA[i - 1]; ij < IA[i]; ij++) {
            if (JA[ij - 1] == i) {
                A_da[i - 1] = ij;
                break;
            }
        }
        if (A_da[i - 1] == 0) {
            printf(" ERROR! Missing diagonal at row %d\n", i);
            exit(1);
        }
    }
}

static int ilu0_host(int N, int M, float* A, int* JA, int* IA, int* A_da)
{
    (void)M;
    std::vector<int> iw(N, 0);

    for (int i = 2; i <= N; i++) {
        const int IA_i     = IA[i - 1];
        const int IA_ip1m1 = IA[i] - 1;

        for (int ij = IA_i; ij <= IA_ip1m1; ij++) {
            iw[JA[ij - 1] - 1] = ij;
        }

        for (int ik = IA_i; ik <= A_da[i - 1] - 1; ik++) {
            const int   k   = JA[ik - 1];
            const float Aik = A[ik - 1] / A[A_da[k - 1] - 1];
            A[ik - 1] = Aik;

            for (int kj = A_da[k - 1] + 1; kj <= IA[k] - 1; kj++) {
                const int ij = iw[JA[kj - 1] - 1];
                if (ij != 0) {
                    A[ij - 1] -= Aik * A[kj - 1];
                }
            }
        }

        if (A[A_da[i - 1] - 1] == 0.0f) {
            return i;
        }

        for (int ij = IA_i; ij <= IA_ip1m1; ij++) {
            iw[JA[ij - 1] - 1] = 0;
        }
    }

    return 0;
}

/* omp_get_mapped_ptr is NULL under OpenMP USM (maptype-modifier=present). */
extern "C" void* omp_get_mapped_ptr(const void* ptr, int device_num);

static float* csr_device_ptr(float* host_ptr)
{
    void* mapped = omp_get_mapped_ptr(host_ptr, omp_get_default_device());
    return mapped ? static_cast<float*>(mapped) : host_ptr;
}

static int* csr_device_ptr(int* host_ptr)
{
    void* mapped = omp_get_mapped_ptr(host_ptr, omp_get_default_device());
    return mapped ? static_cast<int*>(mapped) : host_ptr;
}

static void copy_omp(int dst_device, void* dst, int src_device, const void* src,
                     size_t nbytes)
{
    if (omp_target_memcpy(dst, src, nbytes, 0, 0,
                          dst_device, src_device) != 0) {
        printf(" ERROR! omp_target_memcpy failed (dst_dev=%d src_dev=%d nbytes=%zu)\n",
               dst_device, src_device, nbytes);
        exit(1);
    }
}

/* oneMKL has no SYCL ILU0; pull CSR to host, factorize, push values back. */
static void ilu0_via_host(float* csr_lu, int* csr_i, int* csr_j, int N, int M)
{
    const int host_device = omp_get_initial_device();
    const int dev_device  = omp_get_default_device();
    const bool usm        = (omp_get_mapped_ptr(csr_lu, dev_device) == nullptr);

    std::vector<float> h_vals(M);
    std::vector<int>   h_ia(N + 1), h_ja(M), h_da(N);

    if (usm) {
        /* USM (maptype-modifier=present): host pointers are valid after diacsr. */
        std::copy(csr_lu, csr_lu + M, h_vals.data());
        std::copy(csr_i, csr_i + N + 1, h_ia.data());
        std::copy(csr_j, csr_j + M, h_ja.data());
    } else {
        float* csr_lu_d = csr_device_ptr(csr_lu);
        int*   csr_i_d  = csr_device_ptr(csr_i);
        int*   csr_j_d  = csr_device_ptr(csr_j);
        copy_omp(host_device, h_vals.data(), dev_device, csr_lu_d, M * sizeof(float));
        copy_omp(host_device, h_ia.data(), dev_device, csr_i_d, (N + 1) * sizeof(int));
        copy_omp(host_device, h_ja.data(), dev_device, csr_j_d, M * sizeof(int));
    }

    compute_diag_ptr(N, h_ia.data(), h_ja.data(), h_da.data());
    const int ierr = ilu0_host(N, M, h_vals.data(), h_ja.data(),
                               h_ia.data(), h_da.data());
    if (ierr != 0) {
        printf(" ERROR! ILU0 zero pivot at row %d\n", ierr);
        exit(1);
    }

    if (usm) {
        std::copy(h_vals.data(), h_vals.data() + M, csr_lu);
    } else {
        float* csr_lu_d = csr_device_ptr(csr_lu);
        copy_omp(dev_device, csr_lu_d, host_device, h_vals.data(), M * sizeof(float));
    }
}

#ifdef __cplusplus
extern "C" {
#endif

/* ******************************************************************* */
/* *** Initialize: ILU0 factorization + SpSV analysis on Intel GPU *** */
/* ******************************************************************* */

void load_external_lu_solver(float* CSR_LU, int* CSR_I,
                         int* CSR_J, int N, int M)
{
    N_global = N;

    /* Fortran passes host addresses from C_LOC inside use_device_addr.
       With maptype-modifier=present (USM), omp_get_mapped_ptr returns NULL
       but the same pointer is valid on the GPU; use it directly for MKL. */
    float* csr_lu_d = csr_device_ptr(CSR_LU);
    int*   csr_i_d  = csr_device_ptr(CSR_I);
    int*   csr_j_d  = csr_device_ptr(CSR_J);

    ilu0_via_host(CSR_LU, CSR_I, CSR_J, N, M);

    sycl_q = new sycl::queue(sycl::gpu_selector_v,
                              sycl::property::queue::in_order{});

    x_32 = sycl::malloc_device<float>(N, *sycl_q);
    y_32 = sycl::malloc_device<float>(N, *sycl_q);
    sycl_q->fill<float>(x_32, 0.0f, N).wait();
    sycl_q->fill<float>(y_32, 0.0f, N).wait();

    mkl_sparse::init_matrix_handle(&L_handle);
    mkl_sparse::init_matrix_handle(&U_handle);

    try {
        mkl_sparse::set_csr_data(*sycl_q, L_handle,
            (std::int64_t)N, (std::int64_t)N, (std::int64_t)M,
            mkl::index_base::one,
            csr_i_d, csr_j_d, csr_lu_d, {}).wait();
    } catch (sycl::exception& e) {
        printf(" ERROR! L CSR Matrix Creation Error: %s\n", e.what());
        exit(1);
    }

    try {
        mkl_sparse::set_csr_data(*sycl_q, U_handle,
            (std::int64_t)N, (std::int64_t)N, (std::int64_t)M,
            mkl::index_base::one,
            csr_i_d, csr_j_d, csr_lu_d, {}).wait();
    } catch (sycl::exception& e) {
        printf(" ERROR! U CSR Matrix Creation Error: %s\n", e.what());
        exit(1);
    }

    try {
        mkl_sparse::optimize_trsv(*sycl_q,
            mkl::uplo::lower, mkl::transpose::nontrans, mkl::diag::unit,
            L_handle, {}).wait();
    } catch (sycl::exception& e) {
        printf(" ERROR! L Analysis Error: %s\n", e.what());
        exit(1);
    }

    try {
        mkl_sparse::optimize_trsv(*sycl_q,
            mkl::uplo::upper, mkl::transpose::nontrans, mkl::diag::nonunit,
            U_handle, {}).wait();
    } catch (sycl::exception& e) {
        printf(" ERROR! U Analysis Error: %s\n", e.what());
        exit(1);
    }
}

#ifdef __cplusplus
}
#endif

/* ******************************************************************* */
/* *** Solve: d2f → forward L solve → backward U solve → f2d      *** */
/* ******************************************************************* */

#ifdef __cplusplus
extern "C" {
#endif

void external_lu_solver(double* x_d)
{
    float*  x32  = x_32;
    float*  y32  = y_32;
    const int n  = N_global;

    sycl_q->parallel_for(sycl::range<1>(n), [=](sycl::id<1> idx) {
        x32[idx] = (float) x_d[idx];
    }).wait();

    try {
        mkl_sparse::trsv(*sycl_q,
            mkl::uplo::lower, mkl::transpose::nontrans, mkl::diag::unit,
            alpha, L_handle, x32, y32, {}).wait();
    } catch (sycl::exception& e) {
        printf(" ERROR! Forward Solve Error: %s\n", e.what());
        exit(1);
    }

    try {
        mkl_sparse::trsv(*sycl_q,
            mkl::uplo::upper, mkl::transpose::nontrans, mkl::diag::nonunit,
            alpha, U_handle, y32, x32, {}).wait();
    } catch (sycl::exception& e) {
        printf(" ERROR! Backward Solve Error: %s\n", e.what());
        exit(1);
    }

    sycl_q->parallel_for(sycl::range<1>(n), [=](sycl::id<1> idx) {
        x_d[idx] = (double) x32[idx];
    }).wait();
}

void unload_external_lu_solver()
{
    mkl_sparse::release_matrix_handle(*sycl_q, &L_handle, {}).wait();
    mkl_sparse::release_matrix_handle(*sycl_q, &U_handle, {}).wait();

    sycl::free(x_32, *sycl_q);
    sycl::free(y_32, *sycl_q);

    sycl_q->wait();
    delete sycl_q;
    sycl_q = nullptr;
}

#ifdef __cplusplus
}
#endif
