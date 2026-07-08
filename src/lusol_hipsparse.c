#include <stdlib.h>
#include <stdio.h>
#include <hip/hip_runtime.h>
#include <hipsparse/hipsparse.h>
#include <omp.h>

/* HIP/clang++ may #define restrict; use __restrict__ for pointer qualifiers. */
#ifdef restrict
#undef restrict
#endif

#define POT3D_HIP_THREADS 256

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

#ifdef __cplusplus
extern "C" {
#endif

hipsparseHandle_t hipsparseHandle=NULL;
hipsparseStatus_t hipsparseStatus;

hipsparseSpMatDescr_t        L_mat;
hipsparseSpMatDescr_t        U_mat;
hipsparseDnVecDescr_t        DenseVecX;
hipsparseDnVecDescr_t        DenseVecY;
hipsparseMatDescr_t          M_described=0;
hipsparseSpSVDescr_t         L_described=0;
hipsparseSpSVDescr_t         U_described=0;
csrilu02Info_t               M_analyzed=0;

/* Single shared buffer for L and U (sequential solves — no overlap) */
void*  spSV_buf      = NULL;
size_t spSV_buf_size = 0;

int struct_zero;
int n_zero;
int N_global;

const hipsparseSolvePolicy_t M_policy = HIPSPARSE_SOLVE_POLICY_USE_LEVEL;
const hipsparseOperation_t   L_trans  = HIPSPARSE_OPERATION_NON_TRANSPOSE;
const hipsparseOperation_t   U_trans  = HIPSPARSE_OPERATION_NON_TRANSPOSE;

const float alpha = 1.0f;

float* __restrict__ x_32;
float* __restrict__ y_32;

/* ******************************************************************* */
/* *** Initialize the hipsparse functions in SP: *** */
/* ******************************************************************* */

void load_external_lu_solver(float* __restrict__ CSR_LU, int* __restrict__ CSR_I,
                             int* __restrict__ CSR_J, int N, int M)
{
  N_global = N;

  extern void* omp_get_mapped_ptr(const void* ptr, int device_num);
  int device_id = omp_get_default_device();
  float* csr_lu_d = (float*) omp_get_mapped_ptr(CSR_LU, device_id);
  int*   csr_i_d  = (int*)   omp_get_mapped_ptr(CSR_I,  device_id);
  int*   csr_j_d  = (int*)   omp_get_mapped_ptr(CSR_J,  device_id);

  hipMalloc((void**)&x_32, (size_t)N_global * sizeof(float));
  hipMalloc((void**)&y_32, (size_t)N_global * sizeof(float));
  hipMemset(x_32, 0, (size_t)N_global * sizeof(float));
  hipMemset(y_32, 0, (size_t)N_global * sizeof(float));

  hipsparseCreate(&hipsparseHandle);

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Set up M (ILU0 preconditioner) — all M state is local; freed before return
  //
  ////////////////////////////////////////////////////////////////////////////////////

  void* Mbuffer   = NULL;
  int   Mbuf_size = 0;

  hipsparseCreateMatDescr(&M_described);
  hipsparseSetMatIndexBase(M_described, HIPSPARSE_INDEX_BASE_ONE);
  hipsparseSetMatType(M_described, HIPSPARSE_MATRIX_TYPE_GENERAL);
  hipsparseStatus = hipsparseCreateCsrilu02Info(&M_analyzed);
  if (hipsparseStatus!=HIPSPARSE_STATUS_SUCCESS){
      printf(" ERROR! ILU0 Info Error =       %s \n",
             hipsparseGetErrorString(hipsparseStatus));
      exit(1);
  }

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Get algorithm buffer size for M
  //
  ////////////////////////////////////////////////////////////////////////////////////

  hipsparseStatus = hipsparseScsrilu02_bufferSize(hipsparseHandle, N, M,
                    M_described, csr_lu_d, csr_i_d, csr_j_d,
                    M_analyzed, &Mbuf_size);

  if (hipsparseStatus!=HIPSPARSE_STATUS_SUCCESS){
      printf(" ERROR! ILU0 Buffer Size Error =      %s \n",
             hipsparseGetErrorString(hipsparseStatus));
      exit(1);
  }

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Allocate buffers
  //
  ////////////////////////////////////////////////////////////////////////////////////

  hipMalloc((void**)&Mbuffer, Mbuf_size);

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Analyze M
  //
  ////////////////////////////////////////////////////////////////////////////////////

  hipsparseStatus = hipsparseScsrilu02_analysis(hipsparseHandle,
                    N, M, M_described, csr_lu_d, csr_i_d, csr_j_d,
              M_analyzed, M_policy, Mbuffer);
  if (hipsparseStatus!=HIPSPARSE_STATUS_SUCCESS){
      printf(" ERROR! ILU0 Analysis Error =       %s \n",
             hipsparseGetErrorString(hipsparseStatus));
      exit(1);
  }

  hipsparseStatus = hipsparseXcsrilu02_zeroPivot(hipsparseHandle,
                    M_analyzed, &struct_zero);
  if (HIPSPARSE_STATUS_ZERO_PIVOT == hipsparseStatus){
      printf(" ERROR! A(%d,%d) is missing\n",
           struct_zero, struct_zero);
      exit(1);
  }

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Set preconditioner (M=LU)
  //
  ////////////////////////////////////////////////////////////////////////////////////

  hipsparseStatus = hipsparseScsrilu02(hipsparseHandle, N, M, M_described,
                    csr_lu_d, csr_i_d, csr_j_d, M_analyzed, M_policy, Mbuffer);
  if (hipsparseStatus!=HIPSPARSE_STATUS_SUCCESS){
      printf(" ERROR! ILU0 Formation Error =      %s \n",
             hipsparseGetErrorString(hipsparseStatus));
      exit(1);
  }

  hipsparseStatus = hipsparseXcsrilu02_zeroPivot(hipsparseHandle,
                    M_analyzed, &n_zero);
  if (HIPSPARSE_STATUS_ZERO_PIVOT == hipsparseStatus){
     printf(" ERROR! M(%d,%d) is zero\n", n_zero, n_zero);
     exit(1);
  }

  /* Free all M state immediately — not needed after factorization */
  hipFree(Mbuffer);
  hipsparseDestroyCsrilu02Info(M_analyzed);
  hipsparseDestroyMatDescr(M_described);

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Set up L
  //
  ////////////////////////////////////////////////////////////////////////////////////

  hipsparseStatus = hipsparseCreateCsr(&L_mat, N, N, M,
      csr_i_d, csr_j_d, csr_lu_d, HIPSPARSE_INDEX_32I, HIPSPARSE_INDEX_32I,
      HIPSPARSE_INDEX_BASE_ONE, HIP_R_32F);
  if (hipsparseStatus!=HIPSPARSE_STATUS_SUCCESS){
      printf(" ERROR! L CSR Matrix Creation Error =      %s \n",
             hipsparseGetErrorString(hipsparseStatus));
      exit(1);
  }

  hipsparseSpSV_createDescr(&L_described);

  hipsparseFillMode_t L_fill_mode = HIPSPARSE_FILL_MODE_LOWER;
  hipsparseSpMatSetAttribute(L_mat, HIPSPARSE_SPMAT_FILL_MODE,
      &L_fill_mode, sizeof(L_fill_mode));

  hipsparseDiagType_t L_type_diag = HIPSPARSE_DIAG_TYPE_UNIT;
  hipsparseSpMatSetAttribute(L_mat, HIPSPARSE_SPMAT_DIAG_TYPE,
      &L_type_diag, sizeof(L_type_diag));

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Set up U
  //
  ////////////////////////////////////////////////////////////////////////////////////

  hipsparseStatus = hipsparseCreateCsr(&U_mat, N, N, M,
      csr_i_d, csr_j_d, csr_lu_d, HIPSPARSE_INDEX_32I, HIPSPARSE_INDEX_32I,
      HIPSPARSE_INDEX_BASE_ONE, HIP_R_32F);

  if (hipsparseStatus!=HIPSPARSE_STATUS_SUCCESS){
      printf(" ERROR! U CSR Matrix Creation Error =      %s \n",
             hipsparseGetErrorString(hipsparseStatus));
      exit(1);
  }

  hipsparseSpSV_createDescr(&U_described);

  hipsparseFillMode_t U_fill_mode = HIPSPARSE_FILL_MODE_UPPER;
  hipsparseSpMatSetAttribute(U_mat, HIPSPARSE_SPMAT_FILL_MODE,
      &U_fill_mode, sizeof(U_fill_mode));

  hipsparseDiagType_t U_type_diag = HIPSPARSE_DIAG_TYPE_NON_UNIT;
  hipsparseSpMatSetAttribute(U_mat, HIPSPARSE_SPMAT_DIAG_TYPE,
      &U_type_diag, sizeof(U_type_diag));

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Set up Dense Vectors
  //
  ////////////////////////////////////////////////////////////////////////////////////

  hipsparseCreateDnVec(&DenseVecX, N, x_32, HIP_R_32F);
  hipsparseCreateDnVec(&DenseVecY, N, y_32, HIP_R_32F);

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Get buffer sizes for L and U — allocate once using the larger of the two
  //
  ////////////////////////////////////////////////////////////////////////////////////

  size_t Lbuf_size = 0;
  size_t Ubuf_size = 0;

  hipsparseStatus = hipsparseSpSV_bufferSize(hipsparseHandle, L_trans,
        &alpha, L_mat, DenseVecX, DenseVecY, HIP_R_32F,
        HIPSPARSE_SPSV_ALG_DEFAULT, L_described, &Lbuf_size);

  if (hipsparseStatus!=HIPSPARSE_STATUS_SUCCESS){
      printf(" ERROR! L Buffer Size Error =       %s \n",
             hipsparseGetErrorString(hipsparseStatus));
      exit(1);
  }

  hipsparseStatus = hipsparseSpSV_bufferSize(hipsparseHandle, U_trans,
        &alpha, U_mat, DenseVecY, DenseVecX, HIP_R_32F,
        HIPSPARSE_SPSV_ALG_DEFAULT, U_described, &Ubuf_size);

  if (hipsparseStatus!=HIPSPARSE_STATUS_SUCCESS){
      printf(" ERROR! U Buffer Size Error =       %s \n",
             hipsparseGetErrorString(hipsparseStatus));
      exit(1);
  }

  /* Single allocation sized to the larger requirement */
  spSV_buf_size = (Lbuf_size > Ubuf_size) ? Lbuf_size : Ubuf_size;
  hipMalloc((void**)&spSV_buf, spSV_buf_size);

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Analyze L then U — reuse the same buffer (sequential, no overlap)
  //
  ////////////////////////////////////////////////////////////////////////////////////

  hipsparseStatus = hipsparseSpSV_analysis(hipsparseHandle, L_trans,
        &alpha, L_mat, DenseVecX, DenseVecY, HIP_R_32F,
        HIPSPARSE_SPSV_ALG_DEFAULT, L_described, spSV_buf);
  if (hipsparseStatus!=HIPSPARSE_STATUS_SUCCESS){
      printf(" ERROR! L Analysis Error =       %s \n",
             hipsparseGetErrorString(hipsparseStatus));
      exit(1);
  }

  hipDeviceSynchronize();

  hipsparseStatus = hipsparseSpSV_analysis(hipsparseHandle, U_trans,
        &alpha, U_mat, DenseVecY, DenseVecX, HIP_R_32F,
        HIPSPARSE_SPSV_ALG_DEFAULT, U_described, spSV_buf);
  if (hipsparseStatus!=HIPSPARSE_STATUS_SUCCESS){
      printf(" ERROR! U Analysis Error =       %s \n",
             hipsparseGetErrorString(hipsparseStatus));
      exit(1);
  }

  hipDeviceSynchronize();
}

#ifdef __cplusplus
}
#endif

/* ******************************************************************** */
/* *** HIP kernels: double <-> float conversion entirely on device. *** */
/* ******************************************************************** */

__global__ void d2f_kernel(float* __restrict__ out,
                           const double* __restrict__ in, int n)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) out[i] = (float)in[i];
}

__global__ void f2d_kernel(double* __restrict__ out,
                           const float* __restrict__ in, int n)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) out[i] = (double)in[i];
}

#ifdef __cplusplus
extern "C" {
#endif

/* ******************************************************************** */
/* *** Solve using an OMP/HIP device pointer.  Converts x_dev       *** */
/* *** (double) to x_32 (float) on device, solves L and U, then     *** */
/* *** converts the result back.                                    *** */
/* ******************************************************************** */

void external_lu_solver(void* x_dev)
{
  int blocks = (N_global + POT3D_HIP_THREADS - 1) / POT3D_HIP_THREADS;

  d2f_kernel<<<blocks, POT3D_HIP_THREADS>>>(x_32, (const double*)x_dev, N_global);
  hipDeviceSynchronize();

  hipsparseStatus = hipsparseSpSV_solve(hipsparseHandle, L_trans,
      &alpha, L_mat, DenseVecX, DenseVecY, HIP_R_32F,
      HIPSPARSE_SPSV_ALG_DEFAULT, L_described);
  if (hipsparseStatus!=HIPSPARSE_STATUS_SUCCESS){
      printf(" ERROR! Forward Solve Error =       %s \n",
             hipsparseGetErrorString(hipsparseStatus));
      exit(1);
  }
  hipDeviceSynchronize();

  hipsparseStatus = hipsparseSpSV_solve(hipsparseHandle, U_trans,
      &alpha, U_mat, DenseVecY, DenseVecX, HIP_R_32F,
      HIPSPARSE_SPSV_ALG_DEFAULT, U_described);
  if (hipsparseStatus!=HIPSPARSE_STATUS_SUCCESS){
      printf(" ERROR! Backward Solve Error =      %s \n",
             hipsparseGetErrorString(hipsparseStatus));
      exit(1);
  }
  hipDeviceSynchronize();

  f2d_kernel<<<blocks, POT3D_HIP_THREADS>>>((double*)x_dev, x_32, N_global);
  hipDeviceSynchronize();
}

/* ******************************************************************* */
/* *** Free the Memory used by hipsparse: *** */
/* ******************************************************************* */

void unload_external_lu_solver()
{
  hipsparseSpSV_destroyDescr(L_described);
  hipsparseSpSV_destroyDescr(U_described);
  hipsparseDestroySpMat(L_mat);
  hipsparseDestroySpMat(U_mat);
  hipsparseDestroyDnVec(DenseVecX);
  hipsparseDestroyDnVec(DenseVecY);

  hipFree(spSV_buf);
  hipFree(x_32);
  hipFree(y_32);
  hipsparseDestroy(hipsparseHandle);

  hipDeviceSynchronize();
}

#ifdef __cplusplus
}
#endif
