#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cusparse.h>

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

cusparseHandle_t cusparseHandle=NULL;
cusparseStatus_t cusparseStatus;

cusparseSpMatDescr_t        L_mat;
cusparseSpMatDescr_t        U_mat;
cusparseDnVecDescr_t        DenseVecX;
cusparseDnVecDescr_t        DenseVecY;
cusparseMatDescr_t          M_described=0;
cusparseSpSVDescr_t         L_described=0;
cusparseSpSVDescr_t         U_described=0;
csrilu02Info_t              M_analyzed=0;

void * Mbuffer;
void * Lbuffer;
void * Ubuffer;

int Mbuf_size;
size_t Lbuf_size;
size_t Ubuf_size;

int struct_zero;
int n_zero;
int N_global;

const cusparseSolvePolicy_t M_policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
const cusparseOperation_t   L_trans  = CUSPARSE_OPERATION_NON_TRANSPOSE;
const cusparseOperation_t   U_trans  = CUSPARSE_OPERATION_NON_TRANSPOSE;

const float alpha = 1.0f;

float* restrict x_32;
float* restrict y_32;

/* ******************************************************************* */
/* *** Initialize the cusparse functions in SP: *** */
/* ******************************************************************* */

void load_lusol_cusparse(float* restrict CSR_LU, int* restrict CSR_I, 
                         int* restrict CSR_J,int N, int M)
{
  N_global = N;

  // Allocate global scratch arrays and initialize to 0.
  cudaMalloc((void**)&x_32,N_global*sizeof(float));
  cudaMalloc((void**)&y_32,N_global*sizeof(float));
  cudaMemset((void*) x_32,0,N_global*sizeof(float));
  cudaMemset((void*) y_32,0,N_global*sizeof(float));

  // Setup cusparse.
  cusparseCreate(&cusparseHandle);

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Set up M
  //
  ////////////////////////////////////////////////////////////////////////////////////

  cusparseCreateMatDescr(&M_described);
  cusparseSetMatIndexBase(M_described, CUSPARSE_INDEX_BASE_ONE);
  cusparseSetMatType(M_described, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseStatus = cusparseCreateCsrilu02Info(&M_analyzed);
  if (cusparseStatus!=0){
      printf(" ERROR! ILU0 Info Error =       %s \n",
             cusparseGetErrorString(cusparseStatus));
      exit(1);
  }

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Get algorithm buffer size for M
  //
  ////////////////////////////////////////////////////////////////////////////////////

  cusparseStatus = cusparseScsrilu02_bufferSize(cusparseHandle, N, M,
                    M_described, CSR_LU, CSR_I, CSR_J,
                    M_analyzed, &Mbuf_size);

  if (cusparseStatus!=0){
      printf(" ERROR! ILU0 Buffer Size Error =      %s \n",
             cusparseGetErrorString(cusparseStatus));
      exit(1);
  }

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Allocate buffers
  //
  ////////////////////////////////////////////////////////////////////////////////////

  cudaMalloc((void**)&Mbuffer, Mbuf_size);

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Analyze M
  //
  ////////////////////////////////////////////////////////////////////////////////////

  cusparseStatus = cusparseScsrilu02_analysis(cusparseHandle,
                    N, M, M_described, CSR_LU, CSR_I, CSR_J,
              M_analyzed, M_policy, Mbuffer);
  if (cusparseStatus!=0){
      printf(" ERROR! ILU0 Analysis Error =       %s \n",
             cusparseGetErrorString(cusparseStatus));
      exit(1);
  }

  cusparseStatus = cusparseXcsrilu02_zeroPivot(cusparseHandle,
                    M_analyzed, &struct_zero);
  if (CUSPARSE_STATUS_ZERO_PIVOT == cusparseStatus){
      printf(" ERROR! A(%d,%d) is missing\n",
           struct_zero, struct_zero);
      exit(1);
  }

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Set preconditioner (M=LU)
  //
  ////////////////////////////////////////////////////////////////////////////////////

  cusparseStatus = cusparseScsrilu02(cusparseHandle, N, M, M_described,
                    CSR_LU, CSR_I, CSR_J, M_analyzed, M_policy, Mbuffer);
  if (cusparseStatus!=0){
      printf(" ERROR! ILU0 Formation Error =      %s \n",
             cusparseGetErrorString(cusparseStatus));
      exit(1);
  }

  cusparseStatus = cusparseXcsrilu02_zeroPivot(cusparseHandle,
                    M_analyzed, &n_zero);
  if (CUSPARSE_STATUS_ZERO_PIVOT == cusparseStatus){
     printf(" ERROR! M(%d,%d) is zero\n", n_zero, n_zero);
     exit(1);
  }

  cudaFree(Mbuffer);
  cusparseDestroyCsrilu02Info(M_analyzed);
  cusparseDestroyMatDescr(M_described);

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Set up L
  //
  ////////////////////////////////////////////////////////////////////////////////////

  // Create the sparse matrix
  cusparseStatus = cusparseCreateCsr(&L_mat, N, N, M, CSR_I, CSR_J,
      CSR_LU, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
      CUSPARSE_INDEX_BASE_ONE, CUDA_R_32F);
  if (cusparseStatus!=0){
      printf(" ERROR! L CSR Matrix Creation Error =      %s \n",
             cusparseGetErrorString(cusparseStatus));
      exit(1);
  }

  // Initialize the data structure
  cusparseSpSV_createDescr(&L_described);

  // Set filling mode types
  cusparseFillMode_t L_fill_mode = CUSPARSE_FILL_MODE_LOWER;
  cusparseSpMatSetAttribute(L_mat, CUSPARSE_SPMAT_FILL_MODE,
      &L_fill_mode, sizeof(L_fill_mode));

  // Set diagonal unit types
  cusparseDiagType_t L_type_diag = CUSPARSE_DIAG_TYPE_UNIT;
  cusparseSpMatSetAttribute(L_mat, CUSPARSE_SPMAT_DIAG_TYPE,
      &L_type_diag, sizeof(L_type_diag));

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Set up U
  //
  ////////////////////////////////////////////////////////////////////////////////////

  // Creat the sparse matrix
  cusparseStatus = cusparseCreateCsr(&U_mat, N, N, M, CSR_I, CSR_J,
      CSR_LU, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
      CUSPARSE_INDEX_BASE_ONE, CUDA_R_32F);

  if (cusparseStatus!=0){
      printf(" ERROR! U CSR Matrix Creation Error =      %s \n",
             cusparseGetErrorString(cusparseStatus));
      exit(1);
  }

  // Initialize the data structure
  cusparseSpSV_createDescr(&U_described);

  // Set filling mode types
  cusparseFillMode_t U_fill_mode = CUSPARSE_FILL_MODE_UPPER;
  cusparseSpMatSetAttribute(U_mat, CUSPARSE_SPMAT_FILL_MODE,
      &U_fill_mode, sizeof(U_fill_mode));

  // Set diagonal unit types
  cusparseDiagType_t U_type_diag = CUSPARSE_DIAG_TYPE_NON_UNIT;
  cusparseSpMatSetAttribute(U_mat, CUSPARSE_SPMAT_DIAG_TYPE,
      &U_type_diag, sizeof(U_type_diag));

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Set up Dense Vectors
  //
  ////////////////////////////////////////////////////////////////////////////////////

  cusparseCreateDnVec(&DenseVecX, N, x_32, CUDA_R_32F);
  cusparseCreateDnVec(&DenseVecY, N, y_32, CUDA_R_32F);

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Get algorithm buffer sizes
  //
  ////////////////////////////////////////////////////////////////////////////////////

  cusparseStatus = cusparseSpSV_bufferSize(cusparseHandle, L_trans,
        &alpha, L_mat, DenseVecX, DenseVecY, CUDA_R_32F,
        CUSPARSE_SPSV_ALG_DEFAULT, L_described,&Lbuf_size);

  if (cusparseStatus!=0){
      printf(" ERROR! L Buffer Size Error =       %s \n",
             cusparseGetErrorString(cusparseStatus));
      exit(1);
  }

  cusparseStatus = cusparseSpSV_bufferSize(cusparseHandle, U_trans,
        &alpha, U_mat, DenseVecX, DenseVecY, CUDA_R_32F,
        CUSPARSE_SPSV_ALG_DEFAULT, U_described,&Ubuf_size);

  if (cusparseStatus!=0){
      printf(" ERROR! U Buffer Size Error =       %s \n",
             cusparseGetErrorString(cusparseStatus));
      exit(1);
  }

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Allocate buffers
  //
  ////////////////////////////////////////////////////////////////////////////////////

  cudaMalloc((void**)&Lbuffer, Lbuf_size);
  cudaMalloc((void**)&Ubuffer, Ubuf_size);

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Analyze L
  //
  ////////////////////////////////////////////////////////////////////////////////////

  cusparseStatus = cusparseSpSV_analysis(cusparseHandle, L_trans,
        &alpha, L_mat, DenseVecX, DenseVecY, CUDA_R_32F,
        CUSPARSE_SPSV_ALG_DEFAULT, L_described, Lbuffer);
    if (cusparseStatus!=0){
      printf(" ERROR! L Analysis Error =       %s \n",
             cusparseGetErrorString(cusparseStatus));
      exit(1);
  }

  cudaDeviceSynchronize();

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Analyze U
  //
  ////////////////////////////////////////////////////////////////////////////////////

  cusparseStatus = cusparseSpSV_analysis(cusparseHandle, U_trans,
        &alpha, U_mat, DenseVecX, DenseVecY, CUDA_R_32F,
        CUSPARSE_SPSV_ALG_DEFAULT, U_described, Ubuffer);
    if (cusparseStatus!=0){
      printf(" ERROR! U Analysis Error =       %s \n",
             cusparseGetErrorString(cusparseStatus));
      exit(1);
  }

  cudaDeviceSynchronize();
}


/* ******************************************************************* */
/* *** Do the Solve Phase in SP: *** */
/* ******************************************************************* */

void lusol_cusparse(double* restrict x)
{
  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Update the Dense Vector (already linked to DenseVecX in load)
  //
  ////////////////////////////////////////////////////////////////////////////////////

#pragma omp target teams distribute parallel for is_device_ptr(x,x_32)
  for (int i=0;i<N_global;i++){
    x_32[i] = (float) x[i];
  }

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Forward solve (Ly=x) (DenseVecY already linked to y_32 pointer in load)
  //
  ////////////////////////////////////////////////////////////////////////////////////

  cusparseStatus = cusparseSpSV_solve(cusparseHandle, L_trans,
      &alpha, L_mat, DenseVecX, DenseVecY, CUDA_R_32F,
      CUSPARSE_SPSV_ALG_DEFAULT, L_described);
  if (cusparseStatus!=0){
      printf(" ERROR! Forward Solve Error =       %s \n",
             cusparseGetErrorString(cusparseStatus));
      exit(1);
  }

  cudaDeviceSynchronize();

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Backward solve (Ux=y)
  //
  ////////////////////////////////////////////////////////////////////////////////////

  cusparseStatus = cusparseSpSV_solve(cusparseHandle, U_trans,
      &alpha, U_mat, DenseVecY, DenseVecX, CUDA_R_32F,
      CUSPARSE_SPSV_ALG_DEFAULT, U_described);
  if (cusparseStatus!=0){
      printf(" ERROR! Backward Solve Error =      %s \n",
             cusparseGetErrorString(cusparseStatus));
      exit(1);
  }

  cudaDeviceSynchronize();

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // Convert result to double precision.
  //
  ////////////////////////////////////////////////////////////////////////////////////

#pragma omp target teams distribute parallel for is_device_ptr(x,x_32)
  for (int i=0;i<N_global;i++){
    x[i] = (double) x_32[i];
  }

}

/* ******************************************************************* */
/* *** Free the Memory used by cusparse: *** */
/* ******************************************************************* */


void unload_lusol_cusparse()
{
  cusparseSpSV_destroyDescr(L_described);
  cusparseSpSV_destroyDescr(U_described);
  cusparseDestroySpMat(L_mat);
  cusparseDestroySpMat(U_mat);
  cusparseDestroyDnVec(DenseVecX);
  cusparseDestroyDnVec(DenseVecY);

  cudaFree(Lbuffer);
  cudaFree(Ubuffer);
  cudaFree(x_32);
  cudaFree(y_32);

  cusparseDestroy(cusparseHandle);

  cudaDeviceSynchronize();
}

