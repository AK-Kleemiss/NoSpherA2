#pragma once

// Minimal stand-in for hipBLAS, covering only what itensor_gpu.cu uses. Same caveat as
// the hip_runtime shim: this proves the names resolve, not that AMD's library behaves
// identically. Only a real hipcc build can say that.

#include <cublas_v2.h>

#define hipblasHandle_t cublasHandle_t
#define hipblasStatus_t cublasStatus_t
#define hipblasCreate cublasCreate
#define hipblasDestroy cublasDestroy
#define hipblasSgemm cublasSgemm
#define hipblasDgemm cublasDgemm
#define hipblasOperation_t cublasOperation_t
#define hipblasDgemm cublasDgemm
#define hipblasOperation_t cublasOperation_t
#define HIPBLAS_OP_T CUBLAS_OP_T
#define HIPBLAS_OP_N CUBLAS_OP_N
#define HIPBLAS_STATUS_SUCCESS CUBLAS_STATUS_SUCCESS
