#pragma once

#include "gpu_api.h"

//Device GEMM behind cblas_dgemm's row-major interface, so nos_math's dot_BLAS hands off unseen.
//Size-gated: the operands live on the host, and a small GEMM loses to shipping them.
//False means call BLAS, including whenever there is no device.

NOSPHERA2_GPU_API_BEGIN

bool blas_gpu_available();

//Off unless -gpu_blas
void blas_gpu_set_enabled(bool on);
bool blas_gpu_enabled();

//Smallest GEMM, in flops, that blas_gpu_dgemm ships to the device
double blas_gpu_min_flop();

//Row-major C(m x n) = alpha * op(A) * op(B) + beta * C, arguments and leading dimensions as cblas_dgemm.
bool blas_gpu_dgemm(bool transA, bool transB, int m, int n, int k,
	double alpha, const double* A, int lda, const double* B, int ldb,
	double beta, double* C, int ldc);

NOSPHERA2_GPU_API_END
