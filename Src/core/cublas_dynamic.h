#pragma once

//cuBLAS opened by name at run time (LoadLibrary/dlopen, entry points resolved by hand, nothing linked), so a machine
//without it gets false instead of a loader failure. It differs from the fallbacks in the last digits, so reference
//tests run with -no_gpu_cublas to stay independent of whether a CUDA toolkit is installed.

void cublas_dynamic_set_enabled(bool on);
bool cublas_dynamic_enabled();

//Enabled and loadable
bool cublas_dynamic_available();
bool cublas_dynamic_fast_16f_available();

//Column-major, BLAS conventions
bool cublas_dynamic_gemm(bool transA, bool transB, int m, int n, int k,
	float alpha, const float* A, int lda, const float* B, int ldb,
	float beta, float* C, int ldc);

bool cublas_dynamic_gemm(bool transA, bool transB, int m, int n, int k,
	double alpha, const double* A, int lda, const double* B, int ldb,
	double beta, double* C, int ldc);

//FP32 in and out, FP16 Tensor Core operands with FP32 accumulation; false: use the ordinary SGEMM path
bool cublas_dynamic_gemm_fast_16f(bool transA, bool transB, int m, int n, int k,
	float alpha, const float* A, int lda, const float* B, int ldb,
	float beta, float* C, int ldc);
