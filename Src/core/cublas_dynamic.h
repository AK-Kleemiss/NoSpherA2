#pragma once

//cuBLAS, if the machine happens to have it, loaded by name at run time.
//
//Nothing here is linked and nothing is imported: the library is opened with LoadLibrary or
//dlopen and four entry points are resolved by hand, so a machine without cuBLAS is an
//ordinary false rather than a process the loader refuses to start. That is the whole reason
//for doing it this way rather than linking and delay-loading - the previous arrangement
//needed the import library at build time and shipped half a gigabyte to be useful.
//
//It is off unless asked for. Selecting it automatically would make the same input produce
//different last digits on two machines depending on whether a CUDA toolkit happened to be
//installed, which is the trap the reference logs and the fp32/fp64 tests already document.
//CUTLASS is the default and is always there.

void cublas_dynamic_set_enabled(bool on);
bool cublas_dynamic_enabled();

//True only if it is both asked for and actually loadable.
bool cublas_dynamic_available();

//Column-major, BLAS conventions, matching the GEMM the I tensor performs.
bool cublas_dynamic_gemm(bool transA, bool transB, int m, int n, int k,
	float alpha, const float* A, int lda, const float* B, int ldb,
	float beta, float* C, int ldc);

bool cublas_dynamic_gemm(bool transA, bool transB, int m, int n, int k,
	double alpha, const double* A, int lda, const double* B, int ldb,
	double beta, double* C, int ldc);
