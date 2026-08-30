#pragma once

//cuBLAS, if the machine happens to have it, loaded by name at run time.
//
//Nothing here is linked and nothing is imported: the library is opened with LoadLibrary or
//dlopen and four entry points are resolved by hand, so a machine without cuBLAS is an
//ordinary false rather than a process the loader refuses to start. That is the whole reason
//for doing it this way rather than linking and delay-loading - the previous arrangement
//needed the import library at build time and shipped half a gigabyte to be useful.
//
//Preferred whenever it loads: it is clearly ahead on a datacentre card and within
//measurement noise of the fallbacks on consumer ones, so choosing it costs nothing where it
//does not help.
//
//The consequence has to be lived with rather than wished away. The two disagree in the last
//digits, so the same input gives slightly different output on a machine that has a CUDA
//toolkit and one that does not. That is why a run states which GEMM produced its numbers,
//and why -no_gpu_cublas exists and the reference tests use it: a test left to choose would
//pass or fail on whether a toolkit happened to be installed.

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
