#pragma once

//A device GEMM behind the same row-major interface cblas_dgemm presents, so nos_math's
//dot_BLAS can hand off without the callers knowing. Works against cuBLAS or hipBLAS.
//
//It is size-gated, and that gate is the whole point. With the operands living on the host
//every call pays to ship them, and a small GEMM loses: the SALTED regression GEMM, which
//is m ~ 3400, n ~ 1000, k = 2500, measured 0.515 s on the device against 0.626 s on 16
//CPU threads - 1.22x, not worth a code path. Only shapes well past that are offered to
//the device, and the threshold below was calibrated by measurement, not chosen.
//
//Returns false whenever the caller should just call BLAS, which includes every case where
//no device is present.

bool blas_gpu_available();

//-gpu_blas turns the offload on; off unless asked, like the other GPU paths
void blas_gpu_set_enabled(bool on);
bool blas_gpu_enabled();

//Row-major C(m x n) = alpha * op(A) * op(B) + beta * C, matching cblas_dgemm's arguments.
//lda/ldb/ldc are the row-major leading dimensions cblas would take.
bool blas_gpu_dgemm(bool transA, bool transB, int m, int n, int k,
	double alpha, const double* A, int lda, const double* B, int ldb,
	double beta, double* C, int ldc);
