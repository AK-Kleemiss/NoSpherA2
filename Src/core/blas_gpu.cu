#include "blas_gpu.h"
#include "gpu_backend.h"
#include <cstdio>

//Calibrated on this machine (RTX 4090 Laptop, fp64 at 1/64 rate, 16 CPU threads with MKL):
//the SALTED regression shape, 2*m*n*k = 17 GFLOP, measured 1.22x and was not worth taking.
//The crossover therefore sits above that, and the constant is set an order beyond it so a
//marginal shape stays on the CPU. On a datacentre part, where fp64 is half rate rather
//than 1/64, this is far too conservative and should be re-measured before being trusted.
#define BLAS_GPU_MIN_FLOP 2.0e11

static bool g_blas_gpu = false;
void blas_gpu_set_enabled(bool on) { g_blas_gpu = on; }
bool blas_gpu_enabled() { return g_blas_gpu; }

bool blas_gpu_available()
{
	int n = 0;
	return gpuGetDeviceCount(&n) == gpuSuccess && n > 0;
}

bool blas_gpu_dgemm(const bool transA, const bool transB, const int m, const int n, const int k,
	const double alpha, const double* A, const int lda, const double* B, const int ldb,
	const double beta, double* C, const int ldc)
{
	if (!g_blas_gpu || m <= 0 || n <= 0 || k <= 0) return false;
	if (2.0 * m * n * k < BLAS_GPU_MIN_FLOP) return false;
	if (!blas_gpu_available() || !gpu_blas_runtime_present()) return false;

	static gpublasHandle_t handle = nullptr;
	if (!handle && gpublasCreate(&handle) != GPUBLAS_STATUS_SUCCESS) return false;

	const size_t a_rows = transA ? (size_t)k : (size_t)m;
	const size_t b_rows = transB ? (size_t)n : (size_t)k;
	const size_t a_bytes = sizeof(double) * a_rows * (size_t)lda;
	const size_t b_bytes = sizeof(double) * b_rows * (size_t)ldb;
	const size_t c_bytes = sizeof(double) * (size_t)m * (size_t)ldc;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	if (a_bytes + b_bytes + c_bytes + (1u << 26) > freeb) return false;

	double *dA = nullptr, *dB = nullptr, *dC = nullptr;
	if (gpuMalloc(&dA, a_bytes) != gpuSuccess) return false;
	if (gpuMalloc(&dB, b_bytes) != gpuSuccess) { gpuFree(dA); return false; }
	if (gpuMalloc(&dC, c_bytes) != gpuSuccess) { gpuFree(dA); gpuFree(dB); return false; }

	bool ok = gpuMemcpy(dA, A, a_bytes, gpuMemcpyHostToDevice) == gpuSuccess
	       && gpuMemcpy(dB, B, b_bytes, gpuMemcpyHostToDevice) == gpuSuccess;
	if (ok && beta != 0.0)
		ok = gpuMemcpy(dC, C, c_bytes, gpuMemcpyHostToDevice) == gpuSuccess;

	if (ok) {
		//A row-major matrix read column-major is its own transpose, so computing
		//C^T = op(B)^T * op(A)^T with the operands swapped gives the row-major C.
		const gpublasOperation_t opA = transA ? GPUBLAS_OP_T : GPUBLAS_OP_N;
		const gpublasOperation_t opB = transB ? GPUBLAS_OP_T : GPUBLAS_OP_N;
		ok = gpublasDgemm(handle, opB, opA, n, m, k, &alpha,
			dB, ldb, dA, lda, &beta, dC, ldc) == GPUBLAS_STATUS_SUCCESS;
	}
	if (ok) ok = gpuDeviceSynchronize() == gpuSuccess;
	if (ok) ok = gpuMemcpy(C, dC, c_bytes, gpuMemcpyDeviceToHost) == gpuSuccess;

	gpuFree(dA); gpuFree(dB); gpuFree(dC);
	return ok;
}
