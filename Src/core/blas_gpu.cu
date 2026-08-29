#include "blas_gpu.h"
#include "gpu_backend.h"
#include "sf_gpu.h"   //sf_gpu_fp64_ratio: the device property the threshold below scales on
#include <cstdio>

//Calibrated on this machine (RTX 4090 Laptop, fp64 at 1/64 rate, 16 CPU threads with MKL):
//the SALTED regression shape, 2*m*n*k = 17 GFLOP, measured 1.22x and was not worth taking.
//The crossover therefore sits above that, and the constant is set an order beyond it so a
//marginal shape stays on the CPU.
#define BLAS_GPU_MIN_FLOP_AT_RATIO_64 2.0e11

static bool g_blas_gpu = false;
void blas_gpu_set_enabled(bool on) { g_blas_gpu = on; }
bool blas_gpu_enabled() { return g_blas_gpu; }

//conda-forge packages hipcc without hipBLAS, so a device BLAS is not guaranteed to exist
//just because a device does. Everything below already returns false when there is no device
//and every caller already falls back to the CPU for that, so the no-BLAS build answers the
//same way rather than needing a story of its own.
#ifndef NOSPHERA2_HAVE_GPUBLAS

bool blas_gpu_available() { return false; }

bool blas_gpu_dgemm(const bool, const bool, const int, const int, const int,
	const double, const double*, const int, const double*, const int,
	const double, double*, const int)
{
	return false;
}

#else

//The comment above used to end by saying this would be far too conservative on a datacentre
//part and should be re-measured. It has been: a V100 runs the double-precision GEMM at 6657
//GFLOP/s against 13988 in single, half rate rather than a sixty-fourth. At that speed the
//17 GFLOP shape costs about 2.6 ms of arithmetic against roughly 12 ms of transfers, and
//still beats the host - so the threshold that was right on a consumer card refuses work a
//datacentre card should be taking.
//
//What made the constant right on one machine is the device's double-precision rate, and the
//fp32:fp64 ratio is the one proxy for that available without stopping to run a benchmark.
//Scaling by it keeps this machine's measured value exactly (ratio 64 gives 2.0e11 back) and
//lowers the bar in proportion where doubles are cheap. It is a proxy, not a measurement of
//the crossover on the card in hand: the honest form of that is a calibration GEMM at
//startup, which is worth doing only if this proves too coarse.
static double blas_gpu_min_flop()
{
	const int ratio = sf_gpu_fp64_ratio();
	if (ratio <= 0) return BLAS_GPU_MIN_FLOP_AT_RATIO_64;   //unknown device: stay conservative
	return BLAS_GPU_MIN_FLOP_AT_RATIO_64 * (double)ratio / 64.0;
}

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
	if (2.0 * m * n * k < blas_gpu_min_flop()) return false;
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

#endif //NOSPHERA2_HAVE_GPUBLAS
