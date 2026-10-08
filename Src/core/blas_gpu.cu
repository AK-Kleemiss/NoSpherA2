#include "blas_gpu.h"
#include "tuning.h"
#include "gpu_backend.h"
#include "gemm_gpu.cuh"
#include "sf_gpu.h"
#include <cstdio>
#include <cstdlib>

NOSPHERA2_GPU_API_BEGIN

//Operands cross the bus on every call, so only GEMMs above this many flops go to the device. Set at fp32:fp64
//ratio 64 between break-even and a clear win; -gflops re-derives it, NOSPHERA2_BLAS_GPU_MIN_FLOP overrides it
#define BLAS_GPU_MIN_FLOP_AT_RATIO_64 4.0e9

static bool g_blas_gpu = false;
void blas_gpu_set_enabled(bool on) { g_blas_gpu = on; }
bool blas_gpu_enabled() { return g_blas_gpu; }

//Scaled by the fp32:fp64 ratio, the only proxy for the device rate without benchmarking: a lower bar where doubles are cheap
double blas_gpu_min_flop()
{
	if (const char* env = tuning("NOSPHERA2_BLAS_GPU_MIN_FLOP")) {
		const double v = std::atof(env);
		if (v > 0.0) return v;
	}
	const int ratio = sf_gpu_fp64_ratio();
	if (ratio <= 0) return BLAS_GPU_MIN_FLOP_AT_RATIO_64;   //unknown device: stay conservative
	return BLAS_GPU_MIN_FLOP_AT_RATIO_64 * (double)ratio / 64.0;
}

//Shared with the transform, so a card without code is diagnosed in one place
bool blas_gpu_available() { return sf_gpu_available(); }

bool blas_gpu_dgemm(const bool transA, const bool transB, const int m, const int n, const int k,
	const double alpha, const double* A, const int lda, const double* B, const int ldb,
	const double beta, double* C, const int ldc)
{
	if (!g_blas_gpu || m <= 0 || n <= 0 || k <= 0) return false;
	if (2.0 * m * n * k < blas_gpu_min_flop()) return false;
	if (!blas_gpu_available()) return false;

	//Row-major C read column-major is C^T = op(B)^T op(A)^T, so the column-major GEMM runs with A and B swapped
	const int cm = n, cn = m;
	const int splits = gemm_gpu::split_count(cm, cn, k);
	const size_t p_elems = gemm_gpu::workspace_elems(cm, cn, splits);

	const size_t a_rows = transA ? (size_t)k : (size_t)m;
	const size_t b_rows = transB ? (size_t)n : (size_t)k;
	const size_t a_bytes = sizeof(double) * a_rows * (size_t)lda;
	const size_t b_bytes = sizeof(double) * b_rows * (size_t)ldb;
	const size_t c_bytes = sizeof(double) * (size_t)m * (size_t)ldc;
	const size_t p_bytes = sizeof(double) * p_elems;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	if (a_bytes + b_bytes + c_bytes + p_bytes + (1u << 26) > freeb) return false;

	double *dA = nullptr, *dB = nullptr, *dC = nullptr, *dP = nullptr;
	if (gpuMalloc(&dA, a_bytes) != gpuSuccess) return false;
	if (gpuMalloc(&dB, b_bytes) != gpuSuccess) { gpuFree(dA); return false; }
	if (gpuMalloc(&dC, c_bytes) != gpuSuccess) { gpuFree(dA); gpuFree(dB); return false; }
	if (gpuMalloc(&dP, p_bytes) != gpuSuccess) { gpuFree(dA); gpuFree(dB); gpuFree(dC); return false; }

	bool ok = gpuMemcpy(dA, A, a_bytes, gpuMemcpyHostToDevice) == gpuSuccess
		   && gpuMemcpy(dB, B, b_bytes, gpuMemcpyHostToDevice) == gpuSuccess;
	if (ok && beta != 0.0)
		ok = gpuMemcpy(dC, C, c_bytes, gpuMemcpyHostToDevice) == gpuSuccess;

	if (ok)
		ok = gemm_gpu::launch<double>(transB, transA, cm, cn, k,
			alpha, dB, ldb, dA, lda, beta, dC, ldc, dP);
	if (ok) ok = gpuDeviceSynchronize() == gpuSuccess;
	if (ok) ok = gpuMemcpy(C, dC, c_bytes, gpuMemcpyDeviceToHost) == gpuSuccess;

	gpuFree(dA); gpuFree(dB); gpuFree(dC); gpuFree(dP);
	return ok;
}

NOSPHERA2_GPU_API_END
