#include "blas_gpu.h"
#include "gpu_backend.h"
#include "gemm_gpu.cuh"
#include "sf_gpu.h"   //sf_gpu_fp64_ratio: the device property the threshold below scales on
#include <cstdio>
#include <cstdlib>

//What the device has to earn before it is offered the work. The operands live on the host,
//so every call ships them across and back, and a shape too small to hide that loses. Stated
//in flops so it does not depend on how the caller shaped the matrices.
//
//Measured on square host-resident dgemm, RTX 4090 mobile against 16 Zen4 threads: the
//device draws level near 2e9 flops and is clearly ahead by 8e9. The value sits between the
//two. The predecessor was 2e11 - two orders above the crossover, chosen when there was no
//measurement - which put the whole path out of reach of anything the program actually does.
//Re-derive it elsewhere with -gflops; NOSPHERA2_BLAS_GPU_MIN_FLOP overrides it.
#define BLAS_GPU_MIN_FLOP_AT_RATIO_64 4.0e9

static bool g_blas_gpu = false;
void blas_gpu_set_enabled(bool on) { g_blas_gpu = on; }
bool blas_gpu_enabled() { return g_blas_gpu; }

//The device rate is what made the constant right, and the fp32:fp64 ratio is the only proxy
//for it available without benchmarking. Scaling by it keeps the measured value where it was
//measured and lowers the bar where doubles are cheap.
static double blas_gpu_min_flop()
{
	if (const char* env = std::getenv("NOSPHERA2_BLAS_GPU_MIN_FLOP")) {
		const double v = std::atof(env);
		if (v > 0.0) return v;
	}
	const int ratio = sf_gpu_fp64_ratio();
	if (ratio <= 0) return BLAS_GPU_MIN_FLOP_AT_RATIO_64;   //unknown device: stay conservative
	return BLAS_GPU_MIN_FLOP_AT_RATIO_64 * (double)ratio / 64.0;
}

//Shared with the transform so the "no code for this card" case is diagnosed in one place.
bool blas_gpu_available() { return sf_gpu_available(); }

bool blas_gpu_dgemm(const bool transA, const bool transB, const int m, const int n, const int k,
	const double alpha, const double* A, const int lda, const double* B, const int ldb,
	const double beta, double* C, const int ldc)
{
	if (!g_blas_gpu || m <= 0 || n <= 0 || k <= 0) return false;
	if (2.0 * m * n * k < blas_gpu_min_flop()) return false;
	if (!blas_gpu_available()) return false;

	//A row-major matrix read column-major is its own transpose, so computing
	//C^T = op(B)^T * op(A)^T with the operands swapped gives the row-major C. The GEMM below
	//takes column-major arguments, so the swap is what adapts it to this interface.
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
