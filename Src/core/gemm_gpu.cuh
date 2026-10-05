#pragma once

#include "gpu_backend.h"

//Device GEMM in place of cuBLAS (no static or trimmable build on Windows) and hipBLAS (missing from conda-forge hipcc).
//BLAS conventions: column-major, op(A) m x k, op(B) k x n, C m x n; beta == 0 writes C without reading it.

NOSPHERA2_GPU_API_BEGIN
namespace gemm_gpu {

//4x4 per thread is eight shared loads per sixteen FMAs; the inner loop is bound by that ratio, not by
//occupancy, so one fat tile suits every architecture. Wider wastes half its work at m ~ 100.
struct tile_config { int BM, BN, BK, TM, TN; };
constexpr tile_config TILE{64, 64, 32, 4, 4};

//Transpose flags as template parameters so the indexing folds instead of branching per element.
template <typename T, bool TA>
__device__ inline T a_at(const T* __restrict__ A, const int lda, const int i, const int l)
{
	return TA ? A[(long long)i * lda + l] : A[(long long)l * lda + i];
}

template <typename T, bool TB>
__device__ inline T b_at(const T* __restrict__ B, const int ldb, const int l, const int j)
{
	return TB ? B[(long long)l * ldb + j] : B[(long long)j * ldb + l];
}

//One k-slice of C into its own slot of P; only splitting k fills the device for the I tensor (m, n ~ 100, k ~ 1000s).
template <typename T, bool TA, bool TB, int BM, int BN, int BK, int TM, int TN>
__global__ void gemm_partial_kernel(const int m, const int n, const int k, const int splits,
	const T* __restrict__ A, const int lda, const T* __restrict__ B, const int ldb,
	T* __restrict__ P)
{
	constexpr int NTHREADS = (BM / TM) * (BN / TN);
	__shared__ T As[BK][BM];
	__shared__ T Bs[BK][BN];

	const int tid = threadIdx.y * (BN / TN) + threadIdx.x;
	const int row0 = blockIdx.y * BM + threadIdx.y * TM;
	const int col0 = blockIdx.x * BN + threadIdx.x * TN;

	//Slices are contiguous and cover k exactly; the ends differ by at most one element.
	const int z = blockIdx.z;
	const int k0 = (int)((long long)k * z / splits);
	const int k1 = (int)((long long)k * (z + 1) / splits);

	T acc[TM][TN];
	for (int a = 0; a < TM; a++)
		for (int b = 0; b < TN; b++) acc[a][b] = T(0);

	//Next depth tile prefetched into registers so the global load overlaps the arithmetic; a second
	//shared buffer would make 64 KB in fp64, over the 48 KB static limit.
	constexpr int AREG = (BK * BM) / NTHREADS;
	constexpr int BREG = (BK * BN) / NTHREADS;
	T ra[AREG], rb[BREG];

	const auto fetch = [&](const int kt) {
		for (int r = 0; r < AREG; r++) {
			const int idx = tid + r * NTHREADS;
			const int gi = blockIdx.y * BM + idx % BM;
			const int gl = kt + idx / BM;
			ra[r] = (gi < m && gl < k1) ? a_at<T, TA>(A, lda, gi, gl) : T(0);
		}
		for (int r = 0; r < BREG; r++) {
			const int idx = tid + r * NTHREADS;
			const int gj = blockIdx.x * BN + idx % BN;
			const int gl = kt + idx / BN;
			rb[r] = (gj < n && gl < k1) ? b_at<T, TB>(B, ldb, gl, gj) : T(0);
		}
	};

	fetch(k0);
	for (int kt = k0; kt < k1; kt += BK) {
		//Everyone must be done reading the tile before it is overwritten
		__syncthreads();
		for (int r = 0; r < AREG; r++) {
			const int idx = tid + r * NTHREADS;
			As[idx / BM][idx % BM] = ra[r];
		}
		for (int r = 0; r < BREG; r++) {
			const int idx = tid + r * NTHREADS;
			Bs[idx / BN][idx % BN] = rb[r];
		}
		__syncthreads();

		if (kt + BK < k1) fetch(kt + BK);

		for (int l = 0; l < BK; l++) {
			T av[TM], bv[TN];
			for (int a = 0; a < TM; a++) av[a] = As[l][threadIdx.y * TM + a];
			for (int b = 0; b < TN; b++) bv[b] = Bs[l][threadIdx.x * TN + b];
			for (int a = 0; a < TM; a++)
				for (int b = 0; b < TN; b++) acc[a][b] += av[a] * bv[b];
		}
	}

	T* Pz = P + (long long)z * m * n;
	for (int a = 0; a < TM; a++) {
		const int i = row0 + a;
		if (i >= m) continue;
		for (int b = 0; b < TN; b++) {
			const int j = col0 + b;
			if (j < n) Pz[(long long)j * m + i] = acc[a][b];
		}
	}
}

//Fixed slice order keeps the sum deterministic; atomics would make it differ run to run.
template <typename T>
__global__ void gemm_reduce_kernel(const int m, const int n, const int splits,
	const T alpha, const T* __restrict__ P, const T beta, T* __restrict__ C, const int ldc)
{
	const int i = blockIdx.y * blockDim.y + threadIdx.y;
	const int j = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= m || j >= n) return;
	T s = T(0);
	for (int z = 0; z < splits; z++) s += P[(long long)z * m * n + (long long)j * m + i];
	T* c = C + (long long)j * ldc + i;
	*c = (beta == T(0)) ? alpha * s : alpha * s + beta * *c;
}

//Enough blocks to fill the device, but no slices so thin that reduction and tails outweigh it.
inline int split_count(const int m, const int n, const int k)
{
	constexpr tile_config t = TILE;
	const long long tiles = (long long)((m + t.BM - 1) / t.BM) * ((n + t.BN - 1) / t.BN);
	long long s = 512 / (tiles > 0 ? tiles : 1);
	const long long by_depth = k / (4 * t.BK);
	if (s > by_depth) s = by_depth;
	if (s < 1) s = 1;
	if (s > 64) s = 64;
	return (int)s;
}

inline size_t workspace_elems(const int m, const int n, const int splits)
{
	return (size_t)m * n * splits;
}

template <typename T, int BM, int BN, int BK, int TM, int TN>
inline void launch_tiled(const bool transA, const bool transB,
	const int m, const int n, const int k, const int splits,
	const T* A, const int lda, const T* B, const int ldb, T* P)
{
	const dim3 thr(BN / TN, BM / TM);
	const dim3 blk((n + BN - 1) / BN, (m + BM - 1) / BM, splits);
	if (!transA && !transB)
		gemm_partial_kernel<T, false, false, BM, BN, BK, TM, TN><<<blk, thr>>>(m, n, k, splits, A, lda, B, ldb, P);
	else if (transA && !transB)
		gemm_partial_kernel<T, true, false, BM, BN, BK, TM, TN><<<blk, thr>>>(m, n, k, splits, A, lda, B, ldb, P);
	else if (!transA && transB)
		gemm_partial_kernel<T, false, true, BM, BN, BK, TM, TN><<<blk, thr>>>(m, n, k, splits, A, lda, B, ldb, P);
	else
		gemm_partial_kernel<T, true, true, BM, BN, BK, TM, TN><<<blk, thr>>>(m, n, k, splits, A, lda, B, ldb, P);
}

//P must hold workspace_elems(m, n, split_count(m, n, k)) elements.
template <typename T>
inline bool launch(const bool transA, const bool transB,
	const int m, const int n, const int k,
	const T alpha, const T* A, const int lda, const T* B, const int ldb,
	const T beta, T* C, const int ldc, T* P)
{
	if (m <= 0 || n <= 0 || k <= 0) return false;
	const int splits = split_count(m, n, k);
	launch_tiled<T, TILE.BM, TILE.BN, TILE.BK, TILE.TM, TILE.TN>(
		transA, transB, m, n, k, splits, A, lda, B, ldb, P);

	const dim3 rthr(16, 16);
	const dim3 rblk((n + 15) / 16, (m + 15) / 16);
	gemm_reduce_kernel<T><<<rblk, rthr>>>(m, n, splits, alpha, P, beta, C, ldc);
	return true;
}

} //namespace gemm_gpu
NOSPHERA2_GPU_API_END
