#pragma once

#include "gpu_backend.h"

//Device GEMM, in place of cuBLAS and hipBLAS.
//
//Shipping cuBLAS costs half a gigabyte - cublasLt alone is the larger part of it - and it
//cannot be trimmed: there is no static cuBLAS on Windows, and nvprune refuses a DLL because
//it is not relocatable. Two call sites needed it, one of them a GEMM shape that cuBLAS
//itself runs at well under a fifth of peak. hipBLAS is worse than a size problem: conda-forge
//ships hipcc without it at all, so the AMD build had no BLAS to link against.
//
//BLAS conventions throughout: column-major, op(A) is m x k, op(B) is k x n, C is m x n, and
//beta == 0 means C is written rather than read, so an uninitialised C stays legal.
//
//Templated on the scalar type because the I tensor wants to choose its precision, not
//because two kernels were wanted.

namespace gemm_gpu {

//Tile shape. The narrow alternative spends four shared-memory loads on four FMAs where
//this one spends eight on sixteen, and that ratio is what the inner loop is limited by;
//halving the staged depth gives most of the difference back, which says the same thing.
//The constraint is arithmetic intensity, not occupancy - fewer, fatter blocks win.
//
//Selecting the tile per architecture was tried and removed. Turing has less shared memory
//per SM than Ada, so the narrow tile ought to suit it, and measurement says otherwise: the
//narrow tile is slower there too. One shape, no dispatch.
//
//The depth tile is prefetched, which was the leading explanation for the remaining gap to
//cuBLAS and turned out not to be it - a small consistent gain, not the difference. Untried:
//vectorised shared loads, the one structural difference left. A wider tile is not the
//answer, m being about a hundred here, so a 128-wide tile would discard half its work.
//
//Rates per card are in the vault note rather than here; they age and this file should not.
struct tile_config { int BM, BN, BK, TM, TN; };
constexpr tile_config TILE{64, 64, 32, 4, 4};

//The transpose flags are template parameters rather than arguments so the indexing folds
//away instead of branching once per element of the innermost loop.
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

//One k-slice of C into its own slot of P. Splitting k is what makes the I tensor shape fill
//a device at all: m and n are around a hundred there while k runs to several thousand, so
//the un-split grid is a few dozen blocks whatever the tile size.
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

	//The next depth tile is fetched into registers before the current one is multiplied, so
	//the global load is in flight across the arithmetic instead of being waited on straight
	//after it is issued. Registers rather than a second shared buffer because at this tile
	//size two fp64 buffers come to 64 KB, over the 48 KB a block may declare statically.
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

//Fixed ascending order over the slices, so the sum does not depend on how the device
//scheduled them. An atomic accumulation in the kernel above would have been shorter and
//would have made the result differ run to run.
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

//Enough blocks to occupy a device without splitting k so finely that the reduction and the
//per-slice tails cost more than the parallelism buys. A large GEMM already fills the device
//and comes back with one slice.
inline int split_count(const int m, const int n, const int k)
{
	constexpr tile_config t = TILE;
	const long long tiles = (long long)((m + t.BM - 1) / t.BM) * ((n + t.BN - 1) / t.BN);
	long long s = 512 / (tiles > 0 ? tiles : 1);
	const long long by_depth = k / (4 * t.BK);   //keep each slice worth staging
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
