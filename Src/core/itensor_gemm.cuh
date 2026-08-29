#pragma once

#include "gemm_gpu.cuh"

//The one GEMM the I tensor performs, behind a name so the implementation can change under
//it: C = A^T B, column-major, alpha 1, beta 0, with m and n around a hundred and k in the
//thousands. Nothing else in the program has this shape, and it is the only device GEMM that
//is on a hot path.
//
//CUTLASS where it is available, our own kernel otherwise. Isolated on this shape alone,
//one process so the clocks are common, best split count each, as a fraction of cuBLAS:
//ours 0.55-0.61 and flat in the split count, CUTLASS 0.85-1.12. Our kernel being flat is
//what settled it - the heuristic was already within 2% of the best count available to it,
//so the ceiling was the kernel and no further tuning of the split rule could have moved it.
//
//Whole I tensor stage, which is what the run actually pays, in GFLOP/s:
//
//                          ours   CUTLASS   cuBLAS
//  RTX 4090 mobile         1759      2061     1764
//  RTX 2080 Ti             1181      1663     1794
//
//So CUTLASS is ahead of the cuBLAS build this replaced on Ada and within 7% of it on
//Turing, for a dependency that is headers only. Nothing is linked and nothing ships: it
//compiles into this binary, so the half-gigabyte problem does not come back.
//
//HIP keeps the hand-written kernel: CUTLASS is CUDA-only.

#if defined(NOSPHERA2_HAVE_CUTLASS)
#include "cutlass/cutlass.h"
#include "cutlass/gemm/device/gemm_splitk_parallel.h"
#include "cutlass/layout/matrix.h"
#endif

namespace itensor_gemm {

#if defined(NOSPHERA2_HAVE_CUTLASS)

//op(A) is m x k with element (i,l) at i*lda + l, which is row-major over the column-major
//buffer the caller holds; op(B) and C are column-major as stored.
template <typename T>
using Gemm = cutlass::gemm::device::GemmSplitKParallel<
	T, cutlass::layout::RowMajor,
	T, cutlass::layout::ColumnMajor,
	T, cutlass::layout::ColumnMajor>;

//Slices sized so each carries a few hundred elements of depth. Measured on an RTX 4090
//mobile at k = 4500: 8 slices reaches 0.44 of cuBLAS, 20-32 stays between 0.85 and 1.12,
//and 48 falls back to 0.47. Splitting by depth rather than by SM count keeps that shape
//without asking the device anything.
inline int split_count(const int k)
{
	int s = k / 192;
	if (s < 1) s = 1;
	if (s > 64) s = 64;
	return s;
}

template <typename T>
inline size_t workspace_bytes(const int m, const int n, const int k)
{
	typename Gemm<T>::Arguments args({m, n, k}, {nullptr, k}, {nullptr, k},
		{nullptr, m}, {nullptr, m}, {T(1), T(0)}, split_count(k));
	return Gemm<T>::get_workspace_size(args);
}

template <typename T>
inline bool run(const int m, const int n, const int k,
	const T* A, const int lda, const T* B, const int ldb, T* C, const int ldc, void* ws)
{
	Gemm<T> op;
	typename Gemm<T>::Arguments args({m, n, k}, {A, lda}, {B, ldb},
		{C, ldc}, {C, ldc}, {T(1), T(0)}, split_count(k));
	//initialize + run rather than operator()(args, workspace): in 4.7.1 that overload hands
	//a stream to an initialize() declared with two parameters, and does not compile.
	if (op.initialize(args, ws) != cutlass::Status::kSuccess) return false;
	return op.run() == cutlass::Status::kSuccess;
}

#else

template <typename T>
inline size_t workspace_bytes(const int m, const int n, const int k)
{
	return sizeof(T) * gemm_gpu::workspace_elems(m, n, gemm_gpu::split_count(m, n, k));
}

template <typename T>
inline bool run(const int m, const int n, const int k,
	const T* A, const int lda, const T* B, const int ldb, T* C, const int ldc, void* ws)
{
	return gemm_gpu::launch<T>(true, false, m, n, k,
		T(1), A, lda, B, ldb, T(0), C, ldc, static_cast<T*>(ws));
}

#endif

} //namespace itensor_gemm
