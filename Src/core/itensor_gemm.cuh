#pragma once

#include "gemm_gpu.cuh"
#include "cublas_dynamic.h"
#include <type_traits>

//The one GEMM the I tensor performs, behind a name so the implementation can change under
//it: C = A^T B, column-major, alpha 1, beta 0, with m and n around a hundred and k in the
//thousands. Nothing else in the program has this shape, and it is the only device GEMM that
//is on a hot path.
//
//Three implementations, chosen at run time and per precision because the fastest differs
//by both. cuBLAS leads in single and double alike and is used whenever it loads; CUTLASS
//takes single precision behind it and is beaten by the built-in kernel in double, so
//double falls to the built-in one. -no_gpu_cublas pins whichever of the latter two applies.
//
//Neither fallback is a runtime dependency: CUTLASS is headers compiled into this binary
//and the built-in kernel is ours, so a machine with no CUDA toolkit still runs the device
//path. cuBLAS is opened by name and ignored when absent - nothing links it.
//
//Rates per card and per backend are in the vault note, not here.
//
//HIP keeps the hand-written kernel: CUTLASS is CUDA-only.

#if defined(NOSPHERA2_HAVE_CUTLASS)
#include "cutlass/cutlass.h"
#include "cutlass/gemm/device/gemm_splitk_parallel.h"
#include "cutlass/layout/matrix.h"
#endif

//Named so a log says which GEMM produced its numbers; they differ in the last digits.
#if defined(NOSPHERA2_HAVE_CUTLASS)
#define NOSPHERA2_ITENSOR_GEMM_NAME "CUTLASS"
#else
#define NOSPHERA2_ITENSOR_GEMM_NAME "built-in"
#endif

namespace itensor_gemm {

inline bool& tensor_mode()
{
	static bool on = false;
	return on;
}

inline void set_tensor_mode(const bool on) { tensor_mode() = on; }

#if defined(NOSPHERA2_HAVE_CUTLASS)

//op(A) is m x k with element (i,l) at i*lda + l, which is row-major over the column-major
//buffer the caller holds; op(B) and C are column-major as stored.
template <typename T>
using Gemm = cutlass::gemm::device::GemmSplitKParallel<
	T, cutlass::layout::RowMajor,
	T, cutlass::layout::ColumnMajor,
	T, cutlass::layout::ColumnMajor>;

//Slices sized so each carries a couple of hundred elements of depth. The count matters a
//great deal - too few starves the device and too many drowns it in reduction - and the
//useful band was found by sweeping it on an RTX 4090 mobile at k = 4500. Splitting by depth
//rather than by SM count reproduces that band without asking the device anything, though it
//has not been re-swept on a datacentre part.
inline int split_count(const int k)
{
	int s = k / 192;
	if (s < 1) s = 1;
	if (s > 64) s = 64;
	return s;
}

//CUTLASS was adopted on a single-precision measurement and allowed to take double with it,
//which halved that path until someone ran it. The choice is per precision because the
//answer is per precision.
template <typename T>
constexpr bool cutlass_handles() { return std::is_same<T, float>::value; }

template <typename T>
inline size_t workspace_bytes(const int m, const int n, const int k)
{
	if constexpr (cutlass_handles<T>()) {
		typename Gemm<T>::Arguments args({m, n, k}, {nullptr, k}, {nullptr, k},
			{nullptr, m}, {nullptr, m}, {T(1), T(0)}, split_count(k));
		return Gemm<T>::get_workspace_size(args);
	} else {
		return sizeof(T) * gemm_gpu::workspace_elems(m, n, gemm_gpu::split_count(m, n, k));
	}
}

template <typename T>
inline bool run(const int m, const int n, const int k,
	const T* A, const int lda, const T* B, const int ldb, T* C, const int ldc, void* ws)
{
	//cuBLAS first when it loaded - it is the fastest of the three in both precisions.
	//-no_gpu_cublas pins whichever of the other two handles this one.
	if constexpr (std::is_same<T, float>::value)
		if (tensor_mode() && cublas_dynamic_gemm_fast_16f(true, false, m, n, k,
			T(1), A, lda, B, ldb, T(0), C, ldc)) return true;
	if (cublas_dynamic_gemm(true, false, m, n, k, T(1), A, lda, B, ldb, T(0), C, ldc))
		return true;
	if constexpr (cutlass_handles<T>()) {
		Gemm<T> op;
		typename Gemm<T>::Arguments args({m, n, k}, {A, lda}, {B, ldb},
			{C, ldc}, {C, ldc}, {T(1), T(0)}, split_count(k));
		//initialize + run rather than operator()(args, workspace): in 4.7.1 that overload
		//hands a stream to an initialize() declared with two parameters, and will not compile.
		if (op.initialize(args, ws) != cutlass::Status::kSuccess) return false;
		return op.run() == cutlass::Status::kSuccess;
	} else {
		return gemm_gpu::launch<T>(true, false, m, n, k,
			T(1), A, lda, B, ldb, T(0), C, ldc, static_cast<T*>(ws));
	}
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
	//Same order as the CUTLASS build: cuBLAS when asked for and loadable, ours otherwise.
	//On a HIP build the first call is always false, there being no cuBLAS to open.
	if constexpr (std::is_same<T, float>::value)
		if (tensor_mode() && cublas_dynamic_gemm_fast_16f(true, false, m, n, k,
			T(1), A, lda, B, ldb, T(0), C, ldc)) return true;
	if (cublas_dynamic_gemm(true, false, m, n, k, T(1), A, lda, B, ldb, T(0), C, ldc))
		return true;
	return gemm_gpu::launch<T>(true, false, m, n, k,
		T(1), A, lda, B, ldb, T(0), C, ldc, static_cast<T*>(ws));
}

#endif

} //namespace itensor_gemm
