#pragma once

#include "gemm_gpu.cuh"
#include "cublas_dynamic.h"
#include <type_traits>

//The I tensor's hot-path GEMM: C = A^T B, column-major, alpha 1, beta 0, m and n ~100, k in the thousands.
//Chosen at run time per precision: cuBLAS (opened by name, never linked) whenever it loads, else CUTLASS
//(compiled-in headers, CUDA only) for float and the built-in kernel for double and HIP. -no_gpu_cublas skips cuBLAS

#if defined(NOSPHERA2_HAVE_CUTLASS)
#include "cutlass/cutlass.h"
#include "cutlass/gemm/device/gemm_splitk_parallel.h"
#include "cutlass/layout/matrix.h"
#endif

//Logged because the backends differ in the last digits
#if defined(NOSPHERA2_HAVE_CUTLASS)
#define NOSPHERA2_ITENSOR_GEMM_NAME "CUTLASS"
#else
#define NOSPHERA2_ITENSOR_GEMM_NAME "built-in"
#endif

NOSPHERA2_GPU_API_BEGIN
namespace itensor_gemm {

inline bool& tensor_mode()
{
	static bool on = false;
	return on;
}

inline void set_tensor_mode(const bool on) { tensor_mode() = on; }

#if defined(NOSPHERA2_HAVE_CUTLASS)

//op(A) (m x k, (i,l) at i*lda + l) is row-major over the caller's column-major buffer; op(B) and C column-major
template <typename T>
using Gemm = cutlass::gemm::device::GemmSplitKParallel<
	T, cutlass::layout::RowMajor,
	T, cutlass::layout::ColumnMajor,
	T, cutlass::layout::ColumnMajor>;

//~192 elements of depth per slice: too few slices starve the device, too many drown it in reduction
inline int split_count(const int k)
{
	int s = k / 192;
	if (s < 1) s = 1;
	if (s > 64) s = 64;
	return s;
}

//float only: in double CUTLASS is slower than the built-in kernel
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
	if constexpr (std::is_same<T, float>::value)
		if (tensor_mode() && cublas_dynamic_gemm_fast_16f(true, false, m, n, k,
			T(1), A, lda, B, ldb, T(0), C, ldc)) return true;
	if (cublas_dynamic_gemm(true, false, m, n, k, T(1), A, lda, B, ldb, T(0), C, ldc))
		return true;
	if constexpr (cutlass_handles<T>()) {
		Gemm<T> op;
		typename Gemm<T>::Arguments args({m, n, k}, {A, lda}, {B, ldb},
			{C, ldc}, {C, ldc}, {T(1), T(0)}, split_count(k));
		//Not operator()(args, workspace): in CUTLASS 4.7.1 it passes a stream to a two-parameter initialize() and does not compile
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
	//On HIP there is no cuBLAS to open, so the cuBLAS calls return false
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
NOSPHERA2_GPU_API_END
