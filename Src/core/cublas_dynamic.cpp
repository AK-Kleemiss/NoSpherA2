#include "cublas_dynamic.h"

#include <cstdio>
#include <mutex>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif

//cuBLAS's own headers are deliberately not included. The four entry points below are the
//v2 C API, whose signatures are fixed, so declaring them here keeps the build free of any
//dependency on the toolkit - which is the point: this file must compile on a machine that
//has never seen CUDA.
namespace {

typedef void* cublas_handle;
typedef int cublas_status;          //0 is CUBLAS_STATUS_SUCCESS
enum { OP_N = 0, OP_T = 1 };        //cublasOperation_t
enum { CUDA_R_32F = 0, COMPUTE_32F_FAST_16F = 74, GEMM_DEFAULT_TENSOR_OP = 99 };

typedef cublas_status(*fn_create)(cublas_handle*);
typedef cublas_status(*fn_destroy)(cublas_handle);
typedef cublas_status(*fn_sgemm)(cublas_handle, int, int, int, int, int,
	const float*, const float*, int, const float*, int, const float*, float*, int);
typedef cublas_status(*fn_dgemm)(cublas_handle, int, int, int, int, int,
	const double*, const double*, int, const double*, int, const double*, double*, int);
typedef cublas_status(*fn_gemmex)(cublas_handle, int, int, int, int, int,
	const void*, const void*, int, int, const void*, int, int, const void*, void*, int, int, int, int);

struct Lib {
	bool tried = false;
	bool ok = false;
	fn_create create = nullptr;
	fn_destroy destroy = nullptr;
	fn_sgemm sgemm = nullptr;
	fn_dgemm dgemm = nullptr;
	fn_gemmex gemmex = nullptr;
	cublas_handle handle = nullptr;
};

Lib g_lib;
std::mutex g_mutex;
bool g_enabled = false;

void* open_library()
{
	//Only the major this binary's runtime was built against. Mixing a cuBLAS from one CUDA
	//major with a runtime from another is unsupported, and the failure is not a clean one:
	//the load succeeds and the first call misbehaves.
#ifdef NOSPHERA2_CUDART_MAJOR
	const int major = NOSPHERA2_CUDART_MAJOR;
#else
	const int major = 0;
#endif
	if (major <= 0) return nullptr;
	char name[64];
#ifdef _WIN32
	std::snprintf(name, sizeof(name), "cublas64_%d.dll", major);
	return (void*)LoadLibraryA(name);
#else
	std::snprintf(name, sizeof(name), "libcublas.so.%d", major);
	return dlopen(name, RTLD_NOW | RTLD_LOCAL);
#endif
}

void* symbol(void* lib, const char* name)
{
#ifdef _WIN32
	return (void*)GetProcAddress((HMODULE)lib, name);
#else
	return dlsym(lib, name);
#endif
}

//Caller holds g_mutex
bool load_locked()
{
	if (g_lib.tried) return g_lib.ok;
	g_lib.tried = true;
	void* lib = open_library();
	if (!lib) return false;
	g_lib.create = (fn_create)symbol(lib, "cublasCreate_v2");
	g_lib.destroy = (fn_destroy)symbol(lib, "cublasDestroy_v2");
	g_lib.sgemm = (fn_sgemm)symbol(lib, "cublasSgemm_v2");
	g_lib.dgemm = (fn_dgemm)symbol(lib, "cublasDgemm_v2");
	g_lib.gemmex = (fn_gemmex)symbol(lib, "cublasGemmEx");
	if (!g_lib.create || !g_lib.destroy || !g_lib.sgemm || !g_lib.dgemm) return false;
	if (g_lib.create(&g_lib.handle) != 0) return false;
	g_lib.ok = true;
	return true;
}

} //namespace

void cublas_dynamic_set_enabled(bool on) { g_enabled = on; }
bool cublas_dynamic_enabled() { return g_enabled; }

bool cublas_dynamic_available()
{
	if (!g_enabled) return false;
	std::lock_guard<std::mutex> lock(g_mutex);
	return load_locked();
}

bool cublas_dynamic_fast_16f_available()
{
	if (!g_enabled) return false;
	std::lock_guard<std::mutex> lock(g_mutex);
	return load_locked() && g_lib.gemmex;
}

bool cublas_dynamic_gemm(const bool transA, const bool transB,
	const int m, const int n, const int k,
	const float alpha, const float* A, const int lda, const float* B, const int ldb,
	const float beta, float* C, const int ldc)
{
	if (!g_enabled) return false;
	std::lock_guard<std::mutex> lock(g_mutex);
	if (!load_locked()) return false;
	return g_lib.sgemm(g_lib.handle, transA ? OP_T : OP_N, transB ? OP_T : OP_N,
		m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc) == 0;
}

bool cublas_dynamic_gemm_fast_16f(const bool transA, const bool transB,
	const int m, const int n, const int k,
	const float alpha, const float* A, const int lda, const float* B, const int ldb,
	const float beta, float* C, const int ldc)
{
	if (!g_enabled) return false;
	std::lock_guard<std::mutex> lock(g_mutex);
	if (!load_locked() || !g_lib.gemmex) return false;
	return g_lib.gemmex(g_lib.handle, transA ? OP_T : OP_N, transB ? OP_T : OP_N,
		m, n, k, &alpha, A, CUDA_R_32F, lda, B, CUDA_R_32F, ldb, &beta, C,
		CUDA_R_32F, ldc, COMPUTE_32F_FAST_16F, GEMM_DEFAULT_TENSOR_OP) == 0;
}

bool cublas_dynamic_gemm(const bool transA, const bool transB,
	const int m, const int n, const int k,
	const double alpha, const double* A, const int lda, const double* B, const int ldb,
	const double beta, double* C, const int ldc)
{
	if (!g_enabled) return false;
	std::lock_guard<std::mutex> lock(g_mutex);
	if (!load_locked()) return false;
	return g_lib.dgemm(g_lib.handle, transA ? OP_T : OP_N, transB ? OP_T : OP_N,
		m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc) == 0;
}
