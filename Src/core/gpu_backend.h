#pragma once

//One name for each runtime call the scattering-factor kernel needs, so the same
//source compiles under nvcc and hipcc. HIP is API-compatible here except for the
//fp32:fp64 ratio, which has no equivalent and is inferred from the arch name instead.

#ifdef NOSPHERA2_USE_HIP
#include <hip/hip_runtime.h>
#define gpuError_t hipError_t
#define gpuSuccess hipSuccess
#define gpuGetErrorString hipGetErrorString
#define gpuGetDeviceCount hipGetDeviceCount
#define gpuGetDevice hipGetDevice
#define gpuGetDeviceProperties hipGetDeviceProperties
#define gpuDeviceProp_t hipDeviceProp_t
#define gpuMemGetInfo hipMemGetInfo
#define gpuMalloc hipMalloc
#define gpuFree hipFree
#define gpuMemset hipMemset
#define gpuMemcpy hipMemcpy
#define gpuMemcpyHostToDevice hipMemcpyHostToDevice
#define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
#define gpuGetLastError hipGetLastError
#define gpuDeviceSynchronize hipDeviceSynchronize
#define gpuStream_t hipStream_t
#define gpuStreamCreate hipStreamCreate
#define gpuStreamDestroy hipStreamDestroy
#define gpuStreamSynchronize hipStreamSynchronize
#define gpuMemcpyAsync hipMemcpyAsync
#define gpuHostAlloc hipHostMalloc
#define gpuFreeHost hipHostFree
#include <hipblas/hipblas.h>
#define gpublasHandle_t hipblasHandle_t
#define gpublasStatus_t hipblasStatus_t
#define gpublasCreate hipblasCreate
#define gpublasDestroy hipblasDestroy
#define gpublasSgemm hipblasSgemm
#define gpublasDgemm hipblasDgemm
#define gpublasOperation_t hipblasOperation_t
#define GPUBLAS_OP_T HIPBLAS_OP_T
#define GPUBLAS_OP_N HIPBLAS_OP_N
#define GPUBLAS_STATUS_SUCCESS HIPBLAS_STATUS_SUCCESS
#else
#include <cuda_runtime.h>
#define gpuError_t cudaError_t
#define gpuSuccess cudaSuccess
#define gpuGetErrorString cudaGetErrorString
#define gpuGetDeviceCount cudaGetDeviceCount
#define gpuGetDevice cudaGetDevice
#define gpuGetDeviceProperties cudaGetDeviceProperties
#define gpuDeviceProp_t cudaDeviceProp
#define gpuMemGetInfo cudaMemGetInfo
#define gpuMalloc cudaMalloc
#define gpuFree cudaFree
#define gpuMemset cudaMemset
#define gpuMemcpy cudaMemcpy
#define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
#define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define gpuGetLastError cudaGetLastError
#define gpuDeviceSynchronize cudaDeviceSynchronize
#define gpuStream_t cudaStream_t
#define gpuStreamCreate cudaStreamCreate
#define gpuStreamDestroy cudaStreamDestroy
#define gpuStreamSynchronize cudaStreamSynchronize
#define gpuMemcpyAsync cudaMemcpyAsync
#define gpuHostAlloc(p, n) cudaHostAlloc((p), (n), cudaHostAllocDefault)
#define gpuFreeHost cudaFreeHost
#include <cublas_v2.h>
#define gpublasHandle_t cublasHandle_t
#define gpublasStatus_t cublasStatus_t
#define gpublasCreate cublasCreate
#define gpublasDestroy cublasDestroy
#define gpublasSgemm cublasSgemm
#define gpublasDgemm cublasDgemm
#define gpublasOperation_t cublasOperation_t
#define GPUBLAS_OP_T CUBLAS_OP_T
#define GPUBLAS_OP_N CUBLAS_OP_N
#define GPUBLAS_STATUS_SUCCESS CUBLAS_STATUS_SUCCESS
#endif

//Delay-loading the BLAS DLL stops the loader rejecting the process where no GPU runtime
//exists, but on its own it only moves the failure: the first call then raises the
//delay-load helper exception and kills the process. Measured on the AMD node, exit
//0xC06D007E. So probe for the module before calling in, and treat absence as "no GPU".
#ifdef _WIN32
#include <windows.h>
inline bool gpu_blas_runtime_present()
{
#ifdef NOSPHERA2_USE_HIP
	static const bool ok = (LoadLibraryA("hipblas.dll") != nullptr);
#else
	static const bool ok = (LoadLibraryA("cublas64_13.dll") != nullptr)
	                   || (LoadLibraryA("cublas64_12.dll") != nullptr);
#endif
	return ok;
}
#else
inline bool gpu_blas_runtime_present() { return true; }
#endif
