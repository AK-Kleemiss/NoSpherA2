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
#endif

//No BLAS library is mapped here on purpose. The GEMMs are in gemm_gpu.cuh, which removed
//half a gigabyte of cuBLAS from the package and, on the AMD side, removed a dependency
//conda-forge does not ship at all.
