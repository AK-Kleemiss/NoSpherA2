#pragma once

//One name for each runtime call the scattering-factor kernel needs, so the same
//source compiles under nvcc and hipcc. HIP is API-compatible here except for the
//fp32:fp64 ratio, which has no equivalent and is inferred from the arch name instead.

#ifdef NOSPHERA2_USE_HIP
#include <hip/hip_runtime.h>
//Nothing here links the HIP runtime. Every call below lands in hip_runtime_shim.cpp, which
//opens libamdhip64 / amdhip64_<major>.dll by name the first time it is asked and answers
//hipErrorNoDevice (or a null handle, or nothing) when there is none. So a binary carrying
//AMD kernels starts on a machine without ROCm, on Linux as on Windows, and the kernel
//registration clang runs before main() comes and goes without a runtime.
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
#define gpuMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define gpuGetLastError hipGetLastError
#define gpuDeviceSynchronize hipDeviceSynchronize
#define gpuStream_t hipStream_t
#define gpuStreamCreate hipStreamCreate
#define gpuStreamDestroy hipStreamDestroy
#define gpuStreamSynchronize hipStreamSynchronize
#define gpuStreamCreateNonBlocking(s) hipStreamCreateWithFlags(s, hipStreamNonBlocking)
#define gpuStreamWaitEvent hipStreamWaitEvent
#define gpuEvent_t hipEvent_t
#define gpuEventCreate(e) hipEventCreateWithFlags(e, hipEventDisableTiming)
#define gpuEventDestroy hipEventDestroy
#define gpuEventRecord hipEventRecord
#define gpuMemcpyAsync hipMemcpyAsync
#define gpuHostAlloc hipHostMalloc
#define gpuFreeHost hipHostFree
//The kernels are written for 32-lane warps (lane = threadIdx.x & 31). A gfx9 wavefront is
//64 wide, so the shuffles are given the width explicitly and act within each 32-lane half;
//HIP's _sync variants want a 64-bit mask and add nothing on AMD hardware anyway. The
//streaming load is clang's non-temporal load, which sets the slc bit like __ldcs does on
//NVIDIA.
#define gpuShflDown32(v, o) __shfl_down((v), (o), 32)
#define gpuShflXor32(v, m) __shfl_xor((v), (m), 32)
#define gpuLoadStreaming(p) __builtin_nontemporal_load(p)
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
#define gpuMemcpyDeviceToDevice cudaMemcpyDeviceToDevice
#define gpuGetLastError cudaGetLastError
#define gpuDeviceSynchronize cudaDeviceSynchronize
#define gpuStream_t cudaStream_t
#define gpuStreamCreate cudaStreamCreate
#define gpuStreamDestroy cudaStreamDestroy
#define gpuStreamSynchronize cudaStreamSynchronize
#define gpuStreamCreateNonBlocking(s) cudaStreamCreateWithFlags(s, cudaStreamNonBlocking)
#define gpuStreamWaitEvent cudaStreamWaitEvent
#define gpuEvent_t cudaEvent_t
#define gpuEventCreate(e) cudaEventCreateWithFlags(e, cudaEventDisableTiming)
#define gpuEventDestroy cudaEventDestroy
#define gpuEventRecord cudaEventRecord
#define gpuMemcpyAsync cudaMemcpyAsync
#define gpuHostAlloc(p, n) cudaHostAlloc((p), (n), cudaHostAllocDefault)
#define gpuFreeHost cudaFreeHost
#define gpuShflDown32(v, o) __shfl_down_sync(0xffffffffu, (v), (o))
#define gpuShflXor32(v, m) __shfl_xor_sync(0xffffffffu, (v), (m))
#define gpuLoadStreaming(p) __ldcs(p)
#endif

//No BLAS library is mapped here on purpose. The GEMMs are in gemm_gpu.cuh, which removed
//half a gigabyte of cuBLAS from the package and, on the AMD side, removed a dependency
//conda-forge does not ship at all.
