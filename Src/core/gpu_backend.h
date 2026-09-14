#pragma once

//One name for each runtime call the scattering-factor kernel needs, so the same
//source compiles under nvcc and hipcc. HIP is API-compatible here except for the
//fp32:fp64 ratio, which has no equivalent and is inferred from the arch name instead.

#ifdef NOSPHERA2_USE_HIP
#include <hip/hip_runtime.h>
#define gpuError_t hipError_t
#define gpuSuccess hipSuccess
#define gpuGetErrorString hipGetErrorString
#if defined(_WIN32)
//On Windows the HIP runtime is delay-loaded (NOSPHERA2_HIP_DELAYLOAD in CMakeLists.txt), so
//that a build carrying AMD kernels still starts on a machine without ROCm. Delay-loading
//only postpones the failure, though: the first call into amdhip64 on such a machine raises
//the loader's exception instead of returning an error. Every GPU path begins by counting
//devices, so that one call is answered here, without the runtime, when the DLL is absent.
//(The kernel registration clang runs before main() is the other call that happens without
//being asked for; Src/hip_delayload_hook.cpp takes care of that one.)
//The name is derived from HIP_VERSION_MAJOR exactly as the delay-load flag is, so the two
//cannot disagree. LoadLibraryA is declared by hand rather than through windows.h, whose
//min/max macros this header's includers do not want; the device pass of the compiler sees
//the declaration too and knows no dllimport.
#if defined(__HIP_DEVICE_COMPILE__)
extern "C" void* __stdcall LoadLibraryA(const char*);
#else
extern "C" __declspec(dllimport) void* __stdcall LoadLibraryA(const char*);
#endif
#define NOSPHERA2_GPU_STR_(x) #x
#define NOSPHERA2_GPU_STR(x) NOSPHERA2_GPU_STR_(x)
inline bool nosphera2_hip_runtime_present()
{
	static const bool present =
		LoadLibraryA("amdhip64_" NOSPHERA2_GPU_STR(HIP_VERSION_MAJOR) ".dll") != nullptr;
	return present;
}
inline hipError_t nosphera2_hip_device_count(int* n)
{
	if (!nosphera2_hip_runtime_present()) {
		if (n) *n = 0;
		return hipErrorNoDevice;
	}
	return hipGetDeviceCount(n);
}
#define gpuGetDeviceCount nosphera2_hip_device_count
#else
#define gpuGetDeviceCount hipGetDeviceCount
#endif
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
