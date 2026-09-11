#pragma once

// Minimal stand-in for HIP's runtime header, covering only what sf_gpu.cu uses.
//
// HIP on an NVIDIA host is itself a thin header translation onto CUDA, so compiling
// the HIP branch this way with nvcc checks that the branch's names resolve, that the
// device code is valid, and that gpu_backend.h is symmetric. It does NOT check AMD
// code generation, the AMD device math library, or that any of it runs: there is no
// AMD card here. Treat a pass as "the HIP path is not obviously broken", nothing more.

#include <cuda_runtime.h>
#include <cstring>

#define hipError_t cudaError_t
#define hipSuccess cudaSuccess
#define hipGetErrorString cudaGetErrorString
#define hipGetDeviceCount cudaGetDeviceCount
#define hipGetDevice cudaGetDevice
#define hipMemGetInfo cudaMemGetInfo
#define hipMalloc cudaMalloc
#define hipFree cudaFree
#define hipMemset cudaMemset
#define hipMemcpy cudaMemcpy
#define hipMemcpyHostToDevice cudaMemcpyHostToDevice
#define hipMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define hipGetLastError cudaGetLastError
#define hipDeviceSynchronize cudaDeviceSynchronize
#define hipStream_t cudaStream_t
#define hipStreamCreate cudaStreamCreate
#define hipStreamDestroy cudaStreamDestroy
#define hipStreamSynchronize cudaStreamSynchronize
#define hipMemcpyAsync cudaMemcpyAsync
#define hipHostMalloc(p, n) cudaHostAlloc((p), (n), cudaHostAllocDefault)
#define hipHostFree cudaFreeHost

// HIP carries the architecture name that the fp64-ratio heuristic reads; CUDA has no
// such field, so it is modelled here rather than aliased.
struct hipDeviceProp_t
{
	char gcnArchName[256];
};

inline cudaError_t hipGetDeviceProperties(hipDeviceProp_t* p, int dev)
{
	cudaDeviceProp cp;
	const cudaError_t e = cudaGetDeviceProperties(&cp, dev);
	if (e == cudaSuccess)
		std::strncpy(p->gcnArchName, "gfx000-shim", sizeof(p->gcnArchName) - 1);
	return e;
}
