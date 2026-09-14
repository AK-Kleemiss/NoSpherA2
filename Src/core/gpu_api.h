#pragma once

//The GPU entry points (sf_gpu.h, itensor_gpu.h, grid_gpu.h, salted_gpu.h, aux_density_gpu.h,
//blas_gpu.h) are declared between these two markers, and the kernel sources that define them
//put their whole body between the same two. In an ordinary build - one backend, or none - the
//markers are empty and the functions are global, as they always were.
//
//A build that carries both CUDA and HIP kernels in one binary compiles each kernel source
//twice, once per backend, and the two objects would define the same functions. So those
//compiles are given NOSPHERA2_GPU_BACKEND_NS (nosphera2_cuda or nosphera2_hip) and the
//markers become that namespace: every kernel, every function and every piece of static state
//of one backend lives in its own namespace, and gpu_dispatch.cpp supplies the global names
//by forwarding to whichever backend has a device. Host code never sees the namespaces.
//
//The types the functions take (sf_precision, itensor_gpu_layout, salted_gpu_problem) stay
//global: they are the same on both sides and the forwarder passes them through.
#if defined(NOSPHERA2_GPU_BACKEND_NS)
#define NOSPHERA2_GPU_API_BEGIN namespace NOSPHERA2_GPU_BACKEND_NS {
#define NOSPHERA2_GPU_API_END }
#else
#define NOSPHERA2_GPU_API_BEGIN
#define NOSPHERA2_GPU_API_END
#endif
