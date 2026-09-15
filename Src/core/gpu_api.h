#pragma once

//The GPU entry points (sf_gpu.h, itensor_gpu.h, grid_gpu.h, salted_gpu.h, aux_density_gpu.h,
//blas_gpu.h) are declared between these markers and the kernel sources put their whole body
//between them. In a single-backend build the markers are empty and the functions global. A
//build carrying both CUDA and HIP kernels compiles each kernel source twice, once per backend
//with NOSPHERA2_GPU_BACKEND_NS = nosphera2_cuda or nosphera2_hip, and the markers become that
//namespace; gpu_dispatch.cpp supplies the global names. The argument types (sf_precision,
//itensor_gpu_layout, salted_gpu_problem) stay global and pass through unchanged.
#if defined(NOSPHERA2_GPU_BACKEND_NS)
#define NOSPHERA2_GPU_API_BEGIN namespace NOSPHERA2_GPU_BACKEND_NS {
#define NOSPHERA2_GPU_API_END }
#else
#define NOSPHERA2_GPU_API_BEGIN
#define NOSPHERA2_GPU_API_END
#endif
