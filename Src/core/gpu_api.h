#pragma once

//GPU entry points are declared, and kernel sources defined, between these markers: empty in a
//single-backend build, the namespace nosphera2_cuda or nosphera2_hip when a CUDA+HIP build compiles
//each kernel source once per backend, with gpu_dispatch.cpp supplying the global names.
//Argument types (sf_precision, itensor_gpu_layout, salted_gpu_problem) stay global.
#if defined(NOSPHERA2_GPU_BACKEND_NS)
#define NOSPHERA2_GPU_API_BEGIN namespace NOSPHERA2_GPU_BACKEND_NS {
#define NOSPHERA2_GPU_API_END }
#else
#define NOSPHERA2_GPU_API_BEGIN
#define NOSPHERA2_GPU_API_END
#endif
