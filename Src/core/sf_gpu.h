#pragma once

#include "gpu_api.h"

//Auto: reduced-argument fp32 sincos where fp64 is slow, fp64 on datacentre parts. Tests pin FP32/FP64, since under
//Auto an fp32 test on a V100 (ratio 2) would run the fp64 kernel
enum class sf_precision { Auto, FP32, FP64 };

NOSPHERA2_GPU_API_BEGIN

//calc_SF's non-uniform DFT on CUDA or HIP; false without a device or when it does not fit, the caller runs the CPU loop
bool sf_gpu_available();
//Context creation overlaps the grid work instead of landing inside the measurement
void sf_gpu_warmup_start();
void sf_gpu_warmup_wait();
//"CUDA" or "HIP", for the log line
const char* sf_gpu_backend();
//Single-to-double throughput ratio: 2 on datacentre parts, 32-64 on consumer ones
int sf_gpu_fp64_ratio();
//What Auto resolves to here; the log line uses it so it names the kernel that ran
bool sf_gpu_uses_fp32(const sf_precision prec);
//sf_rows[i]: atom i's 2*smax doubles, interleaved re,im, i.e. the caller's complex<double> storage
bool sf_gpu_run(const int imax, const long long smax,
	const double* k1, const double* k2, const double* k3,
	const double* d1, const double* d2, const double* d3,
	const double* dens, const int* offs, const long long total_points,
	double* const* sf_rows, const sf_precision prec = sf_precision::Auto);

NOSPHERA2_GPU_API_END
