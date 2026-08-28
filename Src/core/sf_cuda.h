#pragma once

//GPU path for calc_SF's non-uniform DFT. Both return false if CUDA is absent,
//the device is missing or the problem will not fit, and the caller falls back to the CPU loop.
bool sf_cuda_available();
//Single-to-double throughput ratio the driver reports: 2 on datacentre parts, 32-64 on consumer
int sf_cuda_fp64_ratio();
bool sf_cuda_run(const int imax, const long long smax,
	const double* k1, const double* k2, const double* k3,
	const double* d1, const double* d2, const double* d3,
	const double* dens, const int* offs, const long long total_points,
	double* sf_re, double* sf_im, const bool force_fp64 = false);
