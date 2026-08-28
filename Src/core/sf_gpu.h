#pragma once

//GPU path for calc_SF's non-uniform DFT, built against CUDA or HIP. All return false
//if no device is present or the problem will not fit, and the caller falls back to the CPU loop.
bool sf_gpu_available();
//Device context creation costs ~95 ms. Start it when grid work begins and wait for it
//just before the transform, so it overlaps instead of landing inside the measurement.
void sf_gpu_warmup_start();
void sf_gpu_warmup_wait();
//"CUDA" or "HIP", for the log line
const char* sf_gpu_backend();
//Single-to-double throughput ratio: 2 on datacentre parts, 32-64 on consumer ones
int sf_gpu_fp64_ratio();
//Which sincos the transform uses. Auto asks the card, which is what a normal run wants: a
//part where fp64 runs at a small fraction of fp32 gets the reduced-argument f32 kernel, a
//datacentre part keeps the double one. The two explicit values exist because a test cannot
//use Auto - on a V100, where the ratio is 2, an "fp32 test" would quietly exercise the fp64
//kernel and pass without ever running the code it is named after.
enum class sf_precision { Auto, FP32, FP64 };
//What Auto resolves to on this card. The log line uses it so the message cannot disagree
//with the kernel that ran.
bool sf_gpu_uses_fp32(const sf_precision prec);
//sf_rows[i] points at atom i's row of 2*smax doubles, interleaved re,im - i.e. straight
//into the caller's complex<double> storage, which the standard guarantees is laid out that way.
bool sf_gpu_run(const int imax, const long long smax,
	const double* k1, const double* k2, const double* k3,
	const double* d1, const double* d2, const double* d3,
	const double* dens, const int* offs, const long long total_points,
	double* const* sf_rows, const sf_precision prec = sf_precision::Auto);
