#include "sf_cuda.h"
#include <cuda_runtime.h>
#include <cstdio>

//Each block owns a tile of k-points for one atom and streams that atom's grid points
//through shared memory. F32 selects the reduced-argument single-precision sincos, which
//consumer cards run 32-64x faster than the double one; the phase itself stays double.
#define SF_TILE_K 128
#define SF_CHUNK 256
#define SF_TWO_PI 6.283185307179586476925286766559
#define SF_INV_TWO_PI 0.15915494309189533576888376337251

template <bool F32>
__global__ void sf_kernel(const int imax, const long long smax,
	const double* __restrict__ k1, const double* __restrict__ k2, const double* __restrict__ k3,
	const double* __restrict__ d1, const double* __restrict__ d2, const double* __restrict__ d3,
	const double* __restrict__ dens, const int* __restrict__ offs,
	double* __restrict__ sf_re, double* __restrict__ sf_im)
{
	__shared__ double s1[SF_CHUNK], s2[SF_CHUNK], s3[SF_CHUNK], sd[SF_CHUNK];
	const int ia = blockIdx.y;
	const long long s = (long long)blockIdx.x * SF_TILE_K + threadIdx.x;
	const bool live = (s < smax);
	const double kx = live ? k1[s] : 0.0;
	const double ky = live ? k2[s] : 0.0;
	const double kz = live ? k3[s] : 0.0;
	double re = 0.0, im = 0.0;
	const int lo = offs[ia], hi = offs[ia + 1];
	for (int base = lo; base < hi; base += SF_CHUNK) {
		const int n = min(SF_CHUNK, hi - base);
		for (int t = threadIdx.x; t < n; t += blockDim.x) {
			s1[t] = d1[base + t];
			s2[t] = d2[base + t];
			s3[t] = d3[base + t];
			sd[t] = dens[base + t];
		}
		__syncthreads();
		if (live) {
			for (int p = 0; p < n; p++) {
				const double w = kx * s1[p] + ky * s2[p] + kz * s3[p];
				const double r = sd[p];
				if (F32) {
					//Reduce in double first: sinf/cosf of a large argument is meaningless
					const double wr = w - SF_TWO_PI * floor(w * SF_INV_TWO_PI + 0.5);
					float sif, cof;
					sincosf((float)wr, &sif, &cof);
					re += r * (double)cof;
					im += r * (double)sif;
				}
				else {
					double si, co;
					sincos(w, &si, &co);
					re += r * co;
					im += r * si;
				}
			}
		}
		__syncthreads();
	}
	if (live) {
		sf_re[(long long)ia * smax + s] = re;
		sf_im[(long long)ia * smax + s] = im;
	}
}

bool sf_cuda_available()
{
	int n = 0;
	return cudaGetDeviceCount(&n) == cudaSuccess && n > 0;
}

//Single-to-double throughput ratio as the driver reports it: 2 on datacentre parts,
//32 or 64 on consumer ones. Above the threshold the double sincos is not worth paying for.
int sf_cuda_fp64_ratio()
{
	int dev = 0, ratio = 0;
	if (cudaGetDevice(&dev) != cudaSuccess) return 0;
	if (cudaDeviceGetAttribute(&ratio, cudaDevAttrSingleToDoublePrecisionPerfRatio, dev) != cudaSuccess) return 0;
	return ratio;
}

#define CUDA_TRY(call) do { const cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
	std::fprintf(stderr, "NoSpherA2 CUDA: %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

bool sf_cuda_run(const int imax, const long long smax,
	const double* k1, const double* k2, const double* k3,
	const double* d1, const double* d2, const double* d3,
	const double* dens, const int* offs, const long long total_points,
	double* sf_re, double* sf_im, const bool force_fp64)
{
	const size_t pts = sizeof(double) * (size_t)total_points;
	const size_t kb = sizeof(double) * (size_t)smax;
	const size_t out = sizeof(double) * (size_t)imax * (size_t)smax;
	size_t freeb = 0, totalb = 0;
	if (cudaMemGetInfo(&freeb, &totalb) != cudaSuccess) return false;
	//Bail to the CPU rather than thrash if the whole problem will not sit in device memory
	if (4 * pts + 3 * kb + 2 * out + (1u << 26) > freeb) return false;
	double *dk1 = nullptr, *dk2 = nullptr, *dk3 = nullptr, *dd1 = nullptr, *dd2 = nullptr, *dd3 = nullptr, *dde = nullptr, *dsr = nullptr, *dsi = nullptr;
	int* dof = nullptr;
	CUDA_TRY(cudaMalloc(&dk1, kb)); CUDA_TRY(cudaMalloc(&dk2, kb)); CUDA_TRY(cudaMalloc(&dk3, kb));
	CUDA_TRY(cudaMalloc(&dd1, pts)); CUDA_TRY(cudaMalloc(&dd2, pts)); CUDA_TRY(cudaMalloc(&dd3, pts));
	CUDA_TRY(cudaMalloc(&dde, pts));
	CUDA_TRY(cudaMalloc(&dof, sizeof(int) * (size_t)(imax + 1)));
	CUDA_TRY(cudaMalloc(&dsr, out)); CUDA_TRY(cudaMalloc(&dsi, out));
	CUDA_TRY(cudaMemcpy(dk1, k1, kb, cudaMemcpyHostToDevice));
	CUDA_TRY(cudaMemcpy(dk2, k2, kb, cudaMemcpyHostToDevice));
	CUDA_TRY(cudaMemcpy(dk3, k3, kb, cudaMemcpyHostToDevice));
	CUDA_TRY(cudaMemcpy(dd1, d1, pts, cudaMemcpyHostToDevice));
	CUDA_TRY(cudaMemcpy(dd2, d2, pts, cudaMemcpyHostToDevice));
	CUDA_TRY(cudaMemcpy(dd3, d3, pts, cudaMemcpyHostToDevice));
	CUDA_TRY(cudaMemcpy(dde, dens, pts, cudaMemcpyHostToDevice));
	CUDA_TRY(cudaMemcpy(dof, offs, sizeof(int) * (size_t)(imax + 1), cudaMemcpyHostToDevice));
	const dim3 grid((unsigned int)((smax + SF_TILE_K - 1) / SF_TILE_K), (unsigned int)imax);
	const int ratio = sf_cuda_fp64_ratio();
	const bool f32 = !force_fp64 && ratio > 4;
	if (f32)
		sf_kernel<true><<<grid, SF_TILE_K>>>(imax, smax, dk1, dk2, dk3, dd1, dd2, dd3, dde, dof, dsr, dsi);
	else
		sf_kernel<false><<<grid, SF_TILE_K>>>(imax, smax, dk1, dk2, dk3, dd1, dd2, dd3, dde, dof, dsr, dsi);
	CUDA_TRY(cudaGetLastError());
	CUDA_TRY(cudaDeviceSynchronize());
	CUDA_TRY(cudaMemcpy(sf_re, dsr, out, cudaMemcpyDeviceToHost));
	CUDA_TRY(cudaMemcpy(sf_im, dsi, out, cudaMemcpyDeviceToHost));
	cudaFree(dk1); cudaFree(dk2); cudaFree(dk3);
	cudaFree(dd1); cudaFree(dd2); cudaFree(dd3); cudaFree(dde);
	cudaFree(dof); cudaFree(dsr); cudaFree(dsi);
	return true;
}
