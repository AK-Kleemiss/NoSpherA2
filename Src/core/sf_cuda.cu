#include "sf_cuda.h"
#include <cuda_runtime.h>
#include <cstdio>
#include <vector>

//One block per k-point; threads in the block split that atom's grid points and reduce
__global__ void sf_kernel(const int imax, const long long smax,
	const double* __restrict__ k1, const double* __restrict__ k2, const double* __restrict__ k3,
	const double* __restrict__ d1, const double* __restrict__ d2, const double* __restrict__ d3,
	const double* __restrict__ dens, const int* __restrict__ offs,
	double* __restrict__ sf_re, double* __restrict__ sf_im)
{
	const long long s = blockIdx.x;
	if (s >= smax) return;
	const double kx = k1[s], ky = k2[s], kz = k3[s];
	extern __shared__ double sh[];
	double* sh_re = sh;
	double* sh_im = sh + blockDim.x;
	for (int i = 0; i < imax; i++) {
		const int lo = offs[i], hi = offs[i + 1];
		double re = 0.0, im = 0.0;
		for (int p = lo + threadIdx.x; p < hi; p += blockDim.x) {
			const double w = kx * d1[p] + ky * d2[p] + kz * d3[p];
			double si, co;
			sincos(w, &si, &co);
			const double r = dens[p];
			re += r * co;
			im += r * si;
		}
		sh_re[threadIdx.x] = re;
		sh_im[threadIdx.x] = im;
		__syncthreads();
		for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
			if (threadIdx.x < stride) {
				sh_re[threadIdx.x] += sh_re[threadIdx.x + stride];
				sh_im[threadIdx.x] += sh_im[threadIdx.x + stride];
			}
			__syncthreads();
		}
		if (threadIdx.x == 0) {
			sf_re[(long long)i * smax + s] = sh_re[0];
			sf_im[(long long)i * smax + s] = sh_im[0];
		}
		__syncthreads();
	}
}

bool sf_cuda_available()
{
	int n = 0;
	return cudaGetDeviceCount(&n) == cudaSuccess && n > 0;
}

#define CUDA_TRY(call) do { const cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
	std::fprintf(stderr, "NoSpherA2 CUDA: %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

bool sf_cuda_run(const int imax, const long long smax,
	const double* k1, const double* k2, const double* k3,
	const double* d1, const double* d2, const double* d3,
	const double* dens, const int* offs, const long long total_points,
	double* sf_re, double* sf_im)
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
	const int threads = 256;
	const size_t shmem = sizeof(double) * 2 * threads;
	sf_kernel<<<(unsigned int)smax, threads, shmem>>>(imax, smax, dk1, dk2, dk3, dd1, dd2, dd3, dde, dof, dsr, dsi);
	CUDA_TRY(cudaGetLastError());
	CUDA_TRY(cudaDeviceSynchronize());
	CUDA_TRY(cudaMemcpy(sf_re, dsr, out, cudaMemcpyDeviceToHost));
	CUDA_TRY(cudaMemcpy(sf_im, dsi, out, cudaMemcpyDeviceToHost));
	cudaFree(dk1); cudaFree(dk2); cudaFree(dk3);
	cudaFree(dd1); cudaFree(dd2); cudaFree(dd3); cudaFree(dde);
	cudaFree(dof); cudaFree(dsr); cudaFree(dsi);
	return true;
}
