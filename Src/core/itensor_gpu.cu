#include "itensor_gpu.h"
#include "gpu_backend.h"
#include "itensor_gemm.cuh"
#include "throughput.h"
#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>

#define SF_TWO_PI 6.283185307179586476925286766559
#define SF_INV_TWO_PI 0.15915494309189533576888376337251

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 I tensor GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

namespace {

//One set of device buffers per scalar type. Only one is ever live, and which one is a
//run-time choice, so both instantiations exist and g_fp64 says which to talk to.
template <typename T>
struct Dev {
	bool ready = false;
	int nmo = 0, packed = 0, n_grids = 0, n_blocks = 0;
	long long n_points = 0;
	//Uploaded once
	T* ao = nullptr;
	int* aos = nullptr;
	unsigned char* skip = nullptr;
	double *d1 = nullptr, *d2 = nullptr, *d3 = nullptr, *w = nullptr;
	//Per reflection scratch
	T *phase_re = nullptr, *phase_im = nullptr;
	//The real and imaginary weighted copies live in one allocation, back to back, and so do
	//the two results. That is what lets the pair of GEMMs below be a single wider one: the
	//imaginary half is simply the next na columns of the same column-major matrix.
	T* wri = nullptr;
	T* cri = nullptr;
	void* gemm_ws = nullptr;   //split-k partials, whichever GEMM is compiled in
	double* I_re = nullptr;
	double* I_im = nullptr;
	//Host copies of the small layout arrays, so the loop stays on the host
	std::vector<int> blk_grid, blk_ps, blk_pc, blk_na, grid_off;
	std::vector<long long> blk_ao_off, blk_aos_off;
};

template <typename T> Dev<T> g;
bool g_fp64 = false;

//The single-precision path keeps the reduced-argument trick the transform uses: the phase
//and its reduction stay in double and only the transcendental drops. In double there is
//nothing to trade, so it takes the argument as it stands.
template <typename T>
__device__ inline void phase_sincos(const double frac, T* s, T* c);

template <>
__device__ inline void phase_sincos<float>(const double frac, float* s, float* c)
{
	sincospif(2.0f * (float)frac, s, c);
}

template <>
__device__ inline void phase_sincos<double>(const double frac, double* s, double* c)
{
	sincospi(2.0 * frac, s, c);
}

//The weight is folded in here so the GEMM operand is exactly what the CPU path multiplies.
template <typename T>
__global__ void phase_kernel(const long long n, const double kx, const double ky, const double kz,
	const double* __restrict__ d1, const double* __restrict__ d2, const double* __restrict__ d3,
	const double* __restrict__ w, T* __restrict__ pre, T* __restrict__ pim)
{
	const long long p = (long long)blockIdx.x * blockDim.x + threadIdx.x;
	if (p >= n) return;
	//kx..kz arrive already divided by 2pi, so t is in turns and sincospi wants 2*frac
	const double t = kx * d1[p] + ky * d2[p] + kz * d3[p];
	const double frac = t - rint(t);
	T s, c;
	phase_sincos<T>(frac, &s, &c);
	const double wp = w[p];
	pre[p] = (T)(wp * (double)c);
	pim[p] = (T)(wp * (double)s);
}

//ao is row-major n_active x np; produce the two phase-weighted copies the GEMM consumes
template <typename T>
__global__ void weight_kernel(const int n_active, const int np, const long long ao_off,
	const int point_base, const T* __restrict__ ao,
	const T* __restrict__ pre, const T* __restrict__ pim,
	T* __restrict__ wre, T* __restrict__ wim)
{
	const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
	const long long total = (long long)n_active * np;
	if (idx >= total) return;
	const int p = (int)(idx % np);
	const T a = ao[ao_off + idx];
	wre[idx] = a * pre[point_base + p];
	wim[idx] = a * pim[point_base + p];
}

//C is symmetric, so only the upper triangle is read back. Distinct (i,j) map to distinct
//packed indices, and launches are ordered on one stream, so no atomics are needed.
template <typename T>
__global__ void accumulate_kernel(const int n_active, const int nmo, const long long aos_off,
	const int* __restrict__ aos, const unsigned char* __restrict__ skip,
	const T* __restrict__ cre, const T* __restrict__ cim,
	const double fre, const double fim,
	double* __restrict__ I_re, double* __restrict__ I_im)
{
	const int i = blockIdx.y * blockDim.y + threadIdx.y;
	const int j = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= n_active || j >= n_active || j < i) return;
	const int mu = aos[aos_off + i];
	const int nu = aos[aos_off + j];
	if (skip[(long long)mu * nmo + nu]) return;
	//The GEMM wrote column-major n_active x n_active; C is symmetric so either index works
	const double re = (double)cre[(long long)j * n_active + i];
	const double im = (double)cim[(long long)j * n_active + i];
	const long long t = (long long)mu * nmo - ((long long)mu * (mu - 1)) / 2 + (nu - mu);
	I_re[t] += re * fre - im * fim;
	I_im[t] += re * fim + im * fre;
}

__global__ void zero_kernel(const int n, double* a, double* b)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n) { a[i] = 0.0; b[i] = 0.0; }
}

//The host holds the AO values in double whatever precision the device will run in. Copying
//through a bounded staging buffer rather than narrowing the whole array first keeps the
//transient allocation independent of the grid size, which on a protein is the difference
//that matters.
template <typename T>
bool upload_ao(T* dst, const double* src, const long long n)
{
	constexpr long long CHUNK = 1 << 22;
	std::vector<T> stage((size_t)std::min<long long>(n, CHUNK));
	for (long long off = 0; off < n; off += CHUNK) {
		const long long m = std::min<long long>(CHUNK, n - off);
		for (long long i = 0; i < m; i++) stage[(size_t)i] = (T)src[off + i];
		if (gpuMemcpy(dst + off, stage.data(), sizeof(T) * (size_t)m,
			gpuMemcpyHostToDevice) != gpuSuccess) return false;
	}
	return true;
}

template <typename T>
void free_impl()
{
	Dev<T>& d = g<T>;
	gpuFree(d.ao); gpuFree(d.aos); gpuFree(d.skip);
	gpuFree(d.d1); gpuFree(d.d2); gpuFree(d.d3); gpuFree(d.w);
	gpuFree(d.phase_re); gpuFree(d.phase_im);
	gpuFree(d.wri); gpuFree(d.cri); gpuFree(d.gemm_ws);
	gpuFree(d.I_re); gpuFree(d.I_im);
	d = Dev<T>{};
}

template <typename T>
bool init_impl(const itensor_gpu_layout& L)
{
	Dev<T>& d = g<T>;
	int max_na = 0;
	long long max_elems = 0;
	size_t max_ws = 0;
	for (int b = 0; b < L.n_blocks; b++) {
		const int na = L.blk_n_active[b];
		const int np = L.blk_point_count[b];
		max_na = std::max(max_na, na);
		//long long: active AOs times points passes 2^31 on a protein-sized grid
		max_elems = std::max(max_elems, (long long)na * np);
		max_ws = std::max(max_ws, itensor_gemm::workspace_bytes<T>(na, 2 * na, np));
	}

	//Block count and shape spread, which is what says whether batching across blocks
	//could ever be one call. Printed only under -gflops.
	if (throughput::enabled()) {
		std::vector<long long> shapes;
		shapes.reserve(L.n_blocks);
		int min_na = L.n_blocks ? L.blk_n_active[0] : 0;
		int min_np = L.n_blocks ? L.blk_point_count[0] : 0, max_np = 0;
		for (int b = 0; b < L.n_blocks; b++) {
			shapes.push_back((long long)L.blk_n_active[b] * 100000 + L.blk_point_count[b]);
			min_na = std::min(min_na, L.blk_n_active[b]);
			min_np = std::min(min_np, L.blk_point_count[b]);
			max_np = std::max(max_np, L.blk_point_count[b]);
		}
		std::sort(shapes.begin(), shapes.end());
		const size_t distinct = std::unique(shapes.begin(), shapes.end()) - shapes.begin();
		//stderr because cout is redirected to the log and moved again mid-run
		std::fprintf(stderr, "I tensor GPU: %d blocks, %zu distinct (n_active, points) shapes,"
		             " n_active %d-%d, points %d-%d\n",
		             L.n_blocks, distinct, min_na, max_na, min_np, max_np);
	}
	const size_t need =
		sizeof(T) * (size_t)L.ao_all_len +
		sizeof(int) * (size_t)L.aos_all_len +
		(size_t)L.nmo * L.nmo +
		sizeof(double) * 4 * (size_t)L.n_points +
		sizeof(T) * 2 * (size_t)L.n_points +
		sizeof(T) * 2 * (size_t)max_elems +
		sizeof(T) * 2 * (size_t)max_na * max_na +
		max_ws +
		sizeof(double) * 2 * (size_t)L.packed;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	if (need + (1u << 26) > freeb) return false;

	GPU_TRY(gpuMalloc(&d.ao, sizeof(T) * (size_t)L.ao_all_len));
	GPU_TRY(gpuMalloc(&d.aos, sizeof(int) * (size_t)L.aos_all_len));
	GPU_TRY(gpuMalloc(&d.skip, (size_t)L.nmo * L.nmo));
	GPU_TRY(gpuMalloc(&d.d1, sizeof(double) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&d.d2, sizeof(double) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&d.d3, sizeof(double) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&d.w, sizeof(double) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&d.phase_re, sizeof(T) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&d.phase_im, sizeof(T) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&d.wri, sizeof(T) * (size_t)max_elems * 2));
	GPU_TRY(gpuMalloc(&d.cri, sizeof(T) * (size_t)max_na * max_na * 2));
	GPU_TRY(gpuMalloc(&d.gemm_ws, max_ws ? max_ws : 1));
	GPU_TRY(gpuMalloc(&d.I_re, sizeof(double) * (size_t)L.packed));
	GPU_TRY(gpuMalloc(&d.I_im, sizeof(double) * (size_t)L.packed));

	if (!upload_ao<T>(d.ao, L.ao_all, L.ao_all_len)) return false;
	GPU_TRY(gpuMemcpy(d.aos, L.aos_all, sizeof(int) * (size_t)L.aos_all_len, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(d.skip, L.skip, (size_t)L.nmo * L.nmo, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(d.d1, L.d1, sizeof(double) * (size_t)L.n_points, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(d.d2, L.d2, sizeof(double) * (size_t)L.n_points, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(d.d3, L.d3, sizeof(double) * (size_t)L.n_points, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(d.w, L.weights, sizeof(double) * (size_t)L.n_points, gpuMemcpyHostToDevice));

	d.nmo = L.nmo; d.packed = L.packed; d.n_grids = L.n_grids; d.n_blocks = L.n_blocks;
	d.n_points = L.n_points;
	d.blk_grid.assign(L.blk_grid, L.blk_grid + L.n_blocks);
	d.blk_ps.assign(L.blk_point_start, L.blk_point_start + L.n_blocks);
	d.blk_pc.assign(L.blk_point_count, L.blk_point_count + L.n_blocks);
	d.blk_na.assign(L.blk_n_active, L.blk_n_active + L.n_blocks);
	d.blk_ao_off.assign(L.blk_ao_off, L.blk_ao_off + L.n_blocks);
	d.blk_aos_off.assign(L.blk_aos_off, L.blk_aos_off + L.n_blocks);
	d.grid_off.assign(L.grid_point_off, L.grid_point_off + L.n_grids + 1);
	d.ready = true;
	return true;
}

template <typename T>
bool reflection_impl(const int num_syms,
	const double* kx, const double* ky, const double* kz,
	const std::complex<double>* factors,
	std::complex<double>* I_r)
{
	Dev<T>& d = g<T>;
	if (!d.ready) return false;
	zero_kernel<<<(d.packed + 255) / 256, 256>>>(d.packed, d.I_re, d.I_im);

	for (int s = 0; s < num_syms; s++) {
		//The CPU path takes sin/cos of k.d directly; scaling k to turns here is what lets
		//the reduction be a rint and the transcendental be sincospi
		phase_kernel<T><<<(unsigned int)((d.n_points + 255) / 256), 256>>>(
			d.n_points, kx[s] * SF_INV_TWO_PI, ky[s] * SF_INV_TWO_PI, kz[s] * SF_INV_TWO_PI,
			d.d1, d.d2, d.d3, d.w, d.phase_re, d.phase_im);
		for (int b = 0; b < d.n_blocks; b++) {
			const int gi = d.blk_grid[b];
			const std::complex<double> f = factors[(size_t)s * d.n_grids + gi];
			if (f.real() == 0.0 && f.imag() == 0.0) continue;
			const int na = d.blk_na[b];
			const int np = d.blk_pc[b];
			const int base = d.grid_off[gi] + d.blk_ps[b];
			const long long total = (long long)na * np;
			//The imaginary half starts exactly where the real one ends, so the two together
			//are one column-major np x 2na matrix and the kernel needs no change
			weight_kernel<T><<<(unsigned int)((total + 255) / 256), 256>>>(
				na, np, d.blk_ao_off[b], base, d.ao, d.phase_re, d.phase_im,
				d.wri, d.wri + total);
			//C = A * W^T, both row-major n_active x np, i.e. column-major np x n_active.
			//One GEMM of width 2na, not two of width na: A is shared and the shape is thin
			//enough to be limited by reading the operands, so widening raises flops per byte.
			itensor_gemm::run<T>(na, 2 * na, np,
				d.ao + d.blk_ao_off[b], np, d.wri, np, d.cri, na, d.gemm_ws);
			const dim3 thr(16, 16);
			const dim3 blk((na + 15) / 16, (na + 15) / 16);
			accumulate_kernel<T><<<blk, thr>>>(na, d.nmo, d.blk_aos_off[b], d.aos, d.skip,
				d.cri, d.cri + (long long)na * na, f.real(), f.imag(), d.I_re, d.I_im);
		}
	}
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuDeviceSynchronize());

	std::vector<double> hr((size_t)d.packed), hi((size_t)d.packed);
	GPU_TRY(gpuMemcpy(hr.data(), d.I_re, sizeof(double) * (size_t)d.packed, gpuMemcpyDeviceToHost));
	GPU_TRY(gpuMemcpy(hi.data(), d.I_im, sizeof(double) * (size_t)d.packed, gpuMemcpyDeviceToHost));
	for (int i = 0; i < d.packed; i++)
		I_r[i] += std::complex<double>(hr[i], hi[i]);
	return true;
}

} //namespace

//Shared with the transform so the "no code for this card" case is diagnosed in one place.
bool itensor_gpu_available() { return sf_gpu_available(); }

bool itensor_gpu_init(const itensor_gpu_layout& L, const sf_precision prec)
{
	itensor_gpu_free();
	if (!itensor_gpu_available()) return false;
	//Auto is not offered here. It would resolve per card, and the I tensor's precision is
	//visible in the reference output, so the same input would produce different logs on
	//different machines. Single precision unless the caller asks for double.
	g_fp64 = (prec == sf_precision::FP64);
	return g_fp64 ? init_impl<double>(L) : init_impl<float>(L);
}

bool itensor_gpu_reflection(const int num_syms,
	const double* kx, const double* ky, const double* kz,
	const std::complex<double>* factors,
	std::complex<double>* I_r)
{
	return g_fp64 ? reflection_impl<double>(num_syms, kx, ky, kz, factors, I_r)
	              : reflection_impl<float>(num_syms, kx, ky, kz, factors, I_r);
}

void itensor_gpu_free()
{
	free_impl<float>();
	free_impl<double>();
}
