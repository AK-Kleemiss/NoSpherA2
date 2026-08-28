#include "itensor_gpu.h"
#include "gpu_backend.h"
#include <cstdio>
#include <vector>
#include <algorithm>

//The contraction is a GEMM and nothing else, so with no device BLAS there is no path here
//worth having. conda-forge packages hipcc without hipBLAS, which is how a build gets here.
//Returning false from init is exactly what this file already does on a machine with no
//device, and XCW logs "the CPU - device unavailable or problem too large" and carries on.
#ifndef NOSPHERA2_HAVE_GPUBLAS

bool itensor_gpu_available() { return false; }
bool itensor_gpu_init(const itensor_gpu_layout&) { return false; }
bool itensor_gpu_reflection(int, const double*, const double*, const double*,
	const std::complex<double>*, std::complex<double>*) { return false; }
void itensor_gpu_free() {}

#else

#define SF_TWO_PI 6.283185307179586476925286766559
#define SF_INV_TWO_PI 0.15915494309189533576888376337251

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 I tensor GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

#define BLAS_TRY(call) do { if ((call) != GPUBLAS_STATUS_SUCCESS) { \
	std::fprintf(stderr, "NoSpherA2 I tensor GPU: BLAS failure at %s:%d\n", __FILE__, __LINE__); \
	return false; } } while (0)

namespace {

struct Dev {
	bool ready = false;
	gpublasHandle_t blas = nullptr;
	int nmo = 0, packed = 0, n_grids = 0, n_blocks = 0;
	long long n_points = 0;
	//Uploaded once
	float* ao = nullptr;
	int* aos = nullptr;
	unsigned char* skip = nullptr;
	double *d1 = nullptr, *d2 = nullptr, *d3 = nullptr, *w = nullptr;
	//Per reflection scratch
	float *phase_re = nullptr, *phase_im = nullptr;
	float *wre = nullptr, *wim = nullptr;
	float *cre = nullptr, *cim = nullptr;
	double* I_re = nullptr;
	double* I_im = nullptr;
	//Host copies of the small layout arrays, so the loop stays on the host
	std::vector<int> blk_grid, blk_ps, blk_pc, blk_na, grid_off;
	std::vector<long long> blk_ao_off, blk_aos_off;
} g;

//Same split the scattering-factor kernel uses: the phase and its reduction stay in
//double, only the transcendental drops to single. The weight is folded in here so the
//GEMM operand is exactly what the CPU path multiplies.
__global__ void phase_kernel(const long long n, const double kx, const double ky, const double kz,
	const double* __restrict__ d1, const double* __restrict__ d2, const double* __restrict__ d3,
	const double* __restrict__ w, float* __restrict__ pre, float* __restrict__ pim)
{
	const long long p = (long long)blockIdx.x * blockDim.x + threadIdx.x;
	if (p >= n) return;
	//kx..kz arrive already divided by 2pi, so t is in turns and sincospi wants 2*frac
	const double t = kx * d1[p] + ky * d2[p] + kz * d3[p];
	const double frac = t - rint(t);
	float s, c;
	sincospif(2.0f * (float)frac, &s, &c);
	const double wp = w[p];
	pre[p] = (float)(wp * (double)c);
	pim[p] = (float)(wp * (double)s);
}

//ao is row-major n_active x np; produce the two phase-weighted copies the GEMMs consume
__global__ void weight_kernel(const int n_active, const int np, const long long ao_off,
	const int point_base, const float* __restrict__ ao,
	const float* __restrict__ pre, const float* __restrict__ pim,
	float* __restrict__ wre, float* __restrict__ wim)
{
	const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
	const long long total = (long long)n_active * np;
	if (idx >= total) return;
	const int p = (int)(idx % np);
	const float a = ao[ao_off + idx];
	wre[idx] = a * pre[point_base + p];
	wim[idx] = a * pim[point_base + p];
}

//C is symmetric, so only the upper triangle is read back. Distinct (i,j) map to distinct
//packed indices, and launches are ordered on one stream, so no atomics are needed.
__global__ void accumulate_kernel(const int n_active, const int nmo, const long long aos_off,
	const int* __restrict__ aos, const unsigned char* __restrict__ skip,
	const float* __restrict__ cre, const float* __restrict__ cim,
	const double fre, const double fim,
	double* __restrict__ I_re, double* __restrict__ I_im)
{
	const int i = blockIdx.y * blockDim.y + threadIdx.y;
	const int j = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= n_active || j >= n_active || j < i) return;
	const int mu = aos[aos_off + i];
	const int nu = aos[aos_off + j];
	if (skip[(long long)mu * nmo + nu]) return;
	//cuBLAS wrote column-major n_active x n_active; C is symmetric so either index works
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

} //namespace

bool itensor_gpu_available()
{
	int n = 0;
	return gpuGetDeviceCount(&n) == gpuSuccess && n > 0;
}

bool itensor_gpu_init(const itensor_gpu_layout& L)
{
	itensor_gpu_free();
	if (!itensor_gpu_available() || !gpu_blas_runtime_present()) return false;

	int max_na = 0;
	long long max_elems = 0;
	for (int b = 0; b < L.n_blocks; b++) {
		max_na = std::max(max_na, L.blk_n_active[b]);
		//long long: active AOs times points passes 2^31 on a protein-sized grid
		max_elems = std::max(max_elems, (long long)L.blk_n_active[b] * L.blk_point_count[b]);
	}
	const size_t need =
		sizeof(float) * (size_t)L.ao_all_len +
		sizeof(int) * (size_t)L.aos_all_len +
		(size_t)L.nmo * L.nmo +
		sizeof(double) * 4 * (size_t)L.n_points +
		sizeof(float) * 2 * (size_t)L.n_points +
		sizeof(float) * 2 * (size_t)max_elems +
		sizeof(float) * 2 * (size_t)max_na * max_na +
		sizeof(double) * 2 * (size_t)L.packed;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	if (need + (1u << 26) > freeb) return false;

	GPU_TRY(gpuMalloc(&g.ao, sizeof(float) * (size_t)L.ao_all_len));
	GPU_TRY(gpuMalloc(&g.aos, sizeof(int) * (size_t)L.aos_all_len));
	GPU_TRY(gpuMalloc(&g.skip, (size_t)L.nmo * L.nmo));
	GPU_TRY(gpuMalloc(&g.d1, sizeof(double) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&g.d2, sizeof(double) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&g.d3, sizeof(double) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&g.w, sizeof(double) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&g.phase_re, sizeof(float) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&g.phase_im, sizeof(float) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&g.wre, sizeof(float) * (size_t)max_elems));
	GPU_TRY(gpuMalloc(&g.wim, sizeof(float) * (size_t)max_elems));
	GPU_TRY(gpuMalloc(&g.cre, sizeof(float) * (size_t)max_na * max_na));
	GPU_TRY(gpuMalloc(&g.cim, sizeof(float) * (size_t)max_na * max_na));
	GPU_TRY(gpuMalloc(&g.I_re, sizeof(double) * (size_t)L.packed));
	GPU_TRY(gpuMalloc(&g.I_im, sizeof(double) * (size_t)L.packed));

	GPU_TRY(gpuMemcpy(g.ao, L.ao_all, sizeof(float) * (size_t)L.ao_all_len, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(g.aos, L.aos_all, sizeof(int) * (size_t)L.aos_all_len, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(g.skip, L.skip, (size_t)L.nmo * L.nmo, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(g.d1, L.d1, sizeof(double) * (size_t)L.n_points, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(g.d2, L.d2, sizeof(double) * (size_t)L.n_points, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(g.d3, L.d3, sizeof(double) * (size_t)L.n_points, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(g.w, L.weights, sizeof(double) * (size_t)L.n_points, gpuMemcpyHostToDevice));

	if (gpublasCreate(&g.blas) != GPUBLAS_STATUS_SUCCESS) return false;

	g.nmo = L.nmo; g.packed = L.packed; g.n_grids = L.n_grids; g.n_blocks = L.n_blocks;
	g.n_points = L.n_points;
	g.blk_grid.assign(L.blk_grid, L.blk_grid + L.n_blocks);
	g.blk_ps.assign(L.blk_point_start, L.blk_point_start + L.n_blocks);
	g.blk_pc.assign(L.blk_point_count, L.blk_point_count + L.n_blocks);
	g.blk_na.assign(L.blk_n_active, L.blk_n_active + L.n_blocks);
	g.blk_ao_off.assign(L.blk_ao_off, L.blk_ao_off + L.n_blocks);
	g.blk_aos_off.assign(L.blk_aos_off, L.blk_aos_off + L.n_blocks);
	g.grid_off.assign(L.grid_point_off, L.grid_point_off + L.n_grids + 1);
	g.ready = true;
	return true;
}

bool itensor_gpu_reflection(const int num_syms,
	const double* kx, const double* ky, const double* kz,
	const std::complex<double>* factors,
	std::complex<double>* I_r)
{
	if (!g.ready) return false;
	zero_kernel<<<(g.packed + 255) / 256, 256>>>(g.packed, g.I_re, g.I_im);

	const float one = 1.0f, zero = 0.0f;
	for (int s = 0; s < num_syms; s++) {
		//The CPU path takes sin/cos of k.d directly; scaling k to turns here is what lets
		//the reduction be a rint and the transcendental be sincospi
		phase_kernel<<<(unsigned int)((g.n_points + 255) / 256), 256>>>(
			g.n_points, kx[s] * SF_INV_TWO_PI, ky[s] * SF_INV_TWO_PI, kz[s] * SF_INV_TWO_PI,
			g.d1, g.d2, g.d3, g.w, g.phase_re, g.phase_im);
		for (int b = 0; b < g.n_blocks; b++) {
			const int gi = g.blk_grid[b];
			const std::complex<double> f = factors[(size_t)s * g.n_grids + gi];
			if (f.real() == 0.0 && f.imag() == 0.0) continue;
			const int na = g.blk_na[b];
			const int np = g.blk_pc[b];
			const int base = g.grid_off[gi] + g.blk_ps[b];
			const long long total = (long long)na * np;
			weight_kernel<<<(unsigned int)((total + 255) / 256), 256>>>(
				na, np, g.blk_ao_off[b], base, g.ao, g.phase_re, g.phase_im, g.wre, g.wim);
			//C = A * W^T with both stored row-major n_active x np, i.e. column-major np x n_active
			BLAS_TRY(gpublasSgemm(g.blas, GPUBLAS_OP_T, GPUBLAS_OP_N, na, na, np,
				&one, g.ao + g.blk_ao_off[b], np, g.wre, np, &zero, g.cre, na));
			BLAS_TRY(gpublasSgemm(g.blas, GPUBLAS_OP_T, GPUBLAS_OP_N, na, na, np,
				&one, g.ao + g.blk_ao_off[b], np, g.wim, np, &zero, g.cim, na));
			const dim3 thr(16, 16);
			const dim3 blk((na + 15) / 16, (na + 15) / 16);
			accumulate_kernel<<<blk, thr>>>(na, g.nmo, g.blk_aos_off[b], g.aos, g.skip,
				g.cre, g.cim, f.real(), f.imag(), g.I_re, g.I_im);
		}
	}
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuDeviceSynchronize());

	std::vector<double> hr((size_t)g.packed), hi((size_t)g.packed);
	GPU_TRY(gpuMemcpy(hr.data(), g.I_re, sizeof(double) * (size_t)g.packed, gpuMemcpyDeviceToHost));
	GPU_TRY(gpuMemcpy(hi.data(), g.I_im, sizeof(double) * (size_t)g.packed, gpuMemcpyDeviceToHost));
	for (int i = 0; i < g.packed; i++)
		I_r[i] += std::complex<double>(hr[i], hi[i]);
	return true;
}

void itensor_gpu_free()
{
	if (g.blas) { gpublasDestroy(g.blas); g.blas = nullptr; }
	gpuFree(g.ao); gpuFree(g.aos); gpuFree(g.skip);
	gpuFree(g.d1); gpuFree(g.d2); gpuFree(g.d3); gpuFree(g.w);
	gpuFree(g.phase_re); gpuFree(g.phase_im);
	gpuFree(g.wre); gpuFree(g.wim); gpuFree(g.cre); gpuFree(g.cim);
	gpuFree(g.I_re); gpuFree(g.I_im);
	g = Dev{};
}

#endif //NOSPHERA2_HAVE_GPUBLAS
