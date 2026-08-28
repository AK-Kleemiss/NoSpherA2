#include "sf_gpu.h"
#include "gpu_backend.h"
#include <cstdio>
#include <cstring>
#include <future>

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
	double2* __restrict__ sf_out)
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
		double2 v;
		v.x = re;
		v.y = im;
		sf_out[(long long)ia * smax + s] = v;
	}
}

bool sf_gpu_available()
{
	int n = 0;
	return gpuGetDeviceCount(&n) == gpuSuccess && n > 0;
}

//Creating the device context costs ~95 ms and would otherwise happen on the first
//allocation, inside the transform. Started when the grids begin and waited for just
//before the transform, it overlaps CPU work that has to happen anyway.
static std::future<void> g_warmup;

void sf_gpu_warmup_start()
{
	if (!g_warmup.valid())
		g_warmup = std::async(std::launch::async, []() { gpuFree(0); });
}

void sf_gpu_warmup_wait()
{
	if (g_warmup.valid())
		g_warmup.get();
}

const char* sf_gpu_backend()
{
#ifdef NOSPHERA2_USE_HIP
	return "HIP";
#else
	return "CUDA";
#endif
}

//Single-to-double throughput ratio: 2 on datacentre parts, 32 or 64 on consumer ones.
//Above the threshold the double sincos is not worth paying for.
int sf_gpu_fp64_ratio()
{
	int dev = 0;
	if (gpuGetDevice(&dev) != gpuSuccess) return 0;
#ifdef NOSPHERA2_USE_HIP
	//HIP exposes no equivalent attribute, so the arch name has to stand in for it.
	//CDNA parts carry the wide fp64 units; every RDNA and APU part does not.
	gpuDeviceProp_t prop;
	if (gpuGetDeviceProperties(&prop, dev) != gpuSuccess) return 0;
	static const char* const cdna[] = { "gfx906", "gfx908", "gfx90a", "gfx940", "gfx941", "gfx942", "gfx950" };
	for (const char* a : cdna)
		if (std::strncmp(prop.gcnArchName, a, std::strlen(a)) == 0) return 2;
	return 32;
#else
	int ratio = 0;
	if (cudaDeviceGetAttribute(&ratio, cudaDevAttrSingleToDoublePrecisionPerfRatio, dev) != cudaSuccess) return 0;
	return ratio;
#endif
}

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

bool sf_gpu_run(const int imax, const long long smax,
	const double* k1, const double* k2, const double* k3,
	const double* d1, const double* d2, const double* d3,
	const double* dens, const int* offs, const long long total_points,
	double* const* sf_rows, const bool force_fp64)
{
	const size_t pts = sizeof(double) * (size_t)total_points;
	const size_t kb = sizeof(double) * (size_t)smax;
	const size_t row = sizeof(double2) * (size_t)smax;
	const size_t out = row * (size_t)imax;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	//Bail to the CPU rather than thrash if the whole problem will not sit in device memory
	if (4 * pts + 3 * kb + out + (1u << 26) > freeb) return false;
	double *dk1 = nullptr, *dk2 = nullptr, *dk3 = nullptr, *dd1 = nullptr, *dd2 = nullptr, *dd3 = nullptr, *dde = nullptr;
	double2* dout = nullptr;
	int* dof = nullptr;
	GPU_TRY(gpuMalloc(&dk1, kb)); GPU_TRY(gpuMalloc(&dk2, kb)); GPU_TRY(gpuMalloc(&dk3, kb));
	GPU_TRY(gpuMalloc(&dd1, pts)); GPU_TRY(gpuMalloc(&dd2, pts)); GPU_TRY(gpuMalloc(&dd3, pts));
	GPU_TRY(gpuMalloc(&dde, pts));
	GPU_TRY(gpuMalloc(&dof, sizeof(int) * (size_t)(imax + 1)));
	GPU_TRY(gpuMalloc(&dout, out));
	GPU_TRY(gpuMemcpy(dk1, k1, kb, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dk2, k2, kb, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dk3, k3, kb, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dd1, d1, pts, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dd2, d2, pts, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dd3, d3, pts, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dde, dens, pts, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dof, offs, sizeof(int) * (size_t)(imax + 1), gpuMemcpyHostToDevice));
	const dim3 grid((unsigned int)((smax + SF_TILE_K - 1) / SF_TILE_K), (unsigned int)imax);
	const int ratio = sf_gpu_fp64_ratio();
	const bool f32 = !force_fp64 && ratio > 4;
	if (f32)
		sf_kernel<true><<<grid, SF_TILE_K>>>(imax, smax, dk1, dk2, dk3, dd1, dd2, dd3, dde, dof, dout);
	else
		sf_kernel<false><<<grid, SF_TILE_K>>>(imax, smax, dk1, dk2, dk3, dd1, dd2, dd3, dde, dof, dout);
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuDeviceSynchronize());
	//Straight into the caller's complex storage, one row per atom: no staging buffer to
	//allocate, zero and scatter, which cost more than the transfer itself did.
	for (int i = 0; i < imax; i++)
		GPU_TRY(gpuMemcpy(sf_rows[i], dout + (size_t)i * (size_t)smax, row, gpuMemcpyDeviceToHost));
	gpuFree(dk1); gpuFree(dk2); gpuFree(dk3);
	gpuFree(dd1); gpuFree(dd2); gpuFree(dd3); gpuFree(dde);
	gpuFree(dof); gpuFree(dout);
	return true;
}
