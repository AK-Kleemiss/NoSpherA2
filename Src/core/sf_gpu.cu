#include "sf_gpu.h"
#include "tuning.h"
#include "gpu_backend.h"
#include <cstdio>
#include <cstring>
#include <future>
#include <vector>
#include <algorithm>
#include <cstdlib>

NOSPHERA2_GPU_API_BEGIN

//Each block owns a tile of k-points for one atom and streams that atom's grid points through
//shared memory.
#define SF_TILE_K 128
#define SF_CHUNK 256
#define SF_TWO_PI 6.283185307179586476925286766559
#define SF_INV_TWO_PI 0.15915494309189533576888376337251

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
		if (live)
			for (int p = 0; p < n; p++) {
				double si, co;
				sincos(kx * s1[p] + ky * s2[p] + kz * s3[p], &si, &co);
				re += sd[p] * co;
				im += sd[p] * si;
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

//No double in the point loop, which a consumer card runs at 1/32 or 1/64 rate: the points arrive as floats relative
//to their atom's centre, so k.r splits into the centre's phase, reduced once per thread in double, and a k.delta that
//floats carry to ~1e-7 turns per turn of |k.delta|. sincospif reduces any argument exactly. Each chunk's sum is
//Kahan-compensated in float and folded into a double, so the error does not grow with the point count.
__global__ void sf_kernel_f32(const int imax, const long long smax,
	const double* __restrict__ k1, const double* __restrict__ k2, const double* __restrict__ k3,
	const double* __restrict__ centre,
	const float* __restrict__ d1, const float* __restrict__ d2, const float* __restrict__ d3,
	const float* __restrict__ dens, const int* __restrict__ offs,
	double2* __restrict__ sf_out)
{
	__shared__ float s1[SF_CHUNK], s2[SF_CHUNK], s3[SF_CHUNK], sd[SF_CHUNK];
	const int ia = blockIdx.y;
	const long long s = (long long)blockIdx.x * SF_TILE_K + threadIdx.x;
	const bool live = (s < smax);
	//In turns, so the reduction is a rint and a subtract
	const double kx = live ? k1[s] * SF_INV_TWO_PI : 0.0;
	const double ky = live ? k2[s] * SF_INV_TWO_PI : 0.0;
	const double kz = live ? k3[s] * SF_INV_TWO_PI : 0.0;
	const double w0 = kx * centre[3 * ia] + ky * centre[3 * ia + 1] + kz * centre[3 * ia + 2];
	const float p0 = (float)(w0 - rint(w0));
	const float fx = (float)kx, fy = (float)ky, fz = (float)kz;
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
			float cre = 0.0f, cim = 0.0f, kre = 0.0f, kim = 0.0f;
			for (int p = 0; p < n; p++) {
				float sif, cof;
				sincospif(2.0f * fmaf(fx, s1[p], fmaf(fy, s2[p], fmaf(fz, s3[p], p0))), &sif, &cof);
				float y = sd[p] * cof - kre;
				float t = cre + y;
				kre = (t - cre) - y;
				cre = t;
				y = sd[p] * sif - kim;
				t = cim + y;
				kim = (t - cim) - y;
				cim = t;
			}
			re += (double)cre;
			im += (double)cim;
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

namespace {
__global__ void probe_kernel(int* p) { if (p) *p = 1; }
}

//A build pinned to other architectures has no code for a present card, and every launch then fails
//as if there were no GPU; one empty launch is the only way to find out, so say so plainly.
bool sf_gpu_available()
{
	int n = 0;
	if (gpuGetDeviceCount(&n) != gpuSuccess || n <= 0) return false;

	static const bool usable = []() {
		probe_kernel<<<1, 1>>>(nullptr);
		const gpuError_t e = gpuDeviceSynchronize();
		//Clear the sticky error either way, so a later launch is judged on its own merits
		(void)gpuGetLastError();
		if (e == gpuSuccess) return true;
		//Compute capability has no HIP equivalent, and NOSPHERA2_CUDA_PORTABLE is wrong advice on AMD.
#ifdef NOSPHERA2_USE_HIP
		std::fprintf(stderr, "NoSpherA2: a GPU is present but this build contains no code "
					 "for it (%s), so every GPU path will use the CPU. Rebuild with this "
					 "card's architecture in CMAKE_HIP_ARCHITECTURES.\n", gpuGetErrorString(e));
#else
		gpuDeviceProp_t prop{};
		int dev = 0;
		if (gpuGetDevice(&dev) == gpuSuccess && gpuGetDeviceProperties(&prop, dev) == gpuSuccess)
			std::fprintf(stderr, "NoSpherA2: a GPU is present (compute %d.%d) but this build "
						 "contains no code for it, so every GPU path will use the CPU. "
						 "Rebuild with -DNOSPHERA2_CUDA_PORTABLE=ON.\n", prop.major, prop.minor);
		else
			std::fprintf(stderr, "NoSpherA2: a GPU is present but unusable (%s); "
						 "every GPU path will use the CPU.\n", gpuGetErrorString(e));
#endif
		return false;
	}();
	return usable;
}

//Creates the context while the grids are built, not on the first allocation inside the transform.
static std::future<void> g_warmup;

void sf_gpu_warmup_start()
{
	int n = 0;
	if (gpuGetDeviceCount(&n) != gpuSuccess || n <= 0) return;
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

//Single-to-double throughput ratio, 2 on datacentre and 32 or 64 on consumer parts.
int sf_gpu_fp64_ratio()
{
	int dev = 0;
	if (gpuGetDevice(&dev) != gpuSuccess) return 0;
#ifdef NOSPHERA2_USE_HIP
	//No HIP attribute: CDNA parts have the wide fp64 units, RDNA and APU parts do not.
	gpuDeviceProp_t prop;
	if (gpuGetDeviceProperties(&prop, dev) != gpuSuccess) return 0;
	static const char* const cdna[] = { "gfx906", "gfx908", "gfx90a", "gfx940", "gfx941", "gfx942", "gfx950" };
	for (int i = 0; i < 7; i++)
		if (std::strncmp(prop.gcnArchName, cdna[i], std::strlen(cdna[i])) == 0) return 2;
	return 32;
#else
	int ratio = 0;
	if (cudaDeviceGetAttribute(&ratio, cudaDevAttrSingleToDoublePrecisionPerfRatio, dev) != cudaSuccess) return 0;
	return ratio;
#endif
}

//Resolved here so the log line cannot name a precision the transform did not use.
bool sf_gpu_uses_fp32(const sf_precision prec)
{
	if (prec == sf_precision::FP32) return true;
	if (prec == sf_precision::FP64) return false;
	return sf_gpu_fp64_ratio() > 4;
}

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

bool sf_gpu_run(const int imax, const long long smax,
	const double* k1, const double* k2, const double* k3,
	const double* d1, const double* d2, const double* d3,
	const double* dens, const int* offs, const long long total_points,
	double* const* sf_rows, const sf_precision prec)
{
	//offs is int-indexed upstream, so a problem past INT_MAX points is already out of range
	if (imax <= 0 || (long long)offs[imax] != total_points) return false;
	const size_t kb = sizeof(double) * (size_t)smax;
	const size_t row = sizeof(double2) * (size_t)smax;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	if (kb * 3 + (1u << 26) >= freeb) return false;
	//Atoms are independent, so oversized problems are batched; a batch costs its grid points four
	//times over plus one output row per atom.
	const size_t budget = freeb - kb * 3 - (1u << 26);
	int batch = imax;
	{
		size_t widest = 0;
		for (int i = 0; i < imax; i++)
			widest = std::max(widest, (size_t)(offs[i + 1] - offs[i]));
		//A single atom that will not fit cannot be split without partial sums
		if (4 * sizeof(double) * widest + row > budget) return false;
		while (batch > 1) {
			size_t need = row * (size_t)batch;
			size_t worst = 0;
			for (int a = 0; a < imax; a += batch)
				worst = std::max(worst, (size_t)(offs[std::min(a + batch, imax)] - offs[a]));
			need += 4 * sizeof(double) * worst;
			if (need <= budget) break;
			batch = (batch + 1) / 2;
		}
	}
	//Forces batching, which a card large enough for the whole problem never exercises.
	if (const char* cap = tuning("NOSPHERA2_GPU_BATCH")) {
		const int c = std::atoi(cap);
		if (c > 0 && c < batch) batch = c;
	}
	size_t widest_pts = 0;
	for (int a = 0; a < imax; a += batch)
		widest_pts = std::max(widest_pts, (size_t)(offs[std::min(a + batch, imax)] - offs[a]));
	const size_t pts = sizeof(double) * widest_pts;
	double *dk1 = nullptr, *dk2 = nullptr, *dk3 = nullptr, *dd1 = nullptr, *dd2 = nullptr, *dd3 = nullptr, *dde = nullptr;
	double2* dout = nullptr;
	int* dof = nullptr;
	GPU_TRY(gpuMalloc(&dk1, kb)); GPU_TRY(gpuMalloc(&dk2, kb)); GPU_TRY(gpuMalloc(&dk3, kb));
	GPU_TRY(gpuMalloc(&dd1, pts)); GPU_TRY(gpuMalloc(&dd2, pts)); GPU_TRY(gpuMalloc(&dd3, pts));
	GPU_TRY(gpuMalloc(&dde, pts));
	GPU_TRY(gpuMalloc(&dof, sizeof(int) * (size_t)(batch + 1)));
	GPU_TRY(gpuMalloc(&dout, row * (size_t)batch));
	GPU_TRY(gpuMemcpy(dk1, k1, kb, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dk2, k2, kb, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dk3, k3, kb, gpuMemcpyHostToDevice));
	const bool f32 = sf_gpu_uses_fp32(prec);
	double* dce = nullptr;
	if (f32) GPU_TRY(gpuMalloc(&dce, sizeof(double) * 3 * (size_t)batch));
	std::vector<int> rel(batch + 1);
	std::vector<float> h1, h2, h3, hd;
	std::vector<double> cen;
	for (int a0 = 0; a0 < imax; a0 += batch) {
		const int na = std::min(batch, imax - a0);
		const int p0 = offs[a0];
		const size_t np = (size_t)(offs[a0 + na] - p0);
		for (int i = 0; i <= na; i++)
			rel[i] = offs[a0 + i] - p0;
		GPU_TRY(gpuMemcpy(dof, rel.data(), sizeof(int) * (size_t)(na + 1), gpuMemcpyHostToDevice));
		const dim3 grid((unsigned int)((smax + SF_TILE_K - 1) / SF_TILE_K), (unsigned int)na);
		if (f32) {
			//Each atom's points about their mean, the centre that keeps |delta| and with it the float error smallest
			h1.resize(np); h2.resize(np); h3.resize(np); hd.resize(np);
			cen.assign(3 * (size_t)na, 0.0);
			for (int i = 0; i < na; i++) {
				const int lo = p0 + rel[i], hi = p0 + rel[i + 1];
				double* c = cen.data() + 3 * i;
				for (int p = lo; p < hi; p++) { c[0] += d1[p]; c[1] += d2[p]; c[2] += d3[p]; }
				for (int x = 0; x < 3; x++) c[x] /= std::max(1, hi - lo);
				for (int p = lo; p < hi; p++) {
					h1[p - p0] = (float)(d1[p] - c[0]);
					h2[p - p0] = (float)(d2[p] - c[1]);
					h3[p - p0] = (float)(d3[p] - c[2]);
					hd[p - p0] = (float)dens[p];
				}
			}
			//The double buffers are reused, half filled
			GPU_TRY(gpuMemcpy(dd1, h1.data(), sizeof(float) * np, gpuMemcpyHostToDevice));
			GPU_TRY(gpuMemcpy(dd2, h2.data(), sizeof(float) * np, gpuMemcpyHostToDevice));
			GPU_TRY(gpuMemcpy(dd3, h3.data(), sizeof(float) * np, gpuMemcpyHostToDevice));
			GPU_TRY(gpuMemcpy(dde, hd.data(), sizeof(float) * np, gpuMemcpyHostToDevice));
			GPU_TRY(gpuMemcpy(dce, cen.data(), sizeof(double) * cen.size(), gpuMemcpyHostToDevice));
			sf_kernel_f32<<<grid, SF_TILE_K>>>(na, smax, dk1, dk2, dk3, dce, reinterpret_cast<float*>(dd1),
				reinterpret_cast<float*>(dd2), reinterpret_cast<float*>(dd3), reinterpret_cast<float*>(dde), dof, dout);
		}
		else {
			GPU_TRY(gpuMemcpy(dd1, d1 + p0, sizeof(double) * np, gpuMemcpyHostToDevice));
			GPU_TRY(gpuMemcpy(dd2, d2 + p0, sizeof(double) * np, gpuMemcpyHostToDevice));
			GPU_TRY(gpuMemcpy(dd3, d3 + p0, sizeof(double) * np, gpuMemcpyHostToDevice));
			GPU_TRY(gpuMemcpy(dde, dens + p0, sizeof(double) * np, gpuMemcpyHostToDevice));
			sf_kernel<<<grid, SF_TILE_K>>>(na, smax, dk1, dk2, dk3, dd1, dd2, dd3, dde, dof, dout);
		}
		GPU_TRY(gpuGetLastError());
		GPU_TRY(gpuDeviceSynchronize());
		//Straight into the caller's complex rows, no staging buffer.
		for (int i = 0; i < na; i++)
			GPU_TRY(gpuMemcpy(sf_rows[a0 + i], dout + (size_t)i * (size_t)smax, row, gpuMemcpyDeviceToHost));
	}
	gpuFree(dk1); gpuFree(dk2); gpuFree(dk3);
	gpuFree(dd1); gpuFree(dd2); gpuFree(dd3); gpuFree(dde);
	gpuFree(dof); gpuFree(dout); gpuFree(dce);
	return true;
}

NOSPHERA2_GPU_API_END
