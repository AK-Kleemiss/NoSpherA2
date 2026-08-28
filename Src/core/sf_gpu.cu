#include "sf_gpu.h"
#include "gpu_backend.h"
#include <cstdio>
#include <cstring>
#include <future>
#include <vector>
#include <algorithm>
#include <cstdlib>

//Each block owns a tile of k-points for one atom and streams that atom's grid points
//through shared memory. F32 keeps the phase and its reduction in double and drops only
//the transcendental and the running sum to single, which consumer parts run 32-64x
//faster. 64, 128 and 256 k-points per tile measure the same, so 128 is not tuned.
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
	//Scaled to turns once per thread, so the phase reduction below is a rint and a
	//subtract instead of a multiply, an add, a floor and an fma. Costs three fp64
	//multiplies against the ~2000 grid points each thread then walks.
	const double kx = live ? k1[s] * (F32 ? SF_INV_TWO_PI : 1.0) : 0.0;
	const double ky = live ? k2[s] * (F32 ? SF_INV_TWO_PI : 1.0) : 0.0;
	const double kz = live ? k3[s] * (F32 ? SF_INV_TWO_PI : 1.0) : 0.0;
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
				const double w = kx * s1[p] + ky * s2[p] + kz * s3[p];
				const double r = sd[p];
				if (F32) {
					//w is in turns here. The reduction stays in double - sincospif of an
					//unreduced argument would be meaningless - and only the transcendental
					//and the running sum drop to fp32.
					const double wr = w - rint(w);
					float sif, cof;
					sincospif(2.0f * (float)wr, &sif, &cof);
					//Compensated, so the error does not grow with the term count. Ten fp32
					//operations still cost far less than two fp64 fmas on a consumer part.
					const float pr = (float)r * cof;
					const float pi = (float)r * sif;
					float y = pr - kre;
					float t = cre + y;
					kre = (t - cre) - y;
					cre = t;
					y = pi - kim;
					t = cim + y;
					kim = (t - cim) - y;
					cim = t;
				}
				else {
					double si, co;
					sincos(w, &si, &co);
					re += r * co;
					im += r * si;
				}
			}
			if (F32) {
				re += (double)cre;
				im += (double)cim;
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

//Resolving Auto lives here rather than at the call site because the log line has to name
//the kernel that actually ran. Repeating the ratio test where the message is written is how
//the two drift apart and the log starts claiming a precision the transform did not use.
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
	//Atoms are independent, so oversized problems go in batches rather than to the CPU.
	//The budget is what is left once the k-points and a margin are accounted for; a batch
	//costs its own grid points four times over plus one output row per atom.
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
	//Every card here has room for the whole problem, so the batching would otherwise never
	//run until it met a protein on someone else's machine.
	if (const char* cap = std::getenv("NOSPHERA2_GPU_BATCH")) {
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
	std::vector<int> rel(batch + 1);
	for (int a0 = 0; a0 < imax; a0 += batch) {
		const int na = std::min(batch, imax - a0);
		const int p0 = offs[a0];
		const size_t np = (size_t)(offs[a0 + na] - p0);
		for (int i = 0; i <= na; i++)
			rel[i] = offs[a0 + i] - p0;
		GPU_TRY(gpuMemcpy(dd1, d1 + p0, sizeof(double) * np, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(dd2, d2 + p0, sizeof(double) * np, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(dd3, d3 + p0, sizeof(double) * np, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(dde, dens + p0, sizeof(double) * np, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(dof, rel.data(), sizeof(int) * (size_t)(na + 1), gpuMemcpyHostToDevice));
		const dim3 grid((unsigned int)((smax + SF_TILE_K - 1) / SF_TILE_K), (unsigned int)na);
		if (f32)
			sf_kernel<true><<<grid, SF_TILE_K>>>(na, smax, dk1, dk2, dk3, dd1, dd2, dd3, dde, dof, dout);
		else
			sf_kernel<false><<<grid, SF_TILE_K>>>(na, smax, dk1, dk2, dk3, dd1, dd2, dd3, dde, dof, dout);
		GPU_TRY(gpuGetLastError());
		GPU_TRY(gpuDeviceSynchronize());
		//Straight into the caller's complex storage, one row per atom: no staging buffer to
		//allocate, zero and scatter, which cost more than the transfer itself did.
		for (int i = 0; i < na; i++)
			GPU_TRY(gpuMemcpy(sf_rows[a0 + i], dout + (size_t)i * (size_t)smax, row, gpuMemcpyDeviceToHost));
	}
	gpuFree(dk1); gpuFree(dk2); gpuFree(dk3);
	gpuFree(dd1); gpuFree(dd2); gpuFree(dd3); gpuFree(dde);
	gpuFree(dof); gpuFree(dout);
	return true;
}
