#include "spherical_density_gpu.h"
#include "aux_density_gpu.h"
#include "gpu_backend.h"
#include <cstdio>

NOSPHERA2_GPU_API_BEGIN

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 spherical density GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

//Below this many point-atom pairs the host loop finishes before the copies do
#define SPH_MIN_WORK (1 << 22)
#define SPH_BLOCK 128

namespace {
//Thakkar::get_interpolated_density on a table slice, including the +-1 index correction of
//log_spline_index against the stored radii, so the interval is the one the CPU picks
__device__ double lookup(const double* r, const double* rho, const int n, const double lincr, const double start, const double dist)
{
	if (dist > r[n - 1]) return 0.0;
	if (dist < r[0]) return rho[0];
	const int max_idx = n - 2;
	int nr = (int)floor(log(dist / start) / lincr);
	if (nr < 0) nr = 0;
	if (nr > max_idx) nr = max_idx;
	if (nr < max_idx && dist >= r[nr + 1]) ++nr;
	else if (nr > 0 && dist < r[nr]) --nr;
	const double v = rho[nr] + (rho[nr + 1] - rho[nr]) / (r[nr + 1] - r[nr]) * (dist - r[nr]);
	return v < 1E-16 ? 0.0 : v;
}

__global__ void spherical_density_kernel(const int nx, const int ny, const int nz, const double* origin, const double* vectors,
	const int n_at, const double* ax, const double* ay, const double* az, const int* at_tab,
	const int* tab_off, const double* r_tab, const double* rho_tab,
	const double lincr, const double start, const double radius_bohr, double* out)
{
	const long long p = (long long)blockIdx.x * blockDim.x + threadIdx.x;
	const long long np = (long long)nx * ny * nz;
	if (p >= np) return;
	const int z = (int)(p % nz), y = (int)((p / nz) % ny), x = (int)(p / ((long long)ny * nz));
	const double px = x * vectors[0] + y * vectors[1] + z * vectors[2] + origin[0];
	const double py = x * vectors[3] + y * vectors[4] + z * vectors[5] + origin[1];
	const double pz = x * vectors[6] + y * vectors[7] + z * vectors[8] + origin[2];
	double dens = 0.0;
	bool inside = false;
	for (int a = 0; a < n_at; a++) {
		const double d = norm3d(px - ax[a], py - ay[a], pz - az[a]);
		inside |= d < radius_bohr;
		const int t = at_tab[a], o = tab_off[t];
		dens += lookup(r_tab + o, rho_tab + o, tab_off[t + 1] - o, lincr, start, d);
	}
	out[p] = inside ? dens : 0.0;
}
}

bool spherical_density_gpu_eval(
	const int nx, const int ny, const int nz, const double* origin, const double* vectors,
	const int n_at, const double* ax, const double* ay, const double* az, const int* at_tab,
	const int n_tab, const int* tab_off, const double* r_tab, const double* rho_tab,
	const double lincr, const double start, const double radius_bohr, double* out)
{
	const long long np = (long long)nx * ny * nz;
	if (np <= 0 || n_at <= 0 || n_tab <= 0) return false;
	if (!aux_density_gpu_enabled() || !aux_density_gpu_available()) return false;
	if (np * n_at < SPH_MIN_WORK) return false;
	const size_t pts = sizeof(double) * (size_t)np, at = sizeof(double) * (size_t)n_at;
	const size_t tab = sizeof(double) * (size_t)tab_off[n_tab];
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	if (pts + 3 * at + 2 * tab + (1u << 26) > freeb) return false;
	double *dout = nullptr, *dax = nullptr, *day = nullptr, *daz = nullptr, *dr = nullptr, *drho = nullptr, *dorig = nullptr, *dvec = nullptr;
	int *dtab = nullptr, *doff = nullptr;
	GPU_TRY(gpuMalloc(&dout, pts));
	GPU_TRY(gpuMalloc(&dax, at)); GPU_TRY(gpuMalloc(&day, at)); GPU_TRY(gpuMalloc(&daz, at));
	GPU_TRY(gpuMalloc(&dtab, sizeof(int) * (size_t)n_at)); GPU_TRY(gpuMalloc(&doff, sizeof(int) * (size_t)(n_tab + 1)));
	GPU_TRY(gpuMalloc(&dr, tab)); GPU_TRY(gpuMalloc(&drho, tab));
	GPU_TRY(gpuMalloc(&dorig, sizeof(double) * 3)); GPU_TRY(gpuMalloc(&dvec, sizeof(double) * 9));
	GPU_TRY(gpuMemcpy(dax, ax, at, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(day, ay, at, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(daz, az, at, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dtab, at_tab, sizeof(int) * (size_t)n_at, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(doff, tab_off, sizeof(int) * (size_t)(n_tab + 1), gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dr, r_tab, tab, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(drho, rho_tab, tab, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dorig, origin, sizeof(double) * 3, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dvec, vectors, sizeof(double) * 9, gpuMemcpyHostToDevice));
	spherical_density_kernel<<<(unsigned)((np + SPH_BLOCK - 1) / SPH_BLOCK), SPH_BLOCK>>>(nx, ny, nz, dorig, dvec, n_at, dax, day, daz, dtab, doff, dr, drho, lincr, start, radius_bohr, dout);
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuDeviceSynchronize());
	GPU_TRY(gpuMemcpy(out, dout, pts, gpuMemcpyDeviceToHost));
	gpuFree(dout); gpuFree(dax); gpuFree(day); gpuFree(daz); gpuFree(dtab); gpuFree(doff);
	gpuFree(dr); gpuFree(drho); gpuFree(dorig); gpuFree(dvec);
	return true;
}

NOSPHERA2_GPU_API_END
