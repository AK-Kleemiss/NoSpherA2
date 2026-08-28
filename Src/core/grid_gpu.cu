#include "grid_gpu.h"
#include "gpu_backend.h"
#include <cstdio>
#include <algorithm>

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 grid GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

//Each thread needs two scratch values per centre, held in shared memory. 48 KB of shared
//per block divided by that is the block size, and below 32 threads the launch is not worth
//making, which puts the ceiling at 384 centres.
#define GRID_SHARED_BYTES (48 * 1024)
#define GRID_MIN_BLOCK 32

static bool g_grid_use_gpu = false;
void grid_gpu_set_enabled(bool on) { g_grid_use_gpu = on; }
bool grid_gpu_enabled() { return g_grid_use_gpu; }

namespace {

//The Becke iterated switching polynomial, three and four iterations
__device__ __forceinline__ double f3d(double x)
{
	double f = x;
	f *= (1.5 - 0.5 * f * f);
	f *= (1.5 - 0.5 * f * f);
	f *= (1.5 - 0.5 * f * f);
	return f;
}

__device__ __forceinline__ double f4d(double x)
{
	double f = x;
	f *= (1.5 - 0.5 * f * f);
	f *= (1.5 - 0.5 * f * f);
	f *= (1.5 - 0.5 * f * f);
	f *= (1.5 - 0.5 * f * f);
	return f;
}

//A verbatim transcription of the chi-present branch of get_integration_weights. The pair
//loop updates both a and b, so it stays sequential per point exactly as on the CPU; the
//parallelism is over points, which are independent.
__global__ void becke_kernel(const int np, const int nc, const int center_index,
	const double* __restrict__ px, const double* __restrict__ py,
	const double* __restrict__ pz, const double* __restrict__ pw,
	const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
	const double* __restrict__ R_v, const double* __restrict__ chi,
	const double far_away, const double cut,
	double* __restrict__ out_x, double* __restrict__ out_y, double* __restrict__ out_z,
	double* __restrict__ out_aw, double* __restrict__ out_becke, double* __restrict__ out_tfvc)
{
	extern __shared__ double sh[];
	double* pa_b = sh + (size_t)threadIdx.x * 2 * nc;
	double* pa_tv = pa_b + nc;

	const int ipoint = blockIdx.x * blockDim.x + threadIdx.x;
	if (ipoint >= np) return;

	const double x = px[ipoint] + cx[center_index];
	const double y = py[ipoint] + cy[center_index];
	const double z = pz[ipoint] + cz[center_index];
	out_x[ipoint] = x;
	out_y[ipoint] = y;
	out_z[ipoint] = z;
	const double temp = pw[ipoint];

	for (int a = 0; a < nc; a++) { pa_b[a] = 1.0; pa_tv[a] = 1.0; }

	for (int a = 0; a < nc; a++) {
		double vx = cx[a] - x, vy = cy[a] - y, vz = cz[a] - z;
		const double dist_a = sqrt(vx * vx + vy * vy + vz * vz);
		if (dist_a > far_away) { pa_b[a] = 0.0; pa_tv[a] = 0.0; continue; }
		const double R_a = R_v[a];
		const double* chi_off = chi + (size_t)a * nc;

		for (int b = a + 1; b < nc; b++) {
			vx = cx[b] - cx[a]; vy = cy[b] - cy[a]; vz = cz[b] - cz[a];
			const double dist_ab = sqrt(vx * vx + vy * vy + vz * vz);
			vx = cx[b] - x; vy = cy[b] - y; vz = cz[b] - z;
			const double dist_b = sqrt(vx * vx + vy * vy + vz * vz);
			const double R_b = R_v[b];

			const double mu_ab = (dist_a - dist_b) / dist_ab;

			//TFVC: JCP 139, 071103 (2013), eq. 7
			const double chi_mod = chi_off[b] * (1.0 - mu_ab);
			double nu_ab = 1.0 + mu_ab;
			nu_ab = (nu_ab - chi_mod) / (nu_ab + chi_mod);
			double f = 1.0 - f4d(nu_ab);
			if (fabs(f) < cut) pa_tv[a] = 0.0;
			else {
				if (pa_tv[a] > 1E-250 || pa_tv[a] < -1E-250) pa_tv[a] *= 0.5 * f;
				else pa_tv[a] = 0.0;
				if (pa_tv[b] > 1E-250 || pa_tv[b] < -1E-250) pa_tv[b] *= 0.5 * (2.0 - f);
				else pa_tv[b] = 0.0;
			}

			//Becke: JCP 88, 2547 (1988), eqs. 11 and A3
			if (fabs(R_a - R_b) > cut) {
				const double chi_becke = R_a / R_b;
				double u_ab = (chi_becke - 1.0) / (chi_becke + 1.0);
				u_ab = u_ab / (u_ab * u_ab - 1.0);
				if (u_ab > 0.5) u_ab = 0.5;
				else if (u_ab < -0.5) u_ab = -0.5;
				nu_ab = mu_ab + u_ab * (1.0 - mu_ab * mu_ab);
			}
			else nu_ab = mu_ab;

			f = 1.0 - f3d(nu_ab);
			if (fabs(f) < cut) pa_b[a] = 0.0;
			else {
				if (pa_b[a] > 1E-250 || pa_b[a] < -1E-250) pa_b[a] *= 0.5 * f;
				else pa_b[a] = 0.0;
				if (pa_b[b] > 1E-250 || pa_b[b] < -1E-250) pa_b[b] *= 0.5 * (2.0 - f);
				else pa_b[b] = 0.0;
			}
		}
	}

	double w_becke = 0.0, w_tfvc = 0.0;
	for (int a = 0; a < nc; a++) { w_becke += pa_b[a]; w_tfvc += pa_tv[a]; }

	const double rb = (fabs(w_becke) > cut) ? pa_b[center_index] / w_becke : 1.0;
	const double rt = (fabs(w_tfvc) > cut) ? pa_tv[center_index] / w_tfvc : 1.0;
	out_becke[ipoint] = temp * rb;
	out_tfvc[ipoint] = temp * rt;
	out_aw[ipoint] = temp;
}

} //namespace

bool grid_gpu_available()
{
	int n = 0;
	return gpuGetDeviceCount(&n) == gpuSuccess && n > 0;
}

bool grid_gpu_becke_weights(const int np, const int nc, const int center_index,
	const double* proto_x, const double* proto_y, const double* proto_z, const double* proto_w,
	const double* cx, const double* cy, const double* cz,
	const double* R_v, const double* chi,
	const double far_away, const double cutoff,
	double* out_x, double* out_y, double* out_z,
	double* out_aw, double* out_becke, double* out_tfvc)
{
	if (np <= 0 || nc <= 1) return false;
	if (!g_grid_use_gpu || !grid_gpu_available()) return false;
	//Two doubles per centre per thread, so a wide molecule shrinks the block
	const int block = (int)(GRID_SHARED_BYTES / (2 * sizeof(double) * (size_t)nc));
	if (block < GRID_MIN_BLOCK) return false;

	const size_t pts = sizeof(double) * (size_t)np;
	const size_t cen = sizeof(double) * (size_t)nc;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	if (10 * pts + 4 * cen + sizeof(double) * (size_t)nc * nc + (1u << 26) > freeb) return false;

	double *dpx = nullptr, *dpy = nullptr, *dpz = nullptr, *dpw = nullptr;
	double *dcx = nullptr, *dcy = nullptr, *dcz = nullptr, *dR = nullptr, *dchi = nullptr;
	double *ox = nullptr, *oy = nullptr, *oz = nullptr, *oaw = nullptr, *ob = nullptr, *ot = nullptr;
	GPU_TRY(gpuMalloc(&dpx, pts)); GPU_TRY(gpuMalloc(&dpy, pts));
	GPU_TRY(gpuMalloc(&dpz, pts)); GPU_TRY(gpuMalloc(&dpw, pts));
	GPU_TRY(gpuMalloc(&dcx, cen)); GPU_TRY(gpuMalloc(&dcy, cen)); GPU_TRY(gpuMalloc(&dcz, cen));
	GPU_TRY(gpuMalloc(&dR, cen));
	GPU_TRY(gpuMalloc(&dchi, sizeof(double) * (size_t)nc * nc));
	GPU_TRY(gpuMalloc(&ox, pts)); GPU_TRY(gpuMalloc(&oy, pts)); GPU_TRY(gpuMalloc(&oz, pts));
	GPU_TRY(gpuMalloc(&oaw, pts)); GPU_TRY(gpuMalloc(&ob, pts)); GPU_TRY(gpuMalloc(&ot, pts));

	GPU_TRY(gpuMemcpy(dpx, proto_x, pts, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dpy, proto_y, pts, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dpz, proto_z, pts, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dpw, proto_w, pts, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dcx, cx, cen, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dcy, cy, cen, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dcz, cz, cen, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dR, R_v, cen, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dchi, chi, sizeof(double) * (size_t)nc * nc, gpuMemcpyHostToDevice));

	const size_t shmem = 2 * sizeof(double) * (size_t)nc * block;
	becke_kernel<<<(np + block - 1) / block, block, shmem>>>(np, nc, center_index,
		dpx, dpy, dpz, dpw, dcx, dcy, dcz, dR, dchi, far_away, cutoff,
		ox, oy, oz, oaw, ob, ot);
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuDeviceSynchronize());

	GPU_TRY(gpuMemcpy(out_x, ox, pts, gpuMemcpyDeviceToHost));
	GPU_TRY(gpuMemcpy(out_y, oy, pts, gpuMemcpyDeviceToHost));
	GPU_TRY(gpuMemcpy(out_z, oz, pts, gpuMemcpyDeviceToHost));
	GPU_TRY(gpuMemcpy(out_aw, oaw, pts, gpuMemcpyDeviceToHost));
	GPU_TRY(gpuMemcpy(out_becke, ob, pts, gpuMemcpyDeviceToHost));
	GPU_TRY(gpuMemcpy(out_tfvc, ot, pts, gpuMemcpyDeviceToHost));

	gpuFree(dpx); gpuFree(dpy); gpuFree(dpz); gpuFree(dpw);
	gpuFree(dcx); gpuFree(dcy); gpuFree(dcz); gpuFree(dR); gpuFree(dchi);
	gpuFree(ox); gpuFree(oy); gpuFree(oz); gpuFree(oaw); gpuFree(ob); gpuFree(ot);
	return true;
}
