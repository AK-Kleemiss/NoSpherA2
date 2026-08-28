#include "grid_gpu.h"
#include "gpu_backend.h"
#include <cstdio>
#include <algorithm>
#include <cstdlib>

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 grid GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

//Each thread needs two scratch values per centre, held in shared memory. 48 KB of shared
//per block divided by that is the block size, and below 32 threads the launch is not worth
//making, which puts the ceiling at 384 centres.
#define GRID_SHARED_BYTES (48 * 1024)
#define GRID_MIN_BLOCK 32
//Neighbours within far_away (10 bohr) of one point; a protein has a few dozen
#define GRID_MAX_NEAR 80

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
__device__ __forceinline__ void becke_point(const int ipoint, const int nc, const int center_index,
	const double* __restrict__ px, const double* __restrict__ py,
	const double* __restrict__ pz, const double* __restrict__ pw,
	const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
	const double* __restrict__ R_v, const double* __restrict__ chi,
	const double far_away, const double cut,
	double* pa_b, double* pa_tv,
	double* __restrict__ out_x, double* __restrict__ out_y, double* __restrict__ out_z,
	double* __restrict__ out_aw, double* __restrict__ out_becke, double* __restrict__ out_tfvc)
{
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

//Wide molecules: the scratch only has to cover the centres that are actually near the
//point. far_away is 10 bohr, so on a protein almost every centre is skipped by the same
//test the CPU applies, its pa stays zero and contributes nothing to either sum - and the
//pair updates written to a far centre are discarded when its own iteration zeroes it.
//Keeping pa only for the near centres is therefore exactly the CPU result, with the
//scratch set by the local environment rather than by the size of the molecule.
//
//A point with more than GRID_MAX_NEAR neighbours raises the flag and the caller falls
//back, rather than silently truncating a partitioning weight.
__global__ void becke_kernel_local(const int np, const int nc, const int center_index,
	const double* __restrict__ px, const double* __restrict__ py,
	const double* __restrict__ pz, const double* __restrict__ pw,
	const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
	const double* __restrict__ R_v, const double* __restrict__ chi,
	const double far_away, const double cut, int* __restrict__ overflow,
	double* __restrict__ out_x, double* __restrict__ out_y, double* __restrict__ out_z,
	double* __restrict__ out_aw, double* __restrict__ out_becke, double* __restrict__ out_tfvc)
{
	extern __shared__ double sh[];
	double* pa_b = sh + (size_t)threadIdx.x * 2 * GRID_MAX_NEAR;
	double* pa_tv = pa_b + GRID_MAX_NEAR;

	const int ipoint = blockIdx.x * blockDim.x + threadIdx.x;
	if (ipoint >= np) return;

	const double x = px[ipoint] + cx[center_index];
	const double y = py[ipoint] + cy[center_index];
	const double z = pz[ipoint] + cz[center_index];
	out_x[ipoint] = x;
	out_y[ipoint] = y;
	out_z[ipoint] = z;
	const double temp = pw[ipoint];

	//Ascending, so a walking index finds a near b without a search
	int near_idx[GRID_MAX_NEAR];
	double near_dist[GRID_MAX_NEAR];
	int n_near = 0;
	for (int a = 0; a < nc; a++) {
		const double vx = cx[a] - x, vy = cy[a] - y, vz = cz[a] - z;
		const double d = sqrt(vx * vx + vy * vy + vz * vz);
		if (d > far_away) continue;
		if (n_near >= GRID_MAX_NEAR) { *overflow = 1; return; }
		near_idx[n_near] = a;
		near_dist[n_near] = d;
		n_near++;
	}

	for (int i = 0; i < n_near; i++) { pa_b[i] = 1.0; pa_tv[i] = 1.0; }

	for (int ia = 0; ia < n_near; ia++) {
		const int a = near_idx[ia];
		const double dist_a = near_dist[ia];
		const double R_a = R_v[a];
		const double* chi_off = chi + (size_t)a * nc;
		int jb = ia + 1;

		for (int b = a + 1; b < nc; b++) {
			//A far b still contributes its factor to a; only its own pa is discarded
			const bool b_near = (jb < n_near && near_idx[jb] == b);
			double vx = cx[b] - cx[a], vy = cy[b] - cy[a], vz = cz[b] - cz[a];
			const double dist_ab = sqrt(vx * vx + vy * vy + vz * vz);
			vx = cx[b] - x; vy = cy[b] - y; vz = cz[b] - z;
			const double dist_b = sqrt(vx * vx + vy * vy + vz * vz);
			const double R_b = R_v[b];

			const double mu_ab = (dist_a - dist_b) / dist_ab;

			const double chi_mod = chi_off[b] * (1.0 - mu_ab);
			double nu_ab = 1.0 + mu_ab;
			nu_ab = (nu_ab - chi_mod) / (nu_ab + chi_mod);
			double f = 1.0 - f4d(nu_ab);
			if (fabs(f) < cut) pa_tv[ia] = 0.0;
			else {
				if (pa_tv[ia] > 1E-250 || pa_tv[ia] < -1E-250) pa_tv[ia] *= 0.5 * f;
				else pa_tv[ia] = 0.0;
				if (b_near) {
					if (pa_tv[jb] > 1E-250 || pa_tv[jb] < -1E-250) pa_tv[jb] *= 0.5 * (2.0 - f);
					else pa_tv[jb] = 0.0;
				}
			}

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
			if (fabs(f) < cut) pa_b[ia] = 0.0;
			else {
				if (pa_b[ia] > 1E-250 || pa_b[ia] < -1E-250) pa_b[ia] *= 0.5 * f;
				else pa_b[ia] = 0.0;
				if (b_near) {
					if (pa_b[jb] > 1E-250 || pa_b[jb] < -1E-250) pa_b[jb] *= 0.5 * (2.0 - f);
					else pa_b[jb] = 0.0;
				}
			}
			if (b_near) jb++;
		}
	}

	double w_becke = 0.0, w_tfvc = 0.0;
	int ic = -1;
	for (int i = 0; i < n_near; i++) {
		w_becke += pa_b[i];
		w_tfvc += pa_tv[i];
		if (near_idx[i] == center_index) ic = i;
	}
	//The owning centre is always within far_away of its own grid points
	const double cb = (ic >= 0) ? pa_b[ic] : 0.0;
	const double ct = (ic >= 0) ? pa_tv[ic] : 0.0;
	out_becke[ipoint] = temp * ((fabs(w_becke) > cut) ? cb / w_becke : 1.0);
	out_tfvc[ipoint] = temp * ((fabs(w_tfvc) > cut) ? ct / w_tfvc : 1.0);
	out_aw[ipoint] = temp;
}

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
	//Narrow molecules keep pa for every centre; wide ones keep it only for the centres
	//near the point, which is what makes a protein fit at all
	int block = (int)(GRID_SHARED_BYTES / (2 * sizeof(double) * (size_t)nc));
	//Only a molecule wider than a few hundred centres reaches the neighbour-local kernel,
	//so without this it would ship untested on every machine here
	const char* force = std::getenv("NOSPHERA2_GRID_LOCAL");
	const bool local_mode = (block < GRID_MIN_BLOCK) || (force && force[0] == '1');
	if (local_mode)
		block = (int)(GRID_SHARED_BYTES / (2 * sizeof(double) * (size_t)GRID_MAX_NEAR));
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

	int* d_over = nullptr;
	int host_over = 0;
	if (local_mode) {
		GPU_TRY(gpuMalloc(&d_over, sizeof(int)));
		GPU_TRY(gpuMemset(d_over, 0, sizeof(int)));
	}
	const size_t shmem = 2 * sizeof(double) * (size_t)(local_mode ? GRID_MAX_NEAR : nc) * block;
	if (local_mode)
		becke_kernel_local<<<(np + block - 1) / block, block, shmem>>>(np, nc, center_index,
			dpx, dpy, dpz, dpw, dcx, dcy, dcz, dR, dchi, far_away, cutoff, d_over,
			ox, oy, oz, oaw, ob, ot);
	else
		becke_kernel<<<(np + block - 1) / block, block, shmem>>>(np, nc, center_index,
			dpx, dpy, dpz, dpw, dcx, dcy, dcz, dR, dchi, far_away, cutoff,
			ox, oy, oz, oaw, ob, ot);
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuDeviceSynchronize());
	if (local_mode) {
		GPU_TRY(gpuMemcpy(&host_over, d_over, sizeof(int), gpuMemcpyDeviceToHost));
		gpuFree(d_over);
		//A point with more neighbours than the scratch holds: hand the atom back rather
		//than return a truncated weight
		if (host_over) {
			gpuFree(dpx); gpuFree(dpy); gpuFree(dpz); gpuFree(dpw);
			gpuFree(dcx); gpuFree(dcy); gpuFree(dcz); gpuFree(dR); gpuFree(dchi);
			gpuFree(ox); gpuFree(oy); gpuFree(oz); gpuFree(oaw); gpuFree(ob); gpuFree(ot);
			return false;
		}
	}

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
