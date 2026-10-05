#include "basin_field_gpu.h"
#include "aux_density_gpu.h"
#include "gpu_backend.h"
#include "gemm_gpu.cuh"
#include <cstdio>

NOSPHERA2_GPU_API_BEGIN

#define BF_BLOCK 128

namespace {

//One thread per (point, function), point fastest so the stores coalesce. Writes the function's K
//components column-major for the GEMM: chi[a * K * n + k * n + p]. A primitive whose centre or
//exponent is below the cutoff is skipped exactly as WFN::exp_table zeroes it on the host.
template <int K>
__global__ void bf_ao_kernel(const int n, const int nao, const double* __restrict__ pts,
	const double* __restrict__ cxyz, const double* __restrict__ cmin, const double cutoff,
	const int* __restrict__ ao_start, const int* __restrict__ pc, const int* __restrict__ pl,
	const double* __restrict__ pe, const double* __restrict__ ps, double* __restrict__ chi)
{
	const long long g = (long long)blockIdx.x * blockDim.x + threadIdx.x;
	if (g >= (long long)n * nao) return;
	const int p = (int)(g % n), a = (int)(g / n);
	const double px = pts[3 * p], py = pts[3 * p + 1], pz = pts[3 * p + 2];
	double t[K];
#pragma unroll
	for (int k = 0; k < K; k++) t[k] = 0.0;
	for (int j = ao_start[a]; j < ao_start[a + 1]; j++)
	{
		const int c = pc[j];
		const double d[3] = { px - cxyz[3 * c], py - cxyz[3 * c + 1], pz - cxyz[3 * c + 2] };
		const double R = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
		if (-cmin[c] * R < cutoff) continue;
		const double e = -pe[j] * R;
		if (e < cutoff) continue;
		const double ex = exp(e), ex2 = 2 * pe[j], s = ps[j];
		//per axis f = x^l, f' and f'' without the exponential, as in WFN::eli_orbital_pass
		double f[3], f1[3], f2[3];
#pragma unroll
		for (int k = 0; k < 3; k++)
		{
			const int l = pl[3 * j + k];
			double pw[13];
			pw[0] = 1.0;
			for (int m = 1; m <= l + 2; m++) pw[m] = pw[m - 1] * d[k];
			f[k] = pw[l];
			f1[k] = (l ? l * pw[l - 1] : 0.0) - ex2 * pw[l + 1];
			if (K == 10)
				f2[k] = (l > 1 ? l * (l - 1) * pw[l - 2] : 0.0) - ex2 * (2 * l + 1) * pw[l] + ex2 * ex2 * pw[l + 2];
		}
		if (K == 4)
		{
			//WFN::computeGrad's factor order
			t[0] += s * (f[0] * f[1] * f[2] * ex);
			t[1] += s * (f1[0] * f[1] * f[2] * ex);
			t[2] += s * (f1[1] * f[0] * f[2] * ex);
			t[3] += s * (f1[2] * f[0] * f[1] * ex);
		}
		else
		{
			//0 value, 1-3 x y z, 4-6 xx yy zz, 7 xy, 8 xz, 9 yz
			t[0] += s * (f[0] * f[1] * f[2] * ex);
			t[1] += s * (f1[0] * f[1] * f[2] * ex);
			t[2] += s * (f[0] * f1[1] * f[2] * ex);
			t[3] += s * (f[0] * f[1] * f1[2] * ex);
			t[4 % K] += s * (f2[0] * f[1] * f[2] * ex);
			t[5 % K] += s * (f[0] * f2[1] * f[2] * ex);
			t[6 % K] += s * (f[0] * f[1] * f2[2] * ex);
			t[7 % K] += s * (f1[0] * f1[1] * f[2] * ex);
			t[8 % K] += s * (f1[0] * f[1] * f1[2] * ex);
			t[9 % K] += s * (f[0] * f1[1] * f1[2] * ex);
		}
	}
	double* out = chi + (size_t)a * K * n + p;
#pragma unroll
	for (int k = 0; k < K; k++) out[(size_t)k * n] = t[k];
}

//One thread per point: phi[mo * K * n + k * n + p] reduced over the occupied MOs with the arithmetic
//of WFN::computeGrad (K = 4) and WFN::computeELIGrad (K = 10).
template <int K>
__global__ void bf_reduce_kernel(const int n, const int nocc, const double* __restrict__ occ,
	const double* __restrict__ phi, double* __restrict__ val, double* __restrict__ grad, double* __restrict__ rho_out)
{
	const int p = blockIdx.x * blockDim.x + threadIdx.x;
	if (p >= n) return;
	const size_t ld = (size_t)K * n;
	if (K == 4)
	{
		double G[3]{ 0, 0, 0 }, rho = 0;
		for (int mo = 0; mo < nocc; mo++)
		{
			const double* q = phi + mo * ld + p;
			const double o = occ[mo], docc = 2 * o, v = q[0];
			G[0] += docc * v * q[n];
			G[1] += docc * v * q[2 * (size_t)n];
			G[2] += docc * v * q[3 * (size_t)n];
			rho += o * v * v;
		}
		for (int k = 0; k < 3; k++) grad[3 * (size_t)p + k] = G[k];
		if (rho_out) rho_out[p] = rho;
		return;
	}
	const int hidx[3][3] ={ {4, 7, 8}, {7, 5, 9}, {8, 9, 6} };
	double rho = 0, tau = 0, G[3]{ 0, 0, 0 }, T[3]{ 0, 0, 0 }, H[3][3]{ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
	for (int mo = 0; mo < nocc; mo++)
	{
		const double* q = phi + mo * ld + p;
		double v[K];
#pragma unroll
		for (int k = 0; k < K; k++) v[k] = q[(size_t)k * n];
		const double o = occ[mo], docc = 2 * o;
		rho += o * v[0] * v[0];
#pragma unroll
		for (int i = 0; i < 3; i++)
		{
			tau += o * v[1 + i] * v[1 + i];
			G[i] += docc * v[0] * v[1 + i];
#pragma unroll
			for (int k = 0; k < 3; k++)
			{
				H[i][k] += docc * (v[0] * v[hidx[i][k] % K] + v[1 + i] * v[1 + k]);
				T[k] += docc * v[1 + i] * v[hidx[i][k] % K];
			}
		}
	}
	if (rho_out) rho_out[p] = rho;
	const double g = rho * tau - 0.25 * (G[0] * G[0] + G[1] * G[1] + G[2] * G[2]);
	double* gr = grad + 3 * (size_t)p;
	if (!(g > 0))
	{
		if (val) val[p] = 0;
		gr[0] = gr[1] = gr[2] = 0;
		return;
	}
	const double f = pow(48 / g, 0.375);
	if (val) val[p] = 0.5 * rho * f;
	for (int k = 0; k < 3; k++)
	{
		const double dg = G[k] * tau + rho * T[k] - 0.5 * (G[0] * H[0][k] + G[1] * H[1][k] + G[2] * H[2][k]);
		gr[k] = 0.5 * f * (G[k] - 0.375 * rho * dg / g);
	}
}

template <typename T>
bool upload(T*& d, const T* h, const size_t n)
{
	if (gpuMalloc(&d, sizeof(T) * (n ? n : 1)) != gpuSuccess) { d = nullptr; return false; }
	return n == 0 || gpuMemcpy(d, h, sizeof(T) * n, gpuMemcpyHostToDevice) == gpuSuccess;
}

} //namespace

bool basin_field_gpu_eval(
	const int K, const int ncen, const double* cxyz, const double* cmin_exp, const double exp_cutoff,
	const int nao, const int* ao_start, const int* prim_center, const int* prim_l, const double* prim_exp, const double* prim_scale,
	const int nocc, const double* coef, const double* occ,
	const int np, const double* pts, double* val, double* grad, double* rho)
{
	if (np <= 0 || nao <= 0 || nocc <= 0 || (K != 4 && K != 10)) return false;
	if (!aux_density_gpu_enabled() || !aux_density_gpu_available()) return false;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;

	//Points per chunk: chi (K x nao), phi and the GEMM workspace (K x nocc each) per point, within half
	//of what is free. The tables are re-sent every call - a few MB against a call's GEMM.
	// ponytail: tables uploaded per call; hold them like esp_gpu does if calls get small and frequent
	const int nprim = ao_start[nao];
	const size_t per_point = sizeof(double) * ((size_t)K * nao + 2 * (size_t)K * nocc + 8);
	size_t budget = freeb / 2;
	if (budget > ((size_t)4 << 30)) budget = (size_t)4 << 30;
	long long nb = (long long)(budget / per_point);
	if (nb > np) nb = np;
	//the GEMM indexes in int: keep K * nb * nao below 2^31
	const long long cap = ((1LL << 31) - 1) / ((long long)K * (nao > nocc ? nao : nocc));
	if (nb > cap) nb = cap;
	if (nb < 1) return false;
	const int chunk = (int)nb;

	double *d_c = nullptr, *d_cmin = nullptr, *d_pe = nullptr, *d_ps = nullptr, *d_coef = nullptr, *d_occ = nullptr;
	int *d_start = nullptr, *d_pc = nullptr, *d_pl = nullptr;
	double *d_pts = nullptr, *d_chi = nullptr, *d_phi = nullptr, *d_P = nullptr, *d_val = nullptr, *d_grad = nullptr, *d_rho = nullptr;
	const size_t ws = gemm_gpu::workspace_elems(K * chunk, nocc, gemm_gpu::split_count(K * chunk, nocc, nao));
	bool ok = upload(d_c, cxyz, 3 * (size_t)ncen) && upload(d_cmin, cmin_exp, (size_t)ncen)
		&& upload(d_start, ao_start, (size_t)nao + 1) && upload(d_pc, prim_center, (size_t)nprim)
		&& upload(d_pl, prim_l, 3 * (size_t)nprim) && upload(d_pe, prim_exp, (size_t)nprim)
		&& upload(d_ps, prim_scale, (size_t)nprim) && upload(d_coef, coef, (size_t)nao * nocc)
		&& upload(d_occ, occ, (size_t)nocc)
		&& gpuMalloc(&d_pts, sizeof(double) * 3 * (size_t)chunk) == gpuSuccess
		&& gpuMalloc(&d_chi, sizeof(double) * (size_t)K * chunk * nao) == gpuSuccess
		&& gpuMalloc(&d_phi, sizeof(double) * (size_t)K * chunk * nocc) == gpuSuccess
		&& gpuMalloc(&d_P, sizeof(double) * ws) == gpuSuccess
		&& gpuMalloc(&d_val, sizeof(double) * (size_t)chunk) == gpuSuccess
		&& gpuMalloc(&d_grad, sizeof(double) * 3 * (size_t)chunk) == gpuSuccess
		&& gpuMalloc(&d_rho, sizeof(double) * (size_t)chunk) == gpuSuccess;

	for (int first = 0; first < np && ok; first += chunk)
	{
		const int n = np - first < chunk ? np - first : chunk;
		ok = gpuMemcpy(d_pts, pts + 3 * (size_t)first, sizeof(double) * 3 * (size_t)n, gpuMemcpyHostToDevice) == gpuSuccess;
		if (!ok) break;
		const long long work = (long long)n * nao;
		const unsigned blocks = (unsigned)((work + BF_BLOCK - 1) / BF_BLOCK);
		if (K == 4)
			bf_ao_kernel<4><<<blocks, BF_BLOCK>>>(n, nao, d_pts, d_c, d_cmin, exp_cutoff, d_start, d_pc, d_pl, d_pe, d_ps, d_chi);
		else
			bf_ao_kernel<10><<<blocks, BF_BLOCK>>>(n, nao, d_pts, d_c, d_cmin, exp_cutoff, d_start, d_pc, d_pl, d_pe, d_ps, d_chi);
		//phi (K n x nocc, column-major) = chi (K n x nao) * coef^T, coef being nocc x nao column-major
		gemm_gpu::launch<double>(false, true, K * n, nocc, nao, 1.0, d_chi, K * n, d_coef, nocc, 0.0, d_phi, K * n, d_P);
		const unsigned rb = (unsigned)((n + BF_BLOCK - 1) / BF_BLOCK);
		if (K == 4)
			bf_reduce_kernel<4><<<rb, BF_BLOCK>>>(n, nocc, d_occ, d_phi, d_val, d_grad, d_rho);
		else
			bf_reduce_kernel<10><<<rb, BF_BLOCK>>>(n, nocc, d_occ, d_phi, d_val, d_grad, d_rho);
		ok = gpuGetLastError() == gpuSuccess && gpuDeviceSynchronize() == gpuSuccess
			&& gpuMemcpy(grad + 3 * (size_t)first, d_grad, sizeof(double) * 3 * (size_t)n, gpuMemcpyDeviceToHost) == gpuSuccess
			&& (!val || K == 4 || gpuMemcpy(val + first, d_val, sizeof(double) * (size_t)n, gpuMemcpyDeviceToHost) == gpuSuccess)
			&& (!rho || gpuMemcpy(rho + first, d_rho, sizeof(double) * (size_t)n, gpuMemcpyDeviceToHost) == gpuSuccess);
	}
	if (!ok)
		std::fprintf(stderr, "NoSpherA2 basin field GPU: %s, falling back to the host\n", gpuGetErrorString(gpuGetLastError()));
	gpuFree(d_c); gpuFree(d_cmin); gpuFree(d_pe); gpuFree(d_ps); gpuFree(d_coef); gpuFree(d_occ);
	gpuFree(d_start); gpuFree(d_pc); gpuFree(d_pl);
	gpuFree(d_pts); gpuFree(d_chi); gpuFree(d_phi); gpuFree(d_P); gpuFree(d_val); gpuFree(d_grad); gpuFree(d_rho);
	return ok;
}

NOSPHERA2_GPU_API_END
