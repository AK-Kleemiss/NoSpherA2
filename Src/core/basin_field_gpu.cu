#include "basin_field_gpu.h"
#include "aux_density_gpu.h"
#include "gpu_backend.h"
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <vector>

NOSPHERA2_GPU_API_BEGIN

//threads per point (a multiple of 64 so a HIP wavefront is never split) and AOs per shared tile
#define BF_T 128
#define BF_TILE 32

namespace {

//The K components of contracted function a at the point, summed over its primitives. A
//primitive whose centre or exponent is below the cutoff is skipped exactly as WFN::exp_table zeroes
//it on the host; returns whether any primitive was inside.
template <int K>
__device__ bool bf_ao(const int a, const double px, const double py, const double pz,
	const double* __restrict__ cxyz, const double* __restrict__ cmin, const double cutoff,
	const int* __restrict__ ao_start, const int* __restrict__ pc, const int* __restrict__ pl,
	const double* __restrict__ pe, const double* __restrict__ ps, double* t)
{
	bool touched = false;
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
		touched = true;
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
	return touched;
}

//One block per point. Pass 1 screens the functions, BF_T at a time, by the distance to their
//centre against their smallest exponent (the test bf_ao applies to that primitive, so nothing it
//would keep is dropped), evaluates the survivors and packs them in ascending function order into
//this point's slice of the workspace (a block scan, so the order and the sums are deterministic).
//Pass 2 gives each thread MOs, sums the packed functions times the MO coefficients through shared
//tiles, and reduces per point with the arithmetic of WFN::computeGrad (K = 4) and
//WFN::computeELIGrad (K = 10): rho, G[3] (K = 4) or rho, tau, G[3], H[3][3], T[3] (K = 10).
template <int K>
__global__ void __launch_bounds__(BF_T) bf_point_kernel(const int nao, const int nocc, const double* __restrict__ pts,
	const double* __restrict__ cxyz, const double* __restrict__ cmin, const double cutoff,
	const int* __restrict__ ao_start, const int* __restrict__ pc, const int* __restrict__ pl,
	const double* __restrict__ pe, const double* __restrict__ ps,
	const int* __restrict__ ao_c, const double* __restrict__ ao_minpe,
	const double* __restrict__ coef, const double* __restrict__ occ,
	int* __restrict__ idx_ws, double* __restrict__ chi_ws,
	double* __restrict__ val, double* __restrict__ grad, double* __restrict__ rho_out, unsigned long long* __restrict__ count)
{
	constexpr int NR = K == 4 ? 4 : 17;
	__shared__ int s_scan[BF_T];
	__shared__ int s_cnt;
	__shared__ int s_idx[BF_TILE];
	__shared__ double s_chi[BF_TILE * K];
	__shared__ double s_red[NR * BF_T];
	const int p = blockIdx.x, tid = threadIdx.x;
	int* idx = idx_ws + (size_t)p * nao;
	double* chi = chi_ws + (size_t)p * nao * K;
	const double px = pts[3 * p], py = pts[3 * p + 1], pz = pts[3 * p + 2];
	if (tid == 0) s_cnt = 0;
	__syncthreads();
	for (int base = 0; base < nao; base += BF_T)
	{
		const int a = base + tid;
		double t[K];
		bool on = false;
		if (a < nao)
		{
			const int c = ao_c[a];
			bool cand = c < 0;
			if (!cand)
			{
				const double dx = px - cxyz[3 * c], dy = py - cxyz[3 * c + 1], dz = pz - cxyz[3 * c + 2];
				cand = -ao_minpe[a] * (dx * dx + dy * dy + dz * dz) >= cutoff;
			}
			if (cand) on = bf_ao<K>(a, px, py, pz, cxyz, cmin, cutoff, ao_start, pc, pl, pe, ps, t);
		}
		s_scan[tid] = on;
		__syncthreads();
		for (int off = 1; off < BF_T; off <<= 1)
		{
			const int v = tid >= off ? s_scan[tid - off] : 0;
			__syncthreads();
			s_scan[tid] += v;
			__syncthreads();
		}
		if (on)
		{
			const int j = s_cnt + s_scan[tid] - 1;
			idx[j] = a;
#pragma unroll
			for (int k = 0; k < K; k++) chi[(size_t)j * K + k] = t[k];
		}
		__syncthreads();
		if (tid == BF_T - 1) s_cnt += s_scan[tid];
		__syncthreads();
	}
	const int cnt = s_cnt;
	if (count && tid == 0) atomicAdd(count, (unsigned long long)cnt);

	const int hidx[3][3] = { {4, 7, 8}, {7, 5, 9}, {8, 9, 6} };
	double acc[NR];
#pragma unroll
	for (int r = 0; r < NR; r++) acc[r] = 0.0;
	for (int mo0 = 0; mo0 < nocc; mo0 += BF_T)
	{
		const int mo = mo0 + tid;
		double v[K];
#pragma unroll
		for (int k = 0; k < K; k++) v[k] = 0.0;
		for (int j0 = 0; j0 < cnt; j0 += BF_TILE)
		{
			const int m = cnt - j0 < BF_TILE ? cnt - j0 : BF_TILE;
			__syncthreads();
			for (int e = tid; e < m * K; e += BF_T) s_chi[e] = chi[(size_t)j0 * K + e];
			if (tid < m) s_idx[tid] = idx[j0 + tid];
			__syncthreads();
			if (mo < nocc)
				for (int j = 0; j < m; j++)
				{
					const double cf = coef[(size_t)s_idx[j] * nocc + mo];
#pragma unroll
					for (int k = 0; k < K; k++) v[k] += cf * s_chi[j * K + k];
				}
		}
		if (mo >= nocc) continue;
		const double o = occ[mo], docc = 2 * o;
		if constexpr (K == 4)
		{
			acc[0] += o * v[0] * v[0];
			acc[1] += docc * v[0] * v[1];
			acc[2] += docc * v[0] * v[2];
			acc[3] += docc * v[0] * v[3];
		}
		else
		{
			//0 rho, 1 tau, 2-4 G, 5-13 H row-major, 14-16 T
			acc[0] += o * v[0] * v[0];
#pragma unroll
			for (int i = 0; i < 3; i++)
			{
				acc[1] += o * v[1 + i] * v[1 + i];
				acc[2 + i] += docc * v[0] * v[1 + i];
#pragma unroll
				for (int k = 0; k < 3; k++)
				{
					acc[5 + 3 * i + k] += docc * (v[0] * v[hidx[i][k] % K] + v[1 + i] * v[1 + k]);
					acc[14 + k] += docc * v[1 + i] * v[hidx[i][k] % K];
				}
			}
		}
	}
#pragma unroll
	for (int r = 0; r < NR; r++) s_red[r * BF_T + tid] = acc[r];
	__syncthreads();
	for (int h = BF_T / 2; h > 0; h >>= 1)
	{
		if (tid < h)
			for (int r = 0; r < NR; r++) s_red[r * BF_T + tid] += s_red[r * BF_T + tid + h];
		__syncthreads();
	}
	if (tid != 0) return;
	double R[NR];
	for (int r = 0; r < NR; r++) R[r] = s_red[r * BF_T];
	const double rho = R[0];
	if (rho_out) rho_out[p] = rho;
	double* gr = grad + 3 * (size_t)p;
	if constexpr (K == 4)
	{
		gr[0] = R[1]; gr[1] = R[2]; gr[2] = R[3];
	}
	else
	{
	const double tau = R[1], *G = R + 2, *H = R + 5, *T = R + 14;
	const double g = rho * tau - 0.25 * (G[0] * G[0] + G[1] * G[1] + G[2] * G[2]);
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
		const double dg = G[k] * tau + rho * T[k] - 0.5 * (G[0] * H[k] + G[1] * H[3 + k] + G[2] * H[6 + k]);
		gr[k] = 0.5 * f * (G[k] - 0.375 * rho * dg / g);
	}
	}
}

template <typename T>
bool upload(T*& d, const T* h, const size_t n)
{
	if (gpuMalloc(&d, sizeof(T) * (n ? n : 1)) != gpuSuccess) { d = nullptr; return false; }
	return n == 0 || gpuMemcpy(d, h, sizeof(T) * n, gpuMemcpyHostToDevice) == gpuSuccess;
}

struct bf_ctx
{
	int K = 0, nao = 0, nocc = 0, chunk = 0;
	double cutoff = 0;
	gpuStream_t s = nullptr;
	double *c = nullptr, *cmin = nullptr, *pe = nullptr, *ps = nullptr, *coef = nullptr, *occ = nullptr, *minpe = nullptr;
	int *start = nullptr, *pc = nullptr, *pl = nullptr, *aoc = nullptr, *idx = nullptr;
	double *pts = nullptr, *chi = nullptr, *val = nullptr, *grad = nullptr, *rho = nullptr;
	//NOS_BASIN_GPU_PROFILE: a sync after every stage and the seconds per stage, printed at close
	bool prof = false;
	double st[3]{ 0, 0, 0 };
	long long prof_points = 0, prof_runs = 0;
	unsigned long long* count = nullptr;
	~bf_ctx()
	{
		if (prof)
		{
			unsigned long long h = 0;
			if (count) gpuMemcpy(&h, count, sizeof(h), gpuMemcpyDeviceToHost);
			const double np = prof_points > 0 ? (double)prof_points : 1.0;
			std::printf("  basin field GPU (K=%d, nao %d, nocc %d): %lld points in %lld runs; copy-in %.3f s, kernel %.3f s, copy-out %.3f s; per point %.1f active functions\n",
				K, nao, nocc, prof_points, prof_runs, st[0], st[1], st[2], h / np);
		}
		gpuFree(count);
		gpuFree(c); gpuFree(cmin); gpuFree(pe); gpuFree(ps); gpuFree(coef); gpuFree(occ); gpuFree(minpe);
		gpuFree(start); gpuFree(pc); gpuFree(pl); gpuFree(aoc); gpuFree(idx);
		gpuFree(pts); gpuFree(chi); gpuFree(val); gpuFree(grad); gpuFree(rho);
		if (s) gpuStreamDestroy(s);
	}
};

} //namespace

void* basin_field_gpu_open(
	const int K, const int ncen, const double* cxyz, const double* cmin_exp, const double exp_cutoff,
	const int nao, const int* ao_start, const int* prim_center, const int* prim_l, const double* prim_exp, const double* prim_scale,
	const int nocc, const double* coef, const double* occ, const int max_points)
{
	if (max_points <= 0 || nao <= 0 || nocc <= 0 || (K != 4 && K != 10)) return nullptr;
	if (!aux_density_gpu_enabled() || !aux_density_gpu_available()) return nullptr;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return nullptr;

	//Points per chunk: the packed workspace (K values and an index per function) per point, within
	//half of what is free
	const int nprim = ao_start[nao];
	const size_t per_point = sizeof(double) * ((size_t)K * nao + 5) + sizeof(int) * (size_t)nao;
	size_t budget = freeb / 2;
	if (budget > ((size_t)4 << 30)) budget = (size_t)4 << 30;
	long long nb = (long long)(budget / per_point);
	if (nb > max_points) nb = max_points;
	if (nb < 1) return nullptr;

	//The screen per function: its centre (-1 if the primitives sit on several, then it is always
	//evaluated) and its smallest exponent
	std::vector<int> aoc(nao, -1);
	std::vector<double> minpe(nao, 0.0);
	for (int a = 0; a < nao; a++)
	{
		if (ao_start[a] == ao_start[a + 1]) continue;
		aoc[a] = prim_center[ao_start[a]];
		minpe[a] = prim_exp[ao_start[a]];
		for (int j = ao_start[a]; j < ao_start[a + 1]; j++)
		{
			if (prim_center[j] != aoc[a]) aoc[a] = -1;
			minpe[a] = std::min(minpe[a], prim_exp[j]);
		}
	}

	bf_ctx* x = new bf_ctx;
	x->K = K; x->nao = nao; x->nocc = nocc; x->chunk = (int)nb; x->cutoff = exp_cutoff;
	x->prof = std::getenv("NOS_BASIN_GPU_PROFILE") != nullptr;
	const size_t ch = (size_t)x->chunk;
	const bool ok = gpuStreamCreateNonBlocking(&x->s) == gpuSuccess
		&& upload(x->c, cxyz, 3 * (size_t)ncen) && upload(x->cmin, cmin_exp, (size_t)ncen)
		&& upload(x->start, ao_start, (size_t)nao + 1) && upload(x->pc, prim_center, (size_t)nprim)
		&& upload(x->pl, prim_l, 3 * (size_t)nprim) && upload(x->pe, prim_exp, (size_t)nprim)
		&& upload(x->ps, prim_scale, (size_t)nprim) && upload(x->coef, coef, (size_t)nao * nocc)
		&& upload(x->occ, occ, (size_t)nocc) && upload(x->aoc, aoc.data(), (size_t)nao) && upload(x->minpe, minpe.data(), (size_t)nao)
		&& gpuMalloc(&x->pts, sizeof(double) * 3 * ch) == gpuSuccess
		&& gpuMalloc(&x->chi, sizeof(double) * (size_t)K * ch * nao) == gpuSuccess
		&& gpuMalloc(&x->idx, sizeof(int) * ch * nao) == gpuSuccess
		&& gpuMalloc(&x->val, sizeof(double) * ch) == gpuSuccess
		&& gpuMalloc(&x->grad, sizeof(double) * 3 * ch) == gpuSuccess
		&& gpuMalloc(&x->rho, sizeof(double) * ch) == gpuSuccess;
	if (!ok)
	{
		std::fprintf(stderr, "NoSpherA2 basin field GPU: %s while allocating, falling back to the host\n", gpuGetErrorString(gpuGetLastError()));
		delete x;
		return nullptr;
	}
	if (x->prof && gpuMalloc(&x->count, sizeof(unsigned long long)) == gpuSuccess)
		gpuMemset(x->count, 0, sizeof(unsigned long long));
	return x;
}

bool basin_field_gpu_run(void* ctx, const int np, const double* pts, double* val, double* grad, double* rho)
{
	bf_ctx* x = (bf_ctx*)ctx;
	if (!x || np <= 0) return false;
	const int K = x->K, nao = x->nao, nocc = x->nocc, chunk = x->chunk;
	bool ok = true;
	gpuError_t err = gpuSuccess;
	auto t = std::chrono::steady_clock::now();
	auto lap = [&](const int s) {
		if (!x->prof) return;
		gpuStreamSynchronize(x->s);
		const auto now = std::chrono::steady_clock::now();
		x->st[s] += std::chrono::duration<double>(now - t).count();
		t = now;
	};
	x->prof_points += np;
	x->prof_runs++;
	for (int first = 0; first < np && ok; first += chunk)
	{
		const int n = np - first < chunk ? np - first : chunk;
		ok = gpuMemcpyAsync(x->pts, pts + 3 * (size_t)first, sizeof(double) * 3 * (size_t)n, gpuMemcpyHostToDevice, x->s) == gpuSuccess;
		if (!ok) break;
		lap(0);
		if (K == 4)
			bf_point_kernel<4><<<n, BF_T, 0, x->s>>>(nao, nocc, x->pts, x->c, x->cmin, x->cutoff, x->start, x->pc, x->pl, x->pe, x->ps,
				x->aoc, x->minpe, x->coef, x->occ, x->idx, x->chi, x->val, x->grad, x->rho, x->count);
		else
			bf_point_kernel<10><<<n, BF_T, 0, x->s>>>(nao, nocc, x->pts, x->c, x->cmin, x->cutoff, x->start, x->pc, x->pl, x->pe, x->ps,
				x->aoc, x->minpe, x->coef, x->occ, x->idx, x->chi, x->val, x->grad, x->rho, x->count);
		err = gpuGetLastError();
		lap(1);
		ok = err == gpuSuccess
			&& gpuMemcpyAsync(grad + 3 * (size_t)first, x->grad, sizeof(double) * 3 * (size_t)n, gpuMemcpyDeviceToHost, x->s) == gpuSuccess
			&& (!val || K == 4 || gpuMemcpyAsync(val + first, x->val, sizeof(double) * (size_t)n, gpuMemcpyDeviceToHost, x->s) == gpuSuccess)
			&& (!rho || gpuMemcpyAsync(rho + first, x->rho, sizeof(double) * (size_t)n, gpuMemcpyDeviceToHost, x->s) == gpuSuccess);
		if (ok)
		{
			err = gpuStreamSynchronize(x->s);
			ok = err == gpuSuccess;
		}
		lap(2);
	}
	if (!ok)
		std::fprintf(stderr, "NoSpherA2 basin field GPU: %s, falling back to the host\n", gpuGetErrorString(err != gpuSuccess ? err : gpuGetLastError()));
	return ok;
}

void basin_field_gpu_close(void* ctx)
{
	delete (bf_ctx*)ctx;
}

NOSPHERA2_GPU_API_END
