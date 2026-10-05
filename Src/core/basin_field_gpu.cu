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

namespace {

struct bf_ctx
{
	int K = 0, nao = 0, nocc = 0, chunk = 0;
	double cutoff = 0;
	double *c = nullptr, *cmin = nullptr, *pe = nullptr, *ps = nullptr, *coef = nullptr, *occ = nullptr;
	int *start = nullptr, *pc = nullptr, *pl = nullptr;
	double *pts = nullptr, *chi = nullptr, *phi = nullptr, *P = nullptr, *val = nullptr, *grad = nullptr, *rho = nullptr;
	~bf_ctx()
	{
		gpuFree(c); gpuFree(cmin); gpuFree(pe); gpuFree(ps); gpuFree(coef); gpuFree(occ);
		gpuFree(start); gpuFree(pc); gpuFree(pl);
		gpuFree(pts); gpuFree(chi); gpuFree(phi); gpuFree(P); gpuFree(val); gpuFree(grad); gpuFree(rho);
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

	//Points per chunk: chi (K x nao), phi and the GEMM workspace (K x nocc each) per point, within half
	//of what is free.
	const int nprim = ao_start[nao];
	const size_t per_point = sizeof(double) * ((size_t)K * nao + 2 * (size_t)K * nocc + 8);
	size_t budget = freeb / 2;
	if (budget > ((size_t)4 << 30)) budget = (size_t)4 << 30;
	long long nb = (long long)(budget / per_point);
	if (nb > max_points) nb = max_points;
	//the GEMM indexes in int: keep K * nb * nao below 2^31
	const long long cap = ((1LL << 31) - 1) / ((long long)K * (nao > nocc ? nao : nocc));
	if (nb > cap) nb = cap;
	//gemm_gpu's reduce pass puts the K * nb rows on grid y in blocks of 16, which stops at 65535
	if (nb > 65535LL * 16 / K) nb = 65535LL * 16 / K;
	if (nb < 1) return nullptr;

	bf_ctx* x = new bf_ctx;
	x->K = K; x->nao = nao; x->nocc = nocc; x->chunk = (int)nb; x->cutoff = exp_cutoff;
	//split_count grows as the row count shrinks, so a short run can need more workspace than a full
	//chunk: take the largest m * splits over every row-tile count a run can have
	const long long mfull = (long long)K * x->chunk;
	size_t ws = 0;
	for (long long mt = 1; (mt - 1) * 64 < mfull; mt++)
	{
		const int m = (int)(mt * 64 < mfull ? mt * 64 : mfull);
		const size_t w = gemm_gpu::workspace_elems(m, nocc, gemm_gpu::split_count(m, nocc, nao));
		if (w > ws) ws = w;
	}
	const size_t ch = (size_t)x->chunk;
	const bool ok = upload(x->c, cxyz, 3 * (size_t)ncen) && upload(x->cmin, cmin_exp, (size_t)ncen)
		&& upload(x->start, ao_start, (size_t)nao + 1) && upload(x->pc, prim_center, (size_t)nprim)
		&& upload(x->pl, prim_l, 3 * (size_t)nprim) && upload(x->pe, prim_exp, (size_t)nprim)
		&& upload(x->ps, prim_scale, (size_t)nprim) && upload(x->coef, coef, (size_t)nao * nocc)
		&& upload(x->occ, occ, (size_t)nocc)
		&& gpuMalloc(&x->pts, sizeof(double) * 3 * ch) == gpuSuccess
		&& gpuMalloc(&x->chi, sizeof(double) * (size_t)K * ch * nao) == gpuSuccess
		&& gpuMalloc(&x->phi, sizeof(double) * (size_t)K * ch * nocc) == gpuSuccess
		&& gpuMalloc(&x->P, sizeof(double) * ws) == gpuSuccess
		&& gpuMalloc(&x->val, sizeof(double) * ch) == gpuSuccess
		&& gpuMalloc(&x->grad, sizeof(double) * 3 * ch) == gpuSuccess
		&& gpuMalloc(&x->rho, sizeof(double) * ch) == gpuSuccess;
	if (!ok)
	{
		std::fprintf(stderr, "NoSpherA2 basin field GPU: %s while allocating, falling back to the host\n", gpuGetErrorString(gpuGetLastError()));
		delete x;
		return nullptr;
	}
	return x;
}

bool basin_field_gpu_run(void* ctx, const int np, const double* pts, double* val, double* grad, double* rho)
{
	bf_ctx* x = (bf_ctx*)ctx;
	if (!x || np <= 0) return false;
	const int K = x->K, nao = x->nao, nocc = x->nocc, chunk = x->chunk;
	bool ok = true;
	gpuError_t err = gpuSuccess;
	for (int first = 0; first < np && ok; first += chunk)
	{
		const int n = np - first < chunk ? np - first : chunk;
		ok = gpuMemcpy(x->pts, pts + 3 * (size_t)first, sizeof(double) * 3 * (size_t)n, gpuMemcpyHostToDevice) == gpuSuccess;
		if (!ok) break;
		const long long work = (long long)n * nao;
		const unsigned blocks = (unsigned)((work + BF_BLOCK - 1) / BF_BLOCK);
		if (K == 4)
			bf_ao_kernel<4><<<blocks, BF_BLOCK>>>(n, nao, x->pts, x->c, x->cmin, x->cutoff, x->start, x->pc, x->pl, x->pe, x->ps, x->chi);
		else
			bf_ao_kernel<10><<<blocks, BF_BLOCK>>>(n, nao, x->pts, x->c, x->cmin, x->cutoff, x->start, x->pc, x->pl, x->pe, x->ps, x->chi);
		//phi (K n x nocc, column-major) = chi (K n x nao) * coef^T, coef being nocc x nao column-major
		gemm_gpu::launch<double>(false, true, K * n, nocc, nao, 1.0, x->chi, K * n, x->coef, nocc, 0.0, x->phi, K * n, x->P);
		const unsigned rb = (unsigned)((n + BF_BLOCK - 1) / BF_BLOCK);
		if (K == 4)
			bf_reduce_kernel<4><<<rb, BF_BLOCK>>>(n, nocc, x->occ, x->phi, x->val, x->grad, x->rho);
		else
			bf_reduce_kernel<10><<<rb, BF_BLOCK>>>(n, nocc, x->occ, x->phi, x->val, x->grad, x->rho);
		err = gpuGetLastError();
		if (err == gpuSuccess) err = gpuDeviceSynchronize();
		ok = err == gpuSuccess
			&& gpuMemcpy(grad + 3 * (size_t)first, x->grad, sizeof(double) * 3 * (size_t)n, gpuMemcpyDeviceToHost) == gpuSuccess
			&& (!val || K == 4 || gpuMemcpy(val + first, x->val, sizeof(double) * (size_t)n, gpuMemcpyDeviceToHost) == gpuSuccess)
			&& (!rho || gpuMemcpy(rho + first, x->rho, sizeof(double) * (size_t)n, gpuMemcpyDeviceToHost) == gpuSuccess);
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
