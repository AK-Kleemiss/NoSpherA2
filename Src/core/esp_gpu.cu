#include "esp_gpu.h"
#include "aux_density_gpu.h"
#include "gpu_backend.h"
#include <cstdio>

NOSPHERA2_GPU_API_BEGIN

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 ESP GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	esp_gpu_free(); return false; } } while (0)

#define ESP_BLOCK 128
#define ESP_CHUNK (1 << 20)
//constants::PI (spherical_harmonic.h), which the host branch of boys() uses
#define ESP_PI 3.1415926535897932384626433832795028

namespace {

//esp_axis_terms in wfn_density.cpp: the (l, r, s) count of one axis with l_i + l_j = L
__device__ int axis_terms(const int L)
{
	int n = 0;
	for (int l = 0; l <= L; l++)
		for (int r = 0; r <= l / 2; r++)
			n += (l - 2 * r) / 2 + 1;
	return n;
}

//the boys() of wfn_density.cpp, same table and same branches
__device__ double boys_dev(const int m, const double T, const double expn,
	const double* tab, const int nT, const int stride, const double step)
{
	if (T >= (nT - 1) * step)
	{
		double f = 0.5 * sqrt(ESP_PI / T);
		for (int n = 1; n <= m; n++)
			f = ((2 * n - 1) * f - expn) / (2 * T);
		return f;
	}
	const int i = (int)(T / step + 0.5);
	const double dT = i * step - T;
	const double* row = tab + (size_t)i * stride + m;
	double f = 0, pw = 1;
	for (int k = 0; k <= 5; k++, pw *= dT / k)
		f += row[k] * pw;
	return f;
}

__global__ void esp_kernel(
	const int n_at, const double* ax, const double* ay, const double* az, const double* q,
	const int npairs, const double* ex_sum, const double* weight, const double* P, const int* L,
	const int* off, const double* coef, const unsigned char* pc_pow, const unsigned char* fn_idx,
	const int nT, const int stride, const double step, const double* tab,
	const long long np, const double* points, double* out)
{
	const long long g = (long long)blockIdx.x * blockDim.x + threadIdx.x;
	if (g >= np) return;
	const double gx = points[3 * g], gy = points[3 * g + 1], gz = points[3 * g + 2];

	double ESP = 0;
	for (int a = 0; a < n_at; a++)
	{
		const double dx = gx - ax[a], dy = gy - ay[a], dz = gz - az[a];
		ESP += q[a] / sqrt(dx * dx + dy * dy + dz * dz);
	}

	double Fn[25], pcp[3][9];
	for (int p = 0; p < npairs; p++)
	{
		const double ex = ex_sum[p];
		const int L0 = L[3 * p], L1 = L[3 * p + 1], L2 = L[3 * p + 2];
		const int MaxFn = L0 + L1 + L2;
		const double PCx = P[3 * p] - gx, PCy = P[3 * p + 1] - gy, PCz = P[3 * p + 2] - gz;
		const double sqpc = PCx * PCx + PCy * PCy + PCz * PCz;
		const int c = off[p];
		if (MaxFn == 0)
		{
			//s-s: one term per axis, every PC power 0, and F_0 never touches exp(-T)
			ESP -= weight[p] * ((coef[c] * coef[c + 1]) * coef[c + 2] * boys_dev(0, ex * sqpc, 0.0, tab, nT, stride, step));
			continue;
		}
		pcp[0][0] = pcp[1][0] = pcp[2][0] = 1.0;
		for (int n = 1; n <= L0; n++) pcp[0][n] = pcp[0][n - 1] * PCx;
		for (int n = 1; n <= L1; n++) pcp[1][n] = pcp[1][n - 1] * PCy;
		for (int n = 1; n <= L2; n++) pcp[2][n] = pcp[2][n - 1] * PCz;

		const double expc = exp(-ex * sqpc);
		Fn[MaxFn] = boys_dev(MaxFn, ex * sqpc, expc, tab, nT, stride, step);
		const double twoexpc = 2 * ex * sqpc;
		for (int nu = MaxFn - 1; nu >= 0; nu--)
			Fn[nu] = (expc + twoexpc * Fn[nu + 1]) / (2 * (nu + 1) - 1);

		const int nl = axis_terms(L0), nm = axis_terms(L1), nn = axis_terms(L2);
		const int cl = c, cm = c + nl, cn = c + nl + nm;
		double term = 0.0;
		for (int l = 0; l < nl; l++)
		{
			const double Al = coef[cl + l] * pcp[0][pc_pow[cl + l]];
			const int ml = fn_idx[cl + l];
			for (int m = 0; m < nm; m++)
			{
				const double lm = Al * (coef[cm + m] * pcp[1][pc_pow[cm + m]]);
				const int mm = ml + fn_idx[cm + m];
				for (int n = 0; n < nn; n++)
					term += lm * (coef[cn + n] * pcp[2][pc_pow[cn + n]]) * Fn[mm + fn_idx[cn + n]];
			}
		}
		ESP -= weight[p] * term;
	}
	out[g] = ESP;
}

//The pair table is the same for every point set of a run (one wavefunction, one cube or one
//mesh), so it is uploaded once and kept until the table pointer changes or the process ends.
const double* held_coef = nullptr;
int held_npairs = 0, held_ncoef = 0, held_nat = 0;
double *d_ax = nullptr, *d_ay = nullptr, *d_az = nullptr, *d_q = nullptr;
double *d_ex = nullptr, *d_w = nullptr, *d_P = nullptr, *d_coef = nullptr, *d_tab = nullptr;
int *d_L = nullptr, *d_off = nullptr;
unsigned char *d_pcp = nullptr, *d_fni = nullptr;

void esp_gpu_free()
{
	gpuFree(d_ax); gpuFree(d_ay); gpuFree(d_az); gpuFree(d_q);
	gpuFree(d_ex); gpuFree(d_w); gpuFree(d_P); gpuFree(d_coef); gpuFree(d_tab);
	gpuFree(d_L); gpuFree(d_off); gpuFree(d_pcp); gpuFree(d_fni);
	d_ax = d_ay = d_az = d_q = d_ex = d_w = d_P = d_coef = d_tab = nullptr;
	d_L = d_off = nullptr;
	d_pcp = d_fni = nullptr;
	held_coef = nullptr;
	held_npairs = held_ncoef = held_nat = 0;
}
}

bool esp_gpu_eval(
	const int n_at, const double* ax, const double* ay, const double* az, const double* q,
	const int npairs, const double* ex_sum, const double* weight, const double* P, const int* L,
	const int* off, const double* coef, const unsigned char* pc_pow, const unsigned char* fn_idx,
	const int nT, const int stride, const double step, const double* boys_tab,
	const int np, const double* points, double* out)
{
	if (np <= 0 || npairs <= 0 || n_at <= 0) return false;
	if (!aux_density_gpu_enabled() || !aux_density_gpu_available()) return false;

	const int ncoef = off[npairs];
	const size_t pairs_bytes = sizeof(double) * (size_t)npairs * 5 + sizeof(int) * ((size_t)npairs * 3 + npairs + 1)
		+ sizeof(double) * (size_t)ncoef + 2 * (size_t)ncoef + sizeof(double) * (size_t)nT * stride;
	const size_t pt_bytes = sizeof(double) * (size_t)(np < ESP_CHUNK ? np : ESP_CHUNK) * 4;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	const bool reuse = (held_coef == coef && held_npairs == npairs && held_ncoef == ncoef && held_nat == n_at);
	if ((reuse ? 0 : pairs_bytes) + pt_bytes + (1u << 26) > freeb) return false;

	if (!reuse)
	{
		esp_gpu_free();
		const size_t at = sizeof(double) * (size_t)n_at;
		GPU_TRY(gpuMalloc(&d_ax, at)); GPU_TRY(gpuMalloc(&d_ay, at)); GPU_TRY(gpuMalloc(&d_az, at)); GPU_TRY(gpuMalloc(&d_q, at));
		GPU_TRY(gpuMalloc(&d_ex, sizeof(double) * (size_t)npairs));
		GPU_TRY(gpuMalloc(&d_w, sizeof(double) * (size_t)npairs));
		GPU_TRY(gpuMalloc(&d_P, sizeof(double) * (size_t)npairs * 3));
		GPU_TRY(gpuMalloc(&d_L, sizeof(int) * (size_t)npairs * 3));
		GPU_TRY(gpuMalloc(&d_off, sizeof(int) * (size_t)(npairs + 1)));
		GPU_TRY(gpuMalloc(&d_coef, sizeof(double) * (size_t)ncoef));
		GPU_TRY(gpuMalloc(&d_pcp, (size_t)ncoef)); GPU_TRY(gpuMalloc(&d_fni, (size_t)ncoef));
		GPU_TRY(gpuMalloc(&d_tab, sizeof(double) * (size_t)nT * stride));
		GPU_TRY(gpuMemcpy(d_ax, ax, at, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_ay, ay, at, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_az, az, at, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_q, q, at, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_ex, ex_sum, sizeof(double) * (size_t)npairs, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_w, weight, sizeof(double) * (size_t)npairs, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_P, P, sizeof(double) * (size_t)npairs * 3, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_L, L, sizeof(int) * (size_t)npairs * 3, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_off, off, sizeof(int) * (size_t)(npairs + 1), gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_coef, coef, sizeof(double) * (size_t)ncoef, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_pcp, pc_pow, (size_t)ncoef, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_fni, fn_idx, (size_t)ncoef, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_tab, boys_tab, sizeof(double) * (size_t)nT * stride, gpuMemcpyHostToDevice));
		held_coef = coef, held_npairs = npairs, held_ncoef = ncoef, held_nat = n_at;
	}

	//chunked, so the device side stays bounded however big the grid is and no single launch runs
	//long enough to meet a display driver's watchdog
	const int chunk = np < ESP_CHUNK ? np : ESP_CHUNK;
	double *d_pts = nullptr, *d_out = nullptr;
	if (gpuMalloc(&d_pts, sizeof(double) * (size_t)chunk * 3) != gpuSuccess) return false;
	if (gpuMalloc(&d_out, sizeof(double) * (size_t)chunk) != gpuSuccess) { gpuFree(d_pts); return false; }
	bool ok = true;
	for (int first = 0; first < np && ok; first += chunk)
	{
		const int n = np - first < chunk ? np - first : chunk;
		ok = gpuMemcpy(d_pts, points + (size_t)3 * first, sizeof(double) * (size_t)n * 3, gpuMemcpyHostToDevice) == gpuSuccess;
		if (!ok) break;
		esp_kernel<<<(unsigned)((n + ESP_BLOCK - 1) / ESP_BLOCK), ESP_BLOCK>>>(
			n_at, d_ax, d_ay, d_az, d_q, npairs, d_ex, d_w, d_P, d_L, d_off, d_coef, d_pcp, d_fni,
			nT, stride, step, d_tab, n, d_pts, d_out);
		ok = gpuGetLastError() == gpuSuccess && gpuDeviceSynchronize() == gpuSuccess
			&& gpuMemcpy(out + first, d_out, sizeof(double) * (size_t)n, gpuMemcpyDeviceToHost) == gpuSuccess;
	}
	gpuFree(d_pts); gpuFree(d_out);
	if (!ok)
	{
		std::fprintf(stderr, "NoSpherA2 ESP GPU: kernel failed, falling back to the host\n");
		esp_gpu_free();
	}
	return ok;
}

NOSPHERA2_GPU_API_END
