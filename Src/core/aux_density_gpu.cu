#include "aux_density_gpu.h"
#include "aux_density.h"
#include "gpu_backend.h"
#include <cstdio>

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 aux density GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

//Below this many point-shell pairs the host loop finishes before the copies do
#define AUX_MIN_WORK (1 << 22)
#define AUX_BLOCK 128

static bool g_aux_use_gpu = false;
void aux_density_gpu_set_enabled(bool on) { g_aux_use_gpu = on; }
bool aux_density_gpu_enabled() { return g_aux_use_gpu; }

bool aux_density_gpu_available()
{
	int n = 0;
	return gpuGetDeviceCount(&n) == gpuSuccess && n > 0;
}

namespace {
__global__ void aux_density_kernel(const int np, const double* x, const double* y, const double* z, const int n_at,
	const double* cx, const double* cy, const double* cz, const double* r2_max,
	const int* sh_start, const int* sh_l, const int* pr_start, const int* coef_off,
	const double* pr_exp, const double* pr_norm, const double* coefs, double* rho)
{
	const int p = blockIdx.x * blockDim.x + threadIdx.x;
	if (p >= np) return;
	rho[p] = aux_density::at(x[p], y[p], z[p], n_at, cx, cy, cz, r2_max, sh_start, sh_l, pr_start, coef_off, pr_exp, pr_norm, coefs);
}
__global__ void aux_density_grad_kernel(const int np, const double* x, const double* y, const double* z, const int n_at,
	const double* cx, const double* cy, const double* cz, const double* r2_max,
	const int* sh_start, const int* sh_l, const int* pr_start, const int* coef_off,
	const double* pr_exp, const double* pr_norm, const double* coefs, double* rho, double* gx, double* gy, double* gz)
{
	const int p = blockIdx.x * blockDim.x + threadIdx.x;
	if (p >= np) return;
	rho[p] = aux_density::at_grad(x[p], y[p], z[p], n_at, cx, cy, cz, r2_max, sh_start, sh_l, pr_start, coef_off, pr_exp, pr_norm, coefs, gx[p], gy[p], gz[p]);
}
__global__ void aux_density_lap_kernel(const int np, const double* x, const double* y, const double* z, const int n_at,
	const double* cx, const double* cy, const double* cz, const double* r2_max,
	const int* sh_start, const int* sh_l, const int* pr_start, const int* coef_off,
	const double* pr_exp, const double* pr_norm, const double* coefs, double* rho, double* gx, double* gy, double* gz, double* lap)
{
	const int p = blockIdx.x * blockDim.x + threadIdx.x;
	if (p >= np) return;
	rho[p] = aux_density::at_lap(x[p], y[p], z[p], n_at, cx, cy, cz, r2_max, sh_start, sh_l, pr_start, coef_off, pr_exp, pr_norm, coefs, gx[p], gy[p], gz[p], lap[p]);
}
}

bool aux_density_gpu_eval(
	const int n_at, const double* cx, const double* cy, const double* cz, const double* r2_max,
	const int n_sh, const int* sh_start, const int* sh_l, const int* pr_start, const int* coef_off,
	const int n_pr, const double* pr_exp, const double* pr_norm,
	const int n_coef, const double* coefs,
	const int np, const double* x, const double* y, const double* z, double* rho, double* gx, double* gy, double* gz, double* lap)
{
	const bool grad = gx != nullptr, lp = lap != nullptr;
	if (np <= 0 || n_at <= 0 || n_sh <= 0) return false;
	if (!g_aux_use_gpu || !aux_density_gpu_available()) return false;
	if ((long long)np * n_sh < AUX_MIN_WORK) return false;
	const size_t pts = sizeof(double) * (size_t)np, at = sizeof(double) * (size_t)n_at;
	const size_t sh = sizeof(int) * (size_t)n_sh, pr = sizeof(double) * (size_t)n_pr;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	if ((lp ? 8 : grad ? 7 : 4) * pts + 4 * at + 5 * sh + 2 * pr + sizeof(double) * (size_t)n_coef + (1u << 26) > freeb) return false;

	double *dx = nullptr, *dy = nullptr, *dz = nullptr, *drho = nullptr, *dgx = nullptr, *dgy = nullptr, *dgz = nullptr, *dlap = nullptr;
	double *dcx = nullptr, *dcy = nullptr, *dcz = nullptr, *dr2 = nullptr, *dexp = nullptr, *dnorm = nullptr, *dcoef = nullptr;
	int *dss = nullptr, *dsl = nullptr, *dps = nullptr, *dco = nullptr;
	GPU_TRY(gpuMalloc(&dx, pts)); GPU_TRY(gpuMalloc(&dy, pts)); GPU_TRY(gpuMalloc(&dz, pts)); GPU_TRY(gpuMalloc(&drho, pts));
	GPU_TRY(gpuMalloc(&dcx, at)); GPU_TRY(gpuMalloc(&dcy, at)); GPU_TRY(gpuMalloc(&dcz, at)); GPU_TRY(gpuMalloc(&dr2, at));
	GPU_TRY(gpuMalloc(&dss, sh + sizeof(int))); GPU_TRY(gpuMalloc(&dsl, sh)); GPU_TRY(gpuMalloc(&dps, sh + sizeof(int))); GPU_TRY(gpuMalloc(&dco, sh));
	GPU_TRY(gpuMalloc(&dexp, pr)); GPU_TRY(gpuMalloc(&dnorm, pr));
	GPU_TRY(gpuMalloc(&dcoef, sizeof(double) * (size_t)n_coef));
	if (grad) { GPU_TRY(gpuMalloc(&dgx, pts)); GPU_TRY(gpuMalloc(&dgy, pts)); GPU_TRY(gpuMalloc(&dgz, pts)); }
	if (lp) GPU_TRY(gpuMalloc(&dlap, pts));

	GPU_TRY(gpuMemcpy(dx, x, pts, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dy, y, pts, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dz, z, pts, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dcx, cx, at, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dcy, cy, at, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dcz, cz, at, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dr2, r2_max, at, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dss, sh_start, sizeof(int) * (size_t)(n_at + 1), gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dsl, sh_l, sh, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dps, pr_start, sh + sizeof(int), gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dco, coef_off, sh, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dexp, pr_exp, pr, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dnorm, pr_norm, pr, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(dcoef, coefs, sizeof(double) * (size_t)n_coef, gpuMemcpyHostToDevice));

	if (lp) aux_density_lap_kernel<<<(np + AUX_BLOCK - 1) / AUX_BLOCK, AUX_BLOCK>>>(np, dx, dy, dz, n_at, dcx, dcy, dcz, dr2, dss, dsl, dps, dco, dexp, dnorm, dcoef, drho, dgx, dgy, dgz, dlap);
	else if (grad) aux_density_grad_kernel<<<(np + AUX_BLOCK - 1) / AUX_BLOCK, AUX_BLOCK>>>(np, dx, dy, dz, n_at, dcx, dcy, dcz, dr2, dss, dsl, dps, dco, dexp, dnorm, dcoef, drho, dgx, dgy, dgz);
	else aux_density_kernel<<<(np + AUX_BLOCK - 1) / AUX_BLOCK, AUX_BLOCK>>>(np, dx, dy, dz, n_at, dcx, dcy, dcz, dr2, dss, dsl, dps, dco, dexp, dnorm, dcoef, drho);
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuDeviceSynchronize());
	GPU_TRY(gpuMemcpy(rho, drho, pts, gpuMemcpyDeviceToHost));
	if (grad) {
		GPU_TRY(gpuMemcpy(gx, dgx, pts, gpuMemcpyDeviceToHost));
		GPU_TRY(gpuMemcpy(gy, dgy, pts, gpuMemcpyDeviceToHost));
		GPU_TRY(gpuMemcpy(gz, dgz, pts, gpuMemcpyDeviceToHost));
		gpuFree(dgx); gpuFree(dgy); gpuFree(dgz);
	}
	if (lp) { GPU_TRY(gpuMemcpy(lap, dlap, pts, gpuMemcpyDeviceToHost)); gpuFree(dlap); }

	gpuFree(dx); gpuFree(dy); gpuFree(dz); gpuFree(drho);
	gpuFree(dcx); gpuFree(dcy); gpuFree(dcz); gpuFree(dr2);
	gpuFree(dss); gpuFree(dsl); gpuFree(dps); gpuFree(dco);
	gpuFree(dexp); gpuFree(dnorm); gpuFree(dcoef);
	return true;
}
