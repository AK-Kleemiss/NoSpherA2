#include "xcw_grid_gpu.h"
#include "gpu_backend.h"
#include "gemm_gpu.cuh"
#include "sf_gpu.h"
#include <algorithm>
#include <cstdio>
#include <vector>

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 XCW grid GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

NOSPHERA2_GPU_API_BEGIN

namespace {

//Block layout on the host, launched from; everything a kernel reads on the device
struct held_grid {
	bool on = false;
	int nmo = 0, nb = 0, npt = 0;
	std::vector<int> na, np, pt_off, ao_off;
	std::vector<long long> x_off;
	float* X = nullptr, * Df = nullptr, * T = nullptr, * P = nullptr;
	int* aos = nullptr;
	unsigned char* keep = nullptr;
	double* w = nullptr, * Dd = nullptr, * rho = nullptr, * out = nullptr;
};
held_grid G;

//D read on its upper triangle and zeroed where the pair is screened, mirrored: the block GEMM takes it whole
__global__ void masked_D_kernel(const int nmo, const double* __restrict__ D, const unsigned char* __restrict__ keep, float* __restrict__ Df)
{
	const long long e = blockIdx.x * (long long)blockDim.x + threadIdx.x;
	if (e >= (long long)nmo * nmo) return;
	const int i = (int)(e / nmo), j = (int)(e % nmo);
	const long long u = (long long)min(i, j) * nmo + max(i, j);
	Df[e] = keep[u] ? (float)D[u] : 0.0f;
}

__global__ void gather_kernel(const int nmo, const int na, const int* __restrict__ a, const float* __restrict__ Df, float* __restrict__ Db)
{
	const int e = blockIdx.x * blockDim.x + threadIdx.x;
	if (e >= na * na) return;
	Db[e] = Df[(long long)a[e / na] * nmo + a[e % na]];
}

//rho_p = 2 w_p Sum_i X_ip T_ip, summed in double as the host does
__global__ void rho_kernel(const int na, const int np, const float* __restrict__ X, const float* __restrict__ T,
	const double* __restrict__ w, double* __restrict__ rho)
{
	const int p = blockIdx.x * blockDim.x + threadIdx.x;
	if (p >= np) return;
	double s = 0;
	for (int i = 0; i < na; i++) s += (double)X[(long long)i * np + p] * T[(long long)i * np + p];
	rho[p] = 2.0 * w[p] * s;
}

__global__ void scale_kernel(const int na, const int np, const float* __restrict__ X, const double* __restrict__ V, float* __restrict__ Y)
{
	const long long e = blockIdx.x * (long long)blockDim.x + threadIdx.x;
	if (e >= (long long)na * np) return;
	Y[e] = X[e] * (float)V[e % np];
}

//Each (i, j) of a block is its own (mu, nu), and blocks run in stream order: no atomics, the same sum every run
__global__ void scatter_kernel(const int nmo, const int na, const int* __restrict__ a, const unsigned char* __restrict__ keep,
	const float* __restrict__ S, double* __restrict__ out)
{
	const int e = blockIdx.x * blockDim.x + threadIdx.x;
	if (e >= na * na) return;
	const int i = e / na, j = e % na;
	if (j < i) return;
	const long long u = (long long)a[i] * nmo + a[j];
	if (keep[u]) out[u] += (double)S[e];
}

unsigned int grid_of(const long long n, const int b) { return (unsigned int)((n + b - 1) / b); }

} //namespace

void xcw_grid_gpu_release()
{
	for (void* p : { (void*)G.X, (void*)G.Df, (void*)G.T, (void*)G.P, (void*)G.aos, (void*)G.keep, (void*)G.w, (void*)G.Dd, (void*)G.rho, (void*)G.out })
		if (p) gpuFree(p);
	G = held_grid{};
}

bool xcw_grid_gpu_held() { return G.on; }

bool xcw_grid_gpu_hold(const int nmo, const int nb, const int* na, const int* np, const int* pt_off, const int* ao_off, const int* aos,
	const long long* x_off, const float* X, const long long nx, const unsigned char* keep, const double* w, const int npt)
{
	xcw_grid_gpu_release();
	if (nmo <= 0 || nb <= 0 || !sf_gpu_available()) return false;
	G.nmo = nmo; G.nb = nb; G.npt = npt;
	G.na.assign(na, na + nb); G.np.assign(np, np + nb); G.pt_off.assign(pt_off, pt_off + nb);
	G.ao_off.assign(ao_off, ao_off + nb); G.x_off.assign(x_off, x_off + nb);
	//T holds a block's na x np product (then its V-scaled values), P the GEMM slices of either shape
	long long t_max = 1, p_max = 1;
	for (int i = 0; i < nb; i++) {
		t_max = std::max(t_max, (long long)na[i] * std::max(np[i], na[i]));
		const int s1 = gemm_gpu::split_count(np[i], na[i], na[i]), s2 = gemm_gpu::split_count(na[i], na[i], np[i]);
		p_max = std::max({ p_max, (long long)gemm_gpu::workspace_elems(np[i], na[i], s1), (long long)gemm_gpu::workspace_elems(na[i], na[i], s2) });
	}
	const long long nao = (long long)ao_off[nb - 1] + na[nb - 1], nn = (long long)nmo * nmo;
	const size_t bytes = sizeof(float) * (nx + nn + 2 * t_max + p_max) + sizeof(int) * nao + nn + sizeof(double) * (2 * (size_t)npt + 2 * nn);
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess || bytes + (1u << 28) > freeb) return false;
	const auto fail = [] { xcw_grid_gpu_release(); return false; };
	if (gpuMalloc(&G.X, sizeof(float) * nx) != gpuSuccess || gpuMalloc(&G.Df, sizeof(float) * nn) != gpuSuccess
		|| gpuMalloc(&G.T, sizeof(float) * 2 * t_max) != gpuSuccess || gpuMalloc(&G.P, sizeof(float) * p_max) != gpuSuccess
		|| gpuMalloc(&G.aos, sizeof(int) * nao) != gpuSuccess || gpuMalloc(&G.keep, nn) != gpuSuccess
		|| gpuMalloc(&G.w, sizeof(double) * npt) != gpuSuccess || gpuMalloc(&G.Dd, sizeof(double) * nn) != gpuSuccess
		|| gpuMalloc(&G.rho, sizeof(double) * npt) != gpuSuccess || gpuMalloc(&G.out, sizeof(double) * nn) != gpuSuccess)
		return fail();
	if (gpuMemcpy(G.X, X, sizeof(float) * nx, gpuMemcpyHostToDevice) != gpuSuccess
		|| gpuMemcpy(G.aos, aos, sizeof(int) * nao, gpuMemcpyHostToDevice) != gpuSuccess
		|| gpuMemcpy(G.keep, keep, nn, gpuMemcpyHostToDevice) != gpuSuccess
		|| gpuMemcpy(G.w, w, sizeof(double) * npt, gpuMemcpyHostToDevice) != gpuSuccess)
		return fail();
	G.on = true;
	return true;
}

bool xcw_grid_gpu_density(const double* D, double* rho)
{
	if (!G.on) return false;
	const long long nn = (long long)G.nmo * G.nmo;
	GPU_TRY(gpuMemcpy(G.Dd, D, sizeof(double) * nn, gpuMemcpyHostToDevice));
	//points of blocks that kept no AO are no block's, and V was here last
	GPU_TRY(gpuMemset(G.rho, 0, sizeof(double) * G.npt));
	masked_D_kernel<<<grid_of(nn, 256), 256>>>(G.nmo, G.Dd, G.keep, G.Df);
	float* Db = G.T + 0, * T = nullptr;
	for (int i = 0; i < G.nb; i++) {
		const int na = G.na[i], np = G.np[i];
		const int* a = G.aos + G.ao_off[i];
		const float* X = G.X + G.x_off[i];
		//Db first, T after it: na^2 + na np fits the 2 t_max held
		T = Db + (long long)na * na;
		gather_kernel<<<grid_of((long long)na * na, 256), 256>>>(G.nmo, na, a, G.Df, Db);
		//row-major T (na x np) = Db X is column-major T^T (np x na) = X^T Db, Db symmetric
		gemm_gpu::launch<float>(false, false, np, na, na, 1.0f, X, np, Db, na, 0.0f, T, np, G.P);
		rho_kernel<<<grid_of(np, 128), 128>>>(na, np, X, T, G.w + G.pt_off[i], G.rho + G.pt_off[i]);
	}
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuMemcpy(rho, G.rho, sizeof(double) * G.npt, gpuMemcpyDeviceToHost));
	return true;
}

bool xcw_grid_gpu_contract(const double* V, double* out)
{
	if (!G.on) return false;
	const long long nn = (long long)G.nmo * G.nmo;
	//V goes where rho came back from: the two are never live together
	GPU_TRY(gpuMemcpy(G.rho, V, sizeof(double) * G.npt, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemset(G.out, 0, sizeof(double) * nn));
	for (int i = 0; i < G.nb; i++) {
		const int na = G.na[i], np = G.np[i];
		const int* a = G.aos + G.ao_off[i];
		const float* X = G.X + G.x_off[i];
		float* Y = G.T, * S = G.T + (long long)na * np;
		scale_kernel<<<grid_of((long long)na * np, 256), 256>>>(na, np, X, G.rho + G.pt_off[i], Y);
		//row-major S (na x na) = X Y^T, symmetric, so column-major alike: op(A) = X from X^T stored, op(B) = Y^T
		gemm_gpu::launch<float>(true, false, na, na, np, 1.0f, X, np, Y, np, 0.0f, S, na, G.P);
		scatter_kernel<<<grid_of((long long)na * na, 256), 256>>>(G.nmo, na, a, G.keep, S, G.out);
	}
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuMemcpy(out, G.out, sizeof(double) * nn, gpuMemcpyDeviceToHost));
	return true;
}

NOSPHERA2_GPU_API_END
