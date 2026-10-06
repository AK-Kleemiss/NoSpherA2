#include "salted_gpu.h"
#include "gpu_backend.h"
#include <cstdio>
#include <vector>
#include <algorithm>

NOSPHERA2_GPU_API_BEGIN

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 SALTED GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

//Sizes the per-thread transform vectors; l21 = 2*lam+1 is 17 for the v7 model, the margin is for
//larger models, and the entry point refuses more
#define SALTED_MAX_L21 64

namespace {

struct descriptor_cache {
	const double *v1_values = nullptr, *v2_values = nullptr;
	const size_t *v1_offsets = nullptr, *v2_offsets = nullptr;
	long long v1_len = 0, v2_len = 0;
	int v1_noff = 0, v2_noff = 0;
	bool conj = false, ready = false;
	double *d_v1 = nullptr, *d_v2 = nullptr;
	size_t *d_v1off = nullptr, *d_v2off = nullptr;
	void clear()
	{
		if (d_v1) gpuFree(d_v1);
		if (d_v2 && d_v2 != d_v1) gpuFree(d_v2);
		if (d_v1off) gpuFree(d_v1off);
		if (d_v2off && d_v2off != d_v1off) gpuFree(d_v2off);
		v1_values = nullptr; v2_values = nullptr; v1_offsets = nullptr; v2_offsets = nullptr;
		v1_len = 0; v2_len = 0; v1_noff = 0; v2_noff = 0; conj = false; ready = false;
		d_v1 = nullptr; d_v2 = nullptr; d_v1off = nullptr; d_v2off = nullptr;
	}
	bool matches(const salted_gpu_problem& q) const
	{
		return ready && v1_values == q.v1_values && v2_values == q.v2_values
			&& v1_offsets == q.v1_offsets && v2_offsets == q.v2_offsets
			&& v1_len == q.v1_len_doubles && v2_len == q.v2_len_doubles
			&& v1_noff == q.v1_noff && v2_noff == q.v2_noff && conj == q.v2_is_conj_of_v1;
	}
	bool upload(const salted_gpu_problem& q)
	{
		if (matches(q)) return true;
		clear();
		GPU_TRY(gpuMalloc(&d_v1, sizeof(double) * q.v1_len_doubles));
		GPU_TRY(gpuMemcpy(d_v1, q.v1_values, sizeof(double) * q.v1_len_doubles, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMalloc(&d_v1off, sizeof(size_t) * q.v1_noff));
		GPU_TRY(gpuMemcpy(d_v1off, q.v1_offsets, sizeof(size_t) * q.v1_noff, gpuMemcpyHostToDevice));
		if (q.v2_is_conj_of_v1) {
			d_v2 = d_v1;
			d_v2off = d_v1off;
		}
		else {
			GPU_TRY(gpuMalloc(&d_v2, sizeof(double) * q.v2_len_doubles));
			GPU_TRY(gpuMemcpy(d_v2, q.v2_values, sizeof(double) * q.v2_len_doubles, gpuMemcpyHostToDevice));
			GPU_TRY(gpuMalloc(&d_v2off, sizeof(size_t) * q.v2_noff));
			GPU_TRY(gpuMemcpy(d_v2off, q.v2_offsets, sizeof(size_t) * q.v2_noff, gpuMemcpyHostToDevice));
		}
		v1_values = q.v1_values; v2_values = q.v2_values; v1_offsets = q.v1_offsets; v2_offsets = q.v2_offsets;
		v1_len = q.v1_len_doubles; v2_len = q.v2_len_doubles; v1_noff = q.v1_noff; v2_noff = q.v2_noff;
		conj = q.v2_is_conj_of_v1; ready = true;
		return true;
	}
};

descriptor_cache g_descriptor_cache;

//block(atom, channel, l) as SALTEDDescriptors lays it out, in interleaved doubles
__device__ __forceinline__ const double* desc_block(const double* v, const size_t* off,
	const int nchannels, const int atom, const int channel, const int l)
{
	return v + 2 * (off[l] + ((size_t)atom * nchannels + channel) * (2 * (size_t)l + 1));
}

//One thread per (atom, output slot): only the nfps selected features are built, the norm comes
//from the host. Consecutive slots write consecutive p entries for every imu, so the stores coalesce.
__global__ void equicomb_kernel(const int natoms, const int nrad2, const int llmax,
	const int l21, const int shells, const int nfps, const bool conj,
	const double* __restrict__ v1, const size_t* __restrict__ v1_off, const int v1_nch,
	const double* __restrict__ v2, const size_t* __restrict__ v2_off, const int v2_nch,
	const double* __restrict__ w3j, const int* __restrict__ llvec0, const int* __restrict__ llvec1,
	const int* __restrict__ runs, const int* __restrict__ c2r_cols,
	const double* __restrict__ c2r_re, const double* __restrict__ c2r_im,
	const int* __restrict__ c2r_cnt, const int* __restrict__ vfps,
	const double* __restrict__ normfact, double* __restrict__ p)
{
	const int atom = blockIdx.y;
	const int slot = blockIdx.x * blockDim.x + threadIdx.x;
	if (atom >= natoms || slot >= nfps) return;
	double* out = p + (size_t)atom * l21 * nfps + slot;
	const int f = vfps[slot];
	//Features past nrad1*nrad2*llmax stay zero; the host zeroed p
	if (f >= shells) return;
	const int il = f % llmax, n2 = (f / llmax) % nrad2, n1 = f / (llmax * nrad2);
	const double* v1p = desc_block(v1, v1_off, v1_nch, atom, n1, llvec0[il]);
	const double* v2p = desc_block(v2, v2_off, v2_nch, atom, n2, llvec1[il]);
	double pc_re[SALTED_MAX_L21], pc_im[SALTED_MAX_L21];
	for (int imu = 0; imu < l21; imu++) {
		const int* run = runs + 4 * (il * l21 + imu);
		const int im1_begin = run[0], im2_begin = run[1], count = run[2], w_off = run[3];
		double acc_r = 0.0, acc_i = 0.0;
		for (int k = 0; k < count; k++) {
			const double wk = w3j[w_off + k];
			const double ar = wk * v1p[2 * (im1_begin + k)];
			const double ai = wk * v1p[2 * (im1_begin + k) + 1];
			const double br = v2p[2 * (im2_begin + k)];
			const double bi = v2p[2 * (im2_begin + k) + 1];
			if (conj) {
				acc_r += ar * br + ai * bi;
				acc_i += ai * br - ar * bi;
			}
			else {
				acc_r += ar * br - ai * bi;
				acc_i += ar * bi + ai * br;
			}
		}
		pc_re[imu] = acc_r;
		pc_im[imu] = acc_i;
	}
	//Two nonzeros per transform row, ascending columns, so the sum matches the CPU
	const double nf = normfact[atom];
	for (int i = 0; i < l21; i++) {
		double preal = 0.0;
		const int nz = c2r_cnt[i];
		for (int k = 0; k < nz; k++) {
			const int j = c2r_cols[2 * i + k];
			preal += c2r_re[2 * i + k] * pc_re[j] - c2r_im[2 * i + k] * pc_im[j];
		}
		out[(size_t)i * nfps] = preal * nf;
	}
}

}

bool salted_gpu_available()
{
	int n = 0;
	return gpuGetDeviceCount(&n) == gpuSuccess && n > 0;
}

void salted_gpu_clear_cache()
{
	g_descriptor_cache.clear();
}

bool salted_gpu_equicomb(const salted_gpu_problem& q)
{
	if (q.natoms <= 0 || q.shells <= 0 || q.nfps <= 0 || q.l21 <= 0) return false;
	if (q.l21 > SALTED_MAX_L21) return false;
	if (!salted_gpu_available()) return false;

	const size_t p_bytes = sizeof(double) * (size_t)q.natoms * q.nfps * q.l21;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	if (p_bytes + (1u << 27) > freeb) return false;

	double *d_w3j = nullptr, *d_c2r_re = nullptr, *d_c2r_im = nullptr, *d_nf = nullptr, *d_p = nullptr;
	int *d_ll0 = nullptr, *d_ll1 = nullptr, *d_runs = nullptr, *d_cols = nullptr, *d_cnt = nullptr, *d_vfps = nullptr;

	if (!g_descriptor_cache.upload(q)) return false;
	GPU_TRY(gpuMalloc(&d_w3j, sizeof(double) * (size_t)q.w3j_len));
	GPU_TRY(gpuMemcpy(d_w3j, q.w3j, sizeof(double) * (size_t)q.w3j_len, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_ll0, sizeof(int) * q.llmax));
	GPU_TRY(gpuMemcpy(d_ll0, q.llvec0, sizeof(int) * q.llmax, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_ll1, sizeof(int) * q.llmax));
	GPU_TRY(gpuMemcpy(d_ll1, q.llvec1, sizeof(int) * q.llmax, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_runs, sizeof(int) * 4 * (size_t)q.llmax * q.l21));
	GPU_TRY(gpuMemcpy(d_runs, q.runs, sizeof(int) * 4 * (size_t)q.llmax * q.l21, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_cols, sizeof(int) * 2 * (size_t)q.l21));
	GPU_TRY(gpuMemcpy(d_cols, q.c2r_cols, sizeof(int) * 2 * (size_t)q.l21, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_c2r_re, sizeof(double) * 2 * (size_t)q.l21));
	GPU_TRY(gpuMemcpy(d_c2r_re, q.c2r_re, sizeof(double) * 2 * (size_t)q.l21, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_c2r_im, sizeof(double) * 2 * (size_t)q.l21));
	GPU_TRY(gpuMemcpy(d_c2r_im, q.c2r_im, sizeof(double) * 2 * (size_t)q.l21, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_cnt, sizeof(int) * (size_t)q.l21));
	GPU_TRY(gpuMemcpy(d_cnt, q.c2r_cnt, sizeof(int) * (size_t)q.l21, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_vfps, sizeof(int) * (size_t)q.nfps));
	GPU_TRY(gpuMemcpy(d_vfps, q.vfps, sizeof(int) * (size_t)q.nfps, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_nf, sizeof(double) * (size_t)q.natoms));
	GPU_TRY(gpuMemcpy(d_nf, q.normfact, sizeof(double) * (size_t)q.natoms, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_p, p_bytes));
	GPU_TRY(gpuMemset(d_p, 0, p_bytes));

	const dim3 thr(128), grid((q.nfps + 127) / 128, q.natoms);
	equicomb_kernel<<<grid, thr>>>(q.natoms, q.nrad2, q.llmax, q.l21,
		q.shells, q.nfps, q.v2_is_conj_of_v1,
		g_descriptor_cache.d_v1, g_descriptor_cache.d_v1off, q.v1_nchannels,
		g_descriptor_cache.d_v2, g_descriptor_cache.d_v2off, q.v2_nchannels,
		d_w3j, d_ll0, d_ll1, d_runs, d_cols, d_c2r_re, d_c2r_im, d_cnt, d_vfps, d_nf, d_p);
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuDeviceSynchronize());

	GPU_TRY(gpuMemcpy(q.p, d_p, p_bytes, gpuMemcpyDeviceToHost));

	gpuFree(d_w3j); gpuFree(d_ll0); gpuFree(d_ll1); gpuFree(d_runs);
	gpuFree(d_cols); gpuFree(d_c2r_re); gpuFree(d_c2r_im); gpuFree(d_cnt);
	gpuFree(d_vfps); gpuFree(d_nf); gpuFree(d_p);
	return true;
}

NOSPHERA2_GPU_API_END
