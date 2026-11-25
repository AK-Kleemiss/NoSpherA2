#include "salted_gpu.h"
#include "itensor_gpu.h"
#include "gpu_backend.h"
#include <cstdio>
#include <vector>
#include <algorithm>

NOSPHERA2_GPU_API_BEGIN

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 SALTED GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

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
//For real densities pc[-mu] = (-1)^mu conj(pc[mu]), so only mu >= 0 is summed and imu = lam +- mu
//is Re(cpm[imu] pc[mu]); in pc[0] the terms of m1 and -m1 are complex conjugates.
__global__ void equicomb_kernel(const int natoms, const int nrad2, const int llmax,
	const int lam, const int shells, const int nfps, const bool conj,
	const double* __restrict__ v1, const size_t* __restrict__ v1_off, const int v1_nch,
	const double* __restrict__ v2, const size_t* __restrict__ v2_off, const int v2_nch,
	const double* __restrict__ w3j, const int* __restrict__ llvec0, const int* __restrict__ llvec1,
	const int* __restrict__ runs, const double* __restrict__ cpm_re, const double* __restrict__ cpm_im,
	const int* __restrict__ vfps, const double* __restrict__ normfact, double* __restrict__ p)
{
	const int atom = blockIdx.y;
	const int slot = blockIdx.x * blockDim.x + threadIdx.x;
	if (atom >= natoms || slot >= nfps) return;
	const int l21 = 2 * lam + 1;
	double* out = p + (size_t)atom * l21 * nfps + slot;
	const int f = vfps[slot];
	//Features past nrad1*nrad2*llmax stay zero; the device buffer was zeroed
	if (f >= shells) return;
	const int il = f % llmax, n2 = (f / llmax) % nrad2, n1 = f / (llmax * nrad2);
	const double* v1p = desc_block(v1, v1_off, v1_nch, atom, n1, llvec0[il]);
	const double* v2p = desc_block(v2, v2_off, v2_nch, atom, n2, llvec1[il]);
	const double nf = normfact[atom], sg = conj ? 1.0 : -1.0;
	{
		const int* run = runs + 4 * (il * l21 + lam);
		const double* a = v1p + 2 * run[0];
		const double* b = v2p + 2 * run[1];
		const double* w = w3j + run[3];
		const int h = run[2] / 2;
		double x = 0.0;
		for (int k = 0; k < h; k++)
			x += w[k] * (a[2 * k] * b[2 * k] + sg * a[2 * k + 1] * b[2 * k + 1]);
		out[(size_t)lam * nfps] = cpm_re[lam] * nf * (2.0 * x + w[h] * a[2 * h] * b[2 * h]);
	}
	for (int mu = 1; mu <= lam; mu++) {
		const int* run = runs + 4 * (il * l21 + lam + mu);
		const double* a = v1p + 2 * run[0];
		const double* b = v2p + 2 * run[1];
		const double* w = w3j + run[3];
		double sr = 0.0, si = 0.0;
		for (int k = 0; k < run[2]; k++) {
			const double ar = w[k] * a[2 * k], ai = w[k] * a[2 * k + 1];
			sr += ar * b[2 * k] + sg * ai * b[2 * k + 1];
			si += ai * b[2 * k] - sg * ar * b[2 * k + 1];
		}
		out[(size_t)(lam + mu) * nfps] = nf * (cpm_re[lam + mu] * sr - cpm_im[lam + mu] * si);
		out[(size_t)(lam - mu) * nfps] = nf * (cpm_re[lam - mu] * sr - cpm_im[lam - mu] * si);
	}
}

}

bool salted_gpu_available()
{
	int n = 0;
	if (gpuGetDeviceCount(&n) != gpuSuccess || n <= 0) return false;
	//Doubles at 1/32 rate over the host's memory bus: the Radeon 780M took 3.67 s for 720 atoms, its 8700G 0.70 s.
	//Once per process, as the properties query is slow and this is asked per lambda
	static const bool apu = itensor_gpu_integrated() && sf_gpu_fp64_ratio() > 4;
	return !apu;
}

void salted_gpu_clear_cache()
{
	g_descriptor_cache.clear();
}

bool salted_gpu_equicomb(const salted_gpu_problem& q)
{
	if (q.natoms <= 0 || q.shells <= 0 || q.nfps <= 0 || q.l21 <= 0) return false;
	if (!salted_gpu_available()) return false;

	const size_t p_bytes = sizeof(double) * (size_t)q.natoms * q.nfps * q.l21;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	if (p_bytes + (1u << 27) > freeb) return false;

	double *d_w3j = nullptr, *d_cpm = nullptr, *d_nf = nullptr, *d_p = nullptr;
	int *d_ll0 = nullptr, *d_ll1 = nullptr, *d_runs = nullptr, *d_vfps = nullptr;

	if (!g_descriptor_cache.upload(q)) return false;
	GPU_TRY(gpuMalloc(&d_w3j, sizeof(double) * (size_t)q.w3j_len));
	GPU_TRY(gpuMemcpy(d_w3j, q.w3j, sizeof(double) * (size_t)q.w3j_len, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_ll0, sizeof(int) * q.llmax));
	GPU_TRY(gpuMemcpy(d_ll0, q.llvec0, sizeof(int) * q.llmax, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_ll1, sizeof(int) * q.llmax));
	GPU_TRY(gpuMemcpy(d_ll1, q.llvec1, sizeof(int) * q.llmax, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_runs, sizeof(int) * 4 * (size_t)q.llmax * q.l21));
	GPU_TRY(gpuMemcpy(d_runs, q.runs, sizeof(int) * 4 * (size_t)q.llmax * q.l21, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_cpm, sizeof(double) * 2 * (size_t)q.l21));
	GPU_TRY(gpuMemcpy(d_cpm, q.cpm_re, sizeof(double) * (size_t)q.l21, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(d_cpm + q.l21, q.cpm_im, sizeof(double) * (size_t)q.l21, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_vfps, sizeof(int) * (size_t)q.nfps));
	GPU_TRY(gpuMemcpy(d_vfps, q.vfps, sizeof(int) * (size_t)q.nfps, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_nf, sizeof(double) * (size_t)q.natoms));
	GPU_TRY(gpuMalloc(&d_p, p_bytes));
	GPU_TRY(gpuMemset(d_p, 0, p_bytes));

	GPU_TRY(gpuMemcpy(d_nf, q.normfact, sizeof(double) * (size_t)q.natoms, gpuMemcpyHostToDevice));

	const dim3 thr(128), grid((q.nfps + 127) / 128, q.natoms);
	equicomb_kernel<<<grid, thr>>>(q.natoms, q.nrad2, q.llmax, q.lam,
		q.shells, q.nfps, q.v2_is_conj_of_v1,
		g_descriptor_cache.d_v1, g_descriptor_cache.d_v1off, q.v1_nchannels,
		g_descriptor_cache.d_v2, g_descriptor_cache.d_v2off, q.v2_nchannels,
		d_w3j, d_ll0, d_ll1, d_runs, d_cpm, d_cpm + q.l21, d_vfps, d_nf, d_p);
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuDeviceSynchronize());

	//Pin p in place for the copy: a pageable copy runs through the driver's single-threaded staging,
	//bound by its slow host memcpy, and registering plus the pinned copy beats it on Linux nodes.
	//A pinned staging buffer of our own was no faster, its host-side memcpy is the same bottleneck.
	//Not on Windows: under WDDM the pageable copy already runs at the PCIe rate and a copy into
	//registered memory is slower
#ifdef _WIN32
	const bool pinned = false;
#else
	const bool pinned = gpuHostRegister(q.p, p_bytes) == gpuSuccess;
#endif
	if (!pinned) (void)gpuGetLastError();
	const gpuError_t copied = gpuMemcpy(q.p, d_p, p_bytes, gpuMemcpyDeviceToHost);
	if (pinned) gpuHostUnregister(q.p);
	GPU_TRY(copied);

	gpuFree(d_w3j); gpuFree(d_ll0); gpuFree(d_ll1); gpuFree(d_runs);
	gpuFree(d_cpm);
	gpuFree(d_vfps); gpuFree(d_nf); gpuFree(d_p);
	return true;
}

NOSPHERA2_GPU_API_END
