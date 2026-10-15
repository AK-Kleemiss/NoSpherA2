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
	//Density matrices, lambda-independent; d_dm2 is d_dm1 when u's are v1's (share)
	double *d_dm1 = nullptr, *d_dm2 = nullptr;
	int dm_nrad1 = -1, dm_nrad2 = -1, dm_lmax1 = -1, dm_lmax2 = -1;
	void clear_dm()
	{
		if (d_dm2 && d_dm2 != d_dm1) gpuFree(d_dm2);
		if (d_dm1) gpuFree(d_dm1);
		d_dm1 = nullptr; d_dm2 = nullptr;
		dm_nrad1 = -1; dm_nrad2 = -1; dm_lmax1 = -1; dm_lmax2 = -1;
	}
	void clear()
	{
		clear_dm();
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

//Start of the l block in one atom's matrices: sum_{k<l} (2k+1)^2, as on the host
__host__ __device__ __forceinline__ int dm_offset(const int l) { return l * (2 * l - 1) * (2 * l + 1) / 3; }

//A = sum_n x conj(x') and B = sum_n x x' per atom and l, one block per atom, one thread per matrix
//entry; layout [atom][are, aim, bre, bim][l blocks] as build_density on the host
__global__ void density_kernel(const int lmax, const int nch, const double* __restrict__ v,
	const size_t* __restrict__ off, const int v_nch, double* __restrict__ out)
{
	const int atom = blockIdx.x;
	const int dsz = dm_offset(lmax + 1);
	double* m = out + (size_t)atom * 4 * dsz;
	for (int e = threadIdx.x; e < dsz; e += blockDim.x) {
		int l = 0;
		while (dm_offset(l + 1) <= e) ++l;
		const int nm = 2 * l + 1, r = e - dm_offset(l), a = r / nm, b = r % nm;
		double ar = 0.0, ai = 0.0, br = 0.0, bi = 0.0;
		for (int n = 0; n < nch; n++) {
			const double* x = desc_block(v, off, v_nch, atom, n, l);
			const double xar = x[2 * a], xai = x[2 * a + 1], xbr = x[2 * b], xbi = x[2 * b + 1];
			ar += xar * xbr + xai * xbi;
			ai += xai * xbr - xar * xbi;
			br += xar * xbr - xai * xbi;
			bi += xai * xbr + xar * xbi;
		}
		m[e] = ar; m[dsz + e] = ai; m[2 * dsz + e] = br; m[3 * dsz + e] = bi;
	}
}

//sum_i Re(q_i)^2 = Re(sum K P + sum G Q) / 2 per atom, the host's pair sum: one block per atom,
//threads over (il, mu), tree reduction, normfact = 1/sqrt or 0 for an empty environment
#define NORM_THREADS 128
__global__ void norm_kernel(const int llmax, const int l21, const double s2,
	const double* __restrict__ dm1, const int dsz1, const double* __restrict__ dm2, const int dsz2,
	const double* __restrict__ w3j, const int* __restrict__ llvec0, const int* __restrict__ llvec1,
	const int* __restrict__ runs, const double* __restrict__ K_re, const double* __restrict__ K_im,
	const double* __restrict__ G_re, const double* __restrict__ G_im, double* __restrict__ normfact)
{
	__shared__ double red[NORM_THREADS];
	const int atom = blockIdx.x;
	const double* m1 = dm1 + (size_t)atom * 4 * dsz1;
	const double* m2 = dm2 + (size_t)atom * 4 * dsz2;
	double t = 0.0;
	for (int w = threadIdx.x; w < llmax * l21; w += blockDim.x) {
		const int il = w / l21, mu = w % l21;
		const int* r = runs + 4 * w;
		if (r[2] == 0) continue;
		const int l1 = llvec0[il], l2 = llvec1[il], nm1 = 2 * l1 + 1, nm2 = 2 * l2 + 1;
		const int o1 = dm_offset(l1), o2 = dm_offset(l2);
		const double *a1r = m1 + o1, *a1i = m1 + dsz1 + o1, *b1r = m1 + 2 * dsz1 + o1, *b1i = m1 + 3 * dsz1 + o1;
		const double *a2r = m2 + o2, *a2i = m2 + dsz2 + o2, *b2r = m2 + 2 * dsz2 + o2, *b2i = m2 + 3 * dsz2 + o2;
		for (int mu2 = 0; mu2 < l21; mu2++) {
			const int mm = mu * l21 + mu2;
			const int* r2 = runs + 4 * (il * l21 + mu2);
			if ((K_re[mm] == 0.0 && K_im[mm] == 0.0 && G_re[mm] == 0.0 && G_im[mm] == 0.0) || r2[2] == 0) continue;
			double Pr = 0.0, Pi = 0.0, Qr = 0.0, Qi = 0.0;
			for (int x = 0; x < r[2]; x++) {
				const double wx = w3j[r[3] + x];
				const int i1 = (r[0] + x) * nm1 + r2[0];
				const int i2 = (r[1] + x) * nm2 + r2[1];
				for (int y = 0; y < r2[2]; y++) {
					const double ww = wx * w3j[r2[3] + y];
					const double ar = a2r[i2 + y], ai = s2 * a2i[i2 + y];
					const double br = b2r[i2 + y], bi = s2 * b2i[i2 + y];
					Pr += ww * (a1r[i1 + y] * ar - a1i[i1 + y] * ai);
					Pi += ww * (a1r[i1 + y] * ai + a1i[i1 + y] * ar);
					Qr += ww * (b1r[i1 + y] * br - b1i[i1 + y] * bi);
					Qi += ww * (b1r[i1 + y] * bi + b1i[i1 + y] * br);
				}
			}
			t += K_re[mm] * Pr - K_im[mm] * Pi + G_re[mm] * Qr - G_im[mm] * Qi;
		}
	}
	red[threadIdx.x] = t;
	__syncthreads();
	for (int s = NORM_THREADS / 2; s > 0; s >>= 1) {
		if (threadIdx.x < s) red[threadIdx.x] += red[threadIdx.x + s];
		__syncthreads();
	}
	if (threadIdx.x == 0) {
		const double inner = 0.5 * red[0];
		normfact[atom] = inner > 0.0 ? 1.0 / sqrt(inner) : 0.0;
	}
}

//Builds the density matrices once per descriptor set and radial split; every later lambda reuses them
bool ensure_density(descriptor_cache& c, const salted_gpu_problem& q)
{
	const bool share = q.v2_is_conj_of_v1 && q.nrad2 == q.nrad1;
	if (c.d_dm1 && c.dm_nrad1 == q.nrad1 && c.dm_nrad2 == q.nrad2) return true;
	c.clear_dm();
	const int lmax1 = q.v1_noff - 1, lmax2 = share ? lmax1 : q.v2_noff - 1;
	if (lmax1 < 0 || lmax2 < 0) return false;
	const size_t dsz1 = dm_offset(lmax1 + 1), dsz2 = dm_offset(lmax2 + 1);
	GPU_TRY(gpuMalloc(&c.d_dm1, sizeof(double) * 4 * dsz1 * q.natoms));
	density_kernel<<<q.natoms, 128>>>(lmax1, q.nrad1, c.d_v1, c.d_v1off, q.v1_nchannels, c.d_dm1);
	GPU_TRY(gpuGetLastError());
	if (share)
		c.d_dm2 = c.d_dm1;
	else {
		GPU_TRY(gpuMalloc(&c.d_dm2, sizeof(double) * 4 * dsz2 * q.natoms));
		density_kernel<<<q.natoms, 128>>>(lmax2, q.nrad2, c.d_v2, c.d_v2off, q.v2_nchannels, c.d_dm2);
		GPU_TRY(gpuGetLastError());
	}
	c.dm_nrad1 = q.nrad1; c.dm_nrad2 = q.nrad2; c.dm_lmax1 = lmax1; c.dm_lmax2 = lmax2;
	return true;
}

//One thread per (atom, output slot): only the nfps selected features are built, the norm comes
//from norm_kernel. Consecutive slots write consecutive p entries for every imu, so the stores coalesce.
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
	GPU_TRY(gpuMalloc(&d_p, p_bytes));
	GPU_TRY(gpuMemset(d_p, 0, p_bytes));

	if (!ensure_density(g_descriptor_cache, q)) return false;
	{
		const descriptor_cache& c = g_descriptor_cache;
		int lmax1 = 0, lmax2 = 0;
		for (int il = 0; il < q.llmax; il++) {
			lmax1 = std::max(lmax1, q.llvec0[il]);
			lmax2 = std::max(lmax2, q.llvec1[il]);
		}
		if (lmax1 > c.dm_lmax1 || lmax2 > c.dm_lmax2) return false;
		const size_t kg = sizeof(double) * (size_t)q.l21 * q.l21;
		double* d_kg = nullptr;
		GPU_TRY(gpuMalloc(&d_kg, 4 * kg));
		const size_t n = (size_t)q.l21 * q.l21;
		GPU_TRY(gpuMemcpy(d_kg, q.K_re, kg, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_kg + n, q.K_im, kg, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_kg + 2 * n, q.G_re, kg, gpuMemcpyHostToDevice));
		GPU_TRY(gpuMemcpy(d_kg + 3 * n, q.G_im, kg, gpuMemcpyHostToDevice));
		norm_kernel<<<q.natoms, NORM_THREADS>>>(q.llmax, q.l21, q.v2_is_conj_of_v1 ? -1.0 : 1.0,
			c.d_dm1, dm_offset(c.dm_lmax1 + 1), c.d_dm2, dm_offset(c.dm_lmax2 + 1),
			d_w3j, d_ll0, d_ll1, d_runs, d_kg, d_kg + n, d_kg + 2 * n, d_kg + 3 * n, d_nf);
		GPU_TRY(gpuGetLastError());
		GPU_TRY(gpuMemcpy(q.normfact, d_nf, sizeof(double) * (size_t)q.natoms, gpuMemcpyDeviceToHost));
		gpuFree(d_kg);
	}

	const dim3 thr(128), grid((q.nfps + 127) / 128, q.natoms);
	equicomb_kernel<<<grid, thr>>>(q.natoms, q.nrad2, q.llmax, q.l21,
		q.shells, q.nfps, q.v2_is_conj_of_v1,
		g_descriptor_cache.d_v1, g_descriptor_cache.d_v1off, q.v1_nchannels,
		g_descriptor_cache.d_v2, g_descriptor_cache.d_v2off, q.v2_nchannels,
		d_w3j, d_ll0, d_ll1, d_runs, d_cols, d_c2r_re, d_c2r_im, d_cnt, d_vfps, d_nf, d_p);
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
	gpuFree(d_cols); gpuFree(d_c2r_re); gpuFree(d_c2r_im); gpuFree(d_cnt);
	gpuFree(d_vfps); gpuFree(d_nf); gpuFree(d_p);
	return true;
}

NOSPHERA2_GPU_API_END
