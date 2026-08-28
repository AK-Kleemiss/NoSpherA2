#include "salted_gpu.h"
#include "gpu_backend.h"
#include <cstdio>
#include <vector>
#include <algorithm>

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 SALTED GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

//l21 reaches 17 with the v7 model (lam 8). The per-shell complex vector lives in shared
//memory rather than in registers: 34 doubles per thread would spill.
//l21 is 2*lam+1 and reaches 17 with the v7 model. It sizes the per-thread transform
//vectors, so it is a stack-frame cost - but halving it to 32 measured neutral, so it
//stays wide enough for models the v7 one does not cover. The entry point refuses more.
#define SALTED_MAX_L21 64

namespace {

//block(atom, channel, l) as SALTEDDescriptors lays it out, in interleaved doubles
__device__ __forceinline__ const double* desc_block(const double* v, const size_t* off,
	const int nchannels, const int atom, const int channel, const int l)
{
	return v + 2 * (off[l] + ((size_t)atom * nchannels + channel) * (2 * (size_t)l + 1));
}

//One block owns one (atom, n1) pair, which is what makes the Wigner-weighted v1 worth
//building: it is reused by all nrad2 shells below, exactly as the CPU loop reuses it.
//Recomputing it per (n1, n2, il) instead costs nrad2 times the traffic on an inner loop
//that only does 0.2 flops per byte, which is why the first version measured no gain.
__global__ void equicomb_kernel(const int natoms, const int nrad2, const int llmax,
	const int l21, const int featsize, const int nfps, const bool conj,
	const int total_terms,
	const double* __restrict__ v1, const size_t* __restrict__ v1_off, const int v1_nch,
	const double* __restrict__ v2, const size_t* __restrict__ v2_off, const int v2_nch,
	const double* __restrict__ w3j, const int* __restrict__ llvec0, const int* __restrict__ llvec1,
	const int* __restrict__ runs, const int* __restrict__ c2r_cols,
	const double* __restrict__ c2r_re, const double* __restrict__ c2r_im,
	const int* __restrict__ c2r_cnt, const int* __restrict__ sel,
	double* __restrict__ out, double* __restrict__ inner)
{
	extern __shared__ double sh[];
	double* wv1_re = sh;
	double* wv1_im = sh + total_terms;

	const int atom = blockIdx.y;
	const int n1 = blockIdx.x;
	if (atom >= natoms) return;

	//Build w3j * v1 once for this (atom, n1): one thread per (il, imu) run
	for (int idx = threadIdx.x; idx < llmax * l21; idx += blockDim.x) {
		const int il = idx / l21;
		const int* run = runs + 4 * idx;
		const int im1_begin = run[0], count = run[2], w_off = run[3];
		const double* v1p = desc_block(v1, v1_off, v1_nch, atom, n1, llvec0[il]);
		for (int k = 0; k < count; k++) {
			const double wk = w3j[w_off + k];
			wv1_re[w_off + k] = wk * v1p[2 * (im1_begin + k)];
			wv1_im[w_off + k] = wk * v1p[2 * (im1_begin + k) + 1];
		}
	}
	__syncthreads();

	//One thread per (n2, il); it walks imu itself, so nothing needs syncing below
	double pc_re[SALTED_MAX_L21], pc_im[SALTED_MAX_L21];
	double thread_inner = 0.0;
	for (int job = threadIdx.x; job < nrad2 * llmax; job += blockDim.x) {
		//n2 varies fastest so neighbouring threads read neighbouring v2 blocks: those
		//are contiguous in the descriptor, where consecutive il are not
		const int n2 = job % nrad2;
		const int il = job / nrad2;
		const double* v2p = desc_block(v2, v2_off, v2_nch, atom, n2, llvec1[il]);
		for (int imu = 0; imu < l21; imu++) {
			const int* run = runs + 4 * (il * l21 + imu);
			const int im2_begin = run[1], count = run[2], w_off = run[3];
			double acc_r = 0.0, acc_i = 0.0;
			for (int k = 0; k < count; k++) {
				const double ar = wv1_re[w_off + k];
				const double ai = wv1_im[w_off + k];
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
		const int ifeat = (n1 * nrad2 + n2) * llmax + il;
		const int slot = sel[ifeat];
		double local_inner = 0.0;
		for (int i = 0; i < l21; i++) {
			double preal = 0.0;
			const int nz = c2r_cnt[i];
			for (int k = 0; k < nz; k++) {
				const int j = c2r_cols[2 * i + k];
				preal += c2r_re[2 * i + k] * pc_re[j] - c2r_im[2 * i + k] * pc_im[j];
			}
			local_inner += preal * preal;
			if (slot >= 0)
				out[((size_t)atom * nfps + slot) * l21 + i] = preal;
		}
		thread_inner += local_inner;
	}
	//One atomic per block instead of one per job: with 642 counters taking about
	//70k updates each, the contention is on the counter, not the arithmetic
	double* red = sh + 2 * total_terms;
	red[threadIdx.x] = thread_inner;
	__syncthreads();
	for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
		if (threadIdx.x < stride) red[threadIdx.x] += red[threadIdx.x + stride];
		__syncthreads();
	}
	if (threadIdx.x == 0) atomicAdd(&inner[atom], red[0]);
}

//p is laid out atom-major with imu striding by nfps, which is what the caller reads
__global__ void normalise_kernel(const int natoms, const int l21, const int nfps,
	const double* __restrict__ out, const double* __restrict__ inner,
	double* __restrict__ p, int* __restrict__ empty)
{
	const int atom = blockIdx.y;
	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (atom >= natoms || idx >= nfps * l21) return;
	const int i = idx / l21;
	const int imu = idx % l21;
	const double in = inner[atom];
	//An empty environment gives an all-zero descriptor; zero is the meaningful answer
	const double nf = (in > 0.0) ? rsqrt(in) : 0.0;
	if (in <= 0.0 && idx == 0) atomicAdd(empty, 1);
	p[(size_t)atom * l21 * nfps + i + imu * nfps] = out[((size_t)atom * nfps + i) * l21 + imu] * nf;
}

} //namespace

bool salted_gpu_available()
{
	int n = 0;
	return gpuGetDeviceCount(&n) == gpuSuccess && n > 0;
}

bool salted_gpu_equicomb(const salted_gpu_problem& q, int* empty_environments)
{
	if (q.natoms <= 0 || q.featsize <= 0 || q.nfps <= 0 || q.l21 <= 0) return false;
	if (q.l21 > SALTED_MAX_L21) return false;
	if (!salted_gpu_available()) return false;

	const size_t out_bytes = sizeof(double) * (size_t)q.natoms * q.nfps * q.l21;
	const size_t inner_bytes = sizeof(double) * (size_t)q.natoms;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	if (out_bytes + inner_bytes + (1u << 27) > freeb) return false;

	double *d_v1 = nullptr, *d_v2 = nullptr, *d_w3j = nullptr, *d_c2r_re = nullptr, *d_c2r_im = nullptr;
	double *d_out = nullptr, *d_inner = nullptr, *d_p = nullptr;
	size_t *d_v1off = nullptr, *d_v2off = nullptr;
	int *d_ll0 = nullptr, *d_ll1 = nullptr, *d_runs = nullptr, *d_cols = nullptr, *d_cnt = nullptr, *d_sel = nullptr, *d_empty = nullptr;

	const size_t v1_vals = q.v1_len_doubles, v2_vals = q.v2_len_doubles;
	GPU_TRY(gpuMalloc(&d_v1, sizeof(double) * v1_vals));
	GPU_TRY(gpuMemcpy(d_v1, q.v1_values, sizeof(double) * v1_vals, gpuMemcpyHostToDevice));
	if (q.v2_is_conj_of_v1) { d_v2 = d_v1; }
	else {
		GPU_TRY(gpuMalloc(&d_v2, sizeof(double) * v2_vals));
		GPU_TRY(gpuMemcpy(d_v2, q.v2_values, sizeof(double) * v2_vals, gpuMemcpyHostToDevice));
	}
	GPU_TRY(gpuMalloc(&d_v1off, sizeof(size_t) * q.v1_noff));
	GPU_TRY(gpuMemcpy(d_v1off, q.v1_offsets, sizeof(size_t) * q.v1_noff, gpuMemcpyHostToDevice));
	if (q.v2_is_conj_of_v1) { d_v2off = d_v1off; }
	else {
		GPU_TRY(gpuMalloc(&d_v2off, sizeof(size_t) * q.v2_noff));
		GPU_TRY(gpuMemcpy(d_v2off, q.v2_offsets, sizeof(size_t) * q.v2_noff, gpuMemcpyHostToDevice));
	}
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
	GPU_TRY(gpuMalloc(&d_sel, sizeof(int) * (size_t)q.featsize));
	GPU_TRY(gpuMemcpy(d_sel, q.sel, sizeof(int) * (size_t)q.featsize, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMalloc(&d_out, out_bytes));
	GPU_TRY(gpuMalloc(&d_inner, inner_bytes));
	GPU_TRY(gpuMalloc(&d_p, sizeof(double) * (size_t)q.natoms * q.l21 * q.nfps));
	GPU_TRY(gpuMalloc(&d_empty, sizeof(int)));
	GPU_TRY(gpuMemset(d_inner, 0, inner_bytes));
	GPU_TRY(gpuMemset(d_empty, 0, sizeof(int)));

	//total_terms is the length w3j is consumed to; the runs table ends at it
	int total_terms = 0;
	for (int r = 0; r < q.llmax * q.l21; r++)
		total_terms = std::max(total_terms, q.runs[4 * r + 3] + q.runs[4 * r + 2]);
	const size_t shmem = sizeof(double) * (2 * (size_t)total_terms + 256);
	//Shared holds the whole Wigner-weighted v1 for one (atom, n1)
	if (shmem > 96u * 1024u) return false;
	const dim3 thr(256), grid(q.nrad1, q.natoms);
	equicomb_kernel<<<grid, thr, shmem>>>(q.natoms, q.nrad2, q.llmax, q.l21,
		q.featsize, q.nfps, q.v2_is_conj_of_v1, total_terms,
		d_v1, d_v1off, q.v1_nchannels, d_v2, d_v2off, q.v2_nchannels,
		d_w3j, d_ll0, d_ll1, d_runs, d_cols, d_c2r_re, d_c2r_im, d_cnt, d_sel,
		d_out, d_inner);
	GPU_TRY(gpuGetLastError());

	const int nrm = q.nfps * q.l21;
	const dim3 nthr(256), ngrid((nrm + 255) / 256, q.natoms);
	normalise_kernel<<<ngrid, nthr>>>(q.natoms, q.l21, q.nfps, d_out, d_inner, d_p, d_empty);
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuDeviceSynchronize());

	GPU_TRY(gpuMemcpy(q.p, d_p, sizeof(double) * (size_t)q.natoms * q.l21 * q.nfps, gpuMemcpyDeviceToHost));
	int host_empty = 0;
	GPU_TRY(gpuMemcpy(&host_empty, d_empty, sizeof(int), gpuMemcpyDeviceToHost));
	if (empty_environments) *empty_environments += host_empty;

	gpuFree(d_v1); if (!q.v2_is_conj_of_v1) gpuFree(d_v2);
	gpuFree(d_v1off); if (!q.v2_is_conj_of_v1) gpuFree(d_v2off);
	gpuFree(d_w3j); gpuFree(d_ll0); gpuFree(d_ll1); gpuFree(d_runs);
	gpuFree(d_cols); gpuFree(d_c2r_re); gpuFree(d_c2r_im); gpuFree(d_cnt);
	gpuFree(d_sel); gpuFree(d_out); gpuFree(d_inner); gpuFree(d_p); gpuFree(d_empty);
	return true;
}
