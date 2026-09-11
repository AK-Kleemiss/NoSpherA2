#include "itensor_gpu.h"
#include "gpu_backend.h"
#include "itensor_gemm.cuh"
#include "throughput.h"
#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <climits>

#define SF_TWO_PI 6.283185307179586476925286766559
#define SF_INV_TWO_PI 0.15915494309189533576888376337251

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 I tensor GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

namespace {

//One set of device buffers per scalar type. Only one is ever live, and which one is a
//run-time choice, so both instantiations exist and g_fp64 says which to talk to.
//
//The blocks are not visited one at a time. A block's GEMM is a hundred by two hundred by a
//thousand, far too small to occupy a device on its own, and a reflection has a few hundred
//of them: issued singly they cost more in launches than in arithmetic. So the blocks are
//sorted into groups of like n_active, every block in a group padded with zeros to the
//group's row count and to the largest point count, and one strided-batched GEMM covers a
//group. The zero padding contributes nothing to the sums, only to the flop count, which
//is why the groups are as narrow as the batch call allows.
template <typename T>
struct Dev {
	bool ready = false;
	int nmo = 0, packed = 0, n_grids = 0, n_blocks = 0, np_max = 0, na_max = 0, n_groups = 0;
	int fac_cap = 0;
	long long n_points = 0;
	double issued_flops = 0.0;
	//Uploaded once: AO values padded to na_pad x np_max per block, in group order
	T* ao = nullptr;
	//Per block (in group order): owning grid, first point, points, n_active, padded rows,
	//offsets into ao (the weighted copy sits at twice that) and into the GEMM results
	int *q_grid = nullptr, *q_base = nullptr, *q_np = nullptr, *q_napad = nullptr;
	long long *q_ao = nullptr, *q_c = nullptr;
	//Per stored pair, the GEMM result elements that feed it, in block order: a CSR whose
	//cursor is advanced as the blocks are consumed in chunks
	int *acc_ptr = nullptr, *acc_q = nullptr, *acc_pos = nullptr, *acc_cur = nullptr;
	double *d1 = nullptr, *d2 = nullptr, *d3 = nullptr, *w = nullptr;
	//Per reflection scratch
	T *phase_re = nullptr, *phase_im = nullptr;
	//The real and imaginary weighted copies of a block sit back to back as one column-major
	//np_max x 2na_pad matrix, so one GEMM of width 2na_pad produces both results.
	T* wri = nullptr;
	T* cri = nullptr;
	long long w_cap = 0, c_cap = 0;
	void* gemm_ws = nullptr;   //split-k partials, whichever GEMM is compiled in
	double* fac[2] = { nullptr, nullptr };
	double* host_fac[2] = { nullptr, nullptr };
	double* I_re[2] = { nullptr, nullptr };
	double* I_im[2] = { nullptr, nullptr };
	double* host_re[2] = { nullptr, nullptr };
	double* host_im[2] = { nullptr, nullptr };
	gpuEvent_t done[2] = {};
	gpuStream_t copy_stream = nullptr;
	//Host copies of the layout, so the chunking stays on the host
	std::vector<int> grp_first, grp_count, grp_na;
	std::vector<long long> h_q_ao, h_q_c;
};

template <typename T> Dev<T> g;
bool g_fp64 = false;
bool g_tensor = false;

//The single-precision path keeps the reduced-argument trick the transform uses: the phase
//and its reduction stay in double and only the transcendental drops. In double there is
//nothing to trade, so it takes the argument as it stands.
template <typename T>
__device__ inline void phase_sincos(const double frac, T* s, T* c);

template <>
__device__ inline void phase_sincos<float>(const double frac, float* s, float* c)
{
	sincospif(2.0f * (float)frac, s, c);
}

template <>
__device__ inline void phase_sincos<double>(const double frac, double* s, double* c)
{
	sincospi(2.0 * frac, s, c);
}

//The weight is folded in here so the GEMM operand is exactly what the CPU path multiplies.
template <typename T>
__global__ void phase_kernel(const long long n, const double kx, const double ky, const double kz,
	const double* __restrict__ d1, const double* __restrict__ d2, const double* __restrict__ d3,
	const double* __restrict__ w, T* __restrict__ pre, T* __restrict__ pim)
{
	const long long p = (long long)blockIdx.x * blockDim.x + threadIdx.x;
	if (p >= n) return;
	//kx..kz arrive already divided by 2pi, so t is in turns and sincospi wants 2*frac
	const double t = kx * d1[p] + ky * d2[p] + kz * d3[p];
	const double frac = t - rint(t);
	T s, c;
	phase_sincos<T>(frac, &s, &c);
	const double wp = w[p];
	pre[p] = (T)(wp * (double)c);
	pim[p] = (T)(wp * (double)s);
}

//Phase-weighted copies of the blocks q0.. of this chunk: blockIdx.z is the block, blockIdx.y
//the AO row. The padding rows hold zeros in ao and so write zeros here; the padding points
//are given a zero phase so nothing is read past a grid's end.
template <typename T>
__global__ void weight_kernel(const int q0, const int np_max, const long long w_base,
	const int* __restrict__ q_napad, const int* __restrict__ q_np, const int* __restrict__ q_base,
	const long long* __restrict__ q_ao, const T* __restrict__ ao,
	const T* __restrict__ pre, const T* __restrict__ pim, T* __restrict__ wri)
{
	const int q = q0 + blockIdx.z;
	const int row = blockIdx.y;
	const int napad = q_napad[q];
	if (row >= napad) return;
	const int p = blockIdx.x * blockDim.x + threadIdx.x;
	if (p >= np_max) return;
	const long long a_off = q_ao[q];
	const T a = ao[a_off + (long long)row * np_max + p];
	T re = T(0), im = T(0);
	if (p < q_np[q]) {
		const int pp = q_base[q] + p;
		re = pre[pp];
		im = pim[pp];
	}
	T* w = wri + 2 * a_off - w_base;
	w[(long long)row * np_max + p] = a * re;
	w[(long long)(napad + row) * np_max + p] = a * im;
}

//One thread per stored pair walks its entries for the blocks of this chunk, in block order,
//and adds their contribution. Fixed order and no atomics, so the result does not depend on
//how the device scheduled the blocks.
template <typename T>
__global__ void gather_kernel(const int packed, const int q_end, const long long c_base,
	const int* __restrict__ acc_ptr, int* __restrict__ acc_cur,
	const int* __restrict__ acc_q, const int* __restrict__ acc_pos,
	const long long* __restrict__ q_c, const int* __restrict__ q_napad, const int* __restrict__ q_grid,
	const T* __restrict__ cri, const double* __restrict__ fre, const double* __restrict__ fim,
	double* __restrict__ I_re, double* __restrict__ I_im)
{
	const int t = blockIdx.x * blockDim.x + threadIdx.x;
	if (t >= packed) return;
	int e = acc_cur[t];
	const int end = acc_ptr[t + 1];
	double sre = 0.0, sim = 0.0;
	for (; e < end; e++) {
		const int q = acc_q[e];
		if (q >= q_end) break;
		const T* c = cri + q_c[q] - c_base;
		const long long napad = q_napad[q];
		const int pos = acc_pos[e];
		//The GEMM wrote column-major na_pad x 2na_pad, the imaginary half after the real
		const double re = (double)c[pos];
		const double im = (double)c[napad * napad + pos];
		const int gi = q_grid[q];
		const double a = fre[gi], b = fim[gi];
		sre += re * a - im * b;
		sim += re * b + im * a;
	}
	acc_cur[t] = e;
	I_re[t] += sre;
	I_im[t] += sim;
}

__global__ void zero_kernel(const int n, double* a, double* b)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n) { a[i] = 0.0; b[i] = 0.0; }
}

template <typename T>
void free_impl()
{
	Dev<T>& d = g<T>;
	gpuFree(d.ao);
	gpuFree(d.q_grid); gpuFree(d.q_base); gpuFree(d.q_np); gpuFree(d.q_napad);
	gpuFree(d.q_ao); gpuFree(d.q_c);
	gpuFree(d.acc_ptr); gpuFree(d.acc_q); gpuFree(d.acc_pos); gpuFree(d.acc_cur);
	gpuFree(d.d1); gpuFree(d.d2); gpuFree(d.d3); gpuFree(d.w);
	gpuFree(d.phase_re); gpuFree(d.phase_im);
	gpuFree(d.wri); gpuFree(d.cri); gpuFree(d.gemm_ws);
	for (int i = 0; i < 2; i++) {
		gpuFree(d.I_re[i]); gpuFree(d.I_im[i]); gpuFree(d.fac[i]);
		if (d.host_re[i]) gpuFreeHost(d.host_re[i]);
		if (d.host_im[i]) gpuFreeHost(d.host_im[i]);
		if (d.host_fac[i]) gpuFreeHost(d.host_fac[i]);
		if (d.done[i]) gpuEventDestroy(d.done[i]);
	}
	if (d.copy_stream) gpuStreamDestroy(d.copy_stream);
	d = Dev<T>{};
}

template <typename T>
bool upload_vec(T** dst, const std::vector<T>& v)
{
	GPU_TRY(gpuMalloc(dst, sizeof(T) * std::max<size_t>(v.size(), 1)));
	if (!v.empty())
		GPU_TRY(gpuMemcpy(*dst, v.data(), sizeof(T) * v.size(), gpuMemcpyHostToDevice));
	return true;
}

template <typename T>
bool init_impl(const itensor_gpu_layout& L)
{
	Dev<T>& d = g<T>;
	const int nb = L.n_blocks;
	//Rows padded to a multiple of eight: what the Tensor Core path wants of its leading
	//dimensions, and narrow enough that the padding stays a few per cent of the work
	auto pad8 = [](const int n) { return (n + 7) & ~7; };
	std::vector<int> order(nb);
	std::iota(order.begin(), order.end(), 0);
	std::stable_sort(order.begin(), order.end(), [&](const int a, const int b) {
		return pad8(L.blk_n_active[a]) < pad8(L.blk_n_active[b]); });
	d.np_max = 0; d.na_max = 0;
	for (int b = 0; b < nb; b++) {
		d.np_max = std::max(d.np_max, L.blk_point_count[b]);
		d.na_max = std::max(d.na_max, pad8(L.blk_n_active[b]));
	}
	std::vector<int> h_grid(nb), h_base(nb), h_np(nb), h_napad(nb);
	d.h_q_ao.assign(nb + 1, 0); d.h_q_c.assign(nb + 1, 0);
	d.grp_first.clear(); d.grp_count.clear(); d.grp_na.clear();
	long long max_blk_ao = 0, max_blk_c = 0;
	size_t max_ws = 0;
	d.issued_flops = 0.0;
	for (int q = 0; q < nb; q++) {
		const int b = order[q];
		const int napad = pad8(L.blk_n_active[b]);
		h_grid[q] = L.blk_grid[b];
		h_base[q] = L.grid_point_off[L.blk_grid[b]] + L.blk_point_start[b];
		h_np[q] = L.blk_point_count[b];
		h_napad[q] = napad;
		d.h_q_ao[q + 1] = d.h_q_ao[q] + (long long)napad * d.np_max;
		d.h_q_c[q + 1] = d.h_q_c[q] + 2LL * napad * napad;
		max_blk_ao = std::max(max_blk_ao, (long long)napad * d.np_max);
		max_blk_c = std::max(max_blk_c, 2LL * napad * napad);
		if (d.grp_na.empty() || d.grp_na.back() != napad) {
			d.grp_first.push_back(q); d.grp_count.push_back(0); d.grp_na.push_back(napad);
			max_ws = std::max(max_ws, itensor_gemm::workspace_bytes<T>(napad, 2 * napad, d.np_max));
		}
		d.grp_count.back()++;
		d.issued_flops += throughput::flops_gemm(napad, 2.0 * napad, d.np_max);
	}
	d.n_groups = (int)d.grp_na.size();
	//Which result elements feed each stored pair, in block order
	std::vector<int> cnt(L.packed + 1, 0);
	long long n_entries = 0;
	for (int q = 0; q < nb; q++) {
		const int b = order[q];
		const int na = L.blk_n_active[b];
		const int* aos = L.aos_all + L.blk_aos_off[b];
		for (int i = 0; i < na; i++)
			for (int j = i; j < na; j++) {
				const int t = L.compact[(long long)aos[i] * L.nmo + aos[j]];
				if (t >= 0) { cnt[t]++; n_entries++; }
			}
	}
	if (n_entries > INT_MAX) return false;
	std::vector<int> ptr(L.packed + 1, 0);
	for (int t = 0; t < L.packed; t++) ptr[t + 1] = ptr[t] + cnt[t];
	std::vector<int> acc_q((size_t)n_entries), acc_pos((size_t)n_entries), fill(ptr.begin(), ptr.end() - 1);
	for (int q = 0; q < nb; q++) {
		const int b = order[q];
		const int na = L.blk_n_active[b];
		const int napad = h_napad[q];
		const int* aos = L.aos_all + L.blk_aos_off[b];
		for (int i = 0; i < na; i++)
			for (int j = i; j < na; j++) {
				const int t = L.compact[(long long)aos[i] * L.nmo + aos[j]];
				if (t < 0) continue;
				acc_q[fill[t]] = q;
				acc_pos[fill[t]] = j * napad + i;
				fill[t]++;
			}
	}
	//The weighted copies and the results are consumed a chunk of blocks at a time, so
	//neither has to hold every block at once; a chunk is as many whole blocks as fit.
	const long long ao_total = d.h_q_ao[nb], c_total = d.h_q_c[nb];
	d.w_cap = std::min(2 * ao_total, std::max(2 * max_blk_ao, (long long)((256u << 20) / sizeof(T))));
	d.c_cap = std::min(c_total, std::max(max_blk_c, (long long)((64u << 20) / sizeof(T))));
	if (throughput::enabled()) {
		double dense = 0.0;
		for (int b = 0; b < nb; b++)
			dense += throughput::flops_gemm(L.blk_n_active[b], 2.0 * L.blk_n_active[b], L.blk_point_count[b]);
		//stderr because cout is redirected to the log and moved again mid-run
		std::fprintf(stderr, "I tensor GPU: %d blocks in %d groups of like n_active, padded to %d rows x %d points"
			" at most, %.1f%% of the issued GEMM work is padding\n",
			nb, d.n_groups, d.na_max, d.np_max, d.issued_flops > 0.0 ? 100.0 * (1.0 - dense / d.issued_flops) : 0.0);
	}
	const size_t need =
		sizeof(T) * (size_t)ao_total +
		sizeof(T) * (size_t)(d.w_cap + d.c_cap) +
		max_ws +
		sizeof(int) * 4 * (size_t)nb + sizeof(long long) * 2 * (size_t)nb +
		sizeof(int) * (2 * (size_t)n_entries + 2 * (size_t)L.packed + 2) +
		sizeof(double) * 4 * (size_t)L.n_points +
		sizeof(T) * 2 * (size_t)L.n_points +
		sizeof(double) * 4 * (size_t)L.packed;
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	if (need + (1u << 26) > freeb) return false;

	GPU_TRY(gpuMalloc(&d.ao, sizeof(T) * (size_t)ao_total));
	GPU_TRY(gpuMemset(d.ao, 0, sizeof(T) * (size_t)ao_total));
	{
		//Staged a block at a time: the padded rows and points stay zero from the memset
		std::vector<T> stage((size_t)max_blk_ao);
		for (int q = 0; q < nb; q++) {
			const int b = order[q];
			const int na = L.blk_n_active[b], np = L.blk_point_count[b];
			const double* src = L.ao_all + L.blk_ao_off[b];
			for (int i = 0; i < na; i++)
				for (int p = 0; p < np; p++) stage[(size_t)i * d.np_max + p] = (T)src[(size_t)i * np + p];
			for (int i = 0; i < na; i++)
				for (int p = np; p < d.np_max; p++) stage[(size_t)i * d.np_max + p] = T(0);
			GPU_TRY(gpuMemcpy(d.ao + d.h_q_ao[q], stage.data(), sizeof(T) * (size_t)na * d.np_max,
				gpuMemcpyHostToDevice));
		}
	}
	if (!upload_vec(&d.q_grid, h_grid) || !upload_vec(&d.q_base, h_base) || !upload_vec(&d.q_np, h_np)
		|| !upload_vec(&d.q_napad, h_napad)) return false;
	if (!upload_vec(&d.q_ao, std::vector<long long>(d.h_q_ao.begin(), d.h_q_ao.end() - 1))) return false;
	if (!upload_vec(&d.q_c, std::vector<long long>(d.h_q_c.begin(), d.h_q_c.end() - 1))) return false;
	if (!upload_vec(&d.acc_ptr, ptr) || !upload_vec(&d.acc_q, acc_q) || !upload_vec(&d.acc_pos, acc_pos)) return false;
	GPU_TRY(gpuMalloc(&d.acc_cur, sizeof(int) * (size_t)(L.packed + 1)));
	GPU_TRY(gpuMalloc(&d.d1, sizeof(double) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&d.d2, sizeof(double) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&d.d3, sizeof(double) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&d.w, sizeof(double) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&d.phase_re, sizeof(T) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&d.phase_im, sizeof(T) * (size_t)L.n_points));
	GPU_TRY(gpuMalloc(&d.wri, sizeof(T) * (size_t)d.w_cap));
	GPU_TRY(gpuMalloc(&d.cri, sizeof(T) * (size_t)d.c_cap));
	GPU_TRY(gpuMalloc(&d.gemm_ws, max_ws ? max_ws : 1));
	for (int i = 0; i < 2; i++) {
		GPU_TRY(gpuMalloc(&d.I_re[i], sizeof(double) * (size_t)L.packed));
		GPU_TRY(gpuMalloc(&d.I_im[i], sizeof(double) * (size_t)L.packed));
		GPU_TRY(gpuHostAlloc((void**)&d.host_re[i], sizeof(double) * (size_t)L.packed));
		GPU_TRY(gpuHostAlloc((void**)&d.host_im[i], sizeof(double) * (size_t)L.packed));
		GPU_TRY(gpuEventCreate(&d.done[i]));
	}
	GPU_TRY(gpuStreamCreateNonBlocking(&d.copy_stream));
	GPU_TRY(gpuMemcpy(d.d1, L.d1, sizeof(double) * (size_t)L.n_points, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(d.d2, L.d2, sizeof(double) * (size_t)L.n_points, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(d.d3, L.d3, sizeof(double) * (size_t)L.n_points, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(d.w, L.weights, sizeof(double) * (size_t)L.n_points, gpuMemcpyHostToDevice));

	d.nmo = L.nmo; d.packed = L.packed; d.n_grids = L.n_grids; d.n_blocks = nb;
	d.n_points = L.n_points;
	d.ready = true;
	return true;
}

//The per-grid factors travel with the reflection; sized on first use since the symmetry
//count is the caller's
template <typename T>
bool ensure_factors(Dev<T>& d, const int n)
{
	if (n <= d.fac_cap) return true;
	GPU_TRY(gpuDeviceSynchronize());
	for (int i = 0; i < 2; i++) {
		gpuFree(d.fac[i]);
		if (d.host_fac[i]) gpuFreeHost(d.host_fac[i]);
		GPU_TRY(gpuMalloc(&d.fac[i], sizeof(double) * 2 * (size_t)n));
		GPU_TRY(gpuHostAlloc((void**)&d.host_fac[i], sizeof(double) * 2 * (size_t)n));
	}
	d.fac_cap = n;
	return true;
}

template <typename T>
bool submit_impl(const int slot, const int num_syms,
	const double* kx, const double* ky, const double* kz,
	const std::complex<double>* factors)
{
	Dev<T>& d = g<T>;
	if (!d.ready || slot < 0 || slot > 1) return false;
	const int nf = num_syms * d.n_grids;
	if (!ensure_factors(d, nf)) return false;
	for (int i = 0; i < nf; i++) {
		d.host_fac[slot][i] = factors[i].real();
		d.host_fac[slot][nf + i] = factors[i].imag();
	}
	GPU_TRY(gpuMemcpyAsync(d.fac[slot], d.host_fac[slot], sizeof(double) * 2 * (size_t)nf, gpuMemcpyHostToDevice, 0));
	zero_kernel<<<(d.packed + 255) / 256, 256>>>(d.packed, d.I_re[slot], d.I_im[slot]);
	for (int s = 0; s < num_syms; s++) {
		//The CPU path takes sin/cos of k.d directly; scaling k to turns here is what lets
		//the reduction be a rint and the transcendental be sincospi
		phase_kernel<T><<<(unsigned int)((d.n_points + 255) / 256), 256>>>(
			d.n_points, kx[s] * SF_INV_TWO_PI, ky[s] * SF_INV_TWO_PI, kz[s] * SF_INV_TWO_PI,
			d.d1, d.d2, d.d3, d.w, d.phase_re, d.phase_im);
		GPU_TRY(gpuMemcpyAsync(d.acc_cur, d.acc_ptr, sizeof(int) * (size_t)d.packed, gpuMemcpyDeviceToDevice, 0));
		int q0 = 0;
		while (q0 < d.n_blocks) {
			int q1 = q0;
			while (q1 < d.n_blocks && q1 - q0 < 65535 && 2 * (d.h_q_ao[q1 + 1] - d.h_q_ao[q0]) <= d.w_cap
				&& d.h_q_c[q1 + 1] - d.h_q_c[q0] <= d.c_cap) q1++;
			const long long w_base = 2 * d.h_q_ao[q0], c_base = d.h_q_c[q0];
			weight_kernel<T><<<dim3((d.np_max + 255) / 256, d.na_max, q1 - q0), 256>>>(
				q0, d.np_max, w_base, d.q_napad, d.q_np, d.q_base, d.q_ao, d.ao,
				d.phase_re, d.phase_im, d.wri);
			//C = A * W^T per block, A row-major na_pad x np_max, i.e. column-major np_max x
			//na_pad, W column-major np_max x 2na_pad; a group is one call
			for (int gr = 0; gr < d.n_groups; gr++) {
				const int qs = std::max(d.grp_first[gr], q0);
				const int qe = std::min(d.grp_first[gr] + d.grp_count[gr], q1);
				if (qe <= qs) continue;
				const int napad = d.grp_na[gr];
				if (!itensor_gemm::run_batched<T>(napad, 2 * napad, d.np_max,
					d.ao + d.h_q_ao[qs], d.np_max, (long long)napad * d.np_max,
					d.wri + 2 * d.h_q_ao[qs] - w_base, d.np_max, 2LL * napad * d.np_max,
					d.cri + d.h_q_c[qs] - c_base, napad, 2LL * napad * napad, qe - qs, d.gemm_ws))
					return false;
			}
			gather_kernel<T><<<(d.packed + 255) / 256, 256>>>(d.packed, q1, c_base,
				d.acc_ptr, d.acc_cur, d.acc_q, d.acc_pos, d.q_c, d.q_napad, d.q_grid, d.cri,
				d.fac[slot] + (size_t)s * d.n_grids, d.fac[slot] + nf + (size_t)s * d.n_grids,
				d.I_re[slot], d.I_im[slot]);
			q0 = q1;
		}
	}
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuEventRecord(d.done[slot], 0));
	return true;
}

template <typename T>
bool collect_impl(const int slot, std::complex<double>* I_r)
{
	Dev<T>& d = g<T>;
	if (!d.ready || slot < 0 || slot > 1) return false;
	GPU_TRY(gpuStreamWaitEvent(d.copy_stream, d.done[slot], 0));
	GPU_TRY(gpuMemcpyAsync(d.host_re[slot], d.I_re[slot], sizeof(double) * (size_t)d.packed,
		gpuMemcpyDeviceToHost, d.copy_stream));
	GPU_TRY(gpuMemcpyAsync(d.host_im[slot], d.I_im[slot], sizeof(double) * (size_t)d.packed,
		gpuMemcpyDeviceToHost, d.copy_stream));
	GPU_TRY(gpuStreamSynchronize(d.copy_stream));
	for (int i = 0; i < d.packed; i++)
		I_r[i] += std::complex<double>(d.host_re[slot][i], d.host_im[slot][i]);
	return true;
}

} //namespace

//Shared with the transform so the "no code for this card" case is diagnosed in one place.
bool itensor_gpu_available() { return sf_gpu_available(); }

const char* itensor_gpu_gemm_name()
{
	if (g_tensor) return "cuBLAS Tensor Core";
	if (cublas_dynamic_available()) return "cuBLAS";
	//CUTLASS covers single precision only, so the double path names a different kernel.
	return g_fp64 ? "built-in" : NOSPHERA2_ITENSOR_GEMM_NAME;
}

double itensor_gpu_issued_flops()
{
	return g_fp64 ? g<double>.issued_flops : g<float>.issued_flops;
}

bool itensor_gpu_init(const itensor_gpu_layout& L, const sf_precision prec, const bool tensor)
{
	itensor_gpu_free();
	if (!itensor_gpu_available()) return false;
	//Auto is not offered here. It would resolve per card, and the I tensor's precision is
	//visible in the reference output, so the same input would produce different logs on
	//different machines. Single precision unless the caller asks for double.
	g_fp64 = (prec == sf_precision::FP64);
	bool tensor_hardware = false;
#ifndef NOSPHERA2_USE_HIP
	int dev = 0;
	gpuDeviceProp_t prop{};
	if (gpuGetDevice(&dev) == gpuSuccess && gpuGetDeviceProperties(&prop, dev) == gpuSuccess)
		tensor_hardware = prop.major >= 7;
#endif
	g_tensor = !g_fp64 && tensor && tensor_hardware && cublas_dynamic_fast_16f_available();
	itensor_gemm::set_tensor_mode(g_tensor);
	const bool ok = g_fp64 ? init_impl<double>(L) : init_impl<float>(L);
	if (!ok) itensor_gpu_free();
	return ok;
}

bool itensor_gpu_submit(const int slot, const int num_syms,
	const double* kx, const double* ky, const double* kz,
	const std::complex<double>* factors)
{
	return g_fp64 ? submit_impl<double>(slot, num_syms, kx, ky, kz, factors)
	              : submit_impl<float>(slot, num_syms, kx, ky, kz, factors);
}

bool itensor_gpu_collect(const int slot, std::complex<double>* I_r)
{
	return g_fp64 ? collect_impl<double>(slot, I_r) : collect_impl<float>(slot, I_r);
}

void itensor_gpu_free()
{
	g_tensor = false;
	itensor_gemm::set_tensor_mode(false);
	free_impl<float>();
	free_impl<double>();
}
