#include "itensor_gpu.h"
#include "gpu_backend.h"
#include "itensor_gemm.cuh"
#include "throughput.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <climits>

#define SF_INV_TWO_PI 0.15915494309189533576888376337251

#define GPU_TRY(call) do { const gpuError_t e_ = (call); if (e_ != gpuSuccess) { \
	std::fprintf(stderr, "NoSpherA2 I tensor GPU: %s at %s:%d\n", gpuGetErrorString(e_), __FILE__, __LINE__); \
	return false; } } while (0)

namespace {

//One set of device buffers per scalar type. Only one is ever live, and which one is a
//run-time choice, so both instantiations exist and g_fp64 says which to talk to.
//
//The contraction is arranged with the reflections as the wide dimension. For a block, the
//products of two of its AOs at every grid point form a table with a row per stored pair
//and a column per point that does not depend on the reflection, and a reflection is one
//column of weighted phases over the same points. So a batch of reflections is one GEMM per
//block, the block's table against the batch's phase columns: a shape wide enough to occupy
//the device, with the point count as the depth, and no padding beyond a few points. Only
//the pairs the caller stores are tabulated, each once, the symmetric half being enough.
//
//The tables can outgrow the device, so the rows of a block are cut into pieces, the pieces
//gathered into chunks that fit, and a chunk's tables rebuilt from the AO values for every
//batch. When everything fits in one chunk the tables are built once.
template <typename T>
struct Dev {
	bool ready = false, table_ready = false;
	int nmo = 0, packed = 0, n_grids = 0, n_blocks = 0, n_pieces = 0, n_chunks = 0;
	int np8_total = 0, np8_max = 0, ncol_cap = 0, batch_max = 0;
	int n_refl[2] = { 0, 0 };
	long long n_entries = 0;
	double issued_flops = 0.0;
	//Points in block order, each block padded to a multiple of eight: the coordinates and
	//weights, and the AO values row-major n_active x np8 per block. The padding points
	//carry zero weight and zero AOs, so nothing reads past a block's end.
	double *d1 = nullptr, *d2 = nullptr, *d3 = nullptr, *w = nullptr;
	T* ao = nullptr;
	//Per block: first padded point, padded points, first row of ao
	int *q_pp = nullptr, *q_np8 = nullptr;
	long long* q_ao = nullptr;
	//Per table row, in block order: the two AO rows, the owning grid and the piece
	int *ent_i = nullptr, *ent_j = nullptr, *ent_grid = nullptr, *ent_piece = nullptr;
	//Per piece: block, first row, first table element within its chunk
	int *pc_blk = nullptr, *pc_e0 = nullptr;
	long long* pc_tab = nullptr;
	//Per stored pair the rows that feed it, in row order, and where each chunk's part of
	//that list begins: acc_ptr[c * packed + t] .. acc_ptr[(c + 1) * packed + t]
	int *acc_ptr = nullptr, *acc_e = nullptr;
	//Per batch: phases column-major np8_total x ncol, the tables of one chunk, the GEMM
	//results column-major ncol x rows
	T *phase = nullptr, *tab = nullptr, *cres = nullptr;
	long long tab_cap = 0, rows_cap = 0;
	void* gemm_ws = nullptr;
	//Two slots so a batch is read back while the next one runs
	double* kvec[2] = { nullptr, nullptr };
	double* fac[2] = { nullptr, nullptr };
	double* host_kf[2] = { nullptr, nullptr };
	double* I_re[2] = { nullptr, nullptr };
	double* I_im[2] = { nullptr, nullptr };
	double* host_re[2] = { nullptr, nullptr };
	double* host_im[2] = { nullptr, nullptr };
	gpuEvent_t done[2] = {};
	gpuStream_t copy_stream = nullptr;
	//The host side of the plan
	struct Piece { int blk, e0, rows; long long tab; };
	struct Chunk { int p0, p1, e0, e1; long long tab0; };
	std::vector<Piece> pieces;
	std::vector<Chunk> chunks;
	std::vector<long long> h_pp;
	std::vector<int> h_np8;
};

template <typename T> Dev<T> g;
bool g_fp64 = false;
bool g_tensor = false;

//The resident tensor and the per-iteration operands. Row-major as on the host, so the
//upload is one copy; the column walk is cut into row chunks whose partial sums are added
//in a fixed order, so a run repeats itself.
constexpr int hold_chunks = 64;
struct Held {
	void* I = nullptr;
	bool fp64 = false;
	int nr = 0, packed = 0;
	double* w = nullptr;
	double* pre = nullptr;
	double* F0 = nullptr;
	double* F = nullptr;
	double* part = nullptr;
	double* out = nullptr;
};
Held g_held;
struct HeldEri {
	double* V = nullptr;
	double* dp = nullptr;
	double* D = nullptr;
	double* J1 = nullptr;
	double* J2 = nullptr;
	double* part = nullptr;
	double* Ka = nullptr;
	double* Kb = nullptr;
	double* K = nullptr;
	int n = 0, npair = 0;
};
HeldEri g_eri;

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

//One column pair per reflection and symmetry operation, the weight folded in so the GEMM
//operand is exactly what the CPU path multiplies. A thread holds one point and walks a
//run of columns, so the coordinates are read once a run and not once a column, and the
//runs are short enough that the double-precision latency has other threads to hide
//behind. Columns past the batch are zeroed so the padding to a multiple of eight holds
//nothing.
constexpr int phase_run = 32;

template <typename T>
__global__ void phase_kernel(const int np8_total, const int ncomb, const int ncol8, const double* __restrict__ kvec,
	const double* __restrict__ d1, const double* __restrict__ d2, const double* __restrict__ d3,
	const double* __restrict__ w, T* __restrict__ phase)
{
	const int p = blockIdx.x * blockDim.x + threadIdx.x;
	if (p >= np8_total) return;
	const double x = d1[p], y = d2[p], z = d3[p], wp = w[p];
	const int c0 = blockIdx.y * phase_run, c1 = min(c0 + phase_run, ncol8 / 2);
	for (int c = c0; c < c1; c++) {
		T re = T(0), im = T(0);
		if (c < ncomb) {
			//kx..kz arrive already divided by 2pi, so t is in turns and sincospi wants 2*frac
			const double t = kvec[3 * c] * x + kvec[3 * c + 1] * y + kvec[3 * c + 2] * z;
			const double frac = t - rint(t);
			T s, co;
			phase_sincos<T>(frac, &s, &co);
			re = (T)(wp * (double)co);
			im = (T)(wp * (double)s);
		}
		phase[(long long)(2 * c) * np8_total + p] = re;
		phase[(long long)(2 * c + 1) * np8_total + p] = im;
	}
}

//The tables of the rows e0.. of a chunk: blockIdx.y strides over the rows, the threads
//over a block's padded points
template <typename T>
__global__ void table_kernel(const int e0, const int n_rows, const long long tab0,
	const int* __restrict__ ent_i, const int* __restrict__ ent_j, const int* __restrict__ ent_piece,
	const int* __restrict__ pc_blk, const int* __restrict__ pc_e0, const long long* __restrict__ pc_tab,
	const int* __restrict__ q_np8, const long long* __restrict__ q_ao,
	const T* __restrict__ ao, T* __restrict__ tab)
{
	const int p = blockIdx.x * blockDim.x + threadIdx.x;
	for (int r = blockIdx.y; r < n_rows; r += gridDim.y) {
		const int e = e0 + r;
		const int pc = ent_piece[e];
		const int b = pc_blk[pc];
		const int np8 = q_np8[b];
		if (p >= np8) continue;
		const T* a = ao + q_ao[b];
		tab[pc_tab[pc] - tab0 + (long long)(e - pc_e0[pc]) * np8 + p] =
			a[(long long)ent_i[e] * np8 + p] * a[(long long)ent_j[e] * np8 + p];
	}
}

//One thread per stored pair and reflection adds up the pair's rows of this chunk, over the
//symmetry operations, with the per-grid factors. Fixed order and no atomics, so the result
//does not depend on how the device scheduled the blocks. Neighbouring threads take
//neighbouring reflections, whose results sit side by side in a row of cres.
template <typename T>
__global__ void gather_kernel(const int packed, const int n_refl, const int ns, const int n_grids,
	const int ncol8, const int e_base, const int* __restrict__ acc_ptr, const int* __restrict__ acc_e,
	const int* __restrict__ ent_grid, const T* __restrict__ cres,
	const double* __restrict__ fre, const double* __restrict__ fim,
	double* __restrict__ I_re, double* __restrict__ I_im)
{
	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const int r = idx % n_refl, t = idx / n_refl;
	if (t >= packed) return;
	const int start = acc_ptr[t], end = acc_ptr[packed + t];
	double sre = 0.0, sim = 0.0;
	for (int k = start; k < end; k++) {
		const int e = acc_e[k];
		const T* c = cres + (long long)(e - e_base) * ncol8;
		const int gi = ent_grid[e];
		for (int s = 0; s < ns; s++) {
			const int j = r * ns + s;
			const double re = (double)c[2 * j], im = (double)c[2 * j + 1];
			const double a = fre[j * n_grids + gi], b = fim[j * n_grids + gi];
			sre += re * a - im * b;
			sim += re * b + im * a;
		}
	}
	I_re[(long long)r * packed + t] += sre;
	I_im[(long long)r * packed + t] += sim;
}

__global__ void zero_kernel(const long long n, double* a, double* b)
{
	const long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n) { a[i] = 0.0; b[i] = 0.0; }
}

template <typename T>
void free_impl()
{
	Dev<T>& d = g<T>;
	gpuFree(d.d1); gpuFree(d.d2); gpuFree(d.d3); gpuFree(d.w); gpuFree(d.ao);
	gpuFree(d.q_pp); gpuFree(d.q_np8); gpuFree(d.q_ao);
	gpuFree(d.ent_i); gpuFree(d.ent_j); gpuFree(d.ent_grid); gpuFree(d.ent_piece);
	gpuFree(d.pc_blk); gpuFree(d.pc_e0); gpuFree(d.pc_tab);
	gpuFree(d.acc_ptr); gpuFree(d.acc_e);
	gpuFree(d.phase); gpuFree(d.tab); gpuFree(d.cres); gpuFree(d.gemm_ws);
	for (int i = 0; i < 2; i++) {
		gpuFree(d.kvec[i]); gpuFree(d.fac[i]); gpuFree(d.I_re[i]); gpuFree(d.I_im[i]);
		if (d.host_kf[i]) gpuFreeHost(d.host_kf[i]);
		if (d.host_re[i]) gpuFreeHost(d.host_re[i]);
		if (d.host_im[i]) gpuFreeHost(d.host_im[i]);
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
	//Points padded to a multiple of eight per block: the leading dimension of both GEMM
	//operands, which the Tensor Core kernels want that way or cuBLAS quietly runs the plain
	//ones
	std::vector<int> h_pp(nb), h_np8(nb);
	std::vector<long long> h_ao(nb);
	d.np8_total = 0; d.np8_max = 0;
	long long ao_total = 0;
	for (int b = 0; b < nb; b++) {
		const int np8 = (L.blk_point_count[b] + 7) & ~7;
		h_pp[b] = d.np8_total; h_np8[b] = np8; h_ao[b] = ao_total;
		d.np8_total += np8;
		d.np8_max = std::max(d.np8_max, np8);
		ao_total += (long long)L.blk_n_active[b] * np8;
	}
	//The stored pairs of every block, one table row each
	std::vector<int> ent_i, ent_j, ent_t, ent_grid, blk_e0(nb + 1, 0);
	for (int b = 0; b < nb; b++) {
		const int na = L.blk_n_active[b];
		const int* aos = L.aos_all + L.blk_aos_off[b];
		for (int i = 0; i < na; i++)
			for (int j = i; j < na; j++) {
				const int t = L.compact[(long long)aos[i] * L.nmo + aos[j]];
				if (t < 0) continue;
				ent_i.push_back(i); ent_j.push_back(j); ent_t.push_back(t); ent_grid.push_back(L.blk_grid[b]);
			}
		blk_e0[b + 1] = (int)ent_i.size();
	}
	d.n_entries = (long long)ent_i.size();
	if (d.n_entries > INT_MAX) return false;
	d.issued_flops = 0.0;
	for (int b = 0; b < nb; b++)
		d.issued_flops += throughput::flops_gemm(2, blk_e0[b + 1] - blk_e0[b], h_np8[b]);
	//Columns per batch: two per reflection and symmetry operation, as many as the phase
	//buffer and the result rows allow
	const long long budget = 512LL << 20;
	d.ncol_cap = (int)std::min<long long>(512, budget / (sizeof(T) * d.np8_total)) & ~7;
	d.batch_max = (int)std::min<long long>(d.ncol_cap / 2, (256LL << 20) / (16 * L.packed));
	if (d.batch_max < 1) return false;
	d.ncol_cap = std::max(8, std::min(d.ncol_cap, (2 * d.batch_max + 7) & ~7));
	//Pieces of at most piece_cap table elements, chunks of at most tab_cap: the whole table
	//when the device holds it, so it is built once, else the cap halves until the plan fits
	//and the tables are rebuilt every batch
	const long long piece_cap = (128LL << 20) / sizeof(T);
	size_t freeb = 0, totalb = 0;
	if (gpuMemGetInfo(&freeb, &totalb) != gpuSuccess) return false;
	long long tab_total = 0, max_piece = 0;
	for (int b = 0; b < nb; b++) tab_total += (long long)(blk_e0[b + 1] - blk_e0[b]) * h_np8[b];
	d.pieces.clear();
	for (int b = 0; b < nb; b++) {
		const int rows_cap = (int)std::max<long long>(1, piece_cap / h_np8[b]);
		for (int e = blk_e0[b]; e < blk_e0[b + 1]; e += rows_cap) {
			typename Dev<T>::Piece pc;
			pc.blk = b; pc.e0 = e; pc.rows = std::min(rows_cap, blk_e0[b + 1] - e); pc.tab = 0;
			d.pieces.push_back(pc);
			max_piece = std::max(max_piece, (long long)pc.rows * h_np8[b]);
		}
	}
	d.n_pieces = (int)d.pieces.size();
	size_t max_ws = 0, need = 0;
	d.tab_cap = tab_total;
	for (;;) {
		d.chunks.clear();
		d.rows_cap = 0;
		long long tab = 0;
		for (int p = 0; p < d.n_pieces;) {
			typename Dev<T>::Chunk ch;
			ch.p0 = p; ch.tab0 = tab; ch.e0 = d.pieces[p].e0;
			long long used = 0;
			while (p < d.n_pieces && used + (long long)d.pieces[p].rows * h_np8[d.pieces[p].blk] <= d.tab_cap) {
				d.pieces[p].tab = tab;
				const long long sz = (long long)d.pieces[p].rows * h_np8[d.pieces[p].blk];
				used += sz; tab += sz; p++;
			}
			ch.p1 = p; ch.e1 = d.pieces[p - 1].e0 + d.pieces[p - 1].rows;
			d.rows_cap = std::max(d.rows_cap, (long long)(ch.e1 - ch.e0));
			d.chunks.push_back(ch);
		}
		d.n_chunks = (int)d.chunks.size();
		max_ws = 0;
		for (int p = 0; p < d.n_pieces; p++)
			max_ws = std::max(max_ws, itensor_gemm::workspace_bytes<T>(d.ncol_cap, d.pieces[p].rows, h_np8[d.pieces[p].blk]));
		need = sizeof(double) * 4 * (size_t)d.np8_total + sizeof(T) * (size_t)ao_total
			+ sizeof(T) * (size_t)d.np8_total * d.ncol_cap
			+ sizeof(T) * (size_t)d.tab_cap + sizeof(T) * (size_t)d.rows_cap * d.ncol_cap + max_ws
			+ sizeof(int) * (5 * (size_t)d.n_entries + (size_t)(d.n_chunks + 1) * L.packed + 2 * (size_t)nb + 2 * (size_t)d.n_pieces)
			+ sizeof(long long) * ((size_t)nb + d.n_pieces)
			+ 2 * sizeof(double) * ((size_t)d.ncol_cap / 2 * (3 + 2 * (size_t)L.n_grids) + 2 * (size_t)d.batch_max * L.packed);
		if (need + (1u << 28) <= freeb || d.tab_cap <= max_piece) break;
		d.tab_cap = std::max(max_piece, d.tab_cap / 2);
	}
	if (need + (1u << 28) > freeb) return false;
	if (throughput::enabled())
		std::fprintf(stderr, "I tensor GPU: %d blocks, %lld table rows in %d pieces and %d chunks, %d columns a batch, %.0f MB\n",
			nb, d.n_entries, d.n_pieces, d.n_chunks, d.ncol_cap, need / 1048576.0);
	//Which rows feed each stored pair, and where each chunk's share of them starts
	std::vector<int> ent_piece(d.n_entries);
	for (int p = 0; p < d.n_pieces; p++)
		for (int e = d.pieces[p].e0; e < d.pieces[p].e0 + d.pieces[p].rows; e++) ent_piece[e] = p;
	std::vector<int> cnt(L.packed + 1, 0), acc_e(d.n_entries), acc_ptr((size_t)(d.n_chunks + 1) * L.packed);
	for (long long e = 0; e < d.n_entries; e++) cnt[ent_t[e] + 1]++;
	for (int t = 0; t < L.packed; t++) cnt[t + 1] += cnt[t];
	{
		std::vector<int> fill(cnt.begin(), cnt.end() - 1);
		for (int e = 0; e < (int)d.n_entries; e++) acc_e[fill[ent_t[e]]++] = e;
	}
	for (int c = 0; c <= d.n_chunks; c++) {
		const int e_lim = c < d.n_chunks ? d.chunks[c].e0 : (int)d.n_entries;
		for (int t = 0; t < L.packed; t++) {
			int k = c == 0 ? cnt[t] : acc_ptr[(size_t)(c - 1) * L.packed + t];
			while (k < cnt[t + 1] && acc_e[k] < e_lim) k++;
			acc_ptr[(size_t)c * L.packed + t] = k;
		}
	}
	std::vector<int> pc_blk(d.n_pieces), pc_e0(d.n_pieces);
	std::vector<long long> pc_tab(d.n_pieces);
	for (int p = 0; p < d.n_pieces; p++) {
		pc_blk[p] = d.pieces[p].blk; pc_e0[p] = d.pieces[p].e0; pc_tab[p] = d.pieces[p].tab;
	}
	//Points and AO values in the padded layout
	{
		std::vector<double> pd1(d.np8_total, 0.0), pd2(d.np8_total, 0.0), pd3(d.np8_total, 0.0), pw(d.np8_total, 0.0);
		for (int b = 0; b < nb; b++) {
			const int src = L.grid_point_off[L.blk_grid[b]] + L.blk_point_start[b];
			for (int p = 0; p < L.blk_point_count[b]; p++) {
				pd1[h_pp[b] + p] = L.d1[src + p]; pd2[h_pp[b] + p] = L.d2[src + p];
				pd3[h_pp[b] + p] = L.d3[src + p]; pw[h_pp[b] + p] = L.weights[src + p];
			}
		}
		if (!upload_vec(&d.d1, pd1) || !upload_vec(&d.d2, pd2) || !upload_vec(&d.d3, pd3) || !upload_vec(&d.w, pw))
			return false;
	}
	GPU_TRY(gpuMalloc(&d.ao, sizeof(T) * (size_t)std::max<long long>(ao_total, 1)));
	{
		std::vector<T> stage;
		for (int b = 0; b < nb; b++) {
			const int na = L.blk_n_active[b], np = L.blk_point_count[b], np8 = h_np8[b];
			const double* src = L.ao_all + L.blk_ao_off[b];
			stage.assign((size_t)na * np8, T(0));
			for (int i = 0; i < na; i++)
				for (int p = 0; p < np; p++) stage[(size_t)i * np8 + p] = (T)src[(size_t)i * np + p];
			GPU_TRY(gpuMemcpy(d.ao + h_ao[b], stage.data(), sizeof(T) * stage.size(), gpuMemcpyHostToDevice));
		}
	}
	if (!upload_vec(&d.q_pp, h_pp) || !upload_vec(&d.q_np8, h_np8) || !upload_vec(&d.q_ao, h_ao)) return false;
	if (!upload_vec(&d.ent_i, ent_i) || !upload_vec(&d.ent_j, ent_j) || !upload_vec(&d.ent_grid, ent_grid)
		|| !upload_vec(&d.ent_piece, ent_piece)) return false;
	if (!upload_vec(&d.pc_blk, pc_blk) || !upload_vec(&d.pc_e0, pc_e0) || !upload_vec(&d.pc_tab, pc_tab)) return false;
	if (!upload_vec(&d.acc_ptr, acc_ptr) || !upload_vec(&d.acc_e, acc_e)) return false;
	GPU_TRY(gpuMalloc(&d.phase, sizeof(T) * (size_t)d.np8_total * d.ncol_cap));
	GPU_TRY(gpuMalloc(&d.tab, sizeof(T) * (size_t)d.tab_cap));
	GPU_TRY(gpuMalloc(&d.cres, sizeof(T) * (size_t)d.rows_cap * d.ncol_cap));
	GPU_TRY(gpuMalloc(&d.gemm_ws, max_ws ? max_ws : 1));
	const size_t ncomb = (size_t)d.ncol_cap / 2, nI = (size_t)d.batch_max * L.packed;
	for (int i = 0; i < 2; i++) {
		GPU_TRY(gpuMalloc(&d.kvec[i], sizeof(double) * 3 * ncomb));
		GPU_TRY(gpuMalloc(&d.fac[i], sizeof(double) * 2 * ncomb * L.n_grids));
		GPU_TRY(gpuHostAlloc((void**)&d.host_kf[i], sizeof(double) * ncomb * (3 + 2 * (size_t)L.n_grids)));
		GPU_TRY(gpuMalloc(&d.I_re[i], sizeof(double) * nI));
		GPU_TRY(gpuMalloc(&d.I_im[i], sizeof(double) * nI));
		GPU_TRY(gpuHostAlloc((void**)&d.host_re[i], sizeof(double) * nI));
		GPU_TRY(gpuHostAlloc((void**)&d.host_im[i], sizeof(double) * nI));
		GPU_TRY(gpuEventCreate(&d.done[i]));
	}
	GPU_TRY(gpuStreamCreateNonBlocking(&d.copy_stream));
	d.h_pp.assign(h_pp.begin(), h_pp.end());
	d.h_np8 = h_np8;
	d.nmo = L.nmo; d.packed = L.packed; d.n_grids = L.n_grids; d.n_blocks = nb;
	d.table_ready = false;
	d.ready = true;
	return true;
}

template <typename T>
int batch_impl(const int num_syms)
{
	const Dev<T>& d = g<T>;
	if (!d.ready) return 0;
	return std::max(1, std::min(d.batch_max, d.ncol_cap / (2 * num_syms)));
}

template <typename T>
bool submit_impl(const int slot, const int n_refl, const int num_syms,
	const double* kx, const double* ky, const double* kz,
	const std::complex<double>* factors)
{
	Dev<T>& d = g<T>;
	if (!d.ready || slot < 0 || slot > 1 || n_refl < 1 || n_refl > d.batch_max) return false;
	const int ncomb = n_refl * num_syms, ncol8 = (2 * ncomb + 7) & ~7;
	if (ncol8 > d.ncol_cap) return false;
	const size_t nf = (size_t)ncomb * d.n_grids;
	double* const hk = d.host_kf[slot];
	double* const hf = hk + 3 * (size_t)ncomb;
	for (int c = 0; c < ncomb; c++) {
		//The CPU path takes sin/cos of k.d directly; scaling k to turns here is what lets
		//the reduction be a rint and the transcendental be sincospi
		hk[3 * c] = kx[c] * SF_INV_TWO_PI; hk[3 * c + 1] = ky[c] * SF_INV_TWO_PI; hk[3 * c + 2] = kz[c] * SF_INV_TWO_PI;
	}
	for (size_t i = 0; i < nf; i++) {
		hf[i] = factors[i].real();
		hf[nf + i] = factors[i].imag();
	}
	GPU_TRY(gpuMemcpyAsync(d.kvec[slot], hk, sizeof(double) * 3 * (size_t)ncomb, gpuMemcpyHostToDevice, 0));
	GPU_TRY(gpuMemcpyAsync(d.fac[slot], hf, sizeof(double) * 2 * nf, gpuMemcpyHostToDevice, 0));
	const long long nI = (long long)n_refl * d.packed;
	zero_kernel<<<(unsigned int)((nI + 255) / 256), 256>>>(nI, d.I_re[slot], d.I_im[slot]);
	phase_kernel<T><<<dim3((d.np8_total + 255) / 256, (ncol8 / 2 + phase_run - 1) / phase_run), 256>>>(
		d.np8_total, ncomb, ncol8, d.kvec[slot], d.d1, d.d2, d.d3, d.w, d.phase);
	for (int c = 0; c < d.n_chunks; c++) {
		const typename Dev<T>::Chunk& ch = d.chunks[c];
		const int rows = ch.e1 - ch.e0;
		if (d.n_chunks > 1 || !d.table_ready)
			table_kernel<T><<<dim3((d.np8_max + 255) / 256, std::min(rows, 65535)), 256>>>(ch.e0, rows, ch.tab0,
				d.ent_i, d.ent_j, d.ent_piece, d.pc_blk, d.pc_e0, d.pc_tab, d.q_np8, d.q_ao, d.ao, d.tab);
		//C = P^T * table per piece: P column-major np8 x ncol8 from the block's first point,
		//the table column-major np8 x rows, C column-major ncol8 x rows
		for (int p = ch.p0; p < ch.p1; p++) {
			const typename Dev<T>::Piece& pc = d.pieces[p];
			const int np8 = d.h_np8[pc.blk];
			if (!itensor_gemm::run<T>(ncol8, pc.rows, np8,
				d.phase + d.h_pp[pc.blk], d.np8_total,
				d.tab + pc.tab - ch.tab0, np8,
				d.cres + (long long)(pc.e0 - ch.e0) * ncol8, ncol8, d.gemm_ws))
				return false;
		}
		gather_kernel<T><<<(unsigned int)((nI + 255) / 256), 256>>>(d.packed, n_refl, num_syms, d.n_grids,
			ncol8, ch.e0, d.acc_ptr + (size_t)c * d.packed, d.acc_e, d.ent_grid, d.cres,
			d.fac[slot], d.fac[slot] + nf, d.I_re[slot], d.I_im[slot]);
	}
	d.table_ready = true;
	d.n_refl[slot] = n_refl;
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuEventRecord(d.done[slot], 0));
	return true;
}

template <typename T>
bool collect_impl(const int slot, std::complex<double>* I_r, const long long row_stride)
{
	Dev<T>& d = g<T>;
	if (!d.ready || slot < 0 || slot > 1 || d.n_refl[slot] < 1) return false;
	const size_t nI = (size_t)d.n_refl[slot] * d.packed;
	GPU_TRY(gpuStreamWaitEvent(d.copy_stream, d.done[slot], 0));
	GPU_TRY(gpuMemcpyAsync(d.host_re[slot], d.I_re[slot], sizeof(double) * nI, gpuMemcpyDeviceToHost, d.copy_stream));
	GPU_TRY(gpuMemcpyAsync(d.host_im[slot], d.I_im[slot], sizeof(double) * nI, gpuMemcpyDeviceToHost, d.copy_stream));
	GPU_TRY(gpuStreamSynchronize(d.copy_stream));
	for (int r = 0; r < d.n_refl[slot]; r++)
		for (int i = 0; i < d.packed; i++)
			I_r[r * row_stride + i] += std::complex<double>(d.host_re[slot][(size_t)r * d.packed + i],
				d.host_im[slot][(size_t)r * d.packed + i]);
	d.n_refl[slot] = 0;
	return true;
}

} //namespace

template <typename T>
__global__ void hold_rows_kernel(const T* I, const double* w, const double* F0, double* F, const int packed)
{
	const size_t base = (size_t)blockIdx.x * packed * 2;
	double sr = 0.0, si = 0.0;
	for (int k = threadIdx.x; k < packed; k += blockDim.x) {
		sr += (double)I[base + 2 * k] * w[k];
		si += (double)I[base + 2 * k + 1] * w[k];
	}
	for (int o = 16; o > 0; o >>= 1) {
		sr += __shfl_down_sync(0xffffffffu, sr, o);
		si += __shfl_down_sync(0xffffffffu, si, o);
	}
	__shared__ double red[2][8];
	if ((threadIdx.x & 31) == 0) { red[0][threadIdx.x >> 5] = sr; red[1][threadIdx.x >> 5] = si; }
	__syncthreads();
	if (threadIdx.x == 0) {
		sr = F0[2 * blockIdx.x]; si = F0[2 * blockIdx.x + 1];
		for (int i = 0; i < blockDim.x / 32; i++) { sr += red[0][i]; si += red[1][i]; }
		F[2 * blockIdx.x] = sr; F[2 * blockIdx.x + 1] = si;
	}
}

template <typename T>
__global__ void hold_cols_kernel(const T* I, const double* pre, double* part, const int nr, const int packed, const int rchunk)
{
	const int k = blockIdx.x * blockDim.x + threadIdx.x;
	if (k >= packed) return;
	const int r1 = min(nr, (int)(blockIdx.y + 1) * rchunk);
	double sum = 0.0;
	for (int r = blockIdx.y * rchunk; r < r1; r++) {
		const size_t e = ((size_t)r * packed + k) * 2;
		sum += pre[2 * r] * (double)I[e] - pre[2 * r + 1] * (double)I[e + 1];
	}
	part[(size_t)blockIdx.y * packed + k] = sum;
}

__global__ void hold_sum_kernel(const double* part, double* out, const int chunks, const int packed)
{
	const int k = blockIdx.x * blockDim.x + threadIdx.x;
	if (k >= packed) return;
	double sum = 0.0;
	for (int c = 0; c < chunks; c++) sum += part[(size_t)c * packed + k];
	out[k] = sum;
}

template <typename T>
bool hold_impl(const std::complex<T>* I, const int nr, const int packed)
{
	itensor_gpu_release();
	if (!itensor_gpu_available() || nr <= 0 || packed <= 0) return false;
	Held& h = g_held;
	const size_t bytes = sizeof(std::complex<T>) * (size_t)nr * packed;
	size_t free_b = 0, total_b = 0;
	if (gpuMemGetInfo(&free_b, &total_b) != gpuSuccess || bytes + (size_t)64 * 1048576 > free_b) return false;
	if (gpuMalloc(&h.I, bytes) != gpuSuccess) { gpuGetLastError(); return false; }
	h.fp64 = sizeof(T) == sizeof(double);
	h.nr = nr; h.packed = packed;
	const bool ok = gpuMemcpy(h.I, I, bytes, gpuMemcpyHostToDevice) == gpuSuccess
		&& gpuMalloc(&h.w, sizeof(double) * packed) == gpuSuccess
		&& gpuMalloc(&h.pre, sizeof(double) * 2 * nr) == gpuSuccess
		&& gpuMalloc(&h.F0, sizeof(double) * 2 * nr) == gpuSuccess
		&& gpuMalloc(&h.F, sizeof(double) * 2 * nr) == gpuSuccess
		&& gpuMalloc(&h.part, sizeof(double) * (size_t)hold_chunks * packed) == gpuSuccess
		&& gpuMalloc(&h.out, sizeof(double) * packed) == gpuSuccess;
	if (!ok) { gpuGetLastError(); itensor_gpu_release(); }
	return ok;
}

//Shared with the transform so the "no code for this card" case is diagnosed in one place.
bool itensor_gpu_available() { return sf_gpu_available(); }

bool itensor_gpu_hold(const std::complex<float>* I, const int nr, const int packed) { return hold_impl(I, nr, packed); }
bool itensor_gpu_hold(const std::complex<double>* I, const int nr, const int packed) { return hold_impl(I, nr, packed); }
bool itensor_gpu_held() { return g_held.I != nullptr; }

bool itensor_gpu_rows(const double* w, const std::complex<double>* F0, std::complex<double>* F)
{
	Held& h = g_held;
	if (!h.I) return false;
	GPU_TRY(gpuMemcpy(h.w, w, sizeof(double) * h.packed, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(h.F0, F0, sizeof(double) * 2 * h.nr, gpuMemcpyHostToDevice));
	if (h.fp64) hold_rows_kernel<double><<<h.nr, 256>>>((const double*)h.I, h.w, h.F0, h.F, h.packed);
	else hold_rows_kernel<float><<<h.nr, 256>>>((const float*)h.I, h.w, h.F0, h.F, h.packed);
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuMemcpy(F, h.F, sizeof(double) * 2 * h.nr, gpuMemcpyDeviceToHost));
	return true;
}

bool itensor_gpu_cols(const std::complex<double>* pre, double* out)
{
	Held& h = g_held;
	if (!h.I) return false;
	GPU_TRY(gpuMemcpy(h.pre, pre, sizeof(double) * 2 * h.nr, gpuMemcpyHostToDevice));
	const int rchunk = (h.nr + hold_chunks - 1) / hold_chunks;
	const dim3 grid((h.packed + 255) / 256, hold_chunks);
	if (h.fp64) hold_cols_kernel<double><<<grid, 256>>>((const double*)h.I, h.pre, h.part, h.nr, h.packed, rchunk);
	else hold_cols_kernel<float><<<grid, 256>>>((const float*)h.I, h.pre, h.part, h.nr, h.packed, rchunk);
	GPU_TRY(gpuGetLastError());
	hold_sum_kernel<<<(h.packed + 255) / 256, 256>>>(h.part, h.out, hold_chunks, h.packed);
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuMemcpy(out, h.out, sizeof(double) * h.packed, gpuMemcpyDeviceToHost));
	return true;
}

//One block per row ab of the packed integrals, as XCW::eri_JK walks them: J1 is the row's
//half of the symmetric matvec, the four scatters of every integral go into the row's two K
//columns, each warp taking segments c of the row with the lanes over d and its own copy of
//the columns in shared memory, so nothing is atomic. The multiplicities of the pairs are
//folded into four copies of the two density columns. The partials of a row leave as K(:, a)
//and K(:, b); the diagonal cd == ab was counted with the doubled weight and is taken back by
//thread 0.
__global__ void eri_jk_kernel(const double* V, const double* dp, const double* D, double* J1, double* Ka, double* Kb, const int n)
{
	extern __shared__ double sh[];
	__shared__ double red[32];
	const int ab = blockIdx.x, a = (int)((sqrt(8.0 * ab + 1.0) - 1.0) / 2.0), b = ab - a * (a + 1) / 2;
	const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5, nwarp = blockDim.x >> 5;
	double* Da1 = sh;
	double* Db1 = sh + n;
	double* Da2 = sh + 2 * n;
	double* Db2 = sh + 3 * n;
	double* sKa = sh + (4 + 2 * warp) * n;
	double* sKb = sKa + n;
	const double* v = V + (size_t)ab * (ab + 1) / 2;
	const double fab = a == b ? 2.0 : 4.0;
	for (int d = threadIdx.x; d < n; d += blockDim.x) {
		Da1[d] = fab * D[d + a * n];
		Db1[d] = fab * D[d + b * n];
		Da2[d] = 2.0 * Da1[d];
		Db2[d] = 2.0 * Db1[d];
	}
	for (int d = threadIdx.x; d < 2 * nwarp * n; d += blockDim.x) sh[4 * n + d] = 0.0;
	__syncthreads();
	double j = 0.0;
	for (int c = warp; c <= a; c += nwarp) {
		const int m = c < a ? c + 1 : b + 1, off = c * (c + 1) / 2;
		double sa = 0.0, sb = 0.0;
		for (int d = lane; d < m; d += 32) {
			const double t = v[off + d];
			const double* Daw = d < c ? Da2 : Da1;
			const double* Dbw = d < c ? Db2 : Db1;
			j += t * dp[off + d];
			sa += t * Dbw[d];
			sb += t * Daw[d];
			sKa[d] += t * Dbw[c];
			sKb[d] += t * Daw[c];
		}
		for (int o = 16; o > 0; o >>= 1) {
			sa += __shfl_down_sync(0xffffffffu, sa, o);
			sb += __shfl_down_sync(0xffffffffu, sb, o);
		}
		if (lane == 0) {
			sKa[c] += sa;
			sKb[c] += sb;
		}
	}
	for (int o = 16; o > 0; o >>= 1) j += __shfl_down_sync(0xffffffffu, j, o);
	if (lane == 0) red[warp] = j;
	__syncthreads();
	if (threadIdx.x == 0) {
		for (int w = 1; w < nwarp; w++) j += red[w];
		J1[ab] = j;
		const double h = (a == b ? 1.0 : 4.0) * v[ab];
		sh[4 * n + a] -= h * D[b + b * n];
		sh[5 * n + a] -= h * D[b + a * n];
		sh[4 * n + b] -= h * D[b + a * n];
		sh[5 * n + b] -= h * D[a + a * n];
	}
	__syncthreads();
	for (int d = threadIdx.x; d < 2 * n; d += blockDim.x) {
		double sum = 0.0;
		for (int w = 0; w < nwarp; w++) sum += sh[(4 + 2 * w) * n + d];
		if (d < n) Ka[(size_t)ab * n + d] = sum;
		else Kb[(size_t)ab * n + d - n] = sum;
	}
}

//K(d, a) is the sum of the row partials over the rows of a, K(:, a) from the rows ab with
//b <= a and K(:, b) from the rows a'b with a' >= a
__global__ void eri_kred_kernel(const double* Ka, const double* Kb, double* K, const int n)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= n * n) return;
	const int d = i % n, a = i / n;
	double sum = 0.0;
	for (int b = 0; b <= a; b++) sum += Ka[(size_t)(a * (a + 1) / 2 + b) * n + d];
	for (int a2 = a; a2 < n; a2++) sum += Kb[(size_t)(a2 * (a2 + 1) / 2 + a) * n + d];
	K[i] = sum;
}

//The other half of the symmetric matvec, J2[cd] = sum_{ab > cd} V[ab][cd] dp[ab]: consecutive
//threads read consecutive elements of a row, the rows cut into chunks for parallelism
__global__ void eri_jt_kernel(const double* V, const double* dp, double* part, const int npair, const int chunk)
{
	const int cd = blockIdx.x * blockDim.x + threadIdx.x;
	if (cd >= npair) return;
	const int ab0 = max(cd + 1, (int)blockIdx.y * chunk), ab1 = min(npair, (int)(blockIdx.y + 1) * chunk);
	double sum = 0.0;
	for (int ab = ab0; ab < ab1; ab++) sum += V[(size_t)ab * (ab + 1) / 2 + cd] * dp[ab];
	part[(size_t)blockIdx.y * npair + cd] = sum;
}

bool eri_gpu_hold(const double* eri, const int n)
{
	eri_gpu_release();
	if (!itensor_gpu_available() || n <= 0 || 12 * n * sizeof(double) > 48 * 1024) return false;
	HeldEri& h = g_eri;
	const int npair = n * (n + 1) / 2;
	const size_t bytes = sizeof(double) * (size_t)npair * (npair + 1) / 2, kbytes = sizeof(double) * (size_t)npair * n;
	size_t free_b = 0, total_b = 0;
	if (gpuMemGetInfo(&free_b, &total_b) != gpuSuccess || bytes + 2 * kbytes + (size_t)64 * 1048576 > free_b) return false;
	if (gpuMalloc(&h.V, bytes) != gpuSuccess) { gpuGetLastError(); return false; }
	h.n = n; h.npair = npair;
	const bool ok = gpuMemcpy(h.V, eri, bytes, gpuMemcpyHostToDevice) == gpuSuccess
		&& gpuMalloc(&h.dp, sizeof(double) * npair) == gpuSuccess
		&& gpuMalloc(&h.D, sizeof(double) * n * n) == gpuSuccess
		&& gpuMalloc(&h.J1, sizeof(double) * npair) == gpuSuccess
		&& gpuMalloc(&h.J2, sizeof(double) * npair) == gpuSuccess
		&& gpuMalloc(&h.part, sizeof(double) * (size_t)hold_chunks * npair) == gpuSuccess
		&& gpuMalloc(&h.Ka, kbytes) == gpuSuccess
		&& gpuMalloc(&h.Kb, kbytes) == gpuSuccess
		&& gpuMalloc(&h.K, sizeof(double) * n * n) == gpuSuccess;
	if (!ok) { gpuGetLastError(); eri_gpu_release(); }
	return ok;
}

bool eri_gpu_JK(const double* D, double* J, double* K)
{
	HeldEri& h = g_eri;
	if (!h.V) return false;
	const int n = h.n, npair = h.npair;
	std::vector<double> dp(npair), j1(npair), j2(npair);
	for (int a = 0, ab = 0; a < n; a++)
		for (int b = 0; b <= a; b++, ab++) dp[ab] = (a == b ? 1.0 : 2.0) * D[a + b * n];
	GPU_TRY(gpuMemcpy(h.dp, dp.data(), sizeof(double) * npair, gpuMemcpyHostToDevice));
	GPU_TRY(gpuMemcpy(h.D, D, sizeof(double) * n * n, gpuMemcpyHostToDevice));
	eri_jk_kernel<<<npair, 128, sizeof(double) * 12 * n>>>(h.V, h.dp, h.D, h.J1, h.Ka, h.Kb, n);
	GPU_TRY(gpuGetLastError());
	eri_kred_kernel<<<(n * n + 255) / 256, 256>>>(h.Ka, h.Kb, h.K, n);
	GPU_TRY(gpuGetLastError());
	const int chunk = (npair + hold_chunks - 1) / hold_chunks;
	eri_jt_kernel<<<dim3((npair + 255) / 256, hold_chunks), 256>>>(h.V, h.dp, h.part, npair, chunk);
	GPU_TRY(gpuGetLastError());
	hold_sum_kernel<<<(npair + 255) / 256, 256>>>(h.part, h.J2, hold_chunks, npair);
	GPU_TRY(gpuGetLastError());
	GPU_TRY(gpuMemcpy(j1.data(), h.J1, sizeof(double) * npair, gpuMemcpyDeviceToHost));
	GPU_TRY(gpuMemcpy(j2.data(), h.J2, sizeof(double) * npair, gpuMemcpyDeviceToHost));
	GPU_TRY(gpuMemcpy(K, h.K, sizeof(double) * n * n, gpuMemcpyDeviceToHost));
	for (int a = 0, ab = 0; a < n; a++)
		for (int b = 0; b <= a; b++, ab++) J[a + b * n] = J[b + a * n] = j1[ab] + j2[ab];
	for (int a = 0; a < n; a++)
		for (int b = 0; b <= a; b++) K[a + b * n] = K[b + a * n] = 0.125 * (K[a + b * n] + K[b + a * n]);
	return true;
}

void eri_gpu_release()
{
	HeldEri& h = g_eri;
	gpuFree(h.V); gpuFree(h.dp); gpuFree(h.D); gpuFree(h.J1); gpuFree(h.J2); gpuFree(h.part); gpuFree(h.Ka); gpuFree(h.Kb); gpuFree(h.K);
	h = HeldEri{};
}

void itensor_gpu_release()
{
	Held& h = g_held;
	gpuFree(h.I); gpuFree(h.w); gpuFree(h.pre); gpuFree(h.F0); gpuFree(h.F); gpuFree(h.part); gpuFree(h.out);
	h = Held{};
}

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

int itensor_gpu_batch(const int num_syms)
{
	return g_fp64 ? batch_impl<double>(num_syms) : batch_impl<float>(num_syms);
}

bool itensor_gpu_submit(const int slot, const int n_refl, const int num_syms,
	const double* kx, const double* ky, const double* kz,
	const std::complex<double>* factors)
{
	return g_fp64 ? submit_impl<double>(slot, n_refl, num_syms, kx, ky, kz, factors)
	              : submit_impl<float>(slot, n_refl, num_syms, kx, ky, kz, factors);
}

bool itensor_gpu_collect(const int slot, std::complex<double>* I_r, const long long row_stride)
{
	return g_fp64 ? collect_impl<double>(slot, I_r, row_stride) : collect_impl<float>(slot, I_r, row_stride);
}

void itensor_gpu_free()
{
	g_tensor = false;
	itensor_gemm::set_tensor_mode(false);
	free_impl<float>();
	free_impl<double>();
}
