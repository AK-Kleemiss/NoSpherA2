#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "mo_class.h"
#include "cube.h"
#include "constants.h"
#include "fchk.h"
#include "basis_set.h"
#include "nos_math.h"
#include "libCintMain.h"
#include "eli_family.h"
#include "integrator.h"
#include "cell.h"
#include "wfn_kernels.h"
#ifdef NOSPHERA2_USE_GPU
#include "aux_density_gpu.h"
#include "esp_gpu.h"
#include "basin_field_gpu.h"
#endif

#if defined(__aarch64__) || defined(_M_ARM64)
namespace
{
	//Set by compute_dens_batch: compute_dens_cartesian sums the primitives into this scratch and returns before the MOs
	thread_local ao_scratch<1>* batch_ao = nullptr;
}
#endif

const double WFN::compute_dens(
	const d3 &Pos,
	vec2 &d,
	vec &phi) const
{
	//if (d_f_switch)
	//{
	//    err_checkf(d.size() >= 5, "d is too small!", std::cout);
	//    err_checkf(phi.size() >= get_nmo(true), "phi is too small!", std::cout);
	//    return compute_dens_spherical(Pos1, Pos2, Pos3, d, phi);
	//}
	//else
	//{
	err_checkf(d.size() >= ncen, "d is too small!", std::cout);
	//the kernel accumulates every MO, virtuals included; callers size phi by the occupied count
	if (phi.size() < (size_t)nmo) phi.resize(nmo);
	return compute_dens_cartesian(Pos, d, phi);
	//}
};

const double WFN::compute_dens(
	const d3 &Pos)
	const
{
	//Called per point of the basin lookup: the ncen + 1 allocations per call were
	//most of its allocator time, so the buffers persist per thread (compute_dens_cartesian zeroes phi)
	thread_local vec2 d;
	thread_local vec phi;
	if (phi.size() < (size_t)nmo) phi.resize(nmo);
	if (d.size() < (size_t)ncen) d.resize(ncen, vec(16, 0.0));
	return compute_dens_cartesian(Pos, d, phi);
};

const double WFN::compute_spin_dens(
	const d3 &Pos,
	vec2 &d,
	vec &phi) const
{
	//if (d_f_switch)
	//{
	//    err_checkf(d.size() >= 5, "d is too small!", std::cout);
	//    err_checkf(phi.size() >= get_nmo(true), "phi is too small!", std::cout);
	//    err_not_impl_f("Nah.. not yet implemented correctly", std::cout);
	//    return compute_dens_spherical(Pos1, Pos2, Pos3, d, phi);
	//}
	//else
	//{
	err_checkf(d.size() >= ncen, "d is too small!", std::cout);
	if (phi.size() < (size_t)nmo) phi.resize(nmo);
	return compute_spin_dens_cartesian(Pos, d, phi);
	//}
};

const double WFN::compute_spin_dens(
	const d3 &Pos) const
{
	vec2 d;
	vec phi(nmo, 0.0);

	//if (d_f_switch)
	//{
	//    d.resize(5);
	//    for (int i = 0; i < 5; i++)
	//        d[i].resize(ncen, 0.0);
	//    err_not_impl_f("Nah.. not yet implemented correctly", std::cout);
	//    return compute_dens_spherical(Pos1, Pos2, Pos3, d, phi);
	//}
	//else
	//{
	d.resize(ncen);
	for (int i = 0; i < ncen; i++)
		d[i].resize(16, 0.0);
	return compute_spin_dens_cartesian(Pos, d, phi);
	//}
};

const double WFN::compute_dens_cartesian(
	const d3 &Pos,
	vec2 &d,
	vec &phi) const
{
	std::fill(phi.begin(), phi.end(), 0.0);
	double Rho = 0.0;
	int j;
	double ex, *d_;

	// precalculate some distances and powers of distances for faster computation
	for (j = 0; j < ncen; j++)
		cart_powers(d[j].data(), Pos, atoms[j]);

	// Pre-cache frequently accessed data
	const int *centers_data = centers.data();
	const int *types_data = types.data();
	const double *exponents_data = exponents.data();
	const MO *MOs_data = MOs.data();
	double *phi_data = phi.data();
	//Builds the exponent groups and the contraction the loop below reads
	if (!get_coef_primitive_major()) return 0.0;
	thread_local vec exps;
	if (exps.size() < group_exponent.size()) exps.resize(group_exponent.size());
	exp_table([&d](const int c) { return d[c][3]; }, exps.data());
	const int *group = prim_exp_group.data();
#if defined(__aarch64__) || defined(_M_ARM64)
	thread_local ao_scratch<1> ao_own;
	ao_scratch<1> &ao = batch_ao ? *batch_ao : ao_own;
#else
	thread_local ao_scratch<1> ao;
#endif
	ao.begin(coef_ao_major.size() / nmo);

	for (j = 0; j < nex; j++)
	{
		ex = exps[group[j]];
		if (ex == 0.0 || prim_ao[j] < 0)
		{ // corresponds to cutoff of maximum density contribution of 1E-5
			continue;
		}
		d_ = d[centers_data[j] - 1].data();
		ex *= cart_monomial(types_data[j], d_);

		ao.add(prim_ao[j], prim_ao_scale[j], &ex);
	}
#if defined(__aarch64__) || defined(_M_ARM64)
	if (batch_ao) return 0.0; //compute_dens_batch multiplies into the MOs itself
#endif
	ao.to_mo(coef_ao_major.data(), nmo, phi_data, [MOs_data](const int mo) { return MOs_data[mo].get_occ(); });

	// use pointer arithmetic and minimize overhead
	const double *phi_ptr = phi_data;
	const double *phi_end = phi_ptr + nmo;
	const MO *mo_ptr = MOs_data;

	for (; phi_ptr != phi_end; ++phi_ptr, ++mo_ptr)
	{
		const double &phi_val = *phi_ptr;
		Rho += mo_ptr->get_occ() * phi_val * phi_val;
	}

	return Rho;
}

#if defined(__aarch64__) || defined(_M_ARM64)
namespace
{
	//R consecutive union rows c[r], weights w[r][b], into phi[b * nmo + mo] of eight points over the occupied runs. Each
	//(point, MO) still adds the rows in union order, one fma each, so only the phi loads and stores are shared by R FMAs
	template <int R>
	void batch_rows(const double* const* c, const double (*w)[8], const ivec& runs, const int nmo, double* phi)
	{
		float64x2_t wv[R][4];
		for (int r = 0; r < R; r++)
			for (int k = 0; k < 4; k++) wv[r][k] = vld1q_f64(w[r] + 2 * k);
		for (size_t q = 0; q < runs.size(); q += 2) {
			int mo = runs[q];
			for (; mo + 2 <= runs[q + 1]; mo += 2) {
				float64x2_t cv[R];
				for (int r = 0; r < R; r++) cv[r] = vld1q_f64(c[r] + mo);
				for (int k = 0; k < 4; k++) {
					double* p0 = phi + (size_t)(2 * k) * nmo + mo;
					double* p1 = p0 + nmo;
					float64x2_t a0 = vld1q_f64(p0), a1 = vld1q_f64(p1);
					for (int r = 0; r < R; r++) {
						a0 = vfmaq_laneq_f64(a0, cv[r], wv[r][k], 0);
						a1 = vfmaq_laneq_f64(a1, cv[r], wv[r][k], 1);
					}
					vst1q_f64(p0, a0);
					vst1q_f64(p1, a1);
				}
			}
			if (mo < runs[q + 1])
				for (int b = 0; b < 8; b++) {
					double& p = phi[(size_t)b * nmo + mo];
					for (int r = 0; r < R; r++) p = std::fma(w[r][b], c[r][mo], p);
				}
		}
	}
}

//ARM: to_mo reads one coefficient row per touched function and point, which a Cortex-A72 fetches from L2 or memory at
//2-5 ns per FMA. Eight neighbouring points share the rows: the union of their functions is walked once and each row
//feeds all eight. A point that did not touch a function adds fma(0, c, phi) = phi, so every point keeps its own sum
//order and rho is bit-identical to compute_dens, provided each point's functions came in ascending order (they are
//numbered by first primitive, so they usually do); a group where one did not, and a last partial group, multiply per
//point. Microbench on the A72: 1.5x at 300 functions x 91 MOs, 3.6x at 600 x 200. x64 measured no gain, its L2 holds
//the rows.
//Four union rows per sweep (batch_rows<4>) load and store each phi pair once per four FMAs instead of once per FMA:
//1.25-1.52x more on the A72, WFN on sucrose 1.58 -> 1.05 s. The dgemm batch (#else) is 1.73x slower there (1.81 s),
//so ARM keeps this one.
//The groups are neighbours in space, not in index: an atom grid comes shell by shell in Lebedev order, whose
//consecutive points lie all over the sphere, so the union of eight points' functions, and the zero FMAs with it,
//grew far past each point's own list. Morton order over the points' bounding box, as XCW orders its atom grids.
void WFN::compute_dens_batch(const int n, const double* x, const double* y, const double* z, double* rho) const
{
	constexpr int B = 8;
	if (!get_coef_primitive_major()) {
		std::fill(rho, rho + n, 0.0);
		return;
	}
	const int nao = (int)(coef_ao_major.size() / nmo);
	const double* coef = coef_ao_major.data();
	const MO* MOs_data = MOs.data();
	auto occ = [MOs_data](const int mo) { return MOs_data[mo].get_occ(); };
	ivec runs; //[begin, end) of the occupied MOs, as ao_scratch::to_mo finds them
	for (int mo = 0; mo < nmo; mo++)
		if (occ(mo) != 0.0) {
			if (runs.empty() || runs.back() != mo) { runs.push_back(mo); runs.push_back(mo + 1); }
			else runs.back() = mo + 1;
		}
	std::vector<std::pair<uint64_t, int>> ord(n); //Morton key, point
	{
		double lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 };
		for (int p = 0; p < n; p++) {
			const double v[3] = { x[p], y[p], z[p] };
			for (int k = 0; k < 3; k++) lo[k] = std::min(lo[k], v[k]), hi[k] = std::max(hi[k], v[k]);
		}
		auto spread = [](uint64_t v) { //21 bits to every third bit
			v &= 0x1fffff;
			v = (v | v << 32) & 0x1f00000000ffffull;
			v = (v | v << 16) & 0x1f0000ff0000ffull;
			v = (v | v << 8) & 0x100f00f00f00f00full;
			v = (v | v << 4) & 0x10c30c30c30c30c3ull;
			return (v | v << 2) & 0x1249249249249249ull;
		};
		for (int p = 0; p < n; p++) {
			const double v[3] = { x[p], y[p], z[p] };
			uint64_t key = 0;
			for (int k = 0; k < 3; k++) {
				const double s = hi[k] - lo[k] > 1e-12 ? (v[k] - lo[k]) / (hi[k] - lo[k]) : 0.0;
				key |= spread(static_cast<uint64_t>(std::clamp(s, 0.0, 1.0) * 2097151.0)) << k;
			}
			ord[p] = { key, p };
		}
		std::sort(ord.begin(), ord.end());
	}
	const int ng = (n + B - 1) / B;
#pragma omp parallel for schedule(dynamic, 8)
	for (int g = 0; g < ng; g++) {
		thread_local ao_scratch<1> ao[B];
		thread_local vec2 d;
		thread_local vec phi, phi1;
		thread_local std::vector<char> mark;
		thread_local ivec uni;
		if (d.size() < (size_t)ncen) d.resize(ncen, vec(16, 0.0));
		if (mark.size() < (size_t)nao) mark.resize(nao, 0);
		if (phi1.size() < (size_t)nmo) phi1.resize(nmo);
		phi.assign((size_t)B * nmo, 0.0);
		const int p0 = g * B, nb = std::min(B, n - p0);
		bool ascending = nb == B;
		for (int b = 0; b < nb; b++) {
			const int p = ord[p0 + b].second;
			batch_ao = &ao[b];
			compute_dens_cartesian(d3{ x[p], y[p], z[p] }, d, phi1);
			batch_ao = nullptr;
			ascending = ascending && std::is_sorted(ao[b].touched.begin(), ao[b].touched.end());
		}
		if (!ascending)
			for (int b = 0; b < nb; b++) ao[b].to_mo(coef, nmo, phi.data() + (size_t)b * nmo, occ);
		else {
			uni.clear();
			for (int b = 0; b < B; b++)
				for (const int a : ao[b].touched) mark[a] = 1;
			for (int a = 0; a < nao; a++)
				if (mark[a]) { mark[a] = 0; uni.push_back(a); }
			constexpr int R = 4;
			const double* c[R];
			double w[R][B];
			const int nu = (int)uni.size();
			for (int i = 0; i < nu; i += R) {
				const int nr = std::min(R, nu - i);
				for (int r = 0; r < nr; r++) {
					const int a = uni[i + r];
					c[r] = coef + (size_t)a * nmo;
					for (int b = 0; b < B; b++) w[r][b] = ao[b].hit[a] ? ao[b].v[a] : 0.0;
				}
				if (nr == R) batch_rows<R>(c, w, runs, nmo, phi.data());
				else
					for (int r = 0; r < nr; r++) batch_rows<1>(c + r, w + r, runs, nmo, phi.data());
			}
			for (int b = 0; b < B; b++)
				for (const int a : ao[b].touched) ao[b].hit[a] = 0;
		}
		for (int b = 0; b < nb; b++) {
			const double* pb = phi.data() + (size_t)b * nmo;
			double Rho = 0.0;
			for (int mo = 0; mo < nmo; mo++) Rho += MOs_data[mo].get_occ() * pb[mo] * pb[mo];
			rho[ord[p0 + b].second] = Rho;
		}
	}
}
#else
//compute_dens_cartesian at n points, 64 at a time: the contracted functions of a block (the union over its points of
//those the exponent cutoff keeps) go into one column-major table X, and one dgemm X * C into the occupied MOs replaces
//the per-point multiply, which streamed every touched coefficient row once per point (47 % of a sucrose grid run)
void WFN::compute_dens_batch(const int n, const double* x, const double* y, const double* z, double* rho) const
{
	if (n <= 0) return;
	if (!get_coef_primitive_major()) {
		std::fill(rho, rho + n, 0.0);
		return;
	}
	const int nao = (int)(coef_ao_major.size() / nmo);
	ivec occ_mo;
	vec occ;
	for (int mo = 0; mo < nmo; mo++)
		if (MOs[mo].get_occ() != 0.0) occ_mo.push_back(mo), occ.push_back(MOs[mo].get_occ());
	const int nocc = (int)occ_mo.size();
	if (nocc == 0) {
		std::fill(rho, rho + n, 0.0);
		return;
	}
	//The occupied columns of coef_ao_major, [nao * nocc]; a block gathers its rows from here
	vec C((size_t)nao * nocc);
	for (int a = 0; a < nao; a++)
		for (int m = 0; m < nocc; m++) C[(size_t)a * nocc + m] = coef_ao_major[(size_t)a * nmo + occ_mo[m]];
	constexpr int B = 64;
	const int nblk = (n + B - 1) / B;
#pragma omp parallel
	{
		vec d((size_t)ncen * 16), exps(group_exponent.size()), X((size_t)B * nao), Cg, Phi((size_t)B * nocc);
		Cg.reserve((size_t)nao * nocc);
		ivec col(nao, -1), cols;
		const int* group = prim_exp_group.data();
#pragma omp for schedule(dynamic)
		for (int b = 0; b < nblk; b++) {
			const int p0 = b * B, nb = std::min(B, n - p0);
			cols.clear();
			for (int p = 0; p < nb; p++) {
				const d3 Pos{ x[p0 + p], y[p0 + p], z[p0 + p] };
				for (int c = 0; c < ncen; c++) cart_powers(d.data() + 16 * c, Pos, atoms[c]);
				exp_table([&d](const int c) { return d[16 * c + 3]; }, exps.data());
				for (int j = 0; j < nex; j++) {
					double ex = exps[group[j]];
					const int a = prim_ao[j];
					if (ex == 0.0 || a < 0) continue;
					ex *= cart_monomial(types[j], d.data() + 16 * (centers[j] - 1));
					if (col[a] < 0) {
						col[a] = (int)cols.size();
						cols.push_back(a);
						std::fill_n(X.data() + (size_t)B * col[a], B, 0.0);
					}
					X[(size_t)B * col[a] + p] += prim_ao_scale[j] * ex;
				}
			}
			const int ncol = (int)cols.size();
			double* r = rho + p0;
			std::fill(r, r + nb, 0.0);
			if (ncol == 0) continue;
			Cg.resize((size_t)ncol * nocc);
			for (int c = 0; c < ncol; c++) {
				std::copy_n(C.data() + (size_t)cols[c] * nocc, nocc, Cg.data() + (size_t)c * nocc);
				col[cols[c]] = -1;
			}
			//Phi[m * B + p] = sum_c X[c * B + p] Cg[c * nocc + m]
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, nb, nocc, ncol, 1.0, X.data(), B, Cg.data(), nocc, 0.0, Phi.data(), B);
			for (int m = 0; m < nocc; m++) {
				const double* f = Phi.data() + (size_t)m * B;
				for (int p = 0; p < nb; p++) r[p] += occ[m] * f[p] * f[p];
			}
		}
	}
}
#endif

const double WFN::eval_ao(
	std::array<double, 4>& d,
	const std::vector<primitive>& prims,
	const int &m,
	const double& root_d3
) const
{

	// normalize distances for spherical harmonic
	const int type = prims[0].get_type();
	const double scale = 1 / root_d3;
	d[0] *= scale;
	d[1] *= scale;
	d[2] *= scale;
	const double rl = std::pow(root_d3, type);

	double radial = 0;
	const primitive* p = prims.data();
	const primitive* const p_end = p + prims.size();

	for (; p != p_end; p++) {
		radial += p->eval_gaussian_unnormalized(d[3]);
	}
	radial *= rl;

	return radial * constants::spherical_harmonic(type, m, d.data());
	// err_checkf(coef_counter == exp_coefs, "WRONG NUMBER OF COEFFICIENTS! " + std::to_string(coef_counter) + " vs. " + std::to_string(exp_coefs), std::cout);
}

const double WFN::compute_g_cartesian(
	const d3 &Pos,
	vec2 &d,
	vec &phi) const
{
	std::fill(phi.begin(), phi.end(), 0.0);
	double g = 0.0;
	int j;
	double ex, *d_;

	// precalculate some distances and powers of distances for faster computation
	for (j = 0; j < ncen; j++)
	{
		const atom &a = atoms[j];
		d_ = d[j].data();
		d_[0] = Pos[0] - a.get_coordinate(0);
		d_[1] = Pos[1] - a.get_coordinate(1);
		d_[2] = Pos[2] - a.get_coordinate(2);
		d_[4] = d_[0] * d_[0];
		d_[5] = d_[1] * d_[1];
		d_[6] = d_[2] * d_[2];
		d_[3] = d_[4] + d_[5] + d_[6];
		d_[7] = d_[0] * d_[4];
		d_[8] = d_[1] * d_[5];
		d_[9] = d_[2] * d_[6];
		d_[10] = d_[0] * d_[7];
		d_[11] = d_[1] * d_[8];
		d_[12] = d_[2] * d_[9];
		d_[13] = d_[0] * d_[10];
		d_[14] = d_[1] * d_[11];
		d_[15] = d_[2] * d_[12];
	}

	// Pre-cache frequently accessed data
	const int *centers_data = centers.data();
	const int *types_data = types.data();
	const double *exponents_data = exponents.data();
	double *phi_data = phi.data();
	const double exp_cutoff = constants::exp_cutoff;

	for (j = 0; j < nex; j++)
	{
		d_ = d[centers_data[j] - 1].data();
		ex = -exponents_data[j] * d_[3];
		if (ex < exp_cutoff)
		{ // corresponds to cutoff of maximum density contribution of 1E-5
			continue;
		}
		ex = exp(ex);
		switch (types_data[j])
		{
		case 0:  break;
		case 1:  break;// 0, 0, 0,
		case 2:  ex *= d_[0]; break;// 1, 0, 0,
		case 3:  ex *= d_[1]; break;// 0, 1, 0,
		case 4:  ex *= d_[2]; break;// 0, 0, 1,
		case 5:  ex *= d_[4]; break;// 2, 0, 0,
		case 6:  ex *= d_[5]; break;// 0, 2, 0,
		case 7:  ex *= d_[6]; break;// 0, 0, 2,
		case 8:  ex *= d_[0] * d_[1]; break;// 1, 1, 0,
		case 9:  ex *= d_[0] * d_[2]; break;// 1, 0, 1,
		case 10: ex *= d_[1] * d_[2]; break;// 0, 1, 1,
		case 11: ex *= d_[7]; break;// 3, 0, 0,
		case 12: ex *= d_[8]; break;// 0, 3, 0,
		case 13: ex *= d_[9]; break;// 0, 0, 3,
		case 14: ex *= d_[4] * d_[1]; break;// 2, 1, 0,
		case 15: ex *= d_[4] * d_[2]; break;// 2, 0, 1,
		case 16: ex *= d_[5] * d_[2]; break;// 0, 2, 1,
		case 17: ex *= d_[0] * d_[5]; break;// 1, 2, 0,
		case 18: ex *= d_[0] * d_[6]; break;// 1, 0, 2,
		case 19: ex *= d_[1] * d_[6]; break;// 0, 1, 2,
		case 20: ex *= d_[0] * d_[1] * d_[2]; break;// 1, 1, 1,
		case 21: ex *= d_[12]; break;// 0, 0, 4,
		case 22: ex *= d_[1] * d_[9]; break;// 0, 1, 3,
		case 23: ex *= d_[5] * d_[6]; break;// 0, 2, 2,
		case 24: ex *= d_[8] * d_[2]; break;// 0, 3, 1,
		case 25: ex *= d_[11]; break;// 0, 4, 0,
		case 26: ex *= d_[0] * d_[9]; break;// 1, 0, 3,
		case 27: ex *= d_[0] * d_[1] * d_[6]; break;// 1, 1, 2,
		case 28: ex *= d_[0] * d_[5] * d_[2]; break;// 1, 2, 1,
		case 29: ex *= d_[0] * d_[8]; break;// 1, 3, 0,
		case 30: ex *= d_[4] * d_[6]; break;// 2, 0, 2,
		case 31: ex *= d_[4] * d_[1] * d_[2]; break;// 2, 1, 1,
		case 32: ex *= d_[4] * d_[5]; break;// 2, 2, 0,
		case 33: ex *= d_[7] * d_[2]; break;// 3, 0, 1,
		case 34: ex *= d_[7] * d_[1]; break;// 3, 1, 0,
		case 35: ex *= d_[10]; break;// 4, 0, 0,
		case 36: ex *= d_[15]; break;// 0, 0, 5,
		case 37: ex *= d_[1] * d_[12]; break;// 0, 1, 4,
		case 38: ex *= d_[5] * d_[9]; break;// 0, 2, 3,
		case 39: ex *= d_[8] * d_[6]; break;// 0, 3, 2,
		case 40: ex *= d_[11] * d_[2]; break;// 0, 4, 1,
		case 41: ex *= d_[14]; break;// 0, 5, 0,
		case 42: ex *= d_[0] * d_[12]; break;// 1, 0, 4,
		case 43: ex *= d_[0] * d_[1] * d_[9]; break;// 1, 1, 3,
		case 44: ex *= d_[0] * d_[5] * d_[6]; break;// 1, 2, 2,
		case 45: ex *= d_[0] * d_[8] * d_[2]; break;// 1, 3, 1,
		case 46: ex *= d_[0] * d_[11]; break;// 1, 4, 0,
		case 47: ex *= d_[4] * d_[9]; break;// 2, 0, 3,
		case 48: ex *= d_[4] * d_[1] * d_[6]; break;// 2, 1, 2,
		case 49: ex *= d_[4] * d_[5] * d_[2]; break;// 2, 2, 1,
		case 50: ex *= d_[4] * d_[8]; break;// 2, 3, 0,
		case 51: ex *= d_[7] * d_[6]; break;// 3, 0, 2,
		case 52: ex *= d_[7] * d_[1] * d_[2]; break;// 3, 1, 1,
		case 53: ex *= d_[7] * d_[5]; break;// 3, 2, 0,
		case 54: ex *= d_[10] * d_[2]; break;// 4, 0, 1,
		case 55: ex *= d_[10] * d_[1]; break;// 4, 1, 0,
		case 56: ex *= d_[13]; break;// 5, 0, 0,
		case 57: ex *= d_[15] * d_[2]; break;// 0, 0, 6,
		case 58: ex *= d_[1] * d_[15]; break;// 0, 1, 5,
		case 59: ex *= d_[5] * d_[12]; break;// 0, 2, 4,
		case 60: ex *= d_[8] * d_[9]; break;// 0, 3, 3,
		case 61: ex *= d_[11] * d_[6]; break;// 0, 4, 2,
		case 62: ex *= d_[14] * d_[2]; break;// 0, 5, 1,
		case 63: ex *= d_[14] * d_[1]; break;// 0, 6, 0,
		case 64: ex *= d_[0] * d_[15]; break;// 1, 0, 5,
		case 65: ex *= d_[0] * d_[1] * d_[12]; break;// 1, 1, 4,
		case 66: ex *= d_[0] * d_[5] * d_[9]; break;// 1, 2, 3,
		case 67: ex *= d_[0] * d_[8] * d_[6]; break;// 1, 3, 2,
		case 68: ex *= d_[0] * d_[11] * d_[2]; break;// 1, 4, 1,
		case 69: ex *= d_[0] * d_[14]; break;// 1, 5, 0,
		case 70: ex *= d_[4] * d_[12]; break;// 2, 0, 4,
		case 71: ex *= d_[4] * d_[1] * d_[9]; break;// 2, 1, 3,
		case 72: ex *= d_[4] * d_[5] * d_[6]; break;// 2, 2, 2,
		case 73: ex *= d_[4] * d_[8] * d_[2]; break;// 2, 3, 1,
		case 74: ex *= d_[4] * d_[11]; break;// 2, 4, 0,
		case 75: ex *= d_[7] * d_[9]; break;// 3, 0, 3,
		case 76: ex *= d_[7] * d_[1] * d_[6]; break;// 3, 1, 2,
		case 77: ex *= d_[7] * d_[5] * d_[2]; break;// 3, 2, 1,
		case 78: ex *= d_[7] * d_[8]; break;// 3, 3, 0,
		case 79: ex *= d_[10] * d_[6]; break;// 4, 0, 2,
		case 80: ex *= d_[10] * d_[1] * d_[2]; break;// 4, 1, 1,
		case 81: ex *= d_[10] * d_[5]; break;// 4, 2, 0,
		case 82: ex *= d_[13] * d_[2]; break;// 5, 0, 1,
		case 83: ex *= d_[13] * d_[1]; break;// 5, 1, 0,
		case 84: ex *= d_[13] * d_[0]; break;// 6, 0, 0,
		case 85: ex *= d_[15] * d_[6]; break;// 0, 0, 7,
		case 86: ex *= d_[1] * d_[15] * d_[2]; break;// 0, 1, 6,
		case 87: ex *= d_[5] * d_[15]; break;// 0, 2, 5,
		case 88: ex *= d_[8] * d_[12]; break;// 0, 3, 4,
		case 89: ex *= d_[11] * d_[9]; break;// 0, 4, 3,
		case 90: ex *= d_[14] * d_[6]; break;// 0, 5, 2,
		case 91: ex *= d_[14] * d_[1] * d_[2]; break;// 0, 6, 1,
		case 92: ex *= d_[14] * d_[5]; break;// 0, 7, 0,
		case 93: ex *= d_[0] * d_[15] * d_[2]; break;// 1, 0, 6,
		case 94: ex *= d_[0] * d_[1] * d_[15]; break;// 1, 1, 5,
		case 95: ex *= d_[0] * d_[5] * d_[12]; break;// 1, 2, 4,
		case 96: ex *= d_[0] * d_[8] * d_[9]; break;// 1, 3, 3,
		case 97: ex *= d_[0] * d_[11] * d_[6]; break;// 1, 4, 2,
		case 98: ex *= d_[0] * d_[14] * d_[2]; break;// 1, 5, 1,
		case 99: ex *= d_[0] * d_[14] * d_[1]; break;// 1, 6, 0,
		case 100: ex *= d_[4] * d_[15]; break;// 2, 0, 5,
		case 101: ex *= d_[4] * d_[1] * d_[12]; break;// 2, 1, 4,
		case 102: ex *= d_[4] * d_[5] * d_[9]; break;// 2, 2, 3,
		case 103: ex *= d_[4] * d_[8] * d_[6]; break;// 2, 3, 2,
		case 104: ex *= d_[4] * d_[11] * d_[2]; break;// 2, 4, 1,
		case 105: ex *= d_[4] * d_[14]; break;// 2, 5, 0,
		case 106: ex *= d_[7] * d_[12]; break;// 3, 0, 4,
		case 107: ex *= d_[7] * d_[1] * d_[9]; break;// 3, 1, 3,
		case 108: ex *= d_[7] * d_[5] * d_[6]; break;// 3, 2, 2,
		case 109: ex *= d_[7] * d_[8] * d_[2]; break;// 3, 3, 1,
		case 110: ex *= d_[7] * d_[11]; break;// 3, 4, 0,
		case 111: ex *= d_[10] * d_[9]; break;// 4, 0, 3,
		case 112: ex *= d_[10] * d_[1] * d_[6]; break;// 4, 1, 2,
		case 113: ex *= d_[10] * d_[5] * d_[2]; break;// 4, 2, 1,
		case 114: ex *= d_[10] * d_[8]; break;// 4, 3, 0,
		case 115: ex *= d_[13] * d_[6]; break;// 5, 0, 2,
		case 116: ex *= d_[13] * d_[1] * d_[2]; break;// 5, 1, 1,
		case 117: ex *= d_[13] * d_[5]; break;// 5, 2, 0,
		case 118: ex *= d_[13] * d_[0] * d_[2]; break;// 6, 0, 1,
		case 119: ex *= d_[13] * d_[0] * d_[1]; break;// 6, 1, 0,
		case 120: ex *= d_[13] * d_[4]; break;// 7, 0, 0,
		case 121: ex *= d_[15] * d_[9]; break;// 0, 0, 8,
		case 122: ex *= d_[1] * d_[15] * d_[6]; break;// 0, 1, 7,
		case 123: ex *= d_[5] * d_[15] * d_[2]; break;// 0, 2, 6,
		case 124: ex *= d_[8] * d_[15]; break;// 0, 3, 5,
		case 125: ex *= d_[11] * d_[12]; break;// 0, 4, 4,
		case 126: ex *= d_[14] * d_[9]; break;// 0, 5, 3,
		case 127: ex *= d_[14] * d_[1] * d_[6]; break;// 0, 6, 2,
		case 128: ex *= d_[14] * d_[5] * d_[2]; break;// 0, 7, 1,
		case 129: ex *= d_[14] * d_[8]; break;// 0, 8, 0,
		case 130: ex *= d_[0] * d_[15] * d_[6]; break;// 1, 0, 7,
		case 131: ex *= d_[0] * d_[1] * d_[15] * d_[2]; break;// 1, 1, 6,
		case 132: ex *= d_[0] * d_[5] * d_[15]; break;// 1, 2, 5,
		case 133: ex *= d_[0] * d_[8] * d_[12]; break;// 1, 3, 4,
		case 134: ex *= d_[0] * d_[11] * d_[9]; break;// 1, 4, 3,
		case 135: ex *= d_[0] * d_[14] * d_[6]; break;// 1, 5, 2,
		case 136: ex *= d_[0] * d_[14] * d_[1] * d_[2]; break;// 1, 6, 1,
		case 137: ex *= d_[0] * d_[14] * d_[5]; break;// 1, 7, 0,
		case 138: ex *= d_[4] * d_[15] * d_[2]; break;// 2, 0, 6,
		case 139: ex *= d_[4] * d_[1] * d_[15]; break;// 2, 1, 5,
		case 140: ex *= d_[4] * d_[5] * d_[12]; break;// 2, 2, 4,
		case 141: ex *= d_[4] * d_[8] * d_[9]; break;// 2, 3, 3,
		case 142: ex *= d_[4] * d_[11] * d_[6]; break;// 2, 4, 2,
		case 143: ex *= d_[4] * d_[14] * d_[2]; break;// 2, 5, 1,
		case 144: ex *= d_[4] * d_[14] * d_[1]; break;// 2, 6, 0,
		case 145: ex *= d_[7] * d_[15]; break;// 3, 0, 5,
		case 146: ex *= d_[7] * d_[1] * d_[12]; break;// 3, 1, 4,
		case 147: ex *= d_[7] * d_[5] * d_[9]; break;// 3, 2, 3,
		case 148: ex *= d_[7] * d_[8] * d_[6]; break;// 3, 3, 2,
		case 149: ex *= d_[7] * d_[11] * d_[2]; break;// 3, 4, 1,
		case 150: ex *= d_[7] * d_[14]; break;// 3, 5, 0,
		case 151: ex *= d_[10] * d_[12]; break;// 4, 0, 4,
		case 152: ex *= d_[10] * d_[1] * d_[9]; break;// 4, 1, 3,
		case 153: ex *= d_[10] * d_[5] * d_[6]; break;// 4, 2, 2,
		case 154: ex *= d_[10] * d_[8] * d_[2]; break;// 4, 3, 1,
		case 155: ex *= d_[10] * d_[11]; break;// 4, 4, 0,
		case 156: ex *= d_[13] * d_[9]; break;// 5, 0, 3,
		case 157: ex *= d_[13] * d_[1] * d_[6]; break;// 5, 1, 2,
		case 158: ex *= d_[13] * d_[5] * d_[2]; break;// 5, 2, 1,
		case 159: ex *= d_[13] * d_[8]; break;// 5, 3, 0,
		case 160: ex *= d_[13] * d_[0] * d_[6]; break;// 6, 0, 2,
		case 161: ex *= d_[13] * d_[0] * d_[1] * d_[2]; break;// 6, 1, 1,
		case 162: ex *= d_[13] * d_[0] * d_[5]; break;// 6, 2, 0,
		case 163: ex *= d_[13] * d_[4] * d_[2]; break;// 7, 0, 1,
		case 164: ex *= d_[13] * d_[4] * d_[1]; break;// 7, 1, 0,
		case 165: ex *= d_[13] * d_[7]; break;// 8, 0, 0,
		case 166: ex *= d_[15] * d_[12]; break;// 0, 0, 9,
		case 167: ex *= d_[1] * d_[15] * d_[9]; break;// 0, 1, 8,
		case 168: ex *= d_[5] * d_[15] * d_[6]; break;// 0, 2, 7,
		case 169: ex *= d_[8] * d_[15] * d_[2]; break;// 0, 3, 6,
		case 170: ex *= d_[11] * d_[15]; break;// 0, 4, 5,
		case 171: ex *= d_[14] * d_[12]; break;// 0, 5, 4,
		case 172: ex *= d_[14] * d_[1] * d_[9]; break;// 0, 6, 3,
		case 173: ex *= d_[14] * d_[5] * d_[6]; break;// 0, 7, 2,
		case 174: ex *= d_[14] * d_[8] * d_[2]; break;// 0, 8, 1,
		case 175: ex *= d_[14] * d_[11]; break;// 0, 9, 0,
		case 176: ex *= d_[0] * d_[15] * d_[9]; break;// 1, 0, 8,
		case 177: ex *= d_[0] * d_[1] * d_[15] * d_[6]; break;// 1, 1, 7,
		case 178: ex *= d_[0] * d_[5] * d_[15] * d_[2]; break;// 1, 2, 6,
		case 179: ex *= d_[0] * d_[8] * d_[15]; break;// 1, 3, 5,
		case 180: ex *= d_[0] * d_[11] * d_[12]; break;// 1, 4, 4,
		case 181: ex *= d_[0] * d_[14] * d_[9]; break;// 1, 5, 3,
		case 182: ex *= d_[0] * d_[14] * d_[1] * d_[6]; break;// 1, 6, 2,
		case 183: ex *= d_[0] * d_[14] * d_[5] * d_[2]; break;// 1, 7, 1,
		case 184: ex *= d_[0] * d_[14] * d_[8]; break;// 1, 8, 0,
		case 185: ex *= d_[4] * d_[15] * d_[6]; break;// 2, 0, 7,
		case 186: ex *= d_[4] * d_[1] * d_[15] * d_[2]; break;// 2, 1, 6,
		case 187: ex *= d_[4] * d_[5] * d_[15]; break;// 2, 2, 5,
		case 188: ex *= d_[4] * d_[8] * d_[12]; break;// 2, 3, 4,
		case 189: ex *= d_[4] * d_[11] * d_[9]; break;// 2, 4, 3,
		case 190: ex *= d_[4] * d_[14] * d_[6]; break;// 2, 5, 2,
		case 191: ex *= d_[4] * d_[14] * d_[1] * d_[2]; break;// 2, 6, 1,
		case 192: ex *= d_[4] * d_[14] * d_[5]; break;// 2, 7, 0,
		case 193: ex *= d_[7] * d_[15] * d_[2]; break;// 3, 0, 6,
		case 194: ex *= d_[7] * d_[1] * d_[15]; break;// 3, 1, 5,
		case 195: ex *= d_[7] * d_[5] * d_[12]; break;// 3, 2, 4,
		case 196: ex *= d_[7] * d_[8] * d_[9]; break;// 3, 3, 3,
		case 197: ex *= d_[7] * d_[11] * d_[6]; break;// 3, 4, 2,
		case 198: ex *= d_[7] * d_[14] * d_[2]; break;// 3, 5, 1,
		case 199: ex *= d_[7] * d_[14] * d_[1]; break;// 3, 6, 0,
		case 200: ex *= d_[10] * d_[15]; break;// 4, 0, 5,
		case 201: ex *= d_[10] * d_[1] * d_[12]; break;// 4, 1, 4,
		case 202: ex *= d_[10] * d_[5] * d_[9]; break;// 4, 2, 3,
		case 203: ex *= d_[10] * d_[8] * d_[6]; break;// 4, 3, 2,
		case 204: ex *= d_[10] * d_[11] * d_[2]; break;// 4, 4, 1,
		case 205: ex *= d_[10] * d_[14]; break;// 4, 5, 0,
		case 206: ex *= d_[13] * d_[12]; break;// 5, 0, 4,
		case 207: ex *= d_[13] * d_[1] * d_[9]; break;// 5, 1, 3,
		case 208: ex *= d_[13] * d_[5] * d_[6]; break;// 5, 2, 2,
		case 209: ex *= d_[13] * d_[8] * d_[2]; break;// 5, 3, 1,
		case 210: ex *= d_[13] * d_[11]; break;// 5, 4, 0,
		case 211: ex *= d_[13] * d_[0] * d_[9]; break;// 6, 0, 3,
		case 212: ex *= d_[13] * d_[0] * d_[1] * d_[6]; break;// 6, 1, 2,
		case 213: ex *= d_[13] * d_[0] * d_[5] * d_[2]; break;// 6, 2, 1,
		case 214: ex *= d_[13] * d_[0] * d_[8]; break;// 6, 3, 0,
		case 215: ex *= d_[13] * d_[4] * d_[6]; break;// 7, 0, 2,
		case 216: ex *= d_[13] * d_[4] * d_[1] * d_[2]; break;// 7, 1, 1,
		case 217: ex *= d_[13] * d_[4] * d_[5]; break;// 7, 2, 0,
		case 218: ex *= d_[13] * d_[7] * d_[2]; break;// 8, 0, 1,
		case 219: ex *= d_[13] * d_[7] * d_[1]; break;// 8, 1, 0,
		case 220: ex *= d_[13] * d_[10]; break;// 9, 0, 0,
		case 221: ex *= d_[15] * d_[15]; break;// 0, 0, 10,
		case 222: ex *= d_[1] * d_[15] * d_[12]; break;// 0, 1, 9,
		case 223: ex *= d_[5] * d_[15] * d_[9]; break;// 0, 2, 8,
		case 224: ex *= d_[8] * d_[15] * d_[6]; break;// 0, 3, 7,
		case 225: ex *= d_[11] * d_[15] * d_[2]; break;// 0, 4, 6,
		case 226: ex *= d_[14] * d_[15]; break;// 0, 5, 5,
		case 227: ex *= d_[14] * d_[1] * d_[12]; break;// 0, 6, 4,
		case 228: ex *= d_[14] * d_[5] * d_[9]; break;// 0, 7, 3,
		case 229: ex *= d_[14] * d_[8] * d_[6]; break;// 0, 8, 2,
		case 230: ex *= d_[14] * d_[11] * d_[2]; break;// 0, 9, 1,
		case 231: ex *= d_[14] * d_[14]; break;// 0, 10, 0,
		case 232: ex *= d_[0] * d_[15] * d_[12]; break;// 1, 0, 9,
		case 233: ex *= d_[0] * d_[1] * d_[15] * d_[9]; break;// 1, 1, 8,
		case 234: ex *= d_[0] * d_[5] * d_[15] * d_[6]; break;// 1, 2, 7,
		case 235: ex *= d_[0] * d_[8] * d_[15] * d_[2]; break;// 1, 3, 6,
		case 236: ex *= d_[0] * d_[11] * d_[15]; break;// 1, 4, 5,
		case 237: ex *= d_[0] * d_[14] * d_[12]; break;// 1, 5, 4,
		case 238: ex *= d_[0] * d_[14] * d_[1] * d_[9]; break;// 1, 6, 3,
		case 239: ex *= d_[0] * d_[14] * d_[5] * d_[6]; break;// 1, 7, 2,
		case 240: ex *= d_[0] * d_[14] * d_[8] * d_[2]; break;// 1, 8, 1,
		case 241: ex *= d_[0] * d_[14] * d_[11]; break;// 1, 9, 0,
		case 242: ex *= d_[4] * d_[15] * d_[9]; break;// 2, 0, 8,
		case 243: ex *= d_[4] * d_[1] * d_[15] * d_[6]; break;// 2, 1, 7,
		case 244: ex *= d_[4] * d_[5] * d_[15] * d_[2]; break;// 2, 2, 6,
		case 245: ex *= d_[4] * d_[8] * d_[15]; break;// 2, 3, 5,
		case 246: ex *= d_[4] * d_[11] * d_[12]; break;// 2, 4, 4,
		case 247: ex *= d_[4] * d_[14] * d_[9]; break;// 2, 5, 3,
		case 248: ex *= d_[4] * d_[14] * d_[1] * d_[6]; break;// 2, 6, 2,
		case 249: ex *= d_[4] * d_[14] * d_[5] * d_[2]; break;// 2, 7, 1,
		case 250: ex *= d_[4] * d_[14] * d_[8]; break;// 2, 8, 0,
		case 251: ex *= d_[7] * d_[15] * d_[6]; break;// 3, 0, 7,
		case 252: ex *= d_[7] * d_[1] * d_[15] * d_[2]; break;// 3, 1, 6,
		case 253: ex *= d_[7] * d_[5] * d_[15]; break;// 3, 2, 5,
		case 254: ex *= d_[7] * d_[8] * d_[12]; break;// 3, 3, 4,
		case 255: ex *= d_[7] * d_[11] * d_[9]; break;// 3, 4, 3,
		case 256: ex *= d_[7] * d_[14] * d_[6]; break;// 3, 5, 2,
		case 257: ex *= d_[7] * d_[14] * d_[1] * d_[2]; break;// 3, 6, 1,
		case 258: ex *= d_[7] * d_[14] * d_[5]; break;// 3, 7, 0,
		case 259: ex *= d_[10] * d_[15] * d_[2]; break;// 4, 0, 6,
		case 260: ex *= d_[10] * d_[1] * d_[15]; break;// 4, 1, 5,
		case 261: ex *= d_[10] * d_[5] * d_[12]; break;// 4, 2, 4,
		case 262: ex *= d_[10] * d_[8] * d_[9]; break;// 4, 3, 3,
		case 263: ex *= d_[10] * d_[11] * d_[6]; break;// 4, 4, 2,
		case 264: ex *= d_[10] * d_[14] * d_[2]; break;// 4, 5, 1,
		case 265: ex *= d_[10] * d_[14] * d_[1]; break;// 4, 6, 0,
		case 266: ex *= d_[13] * d_[15]; break;// 5, 0, 5,
		case 267: ex *= d_[13] * d_[1] * d_[12]; break;// 5, 1, 4,
		case 268: ex *= d_[13] * d_[5] * d_[9]; break;// 5, 2, 3,
		case 269: ex *= d_[13] * d_[8] * d_[6]; break;// 5, 3, 2,
		case 270: ex *= d_[13] * d_[11] * d_[2]; break;// 5, 4, 1,
		case 271: ex *= d_[13] * d_[14]; break;// 5, 5, 0,
		case 272: ex *= d_[13] * d_[0] * d_[12]; break;// 6, 0, 4,
		case 273: ex *= d_[13] * d_[0] * d_[1] * d_[9]; break;// 6, 1, 3,
		case 274: ex *= d_[13] * d_[0] * d_[5] * d_[6]; break;// 6, 2, 2,
		case 275: ex *= d_[13] * d_[0] * d_[8] * d_[2]; break;// 6, 3, 1,
		case 276: ex *= d_[13] * d_[0] * d_[11]; break;// 6, 4, 0,
		case 277: ex *= d_[13] * d_[4] * d_[9]; break;// 7, 0, 3,
		case 278: ex *= d_[13] * d_[4] * d_[1] * d_[6]; break;// 7, 1, 2,
		case 279: ex *= d_[13] * d_[4] * d_[5] * d_[2]; break;// 7, 2, 1,
		case 280: ex *= d_[13] * d_[4] * d_[8]; break;// 7, 3, 0,
		case 281: ex *= d_[13] * d_[7] * d_[6]; break;// 8, 0, 2,
		case 282: ex *= d_[13] * d_[7] * d_[1] * d_[2]; break;// 8, 1, 1,
		case 283: ex *= d_[13] * d_[7] * d_[5]; break;// 8, 2, 0,
		case 284: ex *= d_[13] * d_[10] * d_[2]; break;// 9, 0, 1,
		case 285: ex *= d_[13] * d_[10] * d_[1]; break;// 9, 1, 0,
		case 286: ex *= d_[13] * d_[13]; break;// 10, 0, 0,
		default: break;
		}

		// use pointer arithmetic and cache coefficient pointer
		// This avoids repeated virtual function calls to get_coefficient_f
		double *phi_ptr = phi_data;
		const double *phi_end = phi_ptr + nmo;

		// Cache the coefficient pointer for this primitive across all MOs
		for (; phi_ptr != phi_end; ++phi_ptr)
		{
			*phi_ptr += ex;
		}
	}

	// use pointer arithmetic and minimize overhead
	const double *phi_ptr = phi_data;
	const double *phi_end = phi_ptr + nmo;

	for (; phi_ptr != phi_end; ++phi_ptr)
	{
		const double &phi_val = *phi_ptr;
		g += phi_val * phi_val;
	}

	return g;
}


const double WFN::compute_spin_dens_cartesian(
	const d3 &Pos,
	vec2 &d,
	vec &phi) const
{
	std::fill(phi.begin(), phi.end(), 0.0);
	double alpha = 0.0, beta = 0.0, ex, *d_;
	int j;

	for (j = 0; j < ncen; j++)
	{
		d_ = d[j].data();
		d_[0] = Pos[0] - atoms[j].get_coordinate(0);
		d_[1] = Pos[1] - atoms[j].get_coordinate(1);
		d_[2] = Pos[2] - atoms[j].get_coordinate(2);
		d_[4] = d_[0] * d_[0];
		d_[5] = d_[1] * d_[1];
		d_[6] = d_[2] * d_[2];
		d_[3] = d_[4] + d_[5] + d_[6];
		d_[7] = d_[0] * d_[4];
		d_[8] = d_[1] * d_[5];
		d_[9] = d_[2] * d_[6];
		d_[10] = d_[0] * d_[7];
		d_[11] = d_[1] * d_[8];
		d_[12] = d_[2] * d_[9];
		d_[13] = d_[0] * d_[10];
		d_[14] = d_[1] * d_[11];
		d_[15] = d_[2] * d_[12];
	}

	const int *centers_data = centers.data();
	const int *types_data = types.data();
	const double *exponents_data = exponents.data();
	const MO *MOs_data = MOs.data();
	double *phi_data = phi.data();
	//Primitive-major coefficients: every MO for one primitive is contiguous here, where
	//MOs keeps them nex apart. Identical arithmetic in identical order - one sequential
	//read per primitive instead of nmo scattered ones, and no per-point heap allocation.
	const double *const coefs = get_coef_primitive_major();

	for (j = 0; j < nex; j++)
	{
		d_ = d[centers_data[j] - 1].data();
		ex = -exponents_data[j] * d_[3];
		if (ex < constants::exp_cutoff)
		{ // corresponds to cutoff of ex ~< 1E-20
			continue;
		}
		ex = exp(ex);
		switch (types_data[j])
		{
		case 0:  break;
		case 1:  break;// 0, 0, 0,
		case 2:  ex *= d_[0]; break;// 1, 0, 0,
		case 3:  ex *= d_[1]; break;// 0, 1, 0,
		case 4:  ex *= d_[2]; break;// 0, 0, 1,
		case 5:  ex *= d_[4]; break;// 2, 0, 0,
		case 6:  ex *= d_[5]; break;// 0, 2, 0,
		case 7:  ex *= d_[6]; break;// 0, 0, 2,
		case 8:  ex *= d_[0] * d_[1]; break;// 1, 1, 0,
		case 9:  ex *= d_[0] * d_[2]; break;// 1, 0, 1,
		case 10: ex *= d_[1] * d_[2]; break;// 0, 1, 1,
		case 11: ex *= d_[7]; break;// 3, 0, 0,
		case 12: ex *= d_[8]; break;// 0, 3, 0,
		case 13: ex *= d_[9]; break;// 0, 0, 3,
		case 14: ex *= d_[4] * d_[1]; break;// 2, 1, 0,
		case 15: ex *= d_[4] * d_[2]; break;// 2, 0, 1,
		case 16: ex *= d_[5] * d_[2]; break;// 0, 2, 1,
		case 17: ex *= d_[0] * d_[5]; break;// 1, 2, 0,
		case 18: ex *= d_[0] * d_[6]; break;// 1, 0, 2,
		case 19: ex *= d_[1] * d_[6]; break;// 0, 1, 2,
		case 20: ex *= d_[0] * d_[1] * d_[2]; break;// 1, 1, 1,
		case 21: ex *= d_[12]; break;// 0, 0, 4,
		case 22: ex *= d_[1] * d_[9]; break;// 0, 1, 3,
		case 23: ex *= d_[5] * d_[6]; break;// 0, 2, 2,
		case 24: ex *= d_[8] * d_[2]; break;// 0, 3, 1,
		case 25: ex *= d_[11]; break;// 0, 4, 0,
		case 26: ex *= d_[0] * d_[9]; break;// 1, 0, 3,
		case 27: ex *= d_[0] * d_[1] * d_[6]; break;// 1, 1, 2,
		case 28: ex *= d_[0] * d_[5] * d_[2]; break;// 1, 2, 1,
		case 29: ex *= d_[0] * d_[8]; break;// 1, 3, 0,
		case 30: ex *= d_[4] * d_[6]; break;// 2, 0, 2,
		case 31: ex *= d_[4] * d_[1] * d_[2]; break;// 2, 1, 1,
		case 32: ex *= d_[4] * d_[5]; break;// 2, 2, 0,
		case 33: ex *= d_[7] * d_[2]; break;// 3, 0, 1,
		case 34: ex *= d_[7] * d_[1]; break;// 3, 1, 0,
		case 35: ex *= d_[10]; break;// 4, 0, 0,
		case 36: ex *= d_[15]; break;// 0, 0, 5,
		case 37: ex *= d_[1] * d_[12]; break;// 0, 1, 4,
		case 38: ex *= d_[5] * d_[9]; break;// 0, 2, 3,
		case 39: ex *= d_[8] * d_[6]; break;// 0, 3, 2,
		case 40: ex *= d_[11] * d_[2]; break;// 0, 4, 1,
		case 41: ex *= d_[14]; break;// 0, 5, 0,
		case 42: ex *= d_[0] * d_[12]; break;// 1, 0, 4,
		case 43: ex *= d_[0] * d_[1] * d_[9]; break;// 1, 1, 3,
		case 44: ex *= d_[0] * d_[5] * d_[6]; break;// 1, 2, 2,
		case 45: ex *= d_[0] * d_[8] * d_[2]; break;// 1, 3, 1,
		case 46: ex *= d_[0] * d_[11]; break;// 1, 4, 0,
		case 47: ex *= d_[4] * d_[9]; break;// 2, 0, 3,
		case 48: ex *= d_[4] * d_[1] * d_[6]; break;// 2, 1, 2,
		case 49: ex *= d_[4] * d_[5] * d_[2]; break;// 2, 2, 1,
		case 50: ex *= d_[4] * d_[8]; break;// 2, 3, 0,
		case 51: ex *= d_[7] * d_[6]; break;// 3, 0, 2,
		case 52: ex *= d_[7] * d_[1] * d_[2]; break;// 3, 1, 1,
		case 53: ex *= d_[7] * d_[5]; break;// 3, 2, 0,
		case 54: ex *= d_[10] * d_[2]; break;// 4, 0, 1,
		case 55: ex *= d_[10] * d_[1]; break;// 4, 1, 0,
		case 56: ex *= d_[13]; break;// 5, 0, 0,
		case 57: ex *= d_[15] * d_[2]; break;// 0, 0, 6,
		case 58: ex *= d_[1] * d_[15]; break;// 0, 1, 5,
		case 59: ex *= d_[5] * d_[12]; break;// 0, 2, 4,
		case 60: ex *= d_[8] * d_[9]; break;// 0, 3, 3,
		case 61: ex *= d_[11] * d_[6]; break;// 0, 4, 2,
		case 62: ex *= d_[14] * d_[2]; break;// 0, 5, 1,
		case 63: ex *= d_[14] * d_[1]; break;// 0, 6, 0,
		case 64: ex *= d_[0] * d_[15]; break;// 1, 0, 5,
		case 65: ex *= d_[0] * d_[1] * d_[12]; break;// 1, 1, 4,
		case 66: ex *= d_[0] * d_[5] * d_[9]; break;// 1, 2, 3,
		case 67: ex *= d_[0] * d_[8] * d_[6]; break;// 1, 3, 2,
		case 68: ex *= d_[0] * d_[11] * d_[2]; break;// 1, 4, 1,
		case 69: ex *= d_[0] * d_[14]; break;// 1, 5, 0,
		case 70: ex *= d_[4] * d_[12]; break;// 2, 0, 4,
		case 71: ex *= d_[4] * d_[1] * d_[9]; break;// 2, 1, 3,
		case 72: ex *= d_[4] * d_[5] * d_[6]; break;// 2, 2, 2,
		case 73: ex *= d_[4] * d_[8] * d_[2]; break;// 2, 3, 1,
		case 74: ex *= d_[4] * d_[11]; break;// 2, 4, 0,
		case 75: ex *= d_[7] * d_[9]; break;// 3, 0, 3,
		case 76: ex *= d_[7] * d_[1] * d_[6]; break;// 3, 1, 2,
		case 77: ex *= d_[7] * d_[5] * d_[2]; break;// 3, 2, 1,
		case 78: ex *= d_[7] * d_[8]; break;// 3, 3, 0,
		case 79: ex *= d_[10] * d_[6]; break;// 4, 0, 2,
		case 80: ex *= d_[10] * d_[1] * d_[2]; break;// 4, 1, 1,
		case 81: ex *= d_[10] * d_[5]; break;// 4, 2, 0,
		case 82: ex *= d_[13] * d_[2]; break;// 5, 0, 1,
		case 83: ex *= d_[13] * d_[1]; break;// 5, 1, 0,
		case 84: ex *= d_[13] * d_[0]; break;// 6, 0, 0,
		case 85: ex *= d_[15] * d_[6]; break;// 0, 0, 7,
		case 86: ex *= d_[1] * d_[15] * d_[2]; break;// 0, 1, 6,
		case 87: ex *= d_[5] * d_[15]; break;// 0, 2, 5,
		case 88: ex *= d_[8] * d_[12]; break;// 0, 3, 4,
		case 89: ex *= d_[11] * d_[9]; break;// 0, 4, 3,
		case 90: ex *= d_[14] * d_[6]; break;// 0, 5, 2,
		case 91: ex *= d_[14] * d_[1] * d_[2]; break;// 0, 6, 1,
		case 92: ex *= d_[14] * d_[5]; break;// 0, 7, 0,
		case 93: ex *= d_[0] * d_[15] * d_[2]; break;// 1, 0, 6,
		case 94: ex *= d_[0] * d_[1] * d_[15]; break;// 1, 1, 5,
		case 95: ex *= d_[0] * d_[5] * d_[12]; break;// 1, 2, 4,
		case 96: ex *= d_[0] * d_[8] * d_[9]; break;// 1, 3, 3,
		case 97: ex *= d_[0] * d_[11] * d_[6]; break;// 1, 4, 2,
		case 98: ex *= d_[0] * d_[14] * d_[2]; break;// 1, 5, 1,
		case 99: ex *= d_[0] * d_[14] * d_[1]; break;// 1, 6, 0,
		case 100: ex *= d_[4] * d_[15]; break;// 2, 0, 5,
		case 101: ex *= d_[4] * d_[1] * d_[12]; break;// 2, 1, 4,
		case 102: ex *= d_[4] * d_[5] * d_[9]; break;// 2, 2, 3,
		case 103: ex *= d_[4] * d_[8] * d_[6]; break;// 2, 3, 2,
		case 104: ex *= d_[4] * d_[11] * d_[2]; break;// 2, 4, 1,
		case 105: ex *= d_[4] * d_[14]; break;// 2, 5, 0,
		case 106: ex *= d_[7] * d_[12]; break;// 3, 0, 4,
		case 107: ex *= d_[7] * d_[1] * d_[9]; break;// 3, 1, 3,
		case 108: ex *= d_[7] * d_[5] * d_[6]; break;// 3, 2, 2,
		case 109: ex *= d_[7] * d_[8] * d_[2]; break;// 3, 3, 1,
		case 110: ex *= d_[7] * d_[11]; break;// 3, 4, 0,
		case 111: ex *= d_[10] * d_[9]; break;// 4, 0, 3,
		case 112: ex *= d_[10] * d_[1] * d_[6]; break;// 4, 1, 2,
		case 113: ex *= d_[10] * d_[5] * d_[2]; break;// 4, 2, 1,
		case 114: ex *= d_[10] * d_[8]; break;// 4, 3, 0,
		case 115: ex *= d_[13] * d_[6]; break;// 5, 0, 2,
		case 116: ex *= d_[13] * d_[1] * d_[2]; break;// 5, 1, 1,
		case 117: ex *= d_[13] * d_[5]; break;// 5, 2, 0,
		case 118: ex *= d_[13] * d_[0] * d_[2]; break;// 6, 0, 1,
		case 119: ex *= d_[13] * d_[0] * d_[1]; break;// 6, 1, 0,
		case 120: ex *= d_[13] * d_[4]; break;// 7, 0, 0,
		case 121: ex *= d_[15] * d_[9]; break;// 0, 0, 8,
		case 122: ex *= d_[1] * d_[15] * d_[6]; break;// 0, 1, 7,
		case 123: ex *= d_[5] * d_[15] * d_[2]; break;// 0, 2, 6,
		case 124: ex *= d_[8] * d_[15]; break;// 0, 3, 5,
		case 125: ex *= d_[11] * d_[12]; break;// 0, 4, 4,
		case 126: ex *= d_[14] * d_[9]; break;// 0, 5, 3,
		case 127: ex *= d_[14] * d_[1] * d_[6]; break;// 0, 6, 2,
		case 128: ex *= d_[14] * d_[5] * d_[2]; break;// 0, 7, 1,
		case 129: ex *= d_[14] * d_[8]; break;// 0, 8, 0,
		case 130: ex *= d_[0] * d_[15] * d_[6]; break;// 1, 0, 7,
		case 131: ex *= d_[0] * d_[1] * d_[15] * d_[2]; break;// 1, 1, 6,
		case 132: ex *= d_[0] * d_[5] * d_[15]; break;// 1, 2, 5,
		case 133: ex *= d_[0] * d_[8] * d_[12]; break;// 1, 3, 4,
		case 134: ex *= d_[0] * d_[11] * d_[9]; break;// 1, 4, 3,
		case 135: ex *= d_[0] * d_[14] * d_[6]; break;// 1, 5, 2,
		case 136: ex *= d_[0] * d_[14] * d_[1] * d_[2]; break;// 1, 6, 1,
		case 137: ex *= d_[0] * d_[14] * d_[5]; break;// 1, 7, 0,
		case 138: ex *= d_[4] * d_[15] * d_[2]; break;// 2, 0, 6,
		case 139: ex *= d_[4] * d_[1] * d_[15]; break;// 2, 1, 5,
		case 140: ex *= d_[4] * d_[5] * d_[12]; break;// 2, 2, 4,
		case 141: ex *= d_[4] * d_[8] * d_[9]; break;// 2, 3, 3,
		case 142: ex *= d_[4] * d_[11] * d_[6]; break;// 2, 4, 2,
		case 143: ex *= d_[4] * d_[14] * d_[2]; break;// 2, 5, 1,
		case 144: ex *= d_[4] * d_[14] * d_[1]; break;// 2, 6, 0,
		case 145: ex *= d_[7] * d_[15]; break;// 3, 0, 5,
		case 146: ex *= d_[7] * d_[1] * d_[12]; break;// 3, 1, 4,
		case 147: ex *= d_[7] * d_[5] * d_[9]; break;// 3, 2, 3,
		case 148: ex *= d_[7] * d_[8] * d_[6]; break;// 3, 3, 2,
		case 149: ex *= d_[7] * d_[11] * d_[2]; break;// 3, 4, 1,
		case 150: ex *= d_[7] * d_[14]; break;// 3, 5, 0,
		case 151: ex *= d_[10] * d_[12]; break;// 4, 0, 4,
		case 152: ex *= d_[10] * d_[1] * d_[9]; break;// 4, 1, 3,
		case 153: ex *= d_[10] * d_[5] * d_[6]; break;// 4, 2, 2,
		case 154: ex *= d_[10] * d_[8] * d_[2]; break;// 4, 3, 1,
		case 155: ex *= d_[10] * d_[11]; break;// 4, 4, 0,
		case 156: ex *= d_[13] * d_[9]; break;// 5, 0, 3,
		case 157: ex *= d_[13] * d_[1] * d_[6]; break;// 5, 1, 2,
		case 158: ex *= d_[13] * d_[5] * d_[2]; break;// 5, 2, 1,
		case 159: ex *= d_[13] * d_[8]; break;// 5, 3, 0,
		case 160: ex *= d_[13] * d_[0] * d_[6]; break;// 6, 0, 2,
		case 161: ex *= d_[13] * d_[0] * d_[1] * d_[2]; break;// 6, 1, 1,
		case 162: ex *= d_[13] * d_[0] * d_[5]; break;// 6, 2, 0,
		case 163: ex *= d_[13] * d_[4] * d_[2]; break;// 7, 0, 1,
		case 164: ex *= d_[13] * d_[4] * d_[1]; break;// 7, 1, 0,
		case 165: ex *= d_[13] * d_[7]; break;// 8, 0, 0,
		case 166: ex *= d_[15] * d_[12]; break;// 0, 0, 9,
		case 167: ex *= d_[1] * d_[15] * d_[9]; break;// 0, 1, 8,
		case 168: ex *= d_[5] * d_[15] * d_[6]; break;// 0, 2, 7,
		case 169: ex *= d_[8] * d_[15] * d_[2]; break;// 0, 3, 6,
		case 170: ex *= d_[11] * d_[15]; break;// 0, 4, 5,
		case 171: ex *= d_[14] * d_[12]; break;// 0, 5, 4,
		case 172: ex *= d_[14] * d_[1] * d_[9]; break;// 0, 6, 3,
		case 173: ex *= d_[14] * d_[5] * d_[6]; break;// 0, 7, 2,
		case 174: ex *= d_[14] * d_[8] * d_[2]; break;// 0, 8, 1,
		case 175: ex *= d_[14] * d_[11]; break;// 0, 9, 0,
		case 176: ex *= d_[0] * d_[15] * d_[9]; break;// 1, 0, 8,
		case 177: ex *= d_[0] * d_[1] * d_[15] * d_[6]; break;// 1, 1, 7,
		case 178: ex *= d_[0] * d_[5] * d_[15] * d_[2]; break;// 1, 2, 6,
		case 179: ex *= d_[0] * d_[8] * d_[15]; break;// 1, 3, 5,
		case 180: ex *= d_[0] * d_[11] * d_[12]; break;// 1, 4, 4,
		case 181: ex *= d_[0] * d_[14] * d_[9]; break;// 1, 5, 3,
		case 182: ex *= d_[0] * d_[14] * d_[1] * d_[6]; break;// 1, 6, 2,
		case 183: ex *= d_[0] * d_[14] * d_[5] * d_[2]; break;// 1, 7, 1,
		case 184: ex *= d_[0] * d_[14] * d_[8]; break;// 1, 8, 0,
		case 185: ex *= d_[4] * d_[15] * d_[6]; break;// 2, 0, 7,
		case 186: ex *= d_[4] * d_[1] * d_[15] * d_[2]; break;// 2, 1, 6,
		case 187: ex *= d_[4] * d_[5] * d_[15]; break;// 2, 2, 5,
		case 188: ex *= d_[4] * d_[8] * d_[12]; break;// 2, 3, 4,
		case 189: ex *= d_[4] * d_[11] * d_[9]; break;// 2, 4, 3,
		case 190: ex *= d_[4] * d_[14] * d_[6]; break;// 2, 5, 2,
		case 191: ex *= d_[4] * d_[14] * d_[1] * d_[2]; break;// 2, 6, 1,
		case 192: ex *= d_[4] * d_[14] * d_[5]; break;// 2, 7, 0,
		case 193: ex *= d_[7] * d_[15] * d_[2]; break;// 3, 0, 6,
		case 194: ex *= d_[7] * d_[1] * d_[15]; break;// 3, 1, 5,
		case 195: ex *= d_[7] * d_[5] * d_[12]; break;// 3, 2, 4,
		case 196: ex *= d_[7] * d_[8] * d_[9]; break;// 3, 3, 3,
		case 197: ex *= d_[7] * d_[11] * d_[6]; break;// 3, 4, 2,
		case 198: ex *= d_[7] * d_[14] * d_[2]; break;// 3, 5, 1,
		case 199: ex *= d_[7] * d_[14] * d_[1]; break;// 3, 6, 0,
		case 200: ex *= d_[10] * d_[15]; break;// 4, 0, 5,
		case 201: ex *= d_[10] * d_[1] * d_[12]; break;// 4, 1, 4,
		case 202: ex *= d_[10] * d_[5] * d_[9]; break;// 4, 2, 3,
		case 203: ex *= d_[10] * d_[8] * d_[6]; break;// 4, 3, 2,
		case 204: ex *= d_[10] * d_[11] * d_[2]; break;// 4, 4, 1,
		case 205: ex *= d_[10] * d_[14]; break;// 4, 5, 0,
		case 206: ex *= d_[13] * d_[12]; break;// 5, 0, 4,
		case 207: ex *= d_[13] * d_[1] * d_[9]; break;// 5, 1, 3,
		case 208: ex *= d_[13] * d_[5] * d_[6]; break;// 5, 2, 2,
		case 209: ex *= d_[13] * d_[8] * d_[2]; break;// 5, 3, 1,
		case 210: ex *= d_[13] * d_[11]; break;// 5, 4, 0,
		case 211: ex *= d_[13] * d_[0] * d_[9]; break;// 6, 0, 3,
		case 212: ex *= d_[13] * d_[0] * d_[1] * d_[6]; break;// 6, 1, 2,
		case 213: ex *= d_[13] * d_[0] * d_[5] * d_[2]; break;// 6, 2, 1,
		case 214: ex *= d_[13] * d_[0] * d_[8]; break;// 6, 3, 0,
		case 215: ex *= d_[13] * d_[4] * d_[6]; break;// 7, 0, 2,
		case 216: ex *= d_[13] * d_[4] * d_[1] * d_[2]; break;// 7, 1, 1,
		case 217: ex *= d_[13] * d_[4] * d_[5]; break;// 7, 2, 0,
		case 218: ex *= d_[13] * d_[7] * d_[2]; break;// 8, 0, 1,
		case 219: ex *= d_[13] * d_[7] * d_[1]; break;// 8, 1, 0,
		case 220: ex *= d_[13] * d_[10]; break;// 9, 0, 0,
		case 221: ex *= d_[15] * d_[15]; break;// 0, 0, 10,
		case 222: ex *= d_[1] * d_[15] * d_[12]; break;// 0, 1, 9,
		case 223: ex *= d_[5] * d_[15] * d_[9]; break;// 0, 2, 8,
		case 224: ex *= d_[8] * d_[15] * d_[6]; break;// 0, 3, 7,
		case 225: ex *= d_[11] * d_[15] * d_[2]; break;// 0, 4, 6,
		case 226: ex *= d_[14] * d_[15]; break;// 0, 5, 5,
		case 227: ex *= d_[14] * d_[1] * d_[12]; break;// 0, 6, 4,
		case 228: ex *= d_[14] * d_[5] * d_[9]; break;// 0, 7, 3,
		case 229: ex *= d_[14] * d_[8] * d_[6]; break;// 0, 8, 2,
		case 230: ex *= d_[14] * d_[11] * d_[2]; break;// 0, 9, 1,
		case 231: ex *= d_[14] * d_[14]; break;// 0, 10, 0,
		case 232: ex *= d_[0] * d_[15] * d_[12]; break;// 1, 0, 9,
		case 233: ex *= d_[0] * d_[1] * d_[15] * d_[9]; break;// 1, 1, 8,
		case 234: ex *= d_[0] * d_[5] * d_[15] * d_[6]; break;// 1, 2, 7,
		case 235: ex *= d_[0] * d_[8] * d_[15] * d_[2]; break;// 1, 3, 6,
		case 236: ex *= d_[0] * d_[11] * d_[15]; break;// 1, 4, 5,
		case 237: ex *= d_[0] * d_[14] * d_[12]; break;// 1, 5, 4,
		case 238: ex *= d_[0] * d_[14] * d_[1] * d_[9]; break;// 1, 6, 3,
		case 239: ex *= d_[0] * d_[14] * d_[5] * d_[6]; break;// 1, 7, 2,
		case 240: ex *= d_[0] * d_[14] * d_[8] * d_[2]; break;// 1, 8, 1,
		case 241: ex *= d_[0] * d_[14] * d_[11]; break;// 1, 9, 0,
		case 242: ex *= d_[4] * d_[15] * d_[9]; break;// 2, 0, 8,
		case 243: ex *= d_[4] * d_[1] * d_[15] * d_[6]; break;// 2, 1, 7,
		case 244: ex *= d_[4] * d_[5] * d_[15] * d_[2]; break;// 2, 2, 6,
		case 245: ex *= d_[4] * d_[8] * d_[15]; break;// 2, 3, 5,
		case 246: ex *= d_[4] * d_[11] * d_[12]; break;// 2, 4, 4,
		case 247: ex *= d_[4] * d_[14] * d_[9]; break;// 2, 5, 3,
		case 248: ex *= d_[4] * d_[14] * d_[1] * d_[6]; break;// 2, 6, 2,
		case 249: ex *= d_[4] * d_[14] * d_[5] * d_[2]; break;// 2, 7, 1,
		case 250: ex *= d_[4] * d_[14] * d_[8]; break;// 2, 8, 0,
		case 251: ex *= d_[7] * d_[15] * d_[6]; break;// 3, 0, 7,
		case 252: ex *= d_[7] * d_[1] * d_[15] * d_[2]; break;// 3, 1, 6,
		case 253: ex *= d_[7] * d_[5] * d_[15]; break;// 3, 2, 5,
		case 254: ex *= d_[7] * d_[8] * d_[12]; break;// 3, 3, 4,
		case 255: ex *= d_[7] * d_[11] * d_[9]; break;// 3, 4, 3,
		case 256: ex *= d_[7] * d_[14] * d_[6]; break;// 3, 5, 2,
		case 257: ex *= d_[7] * d_[14] * d_[1] * d_[2]; break;// 3, 6, 1,
		case 258: ex *= d_[7] * d_[14] * d_[5]; break;// 3, 7, 0,
		case 259: ex *= d_[10] * d_[15] * d_[2]; break;// 4, 0, 6,
		case 260: ex *= d_[10] * d_[1] * d_[15]; break;// 4, 1, 5,
		case 261: ex *= d_[10] * d_[5] * d_[12]; break;// 4, 2, 4,
		case 262: ex *= d_[10] * d_[8] * d_[9]; break;// 4, 3, 3,
		case 263: ex *= d_[10] * d_[11] * d_[6]; break;// 4, 4, 2,
		case 264: ex *= d_[10] * d_[14] * d_[2]; break;// 4, 5, 1,
		case 265: ex *= d_[10] * d_[14] * d_[1]; break;// 4, 6, 0,
		case 266: ex *= d_[13] * d_[15]; break;// 5, 0, 5,
		case 267: ex *= d_[13] * d_[1] * d_[12]; break;// 5, 1, 4,
		case 268: ex *= d_[13] * d_[5] * d_[9]; break;// 5, 2, 3,
		case 269: ex *= d_[13] * d_[8] * d_[6]; break;// 5, 3, 2,
		case 270: ex *= d_[13] * d_[11] * d_[2]; break;// 5, 4, 1,
		case 271: ex *= d_[13] * d_[14]; break;// 5, 5, 0,
		case 272: ex *= d_[13] * d_[0] * d_[12]; break;// 6, 0, 4,
		case 273: ex *= d_[13] * d_[0] * d_[1] * d_[9]; break;// 6, 1, 3,
		case 274: ex *= d_[13] * d_[0] * d_[5] * d_[6]; break;// 6, 2, 2,
		case 275: ex *= d_[13] * d_[0] * d_[8] * d_[2]; break;// 6, 3, 1,
		case 276: ex *= d_[13] * d_[0] * d_[11]; break;// 6, 4, 0,
		case 277: ex *= d_[13] * d_[4] * d_[9]; break;// 7, 0, 3,
		case 278: ex *= d_[13] * d_[4] * d_[1] * d_[6]; break;// 7, 1, 2,
		case 279: ex *= d_[13] * d_[4] * d_[5] * d_[2]; break;// 7, 2, 1,
		case 280: ex *= d_[13] * d_[4] * d_[8]; break;// 7, 3, 0,
		case 281: ex *= d_[13] * d_[7] * d_[6]; break;// 8, 0, 2,
		case 282: ex *= d_[13] * d_[7] * d_[1] * d_[2]; break;// 8, 1, 1,
		case 283: ex *= d_[13] * d_[7] * d_[5]; break;// 8, 2, 0,
		case 284: ex *= d_[13] * d_[10] * d_[2]; break;// 9, 0, 1,
		case 285: ex *= d_[13] * d_[10] * d_[1]; break;// 9, 1, 0,
		case 286: ex *= d_[13] * d_[13]; break;// 10, 0, 0,
		default: break;
		}
		// use pointer arithmetic and cache coefficient pointer
		// This avoids repeated virtual function calls to get_coefficient_f
		const double *c_row = coefs + (size_t)j * nmo;
		double *phi_ptr = phi_data;
		for (int mo = 0; mo < nmo; ++mo, ++phi_ptr)
		{
			*phi_ptr += c_row[mo] * ex;
		}
	}

	// use pointer arithmetic and minimize overhead
	const double *phi_ptr = phi_data;
	const double *phi_end = phi_ptr + nmo;
	const MO *mo_ptr = MOs_data;

	for (; phi_ptr != phi_end; ++phi_ptr, ++mo_ptr)
	{
		const double &phi_val = *phi_ptr;
		if (mo_ptr->get_op())
			beta += mo_ptr->get_occ() * phi_val * phi_val;
		else
			alpha += mo_ptr->get_occ() * phi_val * phi_val;
	}

	return alpha - beta;
}

const double WFN::compute_MO_spherical(
	const d3 &Pos,
	const int &MO) const
{
	err_not_impl_f("This one is not tested an will most likely not work, therefore aborting!", std::cout);
	return 0.0;
	(void)MO;
	//err_checkf(d_f_switch, "Only works for spheriacl wavefunctions!", std::cout);
	//int iat;
	//int l = 0;
	//// ex will carry information about radial function
	//double ex;
	//vec2 d(5);
	//for (int i = 0; i < 5; i++)
	//    d[i].resize(ncen);
	//double phi(0.0);

	//for (iat = 0; iat < ncen; iat++)
	//{
	//    d[0][iat] = Pos1 - atoms[iat].x;
	//    d[1][iat] = Pos2 - atoms[iat].y;
	//    d[2][iat] = Pos3 - atoms[iat].z;
	//    d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
	//    d[4][iat] = sqrt(d[3][iat]);
	//}
	///*Here d[0] = x
	//             d[1] = y
	//             d[2] = z
	//             d[3] = r^2
	//             d[4] = r
	//             */

	//for (int j = 0; j < nex; j++)
	//{
	//    iat = centers[j] - 1;
	//    ex = -exponents[j] * d[3][iat];
	//    if (ex < -46.0517)
	//    { // corresponds to cutoff of ex ~< 1E-20
	//        continue;
	//    }
	//    // apply radial function:
	//    ex = exp(ex);
	//    int lam = 0, m = 0;
	//    if (l > 1)
	//    {
	//        if (l <= 4)
	//        {
	//            ex *= d[4][iat];
	//            lam = 1;
	//            if (l == 2)
	//                m = 0;
	//            else if (l == 3)
	//                m = 1;
	//            else if (l == 4)
	//                m = -1;
	//        }
	//        else if (l <= 9)
	//        {
	//            ex *= d[3][iat];
	//            lam = 2;
	//            if (l == 5)
	//                m = 0;
	//            else if (l == 6)
	//                m = 1;
	//            else if (l == 7)
	//                m = -1;
	//            else if (l == 8)
	//                m = 2;
	//            else if (l == 9)
	//                m = -2;
	//        }
	//        else if (l <= 16)
	//        {
	//            ex *= d[3][iat] * d[4][iat];
	//            lam = 3;
	//            if (l == 10)
	//                m = 0;
	//            else if (l == 11)
	//                m = 1;
	//            else if (l == 12)
	//                m = -1;
	//            else if (l == 13)
	//                m = 2;
	//            else if (l == 14)
	//                m = -2;
	//            else if (l == 15)
	//                m = 3;
	//            else if (l == 16)
	//                m = -3;
	//        }
	//        else if (l <= 25)
	//        {
	//            ex *= d[3][iat] * d[3][iat];
	//            lam = 4;
	//            if (l == 17)
	//                m = 0;
	//            else if (l == 18)
	//                m = 1;
	//            else if (l == 19)
	//                m = -1;
	//            else if (l == 20)
	//                m = 2;
	//            else if (l == 21)
	//                m = -2;
	//            else if (l == 22)
	//                m = 3;
	//            else if (l == 23)
	//                m = -3;
	//            else if (l == 24)
	//                m = 4;
	//            else if (l == 25)
	//                m = -4;
	//        }
	//        else if (l <= 36)
	//        {
	//            ex *= pow(d[4][iat], 5);
	//            lam = 5;
	//            if (l == 26)
	//                m = 0;
	//            else if (l == 27)
	//                m = 1;
	//            else if (l == 28)
	//                m = -1;
	//            else if (l == 29)
	//                m = 2;
	//            else if (l == 30)
	//                m = -2;
	//            else if (l == 31)
	//                m = 3;
	//            else if (l == 32)
	//                m = -3;
	//            else if (l == 33)
	//                m = 4;
	//            else if (l == 34)
	//                m = -4;
	//            else if (l == 35)
	//                m = 5;
	//            else if (l == 36)
	//                m = -5;
	//        }
	//    }
	//    // calc spherical harmonic
	//    double d_t[]{d[0][iat], d[1][iat], d[2][iat], d[3][iat], d[4][iat]};
	//    double SH = spherical_harmonic(lam, m, d_t);
	//    SH *= ex;                                 // multiply radial part with spherical harmonic
	//    phi += MOs[MO].get_coefficient_f(j) * SH; // build MO values at this point
	//}
	//shrink_vector<vec>(d);

	//return phi;
}

//const double WFN::compute_dens_spherical(
//	const d3 &Pos,
//	vec2 &d,
//	vec &phi) const
//{
	//err_checkf(d_f_switch, "Only works for spheriacl wavefunctions!", std::cout);
	//std::fill(phi.begin(), phi.end(), 0.0);
	//double Rho = 0.0;
	//int iat;
	//int l;
	//// ex will carry information about radial function
	//double ex;
	//int mo;

	//for (iat = 0; iat < ncen; iat++)
	//{
	//    d[0][iat] = Pos1 - atoms[iat].x;
	//    d[1][iat] = Pos2 - atoms[iat].y;
	//    d[2][iat] = Pos3 - atoms[iat].z;
	//    d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
	//    d[4][iat] = sqrt(d[3][iat]);
	//}
	///*Here d[0] = x
	//             d[1] = y
	//             d[2] = z
	//             d[3] = r^2
	//             d[4] = r
	//             */
	//for (int j = 0; j < nex; j++)
	//{
	//    iat = centers[j] - 1;
	//    ex = -exponents[j] * d[3][iat];
	//    if (ex < -46.0517)
	//    { // corresponds to cutoff of ex ~< 1E-20
	//        continue;
	//    }
	//    ex = exp(ex);
	//    // apply radial function:
	//    l = types[j];
	//    int lam = 0, m = 0;
	//    if (l > 1)
	//    {
	//        if (l <= 4)
	//        {
	//            ex *= d[4][iat];
	//            lam = 1;
	//            if (l == 2)
	//                m = 0;
	//            else if (l == 3)
	//                m = 1;
	//            else if (l == 4)
	//                m = -1;
	//        }
	//        else if (l <= 9)
	//        {
	//            ex *= d[3][iat];
	//            lam = 2;
	//            if (l == 5)
	//                m = 0;
	//            else if (l == 6)
	//                m = 1;
	//            else if (l == 7)
	//                m = -1;
	//            else if (l == 8)
	//                m = 2;
	//            else if (l == 9)
	//                m = -2;
	//        }
	//        else if (l <= 16)
	//        {
	//            ex *= d[3][iat] * d[4][iat];
	//            lam = 3;
	//            if (l == 10)
	//                m = 0;
	//            else if (l == 11)
	//                m = 1;
	//            else if (l == 12)
	//                m = -1;
	//            else if (l == 13)
	//                m = 2;
	//            else if (l == 14)
	//                m = -2;
	//            else if (l == 15)
	//                m = 3;
	//            else if (l == 16)
	//                m = -3;
	//        }
	//        else if (l <= 25)
	//        {
	//            ex *= d[3][iat] * d[3][iat];
	//            lam = 4;
	//            if (l == 17)
	//                m = 0;
	//            else if (l == 18)
	//                m = 1;
	//            else if (l == 19)
	//                m = -1;
	//            else if (l == 20)
	//                m = 2;
	//            else if (l == 21)
	//                m = -2;
	//            else if (l == 22)
	//                m = 3;
	//            else if (l == 23)
	//                m = -3;
	//            else if (l == 24)
	//                m = 4;
	//            else if (l == 25)
	//                m = -4;
	//        }
	//        else if (l <= 36)
	//        {
	//            ex *= pow(d[4][iat], 5);
	//            lam = 5;
	//            if (l == 26)
	//                m = 0;
	//            else if (l == 27)
	//                m = 1;
	//            else if (l == 28)
	//                m = -1;
	//            else if (l == 29)
	//                m = 2;
	//            else if (l == 30)
	//                m = -2;
	//            else if (l == 31)
	//                m = 3;
	//            else if (l == 32)
	//                m = -3;
	//            else if (l == 33)
	//                m = 4;
	//            else if (l == 34)
	//                m = -4;
	//            else if (l == 35)
	//                m = 5;
	//            else if (l == 36)
	//                m = -5;
	//        }
	//    }
	//    double d_t[]{d[0][iat], d[1][iat], d[2][iat], d[3][iat], d[4][iat]};
	//    double SH = spherical_harmonic(lam, m, d_t);
	//    SH *= ex; // multiply radial part with spherical harmonic
	//    auto run = phi.data();
	//    auto run2 = MOs.data();
	//    for (mo = 0; mo < phi.size(); mo++)
	//    {
	//        *run += (*run2).get_coefficient_f(j) * SH; // build MO values at this point
	//        run++, run2++;
	//    }
	//}

	//auto run = phi.data();
	//auto run2 = MOs.data();
	//for (mo = 0; mo < phi.size(); mo++)
	//{
	//    Rho += (*run2).get_occ() * pow(*run, 2);
	//    run++, run2++;
	//}

	//return Rho;
//}

//Transposed MO coefficients, [primitive * nmo + mo], built once and reused by every point.
//The first caller is usually several threads at once: make_chi parallelises over atom
//pairs, and compute_dens underneath it is what asks for this. Unlocked, two threads both
//find it cold and one reallocates under the other. The valid flag is read through an
//atomic_ref so the warm path takes no lock (it is called per grid point from every
//density evaluator); the member stays a plain bool so WFN stays copyable.
#include <atomic>
#if defined(__cpp_lib_atomic_ref) || defined(_MSC_VER)
using nos_bool_ref = std::atomic_ref<bool>;
#else
//libc++ before 19 (Android NDK r27) has no std::atomic_ref; the compiler builtins are the same thing
struct nos_bool_ref {
	bool& b;
	explicit nos_bool_ref(bool& x) : b(x) {}
	bool load(std::memory_order o) const { return __atomic_load_n(&b, static_cast<int>(o)); }
	void store(bool v, std::memory_order o) { __atomic_store_n(&b, v, static_cast<int>(o)); }
};
#endif
const double* WFN::get_coef_primitive_major() const
{
	const int _nmo = get_nmo(false);
	if (_nmo <= 0 || nex <= 0) return nullptr;
	nos_bool_ref valid(coef_primitive_major_valid);
	if (valid.load(std::memory_order_acquire)
		&& coef_primitive_major.size() == (size_t)nex * (size_t)_nmo)
		return coef_primitive_major.data();
	static std::mutex coef_cache_mutex;
	std::lock_guard<std::mutex> lock(coef_cache_mutex);
	if (valid.load(std::memory_order_acquire)
		&& coef_primitive_major.size() == (size_t)nex * (size_t)_nmo)
		return coef_primitive_major.data();
	valid.store(false, std::memory_order_release);
	coef_primitive_major.assign((size_t)nex * (size_t)_nmo, 0.0);
	for (int mo = 0; mo < _nmo; mo++) {
		const double* src = MOs[mo].get_coefficient_ptr();
		if (!src) continue;
		for (int j = 0; j < nex; j++)
			coef_primitive_major[(size_t)j * _nmo + mo] = src[j];
	}
	build_exp_groups();
	build_ao_contraction();
	valid.store(true, std::memory_order_release);
	return coef_primitive_major.data();
}

void WFN::build_ao_contraction() const
{
	const int _nmo = get_nmo(false);
	const double* const cp = coef_primitive_major.data();
	prim_ao.assign(nex, -1);
	prim_ao_scale.assign(nex, 0.0);
	ivec rep;  //the first primitive of every function; its row is the function's row
	std::map<std::pair<int, int>, ivec> by_shape;  //(centre, type) -> functions
	for (int j = 0; j < nex; j++) {
		const double* r = cp + (size_t)j * _nmo;
		int m = 0;
		for (int mo = 1; mo < _nmo; mo++)
			if (std::abs(r[mo]) > std::abs(r[m])) m = mo;
		if (r[m] == 0.0) continue;
		ivec& fs = by_shape[{ centers[j], types[j] }];
		for (const int a : fs) {
			const double* q = cp + (size_t)rep[a] * _nmo;
			if (q[m] == 0.0) continue;
			const double s = r[m] / q[m];
			bool same = true;
			for (int mo = 0; mo < _nmo && same; mo++) same = std::abs(r[mo] - s * q[mo]) <= 1e-12 * std::abs(r[m]);
			if (same) {
				prim_ao[j] = a;
				prim_ao_scale[j] = s;
				break;
			}
		}
		if (prim_ao[j] < 0) {
			prim_ao[j] = (int)rep.size();
			prim_ao_scale[j] = 1.0;
			fs.push_back(prim_ao[j]);
			rep.push_back(j);
		}
	}
	coef_ao_major.resize(rep.size() * (size_t)_nmo);
	for (size_t a = 0; a < rep.size(); a++)
		std::copy_n(cp + (size_t)rep[a] * _nmo, _nmo, coef_ao_major.data() + a * _nmo);
}


void WFN::build_exp_groups() const
{
	std::vector<vec> per_center(ncen);
	ivec local(nex);
	for (int j = 0; j < nex; j++) {
		vec &e = per_center[centers[j] - 1];
		const size_t g = std::find(e.begin(), e.end(), exponents[j]) - e.begin();
		if (g == e.size()) e.push_back(exponents[j]);
		local[j] = (int)g;
	}
	group_exponent.clear();
	center_group_start.assign((size_t)ncen + 1, 0);
	center_min_exponent.assign(ncen, 0.0);
	for (int c = 0; c < ncen; c++) {
		center_group_start[c] = (int)group_exponent.size();
		group_exponent.insert(group_exponent.end(), per_center[c].begin(), per_center[c].end());
		if (!per_center[c].empty()) center_min_exponent[c] = *std::min_element(per_center[c].begin(), per_center[c].end());
	}
	center_group_start[ncen] = (int)group_exponent.size();
	prim_exp_group.resize(nex);
	for (int j = 0; j < nex; j++) prim_exp_group[j] = center_group_start[centers[j] - 1] + local[j];
}

const void WFN::computeValues(
	const d3 &PosGrid, // [3] vector with current position on te grid
	double &Rho,           // Value of Electron Density
	double &normGrad,      // Gradiant Vector
	double *Hess,          // Hessian Matrix, later used to determine lambda2
	double &Elf,           // Value of the ELF
	double &Eli,           // Value of the ELI
	double &Lap            // Value for the Laplacian
) const
{
	d3 Grad;
	double tau;
	computeValues(PosGrid, Rho, Grad, Hess, tau);
	Elf = 0;
	if (Rho > 0)
	{
		normGrad = constants::alpha_coef * sqrt(Grad[0] * Grad[0] + Grad[1] * Grad[1] + Grad[2] * Grad[2]) / pow(Rho, constants::c_43);
		Elf = 1 / (1 + pow(constants::ctelf * pow(Rho, constants::c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
		Eli = 0.5 * Rho * pow(48 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), constants::c_38);
	}
	Lap = Hess[0] + Hess[4] + Hess[8];
};

void WFN::computeValues(const d3 &PosGrid, double &Rho, d3 &Grad, double *Hess, double &tau) const
{
	const int _nmo = get_nmo(false);
	vec phi(10 * _nmo, 0.0);
	double *phi_temp;
	double chi[10]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double d[4]{ 0, 0, 0, 0 };
	int iat = 0, k, j;
	int l[3]{ 0, 0, 0 };
	double ex = 0;
	double xl[3][3]{ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
	Grad = { 0, 0, 0 };
	tau = 0;

	Rho = 0;
	Hess[0] = 0;
	Hess[1] = 0;
	Hess[2] = 0;
	Hess[8] = 0;
	Hess[4] = 0;
	Hess[5] = 0;

	const MO *MOs_data = MOs.data();
	//Primitive-major coefficients: every MO for one primitive is contiguous here, where
	//MOs keeps them nex apart. Identical arithmetic in identical order - one sequential
	//read per primitive instead of nmo scattered ones, and no per-point heap allocation.
	const double *const coefs = get_coef_primitive_major();

	for (j = 0; j < nex; j++)
	{
		iat = get_center(j) - 1;

		constants::type2vector(get_type(j), l);
		d[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
		d[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
		d[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
		d[3] = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
		double temp = -get_exponent(j) * (d[3]);
		if (temp < constants::exp_cutoff)
			continue;
		ex = exp(temp);
		for (k = 0; k < 3; k++)
		{
			if (l[k] == 0)
			{
				xl[k][0] = 1.0;
				xl[k][1] = 0.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 1)
			{
				xl[k][0] = d[k];
				xl[k][1] = 1.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 2)
			{
				xl[k][0] = d[k] * d[k];
				xl[k][1] = 2 * d[k];
				xl[k][2] = 2;
			}
			else if (l[k] == 3)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d[k];
				xl[k][1] = 3 * d2;
				xl[k][2] = 6 * d[k];
			}
			else if (l[k] == 4)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2;
				xl[k][1] = 4 * d2 * d[k];
				xl[k][2] = 12 * d2;
			}
			else if (l[k] == 5)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d[k];
				xl[k][1] = 5 * d2 * d2;
				xl[k][2] = 20 * d2 * d[k];
			}
			else if (l[k] == 6)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d2;
				xl[k][1] = 6 * d2 * d2 * d[k];
				xl[k][2] = 30 * d2 * d2;
			}
			else
			{
				return;
			}
		}
		const double ex2 = 2 * get_exponent(j);
		chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
		chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
		chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
		chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
		const double temp_ex = pow(ex2, 2);
		chi[4] = (xl[0][2] - ex2 * (2 * l[0] + 1) * xl[0][0] + temp_ex * pow(d[0], l[0] + 2)) * xl[1][0] * xl[2][0] * ex;
		chi[5] = (xl[1][2] - ex2 * (2 * l[1] + 1) * xl[1][0] + temp_ex * pow(d[1], l[1] + 2)) * xl[2][0] * xl[0][0] * ex;
		chi[6] = (xl[2][2] - ex2 * (2 * l[2] + 1) * xl[2][0] + temp_ex * pow(d[2], l[2] + 2)) * xl[0][0] * xl[1][0] * ex;
		chi[7] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[2][0] * ex;
		chi[8] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[1][0] * ex;
		chi[9] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * ex;

		const double *c_row = coefs + (size_t)j * nmo;
		double *phi_ptr = phi.data();
		for (int mo = 0; mo < nmo; ++mo, phi_ptr += 10)
		{
			const double c = c_row[mo];
			for (k = 0; k < 10; k++)
				phi_ptr[k] += c * chi[k];
		}
	}
	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = get_MO_occ(mo);
		const double docc = 2 * occ;
		if (occ != 0)
		{
			phi_temp = &phi[mo * 10];
			Rho += occ * pow(*phi_temp, 2);
			Grad[0] += docc * *phi_temp * phi_temp[1];
			Grad[1] += docc * *phi_temp * phi_temp[2];
			Grad[2] += docc * *phi_temp * phi_temp[3];
			Hess[0] += docc * (*phi_temp * phi_temp[4] + pow(phi_temp[1], 2));
			Hess[4] += docc * (*phi_temp * phi_temp[5] + pow(phi_temp[2], 2));
			Hess[8] += docc * (*phi_temp * phi_temp[6] + pow(phi_temp[3], 2));
			Hess[1] += docc * (*phi_temp * phi_temp[7] + phi_temp[1] * phi_temp[2]);
			Hess[2] += docc * (*phi_temp * phi_temp[8] + phi_temp[1] * phi_temp[3]);
			Hess[5] += docc * (*phi_temp * phi_temp[9] + phi_temp[2] * phi_temp[3]);
			tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
		}
	}

	Hess[3] = Hess[1];
	Hess[6] = Hess[2];
	Hess[7] = Hess[5];
};

const double WFN::computeMO(
	const d3 &PosGrid, // [3] array with current position on the grid
	const int &mo) const
{
	double result = 0.0;
	int l[3]{ 0, 0, 0 };
	double ex = 0, *d_;
	// x, y, z, r^2 and the powers 2..5 laid out as in compute_dens_cartesian; per-thread scratch
	// because this is called per grid point per orbital from parallel loops
	thread_local vec d;
	if (d.size() < (size_t)16 * ncen) d.resize((size_t)16 * ncen);
	for (int iat = 0; iat < ncen; iat++)
	{
		d_ = d.data() + (size_t)16 * iat;
		d_[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
		d_[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
		d_[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
		d_[4] = d_[0] * d_[0];
		d_[5] = d_[1] * d_[1];
		d_[6] = d_[2] * d_[2];
		d_[3] = d_[4] + d_[5] + d_[6];
		for (int k = 7; k < 16; k++)
			d_[k] = d_[k - 3] * d_[(k - 7) % 3];
	}
	const double *c = MOs[mo].get_coefficient_ptr();
	for (int j = 0; j < nex; j++)
	{
		d_ = d.data() + (size_t)16 * (get_center(j) - 1);
		ex = -get_exponent(j) * d_[3];
		if (ex < constants::exp_cutoff)
			continue;
		ex = exp(ex);
		constants::type2vector(get_type(j), l);
		// power p of coordinate k sits at k for p = 1 and at 3p - 2 + k above
		for (int k = 0; k < 3; k++)
			if (l[k] == 6) ex *= d_[13 + k] * d_[k];
			else if (l[k]) ex *= d_[(l[k] == 1 ? 0 : 3 * l[k] - 2) + k];
		result += c[j] * ex;
	}
	return result;
}
