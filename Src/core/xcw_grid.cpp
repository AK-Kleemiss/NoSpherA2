#include "pch.h"
#include "structure_factors.h"
#include "AtomGrid.h"
#include "fourier_bessel.h"
#ifdef NOSPHERA2_USE_GPU
#include "xcw_grid_gpu.h"
#endif
#if defined(__SSE2__) || defined(_M_X64)
#include <xmmintrin.h>
#endif

//XCW structure factors and perturbation on build_I's grid instead of through the I tensor, see
//structure_factors::xcw_grid. Both walk the same blocks with the same screening, so the
//tensor path and this one differ by summation order only.
namespace {
	//C = A B, A symmetric na x na with its upper triangle set, B na x np, row-major
	void block_symm(const int na, const int np, const double* a, const double* b, double* c)
	{
		cblas_dsymm(CblasRowMajor, CblasLeft, CblasUpper, na, np, 1.0, a, na, b, np, 0.0, c, np);
	}
	void block_symm(const int na, const int np, const float* a, const float* b, float* c)
	{
		cblas_ssymm(CblasRowMajor, CblasLeft, CblasUpper, na, np, 1.0f, a, na, b, np, 0.0f, c, np);
	}
	//C = A B^T, both na x np row-major
	void block_abt(const int na, const int np, const double* a, const double* b, double* c)
	{
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, na, na, np, 1.0, a, np, b, np, 0.0, c, na);
	}
	void block_abt(const int na, const int np, const float* a, const float* b, float* c)
	{
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, na, na, np, 1.0f, a, np, b, np, 0.0f, c, na);
	}

	//What build_I's integration sets per thread: one BLAS thread inside the region, and in
	//single precision the subnormal AO tails flushed
	struct block_thread {
		unsigned int csr = 0;
		explicit block_thread(const bool fp32)
		{
#if !defined(__APPLE__) && !defined(NSA2_OPENBLAS)
			mkl_set_num_threads_local(1);
#endif
#if defined(__SSE2__) || defined(_M_X64)
			csr = _mm_getcsr();
			if (fp32) _mm_setcsr(csr | 0x8040);
#endif
		}
		~block_thread()
		{
#if defined(__SSE2__) || defined(_M_X64)
			_mm_setcsr(csr);
#endif
		}
	};

	//NOSPHERA2_XCW_GRID_TIME: wall time of each stage, printed per call
	struct stage_clock {
		const bool on = tuning("NOSPHERA2_XCW_GRID_TIME") != nullptr;
		std::chrono::steady_clock::time_point t = std::chrono::steady_clock::now();
		void lap(const char* what)
		{
			if (!on) return;
			const auto n = std::chrono::steady_clock::now();
			std::cout << "XCW grid " << what << " " << std::chrono::duration<double, std::milli>(n - t).count() << " ms\n";
			t = n;
		}
	};
}

//calc_SF_fourier_bessel's transform (fourier_bessel.h) planned once on build_I's grid and reused every
//lambda, atom by atom, with its exact transpose for grid_contract. f(k) = Sum_p rho_p e^{ik.d_p} by
//shells: projection C_s,lm = Sum_{p in s} rho_p Y_lm(u_p), radial T_lm(q_n) = Sum_s J_l(q_n r_s) C_s,lm,
//assembly f = Sum_m w_m(|k|) Y(k^) . T(q_{n0+m}).
struct xcw_fb {
	struct rule { int n = 0, Lr = 0, re = 0, ro = 0; ivec sh; vec Y; }; //points, largest Ls, parity widths, shells, Y [n][re | ro]
	struct proto {
		int ns = 0, L = 0, ne = 0, no = 0;
		ivec Ls, sl; //cut of each shell, first shell reaching l
		std::vector<rule> rules;
		vec J; //[l][node][shell] 4pi (+-1) j_l(q_n r_s), zero past Ls[s]
		int nl() const { return ne + no; }
		int off(const int l) const { return ((l & 1) ? ne : 0) + l * (l - 1) / 2; } //degree l in a [even | odd] row
	};
	std::vector<proto> protos;
	ivec proto_of; //[atom]
	std::vector<std::vector<std::vector<std::pair<int, int>>>> at; //[atom][rule] (point, its place i nsh + slot in W)
	int ne = 0, no = 0, nn = 0, S = 1; //Y(k^) widths for the largest L, |k| nodes, nodes per entry
	ivec n0; //[entry] first node
	vec wt, Ye, Yo; //[entry] S Lagrange weights, Y(k^) [ne], [no]

	//sf[g nc + c] = f_g(k_c)
	void forward(const vec2& rho, cvec& sf) const
	{
		const int ng = static_cast<int>(proto_of.size()), nc = static_cast<int>(n0.size());
		vec2 C(ng), T(ng);
		std::vector<std::pair<int, int>> tp, tr;
		for (int g = 0; g < ng; g++) {
			const proto& P = protos[proto_of[g]];
			C[g].assign(static_cast<size_t>(P.ns) * P.nl(), 0.0);
			T[g].assign(static_cast<size_t>(nn) * P.nl(), 0.0);
			for (int t = 0; t < static_cast<int>(P.rules.size()); t++) tp.emplace_back(g, t);
			for (int l = P.L; l >= 0; l--) if (P.sl[l] < P.ns) tr.emplace_back(g, l);
		}
#pragma omp parallel
		{
			const block_thread bt(false);
			vec W, Ct;
#pragma omp for schedule(dynamic)
			for (int i = 0; i < static_cast<int>(tp.size()); i++) {
				const int g = tp[i].first;
				const proto& P = protos[proto_of[g]];
				const rule& R = P.rules[tp[i].second];
				const int nsh = static_cast<int>(R.sh.size()), rl = R.re + R.ro;
				W.assign(static_cast<size_t>(R.n) * nsh, 0.0);
				for (const auto& pp : at[g][tp[i].second]) W[pp.second] = rho[g][pp.first];
				Ct.resize(static_cast<size_t>(rl) * nsh);
				cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, rl, nsh, R.n, 1.0, R.Y.data(), rl, W.data(), nsh, 0.0, Ct.data(), nsh);
				for (int q = 0; q < nsh; q++) {
					double* c = C[g].data() + static_cast<size_t>(R.sh[q]) * P.nl();
					for (int j = 0; j < R.re; j++) c[j] = Ct[static_cast<size_t>(j) * nsh + q];
					for (int j = 0; j < R.ro; j++) c[P.ne + j] = Ct[static_cast<size_t>(R.re + j) * nsh + q];
				}
			}
#pragma omp for schedule(dynamic)
			for (int i = 0; i < static_cast<int>(tr.size()); i++) {
				const int g = tr[i].first, l = tr[i].second;
				const proto& P = protos[proto_of[g]];
				const int sl = P.sl[l], nl = P.nl(), o = P.off(l);
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nn, 2 * l + 1, P.ns - sl, 1.0,
					P.J.data() + static_cast<size_t>(l) * nn * P.ns + sl, P.ns, C[g].data() + static_cast<size_t>(sl) * nl + o, nl,
					0.0, T[g].data() + o, nl);
			}
		}
#pragma omp parallel for schedule(static)
		for (int c = 0; c < nc; c++) {
			const double* ye = Ye.data() + static_cast<size_t>(c) * ne, * yo = Yo.data() + static_cast<size_t>(c) * no, * w = wt.data() + static_cast<size_t>(c) * S;
			for (int g = 0; g < ng; g++) {
				const proto& P = protos[proto_of[g]];
				const int nl = P.nl();
				double re = 0, im = 0;
				for (int m = 0; m < S; m++) {
					const double* t = T[g].data() + static_cast<size_t>(n0[c] + m) * nl;
					double se = 0, so = 0;
					for (int j = 0; j < P.ne; j++) se += ye[j] * t[j];
					for (int j = 0; j < P.no; j++) so += yo[j] * t[P.ne + j];
					re += w[m] * se;
					im += w[m] * so;
				}
				sf[static_cast<size_t>(g) * nc + c] = cdouble(re, im);
			}
		}
	}

	//The transpose: V_g(p) = Re Sum_c H[g nc + c] f_g's kernel e^{ik_c.d_p}
	void adjoint(const cvec& H, vec2& V) const
	{
		const int ng = static_cast<int>(proto_of.size()), nc = static_cast<int>(n0.size());
		vec2 gT(ng), gC(ng);
		std::vector<std::pair<int, int>> tp, tr;
		for (int g = 0; g < ng; g++) {
			const proto& P = protos[proto_of[g]];
			gC[g].assign(static_cast<size_t>(P.ns) * P.nl(), 0.0);
			for (int t = 0; t < static_cast<int>(P.rules.size()); t++) tp.emplace_back(g, t);
			for (int l = P.L; l >= 0; l--) if (P.sl[l] < P.ns) tr.emplace_back(g, l);
		}
#pragma omp parallel
		{
			const block_thread bt(false);
			vec gCt, gW;
#pragma omp for schedule(dynamic)
			for (int g = 0; g < ng; g++) {
				const proto& P = protos[proto_of[g]];
				const int nl = P.nl();
				gT[g].assign(static_cast<size_t>(nn) * nl, 0.0);
				for (int c = 0; c < nc; c++) {
					const cdouble h = H[static_cast<size_t>(g) * nc + c];
					if (h == cdouble(0, 0)) continue;
					const double* ye = Ye.data() + static_cast<size_t>(c) * ne, * yo = Yo.data() + static_cast<size_t>(c) * no;
					for (int m = 0; m < S; m++) {
						double* t = gT[g].data() + static_cast<size_t>(n0[c] + m) * nl;
						const double a = wt[static_cast<size_t>(c) * S + m] * h.real(), b = -wt[static_cast<size_t>(c) * S + m] * h.imag();
						for (int j = 0; j < P.ne; j++) t[j] += a * ye[j];
						for (int j = 0; j < P.no; j++) t[P.ne + j] += b * yo[j];
					}
				}
			}
#pragma omp for schedule(dynamic)
			for (int i = 0; i < static_cast<int>(tr.size()); i++) {
				const int g = tr[i].first, l = tr[i].second;
				const proto& P = protos[proto_of[g]];
				const int sl = P.sl[l], nl = P.nl(), o = P.off(l);
				cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, P.ns - sl, 2 * l + 1, nn, 1.0,
					P.J.data() + static_cast<size_t>(l) * nn * P.ns + sl, P.ns, gT[g].data() + o, nl,
					0.0, gC[g].data() + static_cast<size_t>(sl) * nl + o, nl);
			}
#pragma omp for schedule(dynamic)
			for (int i = 0; i < static_cast<int>(tp.size()); i++) {
				const int g = tp[i].first;
				const proto& P = protos[proto_of[g]];
				const rule& R = P.rules[tp[i].second];
				const int nsh = static_cast<int>(R.sh.size()), rl = R.re + R.ro;
				gCt.resize(static_cast<size_t>(rl) * nsh);
				for (int q = 0; q < nsh; q++) {
					const double* c = gC[g].data() + static_cast<size_t>(R.sh[q]) * P.nl();
					for (int j = 0; j < R.re; j++) gCt[static_cast<size_t>(j) * nsh + q] = c[j];
					for (int j = 0; j < R.ro; j++) gCt[static_cast<size_t>(R.re + j) * nsh + q] = c[P.ne + j];
				}
				gW.resize(static_cast<size_t>(R.n) * nsh);
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, R.n, nsh, rl, 1.0, R.Y.data(), rl, gCt.data(), nsh, 0.0, gW.data(), nsh);
				for (const auto& pp : at[g][tp[i].second]) V[g][pp.first] = gW[pp.second];
			}
		}
	}
};

void structure_factors::plan_grid_fb(const sf_shell_layout& layout)
{
	XG.fb.reset();
	const int ng = static_cast<int>(XG.blocks.size()), nc = static_cast<int>(XG.inversion_link.size());
	const auto direct = [](const char* why) { std::cout << "XCW grid transform: direct sum, " << why << std::endl; };
	if (tuning("NOSPHERA2_XCW_GRID_DIRECT")) return direct("asked for");
	if (static_cast<int>(layout.proto.size()) < ng || static_cast<int>(layout.idx.size()) < ng) return direct("no shell layout");
	auto fb = std::make_shared<xcw_fb>();
	xcw_fb& F = *fb;
	F.proto_of.assign(ng, -1);
	ivec pmap(layout.protos.size(), -1);
	for (int g = 0; g < ng; g++) {
		const int P = layout.proto[g];
		const ivec& ix = layout.idx[g];
		if (P < 0 || ix.size() != XG.w[g].size()) return direct("an atom grid has no prototype");
		const sf_shell_layout::prototype& pr = layout.protos[P];
		for (size_t p = 0; p < ix.size(); p++)
			if (std::abs(pr.x[ix[p]] - XG.d1[g][p]) > 1e-8 || std::abs(pr.y[ix[p]] - XG.d2[g][p]) > 1e-8 || std::abs(pr.z[ix[p]] - XG.d3[g][p]) > 1e-8)
				return direct("grid points are not their prototype's");
		if (pmap[P] < 0) { pmap[P] = static_cast<int>(F.protos.size()); F.protos.emplace_back(); }
		F.proto_of[g] = pmap[P];
		xcw_fb::proto& Q = F.protos[pmap[P]];
		const ivec& ss = pr.shell_start;
		if (!ix.empty())
			Q.ns = std::max(Q.ns, static_cast<int>(std::upper_bound(ss.begin(), ss.end(), *std::max_element(ix.begin(), ix.end())) - ss.begin()));
	}
	vec qk(nc);
	double qmax = 0;
	for (int c = 0; c < nc; c++) {
		const int i = XG.inversion_link[c][0];
		qmax = std::max(qmax, qk[c] = std::sqrt(k_pt[0][i] * k_pt[0][i] + k_pt[1][i] * k_pt[1][i] + k_pt[2][i] * k_pt[2][i]));
	}
	//shells, cut and Debye limits per prototype as calc_SF_fourier_bessel; NOSPHERA2_XCW_GRID_FB_H sets
	//the node step h r_max of the 4-point Lagrange table, error ~ (h r)^4
	const char* hk = tuning("NOSPHERA2_XCW_GRID_FB_H");
	const double hfac = hk ? std::atof(hk) : 0.25;
	double h = 0.025;
	int L = 0;
	std::vector<vec> rs(F.protos.size());
	for (size_t P = 0; P < layout.protos.size(); P++) {
		if (pmap[P] < 0) continue;
		xcw_fb::proto& Q = F.protos[pmap[P]];
		const sf_shell_layout::prototype& pr = layout.protos[P];
		const ivec& ss = pr.shell_start;
		vec& r = rs[pmap[P]];
		r.resize(Q.ns);
		int most = 1;
		for (int s = 0; s < Q.ns; s++) {
			const int p = ss[s];
			r[s] = std::sqrt(pr.x[p] * pr.x[p] + pr.y[p] * pr.y[p] + pr.z[p] * pr.z[p]);
			most = std::max(most, ss[s + 1] - ss[s]);
		}
		Q.L = fb_lebedev_degree(most) / 2;
		Q.ne = fb_n_even(Q.L); Q.no = fb_n_odd(Q.L);
		Q.Ls.resize(Q.ns);
		for (int s = 0; s < Q.ns; s++) {
			const double x = qmax * r[s];
			Q.Ls[s] = std::min(Q.L, static_cast<int>(std::ceil(x + 9.4 * std::cbrt(x))) + 6);
		}
		Q.sl.assign(Q.L + 1, Q.ns);
		for (int l = 0; l <= Q.L; l++)
			for (int s = 0; s < Q.ns && Q.sl[l] == Q.ns; s++)
				if (Q.Ls[s] >= l) Q.sl[l] = s;
		if (Q.ns > 0) h = std::min(h, hfac / std::max(r[Q.ns - 1], 1e-3));
		L = std::max(L, Q.L);
		//shells on the same Lebedev rule and rotation share Y_lm(u_p)
		ivec rule_of(Q.ns), slot(Q.ns);
		for (int s = 0; s < Q.ns; s++) {
			const int n = ss[s + 1] - ss[s], p1 = ss[s];
			const double i1 = r[s] > 0.0 ? 1.0 / r[s] : 0.0;
			int R = -1;
			for (int t = 0; t < static_cast<int>(Q.rules.size()) && R < 0; t++) {
				const int s0 = Q.rules[t].sh[0], p0 = ss[s0];
				if (Q.rules[t].n != n) continue;
				const double i0 = r[s0] > 0.0 ? 1.0 / r[s0] : 0.0;
				bool same = true;
				for (int i = 0; i < n && same; i++)
					same = std::abs(pr.x[p0 + i] * i0 - pr.x[p1 + i] * i1) < 1e-10 && std::abs(pr.y[p0 + i] * i0 - pr.y[p1 + i] * i1) < 1e-10
					&& std::abs(pr.z[p0 + i] * i0 - pr.z[p1 + i] * i1) < 1e-10;
				if (same) R = t;
			}
			if (R < 0) { R = static_cast<int>(Q.rules.size()); Q.rules.emplace_back(); Q.rules[R].n = n; }
			slot[s] = static_cast<int>(Q.rules[R].sh.size());
			Q.rules[R].sh.push_back(s);
			Q.rules[R].Lr = std::max(Q.rules[R].Lr, Q.Ls[s]);
			rule_of[s] = R;
		}
		const fb_ylm Y(Q.L);
		for (xcw_fb::rule& R : Q.rules) {
			R.re = fb_n_even(R.Lr); R.ro = fb_n_odd(R.Lr);
			const int rl = R.re + R.ro, p0 = ss[R.sh[0]];
			const double ir = r[R.sh[0]] > 0.0 ? 1.0 / r[R.sh[0]] : 0.0; //at r = 0 only l = 0 survives j_l(0)
			R.Y.resize(static_cast<size_t>(R.n) * rl);
			for (int i = 0; i < R.n; i++)
				Y.eval(pr.x[p0 + i] * ir, pr.y[p0 + i] * ir, pr.z[p0 + i] * ir, R.Y.data() + static_cast<size_t>(i) * rl, R.Y.data() + static_cast<size_t>(i) * rl + R.re, R.Lr);
		}
		//each member atom's points by rule, at their place in that rule's W
		for (int g = 0; g < ng; g++) {
			if (F.proto_of[g] != pmap[P]) continue;
			if (F.at.empty()) F.at.resize(ng);
			F.at[g].resize(Q.rules.size());
			const ivec& ix = layout.idx[g];
			for (size_t p = 0; p < ix.size(); p++) {
				const int s = static_cast<int>(std::upper_bound(ss.begin(), ss.end(), ix[p]) - ss.begin()) - 1, t = rule_of[s];
				F.at[g][t].emplace_back(static_cast<int>(p), (ix[p] - ss[s]) * static_cast<int>(Q.rules[t].sh.size()) + slot[s]);
			}
		}
	}
	//|k| nodes: the Lagrange table when it is shorter than the entries, else every entry its own node
	const int nq = static_cast<int>(qmax / h) + 4;
	const bool table = nq < nc;
	F.nn = table ? nq : nc;
	F.S = table ? 4 : 1;
	F.ne = fb_n_even(L); F.no = fb_n_odd(L);
	F.n0.resize(nc); F.wt.resize(static_cast<size_t>(nc) * F.S);
	F.Ye.resize(static_cast<size_t>(nc) * F.ne); F.Yo.resize(static_cast<size_t>(nc) * F.no);
	const fb_ylm Y(L);
	for (int c = 0; c < nc; c++) {
		const int i = XG.inversion_link[c][0];
		const double q = qk[c], iq = q > 0.0 ? 1.0 / q : 0.0;
		Y.eval(k_pt[0][i] * iq, k_pt[1][i] * iq, q > 0.0 ? k_pt[2][i] * iq : 1.0, F.Ye.data() + static_cast<size_t>(c) * F.ne, F.Yo.data() + static_cast<size_t>(c) * F.no);
		double* w = F.wt.data() + static_cast<size_t>(c) * F.S;
		if (!table) { F.n0[c] = c; w[0] = 1.0; continue; }
		F.n0[c] = std::clamp(static_cast<int>(q / h) - 1, 0, nq - 4);
		const double u = q / h - F.n0[c];
		w[0] = -(u - 1) * (u - 2) * (u - 3) / 6.0; w[1] = u * (u - 2) * (u - 3) / 2.0;
		w[2] = -u * (u - 1) * (u - 3) / 2.0; w[3] = u * (u - 1) * (u - 2) / 6.0;
	}
	for (size_t P = 0; P < F.protos.size(); P++) {
		xcw_fb::proto& Q = F.protos[P];
		Q.J.assign(static_cast<size_t>(Q.L + 1) * F.nn * Q.ns, 0.0);
#pragma omp parallel
		{
			vec jl(Q.L + 1);
#pragma omp for schedule(static)
			for (int n = 0; n < F.nn; n++) {
				const double q = table ? n * h : qk[n];
				for (int s = 0; s < Q.ns; s++) {
					fb_sph_bessel(Q.Ls[s], q * rs[P][s], jl.data());
					for (int l = 0; l <= Q.Ls[s]; l++)
						Q.J[(static_cast<size_t>(l) * F.nn + n) * Q.ns + s] = constants::FOUR_PI * (((l >> 1) & 1) ? -1.0 : 1.0) * jl[l];
				}
			}
		}
	}
	std::cout << "XCW grid transform: Fourier-Bessel, L " << L << ", " << F.nn << (table ? " table" : " entry") << " nodes, "
		<< F.protos.size() << " prototypes" << std::endl;
	XG.fb = std::move(fb);
}

void structure_factors::grid_F(const dMatrix2& D)
{
	const int ng = static_cast<int>(XG.blocks.size()), nc = static_cast<int>(XG.inversion_link.size());
	const cvec sfc = grid_sfc(D);
	stage_clock sc;
	//calc_F_calc(cvec2)'s assembly, the atom transform read at the cleaned entry and conjugated
	//at its inverse: rho is real
#pragma omp parallel for
	for (int r = 0; r < model_data.nr; r++) {
		const ivec& lk = XG.lookup[r];
		cdouble sum = scatter_data.anom_correction[r];
		for (int g = 0; g < ng; g++) {
			cdouble t = 0;
			for (int s = 0; s < static_cast<int>(lk.size()); s++) {
				const int e = lk[s];
				cdouble v = sfc[static_cast<size_t>(g) * nc + XG.clean_of[e]];
				if (XG.is_inverse[e]) v = std::conj(v);
				t += v * phase_facts[g][e] * DW_facts[g][e] * translation_phase_facts[r][s];
			}
			sum += t * asym_atoms[g].sym_fact;
		}
		scatter_data.F_calc[r] = sum;
	}
	sc.lap("F assembly");
}

cvec structure_factors::grid_sfc(const dMatrix2& D)
{
	const int ng = static_cast<int>(XG.blocks.size()), nc = static_cast<int>(XG.inversion_link.size());
	//i_double: nothing the SCF reads is rounded to single, as it means for the tensor
	const bool fp32 = opt->cpu_itensor_fp32 && !I_tens.tensor_double;
	std::vector<std::pair<int, int>> list;
	stage_clock sc;
	for (int g = 0; g < ng; g++)
		for (int b = 0; b < static_cast<int>(XG.blocks[g].size()); b++) list.emplace_back(g, b);
	const int nb = static_cast<int>(list.size()), nmo = static_cast<int>(D.extent(0));
	//NOSPHERA2_XCW_GRID_DTOL: per block bound on the |dF| a kept rho may carry; off by default: XCW moves D on
	//every pair each cycle, so on P1 the bound is 1e-5..1 per block and 1e-4 keeps 14 % of blocks for 2.5e-4 in F
	const char* dt = tuning("NOSPHERA2_XCW_GRID_DTOL");
	const double dtol = dt && !XG.on_device ? std::atof(dt) : 0.0;
	const bool fresh = static_cast<int>(XG.bound.size()) != nb || XG.D_last.rows() != nmo;
	if (fresh) {
		XG.bound.assign(nb, 0.0);
		XG.drift.assign(nb, 0.0);
		XG.rho_last.resize(ng);
		for (int g = 0; g < ng; g++) XG.rho_last[g].assign(XG.w[g].size(), 0.0);
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < nb; i++) {
			const GridBlock& b = XG.blocks[list[i].first][list[i].second];
			const int na = static_cast<int>(b.active_aos.size()), np = b.point_count;
			const double* w = XG.w[list[i].first].data() + b.point_start;
			double s = 0;
			for (int p = 0; p < np; p++) {
				double a = 0;
				for (int ii = 0; ii < na; ii++) a += std::abs(b.ao_values[static_cast<size_t>(ii) * np + p]);
				s += w[p] * a * a;
			}
			XG.bound[i] = 2.0 * s;
		}
	}
	vec2& rho = XG.rho_last;
	std::vector<char> redo(nb, 1);
	vec dd(nb, 0.0);
	if (!fresh && dtol > 0) {
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < nb; i++) {
			const GridBlock& b = XG.blocks[list[i].first][list[i].second];
			const int na = static_cast<int>(b.active_aos.size());
			const int* a = b.active_aos.data();
			double m = 0;
			for (int ii = 0; ii < na; ii++)
				for (int jj = ii; jj < na; jj++)
					if (jj == ii || !XG.skip[a[ii]][a[jj]]) m = std::max(m, std::abs(D(a[ii], a[jj]) - XG.D_last(a[ii], a[jj])));
			const double d = dd[i] = XG.drift[i] + m * XG.bound[i];
			if (d < dtol) { XG.drift[i] = d; redo[i] = 0; }
		}
	}
	const int kept = nb - static_cast<int>(std::count(redo.begin(), redo.end(), 1));
	if (XG.D_last.rows() != nmo) XG.D_last.resize(nmo, nmo);
	for (int i = 0; i < nmo; i++)
		for (int j = 0; j < nmo; j++) XG.D_last(i, j) = D(i, j);
	sc.lap("F change test");
#ifdef NOSPHERA2_USE_GPU
	if (XG.on_device) {
		vec Dh(static_cast<size_t>(nmo) * nmo), r(XG.gpu_off[ng]);
		for (int i = 0; i < nmo; i++)
			for (int j = 0; j < nmo; j++) Dh[static_cast<size_t>(i) * nmo + j] = D(i, j);
		err_checkf(xcw_grid_gpu_density(Dh.data(), r.data()), "XCW grid density on the device failed", std::cout);
		for (int g = 0; g < ng; g++) std::copy(r.begin() + XG.gpu_off[g], r.begin() + XG.gpu_off[g + 1], rho[g].begin());
		std::fill(redo.begin(), redo.end(), 0);
	}
#endif
	//w_p rho(p) = 2 w_p Sum_munu D chi_mu chi_nu with D restricted to the block's AOs and the
	//screened pairs zeroed: the Sum_{mu<=nu} (2|4) D I of calc_F_calc on the tensor
#pragma omp parallel
	{
		const block_thread bt(fp32);
		vec Dd, Td;
		std::vector<float> Df, Tf;
#pragma omp for schedule(dynamic)
		for (int i = 0; i < nb; i++) {
			if (!redo[i]) continue;
			XG.drift[i] = 0;
			const int g = list[i].first;
			const GridBlock& b = XG.blocks[g][list[i].second];
			const int na = static_cast<int>(b.active_aos.size()), np = b.point_count;
			const int* a = b.active_aos.data();
			double* r = rho[g].data() + b.point_start;
			std::fill(r, r + np, 0.0);
			auto run = [&](const auto* X, auto& Db, auto& T) {
				using F = std::remove_reference_t<decltype(Db[0])>;
				Db.assign(static_cast<size_t>(na) * na, F(0));
				T.resize(static_cast<size_t>(na) * np);
				for (int ii = 0; ii < na; ii++)
					for (int jj = ii; jj < na; jj++)
						if (jj == ii || !XG.skip[a[ii]][a[jj]]) Db[static_cast<size_t>(ii) * na + jj] = static_cast<F>(D(a[ii], a[jj]));
				block_symm(na, np, Db.data(), X, T.data());
				for (int ii = 0; ii < na; ii++) {
					const auto* x = X + static_cast<size_t>(ii) * np;
					const auto* t = T.data() + static_cast<size_t>(ii) * np;
					for (int p = 0; p < np; p++) r[p] += static_cast<double>(x[p]) * t[p];
				}
			};
			if (fp32) run(b.ao_values_f.data(), Df, Tf); else run(b.ao_values.data(), Dd, Td);
			const double* w = XG.w[g].data() + b.point_start;
			for (int p = 0; p < np; p++) r[p] *= 2.0 * w[p];
		}
	}
	if (sc.on) std::cout << "XCW grid F kept " << kept << " of " << nb << " blocks, bound " << *std::min_element(dd.begin(), dd.end()) << " .. " << *std::max_element(dd.begin(), dd.end()) << "\n";
	sc.lap("F density");
	//ponytail: direct sum over points x cleaned entries, P * nc sincos; the exact reference the
	//Fourier-Bessel transform replaces once it is checked against this
	cvec sfc(static_cast<size_t>(ng) * nc);
	if (XG.fb) XG.fb->forward(rho, sfc);
	else {
#pragma omp parallel for schedule(dynamic)
	for (int c = 0; c < nc; c++) {
		const int idx = XG.inversion_link[c][0];
		const double kx = k_pt[0][idx], ky = k_pt[1][idx], kz = k_pt[2][idx];
		for (int g = 0; g < ng; g++) {
			const double* x = XG.d1[g].data(), * y = XG.d2[g].data(), * z = XG.d3[g].data(), * r = rho[g].data();
			double sr = 0, si = 0;
			for (size_t p = 0; p < rho[g].size(); p++) {
				if (r[p] == 0.0) continue;
				const double ang = kx * x[p] + ky * y[p] + kz * z[p];
				sr += r[p] * std::cos(ang);
				si += r[p] * std::sin(ang);
			}
			sfc[static_cast<size_t>(g) * nc + c] = cdouble(sr, si);
		}
	}
	}
	sc.lap("F transform");
	return sfc;
}

void structure_factors::grid_contract(occ::Mat& out, const cvec& pre)
{
	const int nmo = model_data.nmo, ng = static_cast<int>(XG.blocks.size()), nc = static_cast<int>(XG.inversion_link.size());
	const bool fp32 = opt->cpu_itensor_fp32 && !I_tens.tensor_double;
	out.setZero(nmo, nmo);
	stage_clock sc;
	//Transpose of grid_F's assembly: A(e) = Sum_{r,s -> e} pre_r t_rs, then per atom the factor
	//of e, an inverse entry's conjugate folded onto its cleaned one
	cvec A(model_data.nr_enlarged, cdouble(0, 0));
	for (int r = 0; r < model_data.nr; r++) {
		if (pre[r] == cdouble(0, 0)) continue;
		for (int s = 0; s < static_cast<int>(XG.lookup[r].size()); s++)
			A[XG.lookup[r][s]] += pre[r] * translation_phase_facts[r][s];
	}
	cvec H(static_cast<size_t>(ng) * nc, cdouble(0, 0));
#pragma omp parallel for
	for (int g = 0; g < ng; g++)
		for (int e = 0; e < model_data.nr_enlarged; e++) {
			if (A[e] == cdouble(0, 0)) continue;
			const cdouble G = A[e] * phase_facts[g][e] * DW_facts[g][e] * asym_atoms[g].sym_fact;
			H[static_cast<size_t>(g) * nc + XG.clean_of[e]] += XG.is_inverse[e] ? std::conj(G) : G;
		}
	sc.lap("P assembly");
	//V(p) = w_p Re Sum_c H(c) exp(i k_c d_p), the transpose of the direct sum
	vec kx(nc), ky(nc), kz(nc);
	for (int c = 0; c < nc; c++) {
		const int idx = XG.inversion_link[c][0];
		kx[c] = k_pt[0][idx]; ky[c] = k_pt[1][idx]; kz[c] = k_pt[2][idx];
	}
	vec2 V(ng);
	if (XG.fb) {
		for (int g = 0; g < ng; g++) V[g].assign(XG.w[g].size(), 0.0);
		XG.fb->adjoint(H, V);
		for (int g = 0; g < ng; g++)
			for (size_t p = 0; p < V[g].size(); p++) V[g][p] *= XG.w[g][p];
	}
	else
	for (int g = 0; g < ng; g++) {
		V[g].assign(XG.w[g].size(), 0.0);
		const cdouble* h = H.data() + static_cast<size_t>(g) * nc;
#pragma omp parallel for schedule(static)
		for (int p = 0; p < static_cast<int>(V[g].size()); p++) {
			const double x = XG.d1[g][p], y = XG.d2[g][p], z = XG.d3[g][p];
			double v = 0;
			for (int c = 0; c < nc; c++) {
				if (h[c] == cdouble(0, 0)) continue;
				const double ang = kx[c] * x + ky[c] * y + kz[c] * z;
				v += h[c].real() * std::cos(ang) - h[c].imag() * std::sin(ang);
			}
			V[g][p] = XG.w[g][p] * v;
		}
	}
	sc.lap("P transform");
#ifdef NOSPHERA2_USE_GPU
	if (XG.on_device) {
		vec Vh(XG.gpu_off[ng]);
		for (int g = 0; g < ng; g++) std::copy(V[g].begin(), V[g].end(), Vh.begin() + XG.gpu_off[g]);
		err_checkf(xcw_grid_gpu_contract(Vh.data(), out.data()), "XCW grid contraction on the device failed", std::cout);
		//row-major upper triangle into the column-major Mat: its lower one, mirrored below
		out.transposeInPlace();
		for (int mu = 0; mu < nmo; mu++)
			for (int nu = mu + 1; nu < nmo; nu++) out(nu, mu) = out(mu, nu);
		sc.lap("P blocks");
		return;
	}
#endif
	//out(mu, nu) = Sum_p chi_mu V chi_nu over the blocks both are active on, screened pairs left out
	std::vector<std::pair<int, int>> list;
	for (int g = 0; g < ng; g++)
		for (int b = 0; b < static_cast<int>(XG.blocks[g].size()); b++) list.emplace_back(g, b);
	std::vector<occ::Mat> parts(omp_get_max_threads());
#pragma omp parallel
	{
		const block_thread bt(fp32);
		occ::Mat local = occ::Mat::Zero(nmo, nmo);
		vec Yd, Sd;
		std::vector<float> Yf, Sf;
#pragma omp for schedule(dynamic)
		for (int i = 0; i < static_cast<int>(list.size()); i++) {
			const int g = list[i].first;
			const GridBlock& b = XG.blocks[g][list[i].second];
			const int na = static_cast<int>(b.active_aos.size()), np = b.point_count;
			const int* a = b.active_aos.data();
			const double* v = V[g].data() + b.point_start;
			auto run = [&](const auto* X, auto& Y, auto& S) {
				using F = std::remove_reference_t<decltype(Y[0])>;
				Y.resize(static_cast<size_t>(na) * np);
				S.resize(static_cast<size_t>(na) * na);
				for (int ii = 0; ii < na; ii++)
					for (int p = 0; p < np; p++)
						Y[static_cast<size_t>(ii) * np + p] = X[static_cast<size_t>(ii) * np + p] * static_cast<F>(v[p]);
				block_abt(na, np, X, Y.data(), S.data());
				for (int ii = 0; ii < na; ii++)
					for (int jj = ii; jj < na; jj++)
						if (jj == ii || !XG.skip[a[ii]][a[jj]]) local(a[ii], a[jj]) += static_cast<double>(S[static_cast<size_t>(ii) * na + jj]);
			};
			if (fp32) run(b.ao_values_f.data(), Yf, Sf); else run(b.ao_values.data(), Yd, Sd);
		}
		parts[omp_get_thread_num()].swap(local);
	}
	for (const occ::Mat& m : parts)
		if (m.size() != 0) out += m;
	for (int mu = 0; mu < nmo; mu++)
		for (int nu = mu + 1; nu < nmo; nu++) out(nu, mu) = out(mu, nu);
	sc.lap("P blocks");
}

bool structure_factors::hold_grid_gpu()
{
	XG.on_device = false;
#ifdef NOSPHERA2_USE_GPU
	const int nmo = model_data.nmo, ng = static_cast<int>(XG.blocks.size());
	XG.gpu_off.assign(ng + 1, 0);
	for (int g = 0; g < ng; g++) XG.gpu_off[g + 1] = XG.gpu_off[g] + static_cast<int>(XG.w[g].size());
	ivec na, np, pt, ao_off, aos;
	std::vector<long long> x_off;
	std::vector<float> X;
	vec w;
	for (int g = 0; g < ng; g++) {
		w.insert(w.end(), XG.w[g].begin(), XG.w[g].end());
		for (const GridBlock& b : XG.blocks[g]) {
			na.push_back(static_cast<int>(b.active_aos.size())); np.push_back(b.point_count);
			pt.push_back(XG.gpu_off[g] + b.point_start);
			ao_off.push_back(static_cast<int>(aos.size())); aos.insert(aos.end(), b.active_aos.begin(), b.active_aos.end());
			x_off.push_back(static_cast<long long>(X.size()));
			if (!b.ao_values_f.empty()) X.insert(X.end(), b.ao_values_f.begin(), b.ao_values_f.end());
			else X.insert(X.end(), b.ao_values.begin(), b.ao_values.end());
		}
	}
	std::vector<unsigned char> keep(static_cast<size_t>(nmo) * nmo, 0);
	for (int i = 0; i < nmo; i++)
		for (int j = i; j < nmo; j++) keep[static_cast<size_t>(i) * nmo + j] = j == i || !XG.skip[i][j];
	XG.on_device = !na.empty() && xcw_grid_gpu_hold(nmo, static_cast<int>(na.size()), na.data(), np.data(), pt.data(), ao_off.data(), aos.data(),
		x_off.data(), X.data(), static_cast<long long>(X.size()), keep.data(), w.data(), XG.gpu_off[ng]);
#endif
	return XG.on_device;
}
