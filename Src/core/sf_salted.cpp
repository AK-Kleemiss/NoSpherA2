#include "pch.h"
#include "scattering_factors.h"

namespace {
//Atoms whose shells match (one element of one basis) share their transformed basis functions
struct salted_sf_group {
	int t = 0;                 //template atom
	ivec members;              //indices into asym_atom_list
	int ne = 0, nf = 0;        //columns: the ne even-l ones first, then the odd
	ivec col_shell, col_lm;    //shell offset from the template's first, and l*l + m
	vec col_sign;              //PI3_2 times the sign of i^l
	vec Ce, Co;                //[members][ne] and [members][nf - ne], the coefficients in column order
};

bool salted_same_shells(const aux_density_table& t, const int a, const int b)
{
	const int n = t.sh_start[a + 1] - t.sh_start[a];
	if (n != t.sh_start[b + 1] - t.sh_start[b]) return false;
	for (int i = 0; i < n; i++) {
		const int sa = t.sh_start[a] + i, sb = t.sh_start[b] + i, np = t.pr_start[sa + 1] - t.pr_start[sa];
		if (t.sh_l[sa] != t.sh_l[sb] || np != t.pr_start[sb + 1] - t.pr_start[sb]) return false;
		for (int p = 0; p < np; p++)
			if (t.pr_uniq[t.pr_start[sa] + p] != t.pr_uniq[t.pr_start[sb] + p] || t.pr_norm[t.pr_start[sa] + p] != t.pr_norm[t.pr_start[sb] + p]) return false;
	}
	return true;
}
}

//The transform as one GEMM per group of atoms with the same shells: B[j][k] = i^l PI3_2 R_s(|k|) Y_lm(k^) depends on the
//basis function alone, so the form factors of a group are its coefficient rows times B, even l into the real part and odd
//l into the imaginary. The functions are atom-centred, so there is no phase. A per-atom loop recomputed Y_lm per shell
//and summed scalar complex terms, ~25 ns per atom and reflection on 8 threads. aux_density_table::fourier_atom is the
//one-atom reference. A streaming caller owns one bar for the whole table and passes it in, else every block draws its own.
void calc_SF_SALTED(
	const vec2& k_pt,
	const vec& coefs,
	const aux_density_table& table,
	const ivec& asym_atom_list,
	cvec2& sf,
	ProgressBar* progress)
{
	const int num_asym_atoms = static_cast<int>(asym_atom_list.size());
	const int nk = static_cast<int>(k_pt[0].size());
	sf.resize(num_asym_atoms);
	//in parallel: a whole-sphere table is GBs, and one thread zeroing it page by page left the transform flat in threads
#pragma omp parallel for schedule(dynamic)
	for (int ia = 0; ia < num_asym_atoms; ++ia)
		sf[ia].assign(nk, constants::cnull);
	std::unique_ptr<ProgressBar> local_pb;
	if (!progress)
		local_pb = std::make_unique<ProgressBar>(nk, 60, "#", " ", "Generating scattering factors...");
	ProgressBar& pb = progress ? *progress : *local_pb;
	const int n_uniq = static_cast<int>(table.uniq_exp.size());

	std::vector<salted_sf_group> groups;
	for (int ia = 0; ia < num_asym_atoms; ++ia) {
		const int a = asym_atom_list[ia];
		auto g = std::find_if(groups.begin(), groups.end(), [&](const salted_sf_group& x) { return salted_same_shells(table, x.t, a); });
		if (g == groups.end()) {
			groups.emplace_back();
			g = groups.end() - 1;
			g->t = a;
		}
		g->members.push_back(ia);
	}
	int lmax = 0, nsh_max = 0, nf_max = 0, nm_max = 0;
	for (salted_sf_group& g : groups) {
		const int s0 = table.sh_start[g.t], s1 = table.sh_start[g.t + 1];
		ivec col_off;
		for (int odd = 0; odd < 2; odd++)
			for (int s = s0; s < s1; s++) {
				const int l = table.sh_l[s];
				if ((l & 1) != odd) continue;
				for (int m = 0; m <= 2 * l; m++) {
					g.col_shell.push_back(s - s0);
					g.col_lm.push_back(l * l + m);
					g.col_sign.push_back((l & 2) ? -constants::PI3_2 : constants::PI3_2);
					col_off.push_back(table.coef_off[s] - table.coef_off[s0] + m);
				}
				lmax = std::max(lmax, l);
				if (!odd) g.ne += 2 * l + 1;
			}
		g.nf = static_cast<int>(col_off.size());
		const int no = g.nf - g.ne;
		g.Ce.resize(g.members.size() * g.ne);
		g.Co.resize(g.members.size() * no);
		for (size_t i = 0; i < g.members.size(); i++) {
			const double* c = coefs.data() + table.coef_off[table.sh_start[asym_atom_list[g.members[i]]]];
			for (int j = 0; j < g.ne; j++)
				g.Ce[i * g.ne + j] = c[col_off[j]];
			for (int j = 0; j < no; j++)
				g.Co[i * no + j] = c[col_off[g.ne + j]];
		}
		nsh_max = std::max(nsh_max, s1 - s0);
		nf_max = std::max(nf_max, g.nf);
		nm_max = std::max(nm_max, static_cast<int>(g.members.size()));
	}
	//reflections per GEMM: a block of B stays in L2 at ~150 functions an atom
	constexpr int KB = 64;
	const int nblk = (nk + KB - 1) / KB;
	//unit + 16 - m is 1 at index m and 0 elsewhere for m <= 16, so the contraction returns Y_lm alone
	static constexpr double unit[33] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };
	//the radial factor (H/2)^l exp(-H^2/4a) / a^(l+3/2) depends on (a, l) and |k| alone, so it is tabulated per k over the
	//distinct pairs; rad and Y hold KB reflections from k0
	auto tabulate = [&](vec& rad, vec& Y, const int k0, const int kb) {
		for (int kk = 0; kk < kb; kk++) {
			const double kx = k_pt[0][k0 + kk], ky = k_pt[1][k0 + kk], kz = k_pt[2][k0 + kk];
			const double H2 = kx * kx + ky * ky + kz * kz, H = std::sqrt(H2);
			const double k[3] = { H > 0.0 ? kx / H : 0.0, H > 0.0 ? ky / H : 0.0, H > 0.0 ? kz / H : 1.0 };
			for (int u = 0; u < n_uniq; u++) {
				double Hl_over_2l = 1.0;
				for (int i = 0; i < table.uniq_l[u]; i++)
					Hl_over_2l *= 0.5 * H;
				rad[u * KB + kk] = Hl_over_2l * std::exp(-H2 / (4.0 * table.uniq_exp[u])) / table.uniq_exp_l32[u];
			}
			for (int l = 0; l <= lmax; l++)
				for (int m = 0; m <= 2 * l; m++)
					Y[(l * l + m) * KB + kk] = constants::spherical_harmonic(l, k[0], k[1], k[2], unit + 16 - m);
		}
	};
	//B[j * ld + kk] of group g from the tables, rs the contracted shells
	auto fill_B = [&](const vec& rad, const vec& Y, vec& rs, const salted_sf_group& g, const int kb, double* B, const size_t ld) {
		const int s0 = table.sh_start[g.t];
		for (int i = 0; i < table.sh_start[g.t + 1] - s0; i++) {
			double* r = &rs[i * KB];
			std::fill(r, r + kb, 0.0);
			for (int p = table.pr_start[s0 + i]; p < table.pr_start[s0 + i + 1]; p++) {
				const double w = table.pr_norm[p], * q = &rad[table.pr_uniq[p] * KB];
				for (int kk = 0; kk < kb; kk++) r[kk] += w * q[kk];
			}
		}
		for (int j = 0; j < g.nf; j++) {
			const double c = g.col_sign[j], * r = &rs[g.col_shell[j] * KB], * y = &Y[g.col_lm[j] * KB];
			double* d = B + j * ld;
			for (int kk = 0; kk < kb; kk++) d[kk] = c * r[kk] * y[kk];
		}
	};
	//re (nm x kb) = Ce B[:ne, :], im = Co B[ne:, :]. Not offloaded: GEMM on the device lost on every card tried,
	//V100 included, since B and the form factors cross the bus (branch salted-ft-gpu-attempt)
	auto products = [&](const salted_sf_group& g, const int kb, const double* B, const int ldb, vec& re, vec& im, const int ldc) {
		const int nm = static_cast<int>(g.members.size()), no = g.nf - g.ne;
		const auto gemm = [&](const int k, const double* A, const double* Bk, vec& C) {
			if (k == 0) { std::fill(C.begin(), C.end(), 0.0); return; }
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nm, kb, k, 1.0, A, k, Bk, ldb, 0.0, C.data(), ldc);
		};
		gemm(g.ne, g.Ce.data(), B, re);
		gemm(no, g.Co.data(), B + static_cast<size_t>(g.ne) * ldb, im);
	};

#pragma omp parallel
	{
#if !defined(__APPLE__) && !defined(NSA2_OPENBLAS)
		const int mkl_before = mkl_set_num_threads_local(1);
#endif
		vec rad(static_cast<size_t>(n_uniq) * KB), Y(static_cast<size_t>(lmax + 1) * (lmax + 1) * KB), rs(static_cast<size_t>(nsh_max) * KB),
			Bt(static_cast<size_t>(nf_max) * KB), re(static_cast<size_t>(nm_max) * KB), im(static_cast<size_t>(nm_max) * KB);
#pragma omp for schedule(dynamic)
		for (int b = 0; b < nblk; b++)
		{
			const int k0 = b * KB, kb = std::min(KB, nk - k0);
			tabulate(rad, Y, k0, kb);
			for (const salted_sf_group& g : groups) {
				fill_B(rad, Y, rs, g, kb, Bt.data(), KB);
				products(g, kb, Bt.data(), KB, re, im, KB);
				for (size_t i = 0; i < g.members.size(); i++) {
					cdouble* out = sf[g.members[i]].data() + k0;
					for (int kk = 0; kk < kb; kk++) out[kk] = cdouble(re[i * KB + kk], im[i * KB + kk]);
				}
			}
			pb.update(kb);
		}
#if !defined(__APPLE__) && !defined(NSA2_OPENBLAS)
		mkl_set_num_threads_local(mkl_before);
#endif
	}
}
