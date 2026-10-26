#include "pch.h"
#include "structure_factors.h"
#include "constants.h"
#include <occ/qm/integral_engine_df.h>
#include <occ/gto/gto.h>

//XCW's per-atom form factors from an RI fit of rho(D) instead of the grid weights, see structure_factors::xcw_ri
namespace {
	//real spherical harmonics as occ's spherical functions carry them, the convention of the integrals it hands out:
	//a one-primitive shell e^(-r^2) at the origin on the unit sphere, times e; rows m, columns the directions
	occ::Mat occ_ylm(const int l, const occ::Mat3N& dirs)
	{
		occ::gto::Shell sh(l, { 1.0 }, { { 1.0 } }, { 0.0, 0.0, 0.0 });
		sh.kind = occ::gto::Shell::Kind::Spherical;
		occ::gto::AOBasis b({ occ::core::Atom{ 1, 0.0, 0.0, 0.0 } }, { sh });
		b.set_pure(true);
		return occ::gto::evaluate_basis(b, dirs, 0).phi.transpose() * std::exp(1.0);
	}

	occ::Mat3N unit(const vec& x, const vec& y, const vec& z, vec& r)
	{
		const int n = static_cast<int>(x.size());
		occ::Mat3N d(3, n);
		r.resize(n);
		for (int p = 0; p < n; p++) {
			r[p] = std::sqrt(x[p] * x[p] + y[p] * y[p] + z[p] * z[p]);
			if (r[p] > 0.0) { d(0, p) = x[p] / r[p]; d(1, p) = y[p] / r[p]; d(2, p) = z[p] / r[p]; }
			else { d(0, p) = 0.0; d(1, p) = 0.0; d(2, p) = 1.0; }
		}
		return d;
	}
}

void structure_factors::setup_ri(const occ::qm::AOBasis& ao, const occ::qm::AOBasis& aux_in, const int lmax, const double strength)
{
	err_checkf(XG.mode == 1, "-xcw_partition ri runs on the grid route: not with NOSPHERA2_XCW_GRID=0 or a saved I tensor", std::cout);
	occ::qm::AOBasis aux = aux_in;
	aux.set_pure(true);
	const int ng = model_data.ncen, nc = static_cast<int>(XG.inversion_link.size()), n = static_cast<int>(ao.nbf()), na = static_cast<int>(aux.nbf());
	const size_t bytes = static_cast<size_t>(n) * n * na * sizeof(double);
	err_checkf(free_memory() == 0 || bytes < free_memory(), "-xcw_partition ri holds (mu nu|P), " + std::to_string(bytes >> 20) + " MB, more than -mem leaves", std::cout);
	XR.nbf = n;
	XR.naux = na;
	XR.lmax = lmax;
	XR.n_mom = (lmax + 1) * (lmax + 1);
	occ::qm::IntegralEngineDF df(ao.atoms(), ao.shells(), aux.shells());
	const occ::Mat I = occ::Mat::Identity(n, n);
	XR.B = df.build_b_direct(I, I);
	occ::Mat H = df.aux_engine().one_electron_operator(occ::qm::cint::Operator::coulomb, false);
	//the aux functions by atom: a basis lists its shells atom by atom
	const ivec& bf_atom = aux.bf_to_atom();
	XR.atom_start.assign(ng + 1, 0);
	for (int P = 0; P < na; P++) {
		err_checkf(P == 0 || bf_atom[P] >= bf_atom[P - 1], "The auxiliary basis does not list its shells atom by atom", std::cout);
		XR.atom_start[bf_atom[P] + 1] = P + 1;
	}
	for (int g = 0; g < ng; g++) XR.atom_start[g + 1] = std::max(XR.atom_start[g + 1], XR.atom_start[g]);
	int lmax_aux = 0;
	for (const auto& sh : aux.shells()) lmax_aux = std::max(lmax_aux, static_cast<int>(sh.l));
	//k at every cleaned entry, its direction and Y_lm there
	vec kx(nc), ky(nc), kz(nc), k;
	for (int c = 0; c < nc; c++) {
		const int idx = XG.inversion_link[c][0];
		kx[c] = k_pt[0][idx]; ky[c] = k_pt[1][idx]; kz[c] = k_pt[2][idx];
	}
	const occ::Mat3N kd = unit(kx, ky, kz, k);
	std::vector<occ::Mat> Yk(lmax_aux + 1);
	for (int l = 0; l <= lmax_aux; l++) Yk[l] = occ_ylm(l, kd);
	//F_P(k) = i^l Y_lm(k^) Sum_p c_p pi^(3/2) (k/2)^l e^(-k^2/4a) / a^(l+3/2), the transform of c_p r^l e^(-a r^2) Y_lm
	XR.F.assign(ng, cvec());
	if (lmax >= 0) XR.R = occ::Mat::Zero(static_cast<size_t>(ng) * XR.n_mom, na);
	XR.s2.assign(static_cast<size_t>(ng) * XR.n_mom, 0.0);
	const cdouble il[4] = { { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } };
	for (int g = 0; g < ng; g++) {
		const int p0 = XR.atom_start[g], np = XR.atom_start[g + 1] - p0;
		XR.F[g].assign(static_cast<size_t>(nc) * np, cdouble(0, 0));
		const double r_cov = constants::ang2bohr(constants::covalent_radii[asym_atoms[g].type]);
		for (int l = 0; l <= lmax; l++)
			for (int m = 0; m <= 2 * l; m++)
				XR.s2[static_cast<size_t>(g) * XR.n_mom + l * l + m] = std::pow(strength * (l == 0 ? std::sqrt(constants::FOUR_PI) : std::pow(r_cov, -l)), 2);
	}
	for (size_t s = 0; s < aux.size(); s++) {
		const auto& sh = aux[s];
		const int l = sh.l, P0 = static_cast<int>(aux.first_bf()[s]), g = bf_atom[P0], j0 = P0 - XR.atom_start[g], np = XR.atom_start[g + 1] - XR.atom_start[g];
		double moment = 0.0;
		for (int p = 0; p < sh.num_primitives(); p++)
			moment += sh.contraction_coefficients(p, 0) * std::tgamma(l + 1.5) / (2.0 * std::pow(sh.exponents(p), l + 1.5));
#pragma omp parallel for
		for (int c = 0; c < nc; c++) {
			double rad = 0.0;
			for (int p = 0; p < sh.num_primitives(); p++) {
				const double a = sh.exponents(p);
				rad += sh.contraction_coefficients(p, 0) * constants::PI3_2 * std::pow(0.5 * k[c], l) * std::exp(-k[c] * k[c] / (4.0 * a)) / std::pow(a, l + 1.5);
			}
			for (int m = 0; m <= 2 * l; m++)
				XR.F[g][static_cast<size_t>(j0 + m) * nc + c] = il[l % 4] * (rad * Yk[l](m, c));
		}
		if (l <= lmax)
			for (int m = 0; m <= 2 * l; m++)
				XR.R(static_cast<size_t>(g) * XR.n_mom + l * l + m, P0 + m) = moment;
		(void)np;
	}
	//the grid partition's moments about each nucleus, r^l Y_lm at the points of the atom's grid, and per l the map from
	//constants::spherical_harmonic's real harmonics to occ's, fitted on 64 directions
	if (lmax >= 0) {
		std::mt19937 gen(3);
		std::normal_distribution<double> nd;
		vec x(64), y(64), z(64), r;
		for (int i = 0; i < 64; i++) { x[i] = nd(gen); y[i] = nd(gen); z[i] = nd(gen); }
		const occ::Mat3N dirs = unit(x, y, z, r);
		XR.A.assign(lmax + 1, occ::Mat());
		for (int l = 0; l <= lmax; l++) {
			occ::Mat Yn(2 * l + 1, 64);
			for (int i = 0; i < 64; i++) {
				const double d[3] = { dirs(0, i), dirs(1, i), dirs(2, i) };
				for (int m = -l; m <= l; m++) Yn(l + m, i) = constants::spherical_harmonic(l, m, d);
			}
			XR.A[l] = Yn.transpose().colPivHouseholderQr().solve(occ_ylm(l, dirs).transpose()).transpose();
		}
		XR.M.assign(ng, vec());
		for (int g = 0; g < ng; g++) {
			vec r;
			const occ::Mat3N d = unit(XG.d1[g], XG.d2[g], XG.d3[g], r);
			const int np = static_cast<int>(r.size());
			XR.M[g].assign(static_cast<size_t>(XR.n_mom) * np, 0.0);
			for (int l = 0; l <= lmax; l++) {
				const occ::Mat Y = occ_ylm(l, d);
				for (int m = 0; m <= 2 * l; m++)
					for (int p = 0; p < np; p++)
						XR.M[g][static_cast<size_t>(l * l + m) * np + p] = std::pow(r[p], l) * Y(m, p);
			}
		}
		occ::Vec s2 = Eigen::Map<const occ::Vec>(XR.s2.data(), XR.s2.size());
		H += XR.R.transpose() * s2.asDiagonal() * XR.R;
	}
	XR.H.compute(H);
	err_checkf(XR.H.info() == Eigen::Success, "-xcw_partition ri: the fit metric is not positive definite", std::cout);
	XR.on = true;
	std::cout << "XCW partition: RI fit, " << na << " auxiliary functions, (mu nu|P) " << (bytes >> 20) << " MB"
		<< (lmax >= 0 ? ", moments restrained to the grid partition through l = " + std::to_string(lmax) : std::string()) << std::endl;
}

cvec structure_factors::ri_sfc(const dMatrix2& D)
{
	const int ng = model_data.ncen, nc = static_cast<int>(XG.inversion_link.size()), n = XR.nbf;
	occ::Vec Dv(static_cast<size_t>(n) * n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) Dv(static_cast<size_t>(i) * n + j) = D(i, j);
	occ::Vec g = 2.0 * (XR.B.transpose() * Dv);
	//the fixed targets are F's constant part, left out of a response to dD as the anomalous part is
	if (XR.lmax >= 0 && XR.frozen && !XR.response)
		g += XR.R.transpose() * Eigen::Map<const occ::Vec>(XR.t0.data(), XR.t0.size());
	else if (XR.lmax >= 0) {
		const vec2& rho = grid_rho(D);
		occ::Vec T = occ::Vec::Zero(XR.R.rows());
		for (int a = 0; a < ng; a++) {
			const int np = static_cast<int>(rho[a].size());
			for (int k = 0; k < XR.n_mom; k++) {
				const double* m = XR.M[a].data() + static_cast<size_t>(k) * np;
				double t = 0;
				for (int p = 0; p < np; p++) t += m[p] * rho[a][p];
				T(static_cast<size_t>(a) * XR.n_mom + k) = XR.s2[static_cast<size_t>(a) * XR.n_mom + k] * t;
			}
		}
		g += XR.R.transpose() * T;
	}
	const occ::Vec c = XR.H.solve(g);
	XR.c.assign(c.data(), c.data() + c.size());
	cvec sfc(static_cast<size_t>(ng) * nc);
#pragma omp parallel for schedule(dynamic)
	for (int a = 0; a < ng; a++) {
		const int p0 = XR.atom_start[a], np = XR.atom_start[a + 1] - p0;
		for (int e = 0; e < nc; e++) {
			cdouble s = 0;
			for (int j = 0; j < np; j++) s += XR.F[a][static_cast<size_t>(j) * nc + e] * c(p0 + j);
			sfc[static_cast<size_t>(a) * nc + e] = s;
		}
	}
	return sfc;
}

void structure_factors::ri_contract(const cvec& H, occ::Mat& out)
{
	const int ng = model_data.ncen, nc = static_cast<int>(XG.inversion_link.size()), n = XR.nbf;
	//u_P = Re Sum_c H(g, c) F_P(c), the transpose of ri_sfc's last step, then back through H^-1 and B
	occ::Vec u(XR.naux);
#pragma omp parallel for schedule(dynamic)
	for (int a = 0; a < ng; a++) {
		const int p0 = XR.atom_start[a], np = XR.atom_start[a + 1] - p0;
		for (int j = 0; j < np; j++) {
			double s = 0;
			for (int e = 0; e < nc; e++) s += std::real(H[static_cast<size_t>(a) * nc + e] * XR.F[a][static_cast<size_t>(j) * nc + e]);
			u(p0 + j) = s;
		}
	}
	const occ::Vec w = XR.H.solve(u);
	const occ::Vec o = XR.B * w;
	out.resize(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) out(i, j) = o(static_cast<size_t>(i) * n + j);
	if (XR.lmax < 0 || XR.frozen)
		return;
	//the restraint targets' share: V_p = w_p Sum_k z_k r^l Y_k, z = S^2 R w, through the grid blocks
	const occ::Vec z = XR.R * w;
	vec2 V(ng);
	for (int a = 0; a < ng; a++) {
		const int np = static_cast<int>(XG.w[a].size());
		V[a].assign(np, 0.0);
		for (int k = 0; k < XR.n_mom; k++) {
			const double zk = XR.s2[static_cast<size_t>(a) * XR.n_mom + k] * z(static_cast<size_t>(a) * XR.n_mom + k);
			if (zk == 0.0) continue;
			const double* m = XR.M[a].data() + static_cast<size_t>(k) * np;
			for (int p = 0; p < np; p++) V[a][p] += zk * m[p];
		}
		for (int p = 0; p < np; p++) V[a][p] *= XG.w[a][p];
	}
	occ::Mat t = occ::Mat::Zero(n, n);
	grid_apply_V(V, t);
	out += t;
}

void structure_factors::set_ri_targets(const vec2& moments)
{
	const int ng = model_data.ncen;
	XR.t0.assign(static_cast<size_t>(ng) * XR.n_mom, 0.0);
	for (int a = 0; a < ng; a++)
		for (int l = 0; l <= XR.lmax; l++)
			for (int m = 0; m <= 2 * l; m++) {
				double t = 0;
				for (int k = 0; k <= 2 * l; k++) t += XR.A[l](m, k) * moments[a][l * l + k];
				const size_t row = static_cast<size_t>(a) * XR.n_mom + l * l + m;
				XR.t0[row] = XR.s2[row] * t;
			}
	XR.frozen = true;
}
