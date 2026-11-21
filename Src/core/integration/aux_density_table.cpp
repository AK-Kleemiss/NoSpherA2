#include "pch.h"
#include "aux_density_table.h"
#include "integration_params.h"
#include "constants.h"
#include "atoms.h"
#include "aux_density.h"
#ifdef NOSPHERA2_USE_GPU
#include "aux_density_gpu.h"
#endif

aux_density_table::aux_density_table(const std::vector<atom>& atoms)
{
	n_at = (int)atoms.size();
	cx.resize(n_at), cy.resize(n_at), cz.resize(n_at), r2_max.resize(n_at), Z.resize(n_at), sh_start.resize(n_at + 1);
	for (int a = 0; a < n_at; a++) {
		cx[a] = atoms[a].get_coordinate(0), cy[a] = atoms[a].get_coordinate(1), cz[a] = atoms[a].get_coordinate(2), Z[a] = atoms[a].get_charge() - atoms[a].get_ECP_electrons();
		const std::vector<unsigned int> sc = atoms[a].get_shellcount();
		double alpha_min = DBL_MAX;
		sh_start[a] = n_sh;
		int prim = 0;
		for (int s = 0; s < (int)sc.size(); s++) {
			const int l = atoms[a].get_basis_set_type(prim);
			err_checkf(l <= 8, "Aux basis shells above l = 8 are not supported on the grid", std::cout);
			sh_l.push_back(l), pr_start.push_back(n_pr), coef_off.push_back(n_coef), sh_atom.push_back(a);

			for (int m = -l; m <= l; ++m)
			{
				coef_shell.push_back(n_sh);
				coef_m.push_back(m);
			}

			vec exponents(sc[s]);
			vec coefficients(sc[s]);
			for (int p = 0; p < sc[s]; ++p, prim++)
			{
				const primitive& pr = atoms[a].get_basis_set_entry(prim).get_primitive();
				exponents[p] = pr.get_exp();
				coefficients[p] = pr.get_coef();

				pr_exp_l32.push_back(std::pow(exponents[p], l + 1.5));
				alpha_min = std::min(alpha_min, exponents[p]);
				pr_exp.push_back(exponents[p]);
				int slot = 0;
				while (slot < (int)uniq_exp.size() && !(uniq_exp[slot] == exponents[p] && uniq_l[slot] == l)) slot++;
				if (slot == (int)uniq_exp.size())
					uniq_exp.push_back(exponents[p]), uniq_l.push_back(l), uniq_exp_l32.push_back(pr_exp_l32.back());
				pr_uniq.push_back(slot);
			}
			coefficients = Int_Params::normalize_gto(coefficients, exponents, l);
			{
				// |r^l R| <= N t^(l/2) exp(-a t), t = r^2, N = sum |n_p|, a the smallest exponent; its log h(t) falls beyond
				// t = l / 2a, so bisect there for h = ln 5e-11, half the kernels' 1e-10 so rounding cannot make a skip differ
				double N = 0.0, a_min = DBL_MAX;
				for (int p = 0; p < sc[s]; ++p) N += std::abs(coefficients[p]), a_min = std::min(a_min, exponents[p]);
				const double target = std::log(5E-11);
				auto h = [&](const double t) { return std::log(N) + (l == 0 ? 0.0 : 0.5 * l * std::log(t)) - a_min * t; };
				double lo = l / (2.0 * a_min), hi = std::max(2.0 * lo, 1.0);
				if (h(lo) < target) hi = lo;
				else {
					while (h(hi) >= target) lo = hi, hi *= 2.0;
					for (int it = 0; it < 100; it++) {
						const double mid = 0.5 * (lo + hi);
						if (mid <= lo || mid >= hi) break;
						(h(mid) < target ? hi : lo) = mid;
					}
				}
				sh_r2.push_back(hi);
			}

			pr_norm.insert(
				pr_norm.end(),
				coefficients.begin(),
				coefficients.end());

			n_pr += sc[s], n_coef += 2 * l + 1, n_sh++;
		}
		r2_max[a] = 46.0517 / alpha_min;
	}
	sh_start[n_at] = n_sh, pr_start.push_back(n_pr);

	err_checkf(
		static_cast<int>(sh_l.size()) == n_sh &&
		static_cast<int>(sh_atom.size()) == n_sh &&
		static_cast<int>(coef_off.size()) == n_sh &&
		static_cast<int>(pr_start.size()) == n_sh + 1,
		"Invalid auxiliary shell table",
		std::cout);

	err_checkf(
		static_cast<int>(pr_exp.size()) == n_pr &&
		static_cast<int>(pr_norm.size()) == n_pr &&
		static_cast<int>(pr_exp_l32.size()) == n_pr,
		"Invalid auxiliary primitive table",
		std::cout);

	err_checkf(
		static_cast<int>(coef_shell.size()) == n_coef &&
		static_cast<int>(coef_m.size()) == n_coef,
		"Invalid auxiliary coefficient table",
		std::cout);
}

double aux_density_table::operator()(const double x, const double y, const double z, const double* coefs) const
{
	return aux_density::at(x, y, z, n_at, cx.data(), cy.data(), cz.data(), r2_max.data(), sh_start.data(), sh_l.data(), pr_start.data(), coef_off.data(), pr_exp.data(), pr_norm.data(), coefs, sh_r2.data());
}

double aux_density_table::operator()(const double x, const double y, const double z, const double* coefs, double& gx, double& gy, double& gz) const
{
	return aux_density::at_grad(x, y, z, n_at, cx.data(), cy.data(), cz.data(), r2_max.data(), sh_start.data(), sh_l.data(), pr_start.data(), coef_off.data(), pr_exp.data(), pr_norm.data(), coefs, gx, gy, gz, sh_r2.data());
}

double aux_density_table::operator()(const double x, const double y, const double z, const double* coefs, double& gx, double& gy, double& gz, double& lap) const
{
	return aux_density::at_lap(x, y, z, n_at, cx.data(), cy.data(), cz.data(), r2_max.data(), sh_start.data(), sh_l.data(), pr_start.data(), coef_off.data(), pr_exp.data(), pr_norm.data(), coefs, gx, gy, gz, lap, sh_r2.data());
}

double aux_density_table::operator()(const double x, const double y, const double z, const double* coefs, double& gx, double& gy, double& gz, double* H) const
{
	return aux_density::at_hess(x, y, z, n_at, cx.data(), cy.data(), cz.data(), r2_max.data(), sh_start.data(), sh_l.data(), pr_start.data(), coef_off.data(), pr_exp.data(), pr_norm.data(), coefs, gx, gy, gz, H, sh_r2.data());
}

static inline cdouble apply_i_to_l(
	const int l,
	const double value)
{
	switch (l & 3)
	{
	case 0:
		return cdouble(value, 0.0);

	case 1:
		return cdouble(0.0, value);

	case 2:
		return cdouble(-value, 0.0);

	case 3:
		return cdouble(0.0, -value);
	}

	return constants::cnull;
}

cdouble aux_density_table::fourier_atom(
	const double kx,
	const double ky,
	const double kz,
	const double* coefs,
	const int atom_idx) const
{
	err_checkf(
		atom_idx >= 0 && atom_idx < n_at,
		"Invalid atom index in aux_density_table::fourier_atom",
		std::cout
	);

	const double H2 =
		kx * kx +
		ky * ky +
		kz * kz;

	const double H = std::sqrt(H2);

	double k[4];

	if (H > 0.0) [[likely]]
	{
		k[0] = kx / H;
		k[1] = ky / H;
		k[2] = kz / H;
	}
	else
	{
		k[0] = 0.0;
		k[1] = 0.0;
		k[2] = 1.0;
	}

	k[3] = H;

	cdouble sf = constants::cnull;

	for (int s = sh_start[atom_idx]; s < sh_start[atom_idx + 1]; ++s)
	{
		const int l = sh_l[s];

		if (H == 0.0 && l > 0)
			continue;

		double Hl_over_2l = 1.0;

		for (int i = 0; i < l; ++i)
			Hl_over_2l *= 0.5 * H;

		double radial = 0.0;

		for (int p = pr_start[s]; p < pr_start[s + 1]; ++p)
		{
			radial += pr_norm[p] * Hl_over_2l * std::exp(-H2 / (4.0 * pr_exp[p])) / pr_exp_l32[p];
		}

		const double angular = constants::spherical_harmonic(l, k[0], k[1], k[2], coefs + coef_off[s]);

		sf += apply_i_to_l(l, constants::PI3_2 * radial * angular);
	}

	return sf;
}

double aux_density_table::lap(const double x, const double y, const double z, const double* coefs) const
{
	double gx, gy, gz, lap;
	(*this)(x, y, z, coefs, gx, gy, gz, lap);
	return lap;
}
double aux_density_table::eli(const double x, const double y, const double z, const double* coefs) const
{
	double gx, gy, gz, lap;
	const double rho = (*this)(x, y, z, coefs, gx, gy, gz, lap);
	return aux_density::eli_from_density(rho, gx * gx + gy * gy + gz * gz, lap);
}
double aux_density_table::esp(const double x, const double y, const double z, const double* coefs) const
{
	return aux_density::esp_at(x, y, z, n_at, cx.data(), cy.data(), cz.data(), Z.data(), sh_start.data(), sh_l.data(), pr_start.data(), coef_off.data(), pr_exp.data(), pr_norm.data(), coefs);
}

void calc_aux_density(const aux_density_table& t, const vec& coefficients, const int np, const double* x, const double* y, const double* z, double* rho, double* gx, double* gy, double* gz, double* lap, double* hess)
{
	err_checkf((int)coefficients.size() == t.n_coef, "Coefficient count does not match the auxiliary basis", std::cout);
	err_checkf(hess == nullptr || gx != nullptr, "The Hessian of the fitted density needs the gradient arrays", std::cout);
#ifdef NOSPHERA2_USE_GPU
	if (aux_density_gpu_enabled() && aux_density_gpu_eval(t.n_at, t.cx.data(), t.cy.data(), t.cz.data(), t.r2_max.data(), t.n_sh, t.sh_start.data(), t.sh_l.data(), t.pr_start.data(), t.coef_off.data(), t.n_pr, t.pr_exp.data(), t.pr_norm.data(), t.n_coef, coefficients.data(), np, x, y, z, rho, gx, gy, gz, lap, hess)) {
		static std::atomic<unsigned> announced{ 0 };
		if (constants::first_this_run(announced) && !constants::hide_gpu_notes)
			std::cout << "GPU in use: fitted density on the grid" << std::endl;
		return;
	}
#endif
	if (gx == nullptr) {
#pragma omp parallel for
		for (int p = 0; p < np; p++) rho[p] = t(x[p], y[p], z[p], coefficients.data());
		return;
	}
	if (hess != nullptr) {
#pragma omp parallel for
		for (int p = 0; p < np; p++) {
			double* H = hess + 9 * (size_t)p;
			rho[p] = t(x[p], y[p], z[p], coefficients.data(), gx[p], gy[p], gz[p], H);
			if (lap != nullptr) lap[p] = H[0] + H[4] + H[8];
		}
		return;
	}
	if (lap == nullptr) {
#pragma omp parallel for
		for (int p = 0; p < np; p++) rho[p] = t(x[p], y[p], z[p], coefficients.data(), gx[p], gy[p], gz[p]);
		return;
	}
#pragma omp parallel for
	for (int p = 0; p < np; p++) rho[p] = t(x[p], y[p], z[p], coefficients.data(), gx[p], gy[p], gz[p], lap[p]);
}


double aux_density_table::shell_radial_moment(
	const int shell) const
{
	const int l = sh_l[shell];

	const double prefactor =
		0.5 * std::tgamma(l + 1.5);

	double integral = 0.0;

	for (int p = pr_start[shell];
		p < pr_start[shell + 1];
		++p)
	{
		// pr_exp_l32[p] = alpha^(l + 3/2)
		integral +=
			pr_norm[p]
			* prefactor
			/ pr_exp_l32[p];
	}

	return integral;
}

double aux_density_table::shell_population_integral(
	const int shell) const
{
	err_checkf(
		sh_l[shell] == 0,
		"Population integral requested for non-s auxiliary shell",
		std::cout
	);

	double integral = 0.0;

	for (int p = pr_start[shell];
		p < pr_start[shell + 1];
		++p)
	{
		// For l=0:
		//
		// ∫ exp(-alpha r²) Y00 d³r
		//
		// = pi / (2 alpha^(3/2))
		integral +=
			pr_norm[p]
			* constants::PI
			/ (2.0 * pr_exp_l32[p]);
	}

	return integral;
}

/**
 * Calculates the atomic density for a given list of atoms and coefficients.
 *
 * @param atoms The list of atoms.
 * @param coefs The coefficients used in the calculation.
 * @return The atomic density for each atom.
 */
vec calc_atomic_density(
	const std::vector<atom>& atoms,
	const vec& coefs)
{
	vec atom_elecs(atoms.size(), 0.0);
	int coef_counter = 0;
	for (int a = 0; a < atoms.size(); ++a)
	{
		int prim = 0;
		for (unsigned int shell = 0; shell < atoms[a].get_shellcount().size(); ++shell)
		{
			const int type = atoms[a].get_basis_set_entry(prim).get_type();
			const unsigned int nprim = atoms[a].get_shellcount()[shell];

			// Only s-functions have a non-zero integral
			// over all space.
			if (type != 0)
			{
				coef_counter += 2 * type + 1; prim += nprim;
				continue;
			}

			vec shell_coefs(nprim);
			vec shell_exps(nprim);

			// Collect RAW contraction coefficients and exponents.
			for (unsigned int e = 0; e < nprim; ++e, ++prim)
			{
				const basis_set_entry& bf = atoms[a].get_basis_set_entry(prim);
				shell_coefs[e] = bf.get_coefficient();
				shell_exps[e] = bf.get_exponent();
			}

			//   primitive normalization, contraction normalization
			const vec normalized_coefs = Int_Params::normalize_gto(shell_coefs, shell_exps, 0);

			double radial_integral = 0.0;
			for (unsigned int e = 0; e < nprim; ++e)
			{
				// Integral:
				//
				// ∫ exp(-alpha*r²) Y_00 d³r
				//
				// with Y_00 = 1 / sqrt(4*pi)
				//
				// = pi / (2 * alpha^(3/2))
				radial_integral +=
					normalized_coefs[e] / (2.0 * std::pow(shell_exps[e], 1.5));
			}

			atom_elecs[a] += radial_integral * coefs[coef_counter] * constants::PI;

			++coef_counter;
		}

		atom_elecs[a] += atoms[a].get_ECP_electrons();
	}

	return atom_elecs;
}
