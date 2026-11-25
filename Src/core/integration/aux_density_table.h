#pragma once
#include "convenience.h"

//Flat copy of an aux basis for evaluating the fitted density on many points, see aux_density.h
struct aux_density_table
{
	int n_at = 0, n_sh = 0, n_pr = 0, n_coef = 0;

	vec cx, cy, cz, r2_max, pr_exp, pr_norm;
	// per shell, the squared distance beyond which |r^l R| < 1e-10 for certain: the point kernels skip it there unevaluated
	vec sh_r2;
	// alpha^(l + 3/2), needed for Fourier-Bessel transform
	vec pr_exp_l32;
	ivec sh_start, sh_atom, sh_l, pr_start, coef_off;
	// nuclear charges less ECP electrons, for the potential
	ivec Z;
	// distinct (exponent, l) pairs and the slot of each primitive: the radial Fourier factor depends on nothing else
	vec uniq_exp, uniq_exp_l32;
	ivec uniq_l, pr_uniq;

	// coefficient index -> shell / m
	ivec coef_shell;
	ivec coef_m;

	explicit aux_density_table(const std::vector<atom>& atoms);
	double operator()(const double x, const double y, const double z, const double* coefs) const;
	double operator()(const double x, const double y, const double z, const double* coefs, double& gx, double& gy, double& gz) const;
	double operator()(const double x, const double y, const double z, const double* coefs, double& gx, double& gy, double& gz, double& lap) const;
	//rho with its gradient and Hessian (row-major 3x3)
	double operator()(const double x, const double y, const double z, const double* coefs, double& gx, double& gy, double& gz, double* H) const;
	//Electrostatic potential of nuclei and fitted density, see aux_density::esp_at
	double esp(const double x, const double y, const double z, const double* coefs) const;
	double lap(const double x, const double y, const double z, const double* coefs) const;
	double eli(const double x, const double y, const double z, const double* coefs) const;

	// Convenience function for one atom
	cdouble fourier_atom(double kx, double ky, double kz, const double* coefs, int atom_idx) const;

	// Integral:
	//
	// ∫_0^inf R_l(r) r^(l+2) dr
	//
	// where
	//
	// R_l(r) = sum_p pr_norm[p] exp(-alpha_p r²)
	//
	double shell_radial_moment(int shell) const;

	// Full 3D integral of an s-shell including Y_00.
	// Only valid for l = 0.
	double shell_population_integral(int shell) const;
};
//The fitted density of an aux_density_table (RI fit or SALTED prediction, the same path) on np points:
//rho, with the gradient when gx/gy/gz are given, the Laplacian with lap and the row-major 3x3 Hessian
//with hess (9 per point, needs the gradient arrays); GPU kernel when one is present and enabled, else OpenMP
void calc_aux_density(const aux_density_table& t, const vec& coefficients, const int np, const double* x, const double* y, const double* z, double* rho, double* gx = nullptr, double* gy = nullptr, double* gz = nullptr, double* lap = nullptr, double* hess = nullptr);

vec calc_atomic_density(const std::vector<atom> &atoms, const vec &coefs);
