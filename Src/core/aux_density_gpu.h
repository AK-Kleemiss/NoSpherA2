#pragma once

//GPU path for the fitted density on a point set. The Gordon-Kim repulsion in
//interaction_energy evaluates rho_A and rho_B of the RI or SALTED coefficients on a Becke
//grid over the dimer, and that walk dominated a crystal energy job: every point visits
//every shell of the molecule, the points are independent, and the work per point is a
//few thousand exponentials and polynomials with no data shared between threads.
//
//The device runs aux_density::at, or at_grad when gradients are asked for, from
//aux_density.h, the same inline text the CPU loop in calc_density_ML runs, on the
//flattened basis it is handed. Small point sets stay on the host, the launch and the
//copies cost more than the work there.
//
//Returns false if no device is present, the set is too small to pay, or the arrays will
//not fit, and the caller keeps the OpenMP loop.

bool aux_density_gpu_available();
//On by default; -no_gpu_density turns this off
void aux_density_gpu_set_enabled(bool on);
bool aux_density_gpu_enabled();

bool aux_density_gpu_eval(
	int n_at, const double* cx, const double* cy, const double* cz, const double* r2_max,
	int n_sh, const int* sh_start, const int* sh_l, const int* pr_start, const int* coef_off,
	int n_pr, const double* pr_exp, const double* pr_norm,
	int n_coef, const double* coefs,
	int np, const double* x, const double* y, const double* z, double* rho,
	double* gx = nullptr, double* gy = nullptr, double* gz = nullptr);
