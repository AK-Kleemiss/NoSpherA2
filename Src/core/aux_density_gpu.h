#pragma once

//Fitted density (and gradient, Laplacian) of a flattened aux basis on a point set, one thread per point running
//the same aux_density::at / at_grad / at_lap text as the CPU loop in calc_density_ML. Returns false when no device
//is present, the set is too small to pay for the copies or the arrays do not fit, and the caller keeps the OpenMP loop.

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
	double* gx = nullptr, double* gy = nullptr, double* gz = nullptr, double* lap = nullptr);
