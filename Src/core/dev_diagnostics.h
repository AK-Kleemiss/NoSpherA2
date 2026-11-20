#pragma once
#include "convenience.h"

//Developer diagnostics, each run by one option of the parser (-spherical_aver_hirsh, -draw_orbits, -NNLS_TEST, ...)
double calc_hirsh_grid_averaged_at_r(const WFN &wavy,
	const int i,
	const double &r,
	const int &min_angular = 60,
	const int &max_angular = 1000,
	const bool print = false);
double calc_fukui_averaged_at_r(const WFN &wavy1,
	const WFN &wavy2,
	const double &r,
	const int &min_angular = 60,
	const int &max_angular = 1000,
	const bool print = false);
double subtract_dens_from_gbw(std::filesystem::path &wfn_name_1,
	std::filesystem::path &wfn_name_2, const double &r, const double &resol);
void calc_partition_densities();
cdouble calc_spherically_averaged_at_k(vec2 &d1, vec2 &d2, vec2 &d3, vec2 &dens,
	const double &k);
void spherically_averaged_density(options &opt, const ivec val_els_alpha, const ivec val_els_beta);
void spherical_harmonic_test();
void cube_from_coef_npy(vec &coefs, WFN &wave);
void get1DGridData(WFN &wavy, std::vector<std::shared_ptr<BasisSet>> &aux_basis, const int atom_idx_1, const int atom_idx_2, const int gridpoints = 1000, const double padding = 2.0);
void draw_orbital(const int lambda, const int m, const double resulution = 0.025, const double radius = 3.5);
void test_NNLS();
void write_spherical_atoms();
