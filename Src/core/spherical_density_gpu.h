#pragma once

#include "gpu_api.h"

//Sum of tabulated spherical atomic densities on a cube grid, one thread per grid point, the
//same log-spaced linear lookup as Thakkar::get_interpolated_density and the same
//"0 outside radius_bohr of every atom" rule as evaluate_cube_in_radius. Returns false when no
//device is present, the grid is too small to pay for the copies or the arrays do not fit, and
//the caller keeps the OpenMP loop. Toggled with aux_density_gpu_set_enabled (-no_gpu_density).

NOSPHERA2_GPU_API_BEGIN

//origin[3], vectors[9] row-major as cube::get_pos uses them; out is x-major nx*ny*nz.
//at_tab[a] indexes the tables of atom a; table t occupies [tab_off[t], tab_off[t+1]) of
//r_tab / rho_tab; lincr and start are shared by every table (make_thakkar_interpolators).
bool spherical_density_gpu_eval(
	int nx, int ny, int nz, const double* origin, const double* vectors,
	int n_at, const double* ax, const double* ay, const double* az, const int* at_tab,
	int n_tab, const int* tab_off, const double* r_tab, const double* rho_tab,
	double lincr, double start, double radius_bohr, double* out);

NOSPHERA2_GPU_API_END
