#pragma once

#include "gpu_api.h"

//The XCW grid runner's two block stages (xcw_grid.cpp) on the device: rho from D and Sum_p chi V chi
//from V, each a per-block fp32 GEMM against AO values uploaded once. The transform between them stays
//on the host. False without a device or room; the caller keeps the CPU loop.

NOSPHERA2_GPU_API_BEGIN

//nb blocks, block i: na[i] AOs aos[ao_off[i] ..], np[i] points from pt_off[i] in the npt concatenated
//atom grids, values X[x_off[i] ..] na x np row-major; keep nmo x nmo, upper triangle: the pair counts
bool xcw_grid_gpu_hold(int nmo, int nb, const int* na, const int* np, const int* pt_off, const int* ao_off, const int* aos,
	const long long* x_off, const float* X, long long nx, const unsigned char* keep, const double* w, int npt);
bool xcw_grid_gpu_held();
//rho[p] = 2 w_p Sum_munu D chi_mu chi_nu, D nmo x nmo row-major read on its upper triangle
bool xcw_grid_gpu_density(const double* D, double* rho);
//out nmo x nmo, upper triangle = Sum_p chi_mu V_p chi_nu over the kept pairs, the rest zero
bool xcw_grid_gpu_contract(const double* V, double* out);
void xcw_grid_gpu_release();

NOSPHERA2_GPU_API_END
