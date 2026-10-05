#pragma once

#include "gpu_api.h"

//Density gradient (K = 4, WFN::computeGrad) or ELI-D value and gradient (K = 10,
//WFN::computeELIGrad) of a wavefunction at np points, xyz interleaved. Three steps per chunk of
//points: one thread per (point, contracted function) sums that function's primitives into its K
//components, a device GEMM multiplies them into the occupied MOs, and one thread per point
//reduces the MOs with the host's arithmetic. The tables are the ones WFN::field_grad_gpu builds:
//prim_* per primitive grouped by function (ao_start is the CSR offset, primitives in ascending wfn
//order inside a function, so each function sums in the host's order), coef is [nao x nocc]
//row-major over the occupied MOs only. val (K = 10) and rho may be null.
//open uploads the tables once and sizes the buffers for up to max_points per chunk; run evaluates
//any np (chunked) and may be called many times, from one thread at a time; close frees.
//open returns null and run false when no device is present, -no_gpu_density turned it off or the
//arrays do not fit; the caller then keeps its host loop.

NOSPHERA2_GPU_API_BEGIN

void* basin_field_gpu_open(
	int K, int ncen, const double* cxyz, const double* cmin_exp, double exp_cutoff,
	int nao, const int* ao_start, const int* prim_center, const int* prim_l, const double* prim_exp, const double* prim_scale,
	int nocc, const double* coef, const double* occ, int max_points);
bool basin_field_gpu_run(void* ctx, int np, const double* pts, double* val, double* grad, double* rho);
void basin_field_gpu_close(void* ctx);

NOSPHERA2_GPU_API_END
