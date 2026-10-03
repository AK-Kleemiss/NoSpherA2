#pragma once

#include "gpu_api.h"

//Electrostatic potential of a primitive-pair table on a point set, one thread per point running
//the same text as WFN::computeESP. The pair table is what WFN::build_ESP_pairs returns, flattened:
//P is 3 doubles and L 3 ints per pair, off/coef/pc_pow/fn_idx are the per-pair (l,r,s) blocks.
//boys_tab is the table esp_boys_table() hands out, nT rows of stride values for T = row * step.
//Returns false when no device is present or the arrays do not fit, and the caller keeps the
//OpenMP loop. Shares the -no_gpu_density toggle with the other density kernels.

NOSPHERA2_GPU_API_BEGIN

bool esp_gpu_eval(
	int n_at, const double* ax, const double* ay, const double* az, const double* q,
	int npairs, const double* ex_sum, const double* weight, const double* P, const int* L,
	const int* off, const double* coef, const unsigned char* pc_pow, const unsigned char* fn_idx,
	int nT, int stride, double step, const double* boys_tab,
	int np, const double* points, double* out);

NOSPHERA2_GPU_API_END
