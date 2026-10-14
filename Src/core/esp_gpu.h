#pragma once

#include "gpu_api.h"

//ESP of a WFN::build_ESP_pairs table on a point set, one thread per point, same arithmetic as WFN::computeESP.
//P: 3 doubles, L: 3 ints per pair; off/coef/pc_pow/fn_idx: per-pair (l,r,s) blocks; boys_tab: esp_boys_table(), nT rows of stride, T = row * step.
//false = no device or no room, caller keeps the OpenMP loop; -no_gpu_density switches it off

NOSPHERA2_GPU_API_BEGIN

bool esp_gpu_eval(
	int n_at, const double* ax, const double* ay, const double* az, const double* q,
	int npairs, const double* ex_sum, const double* weight, const double* P, const int* L,
	const int* off, const double* coef, const unsigned char* pc_pow, const unsigned char* fn_idx,
	int nT, int stride, double step, const double* boys_tab,
	int np, const double* points, double* out);

NOSPHERA2_GPU_API_END
