#pragma once

#include <complex>
#include <cstddef>

//GPU path for the XCW I tensor. eval_I spends 99.4% of its time in cblas_dgemm at
//m = n = active AOs, k = the block's grid points, and that shape runs an order of
//magnitude faster in single precision on the device. The AO values do not change with
//the reflection, so they are uploaded once and every reflection reuses them.
//
//Returns false if no device is present or the problem will not fit, and the caller keeps
//the CPU loop. Init must succeed before evaluate is called.

bool itensor_gpu_available();

struct itensor_gpu_layout {
	int nmo = 0;
	int packed = 0;
	int n_grids = 0;             //atom grids carrying the blocks
	int n_blocks = 0;
	const int* blk_grid = nullptr;        //[n_blocks] owning grid
	const int* blk_point_start = nullptr; //[n_blocks] first point within that grid
	const int* blk_point_count = nullptr; //[n_blocks]
	const int* blk_n_active = nullptr;    //[n_blocks]
	const long long* blk_ao_off = nullptr;  //[n_blocks] into ao_all
	const long long* blk_aos_off = nullptr; //[n_blocks] into aos_all
	const float* ao_all = nullptr;        //row-major n_active x point_count per block
	long long ao_all_len = 0;
	const int* aos_all = nullptr;         //AO index per active row, ascending
	long long aos_all_len = 0;
	const unsigned char* skip = nullptr;  //[nmo*nmo], 1 = pair screened out
	const int* grid_point_off = nullptr;  //[n_grids+1] into the flattened point arrays
	const double* d1 = nullptr;           //flattened atom-centred coordinates and weights
	const double* d2 = nullptr;
	const double* d3 = nullptr;
	const double* weights = nullptr;
	long long n_points = 0;
};

bool itensor_gpu_init(const itensor_gpu_layout& L);

//One reflection: phases for each symmetry mate, then every block, accumulated into I_r.
//factors is [num_syms * n_grids], the same product the CPU path applies per (sym, grid).
bool itensor_gpu_reflection(int num_syms,
	const double* kx, const double* ky, const double* kz,
	const std::complex<double>* factors,
	std::complex<double>* I_r);

void itensor_gpu_free();
