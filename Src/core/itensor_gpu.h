#pragma once

#include "sf_gpu.h"
#include <complex>
#include <cstddef>

//GPU path for the XCW I tensor. eval_I spends 99.4% of its time in cblas_dgemm at
//m = n = active AOs, k = the block's grid points, and that shape runs an order of
//magnitude faster in single precision on the device. The AO values do not change with
//the reflection, so they are uploaded once, and the device contracts a batch of
//reflections at a time against tables of AO pair products built from them.
//
//Returns false if no device is present or the problem will not fit, and the caller keeps
//the CPU loop. Init must succeed before evaluate is called.

bool itensor_gpu_available();

//Which GEMM the device path will actually use, for the log: the three do not agree in the
//last digits, so a run says which one produced its numbers.
const char* itensor_gpu_gemm_name();

//Flops of the GEMMs the device path issues for one reflection and one symmetry operation,
//padding included, so the throughput row counts what ran. Valid after init.
double itensor_gpu_issued_flops();

//How many reflections one submit may carry, at this many symmetry operations. Valid
//after init.
int itensor_gpu_batch(int num_syms);

struct itensor_gpu_layout {
	int nmo = 0;
	int packed = 0;              //stored pairs per reflection
	int n_grids = 0;             //atom grids carrying the blocks
	int n_blocks = 0;
	const int* blk_grid = nullptr;        //[n_blocks] owning grid
	const int* blk_point_start = nullptr; //[n_blocks] first point within that grid
	const int* blk_point_count = nullptr; //[n_blocks]
	const int* blk_n_active = nullptr;    //[n_blocks]
	const long long* blk_ao_off = nullptr;  //[n_blocks] into ao_all
	const long long* blk_aos_off = nullptr; //[n_blocks] into aos_all
	//Double whichever precision the device runs in: narrowing, when it happens, is the
	//upload's business, so the caller does not have to know what was chosen.
	const double* ao_all = nullptr;       //row-major n_active x point_count per block
	long long ao_all_len = 0;
	const int* aos_all = nullptr;         //AO index per active row, ascending
	long long aos_all_len = 0;
	const int* compact = nullptr;         //[nmo*nmo], stored index of (mu, nu), -1 if screened out
	const int* grid_point_off = nullptr;  //[n_grids+1] into the flattened point arrays
	const double* d1 = nullptr;           //flattened atom-centred coordinates and weights
	const double* d2 = nullptr;
	const double* d3 = nullptr;
	const double* weights = nullptr;
	long long n_points = 0;
};

//FP64 runs the whole device path in double, phase and GEMM alike; anything else runs it in
//single. sf_precision::Auto is deliberately not honoured - see the definition.
bool itensor_gpu_init(const itensor_gpu_layout& L, sf_precision prec = sf_precision::FP32,
	bool tensor = false);

//Submit a batch of reflections to one of two result slots: kx..kz hold n_refl * num_syms
//scattering vectors, factors n_refl * num_syms * n_grids per-grid prefactors, reflection
//major. collect() adds the batch into n_refl rows of I_r, row_stride apart, after a later
//batch has started, so the device-to-host copy overlaps that calculation.
bool itensor_gpu_submit(int slot, int n_refl, int num_syms,
	const double* kx, const double* ky, const double* kz,
	const std::complex<double>* factors);
bool itensor_gpu_collect(int slot, std::complex<double>* I_r, long long row_stride);

void itensor_gpu_free();
