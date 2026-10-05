#pragma once

#include "sf_gpu.h"
#include <complex>
#include <cstddef>

//GPU path for the XCW I tensor, a GEMM at m = n = active AOs, k = block grid points. The AO values
//do not depend on the reflection, so they are uploaded once and contracted per reflection batch
//against AO pair-product tables. False without a device or room; the caller keeps the CPU loop.

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
	//Double whatever the device precision; the upload narrows.
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

NOSPHERA2_GPU_API_BEGIN

bool itensor_gpu_available();

//For the log: the GEMM backends differ in the last digits.
const char* itensor_gpu_gemm_name();

//Per reflection and symmetry operation, padding included. Valid after init.
double itensor_gpu_issued_flops();

//Reflections per submit. Valid after init.
int itensor_gpu_batch(int num_syms);

//FP64 runs phase and GEMM in double, anything else in single; Auto is deliberately not honoured.
bool itensor_gpu_init(const itensor_gpu_layout& L, sf_precision prec = sf_precision::FP32,
	bool tensor = false);

//kx..kz: n_refl * num_syms scattering vectors; factors: n_refl * num_syms * n_grids per-grid
//prefactors, reflection major. Two slots, so collect() into n_refl rows of I_r, row_stride apart,
//overlaps the next batch.
bool itensor_gpu_submit(int slot, int n_refl, int num_syms,
	const double* kx, const double* ky, const double* kz,
	const std::complex<double>* factors);
bool itensor_gpu_collect(int slot, std::complex<double>* I_r, long long row_stride);

void itensor_gpu_free();

//Finished tensor, row-major nr x packed, held on the device for the memory-bound SCF walks; false
//when it does not fit. Both contractions accumulate in double whatever the stored precision.
bool itensor_gpu_hold(const std::complex<float>* I, int nr, int packed);
bool itensor_gpu_hold(const std::complex<double>* I, int nr, int packed);
bool itensor_gpu_held();
//F[r] = F0[r] + sum_k I[r][k] w[k]
bool itensor_gpu_rows(const double* w, const std::complex<double>* F0, std::complex<double>* F);
//out[k] = Re sum_r pre[r] I[r][k]
bool itensor_gpu_cols(const std::complex<double>* pre, double* out);
void itensor_gpu_release();

//Packed ERIs in the XCW::store_ERIs layout; J and K come back n x n column-major, scaled as
//XCW::eri_JK, from D in the same layout.
bool eri_gpu_hold(const double* eri, int n, int npair, const int* pa, const int* pb, const int* first, const int* idx);
bool eri_gpu_JK(const double* D, double* J, double* K);
void eri_gpu_release();

NOSPHERA2_GPU_API_END
