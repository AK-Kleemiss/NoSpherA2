#pragma once

#include "gpu_api.h"
#include <cstddef>
#include <cstdint>

//GPU SALTED equicomb, the CPU walk exactly: per atom and selected shell triple (n1, n2, il) the
//Wigner-weighted v1 block is contracted with v2 over the surviving m pairs and made real. Only the
//nfps sparsified features are built; the per-atom norm comes from the host.

struct salted_gpu_problem {
	int natoms = 0;
	int nrad1 = 0;            //nspe1 * nrad1, the channel count of v1
	int nrad2 = 0;            //nspe2 * nrad2, the channel count of v2
	int llmax = 0;
	int lam = 0;
	int l21 = 0;              //2 * lam + 1
	int shells = 0;           //nrad1 * nrad2 * llmax; vfps entries past it are zero features
	int nfps = 0;             //sparsified output features
	bool v2_is_conj_of_v1 = false;

	//Flat descriptor storage, as SALTEDDescriptors lays it out
	const double* v1_values = nullptr;   //interleaved re,im
	const size_t* v1_offsets = nullptr;
	int v1_nchannels = 0;
	int v1_noff = 0;
	long long v1_len_doubles = 0;   //2 * complex count
	const double* v2_values = nullptr;
	const size_t* v2_offsets = nullptr;
	int v2_nchannels = 0;
	int v2_noff = 0;
	long long v2_len_doubles = 0;

	const double* w3j = nullptr;
	long long w3j_len = 0;
	const int* llvec0 = nullptr;         //[llmax] l1 per shell
	const int* llvec1 = nullptr;         //[llmax] l2 per shell
	//runs[il * l21 + imu] = {im1_begin, im2_begin, count, w_off}, flattened as 4 ints
	const int* runs = nullptr;
	//complex-to-real rows, two nonzeros each: cols[2i], cols[2i+1] and re/im pairs
	const int* c2r_cols = nullptr;
	const double* c2r_re = nullptr;
	const double* c2r_im = nullptr;
	const int* c2r_cnt = nullptr;
	const int* vfps = nullptr;           //[nfps] shell triple (n1*nrad2+n2)*llmax+il per output slot
	const double* normfact = nullptr;    //[natoms] 1/|feature vector| over all shells, 0 if empty

	double* p = nullptr;                 //[natoms * l21 * nfps], the caller's buffer
};

NOSPHERA2_GPU_API_BEGIN

bool salted_gpu_available();
void salted_gpu_clear_cache();

//Whole lambda block; false means fall back to the CPU.
bool salted_gpu_equicomb(const salted_gpu_problem& prob);

NOSPHERA2_GPU_API_END
