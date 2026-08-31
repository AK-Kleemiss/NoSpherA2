#pragma once

#include <cstddef>
#include <cstdint>

//GPU path for the SALTED descriptor combination. equicomb is 73% of predict() on a
//642-atom protein with the v7 model, and it was measured on 21 Aug to be bound by the
//multiply count rather than by gathering, which is the property that makes a port pay.
//
//The device reproduces the CPU walk exactly: for each atom and each shell triple
//(n1, n2, il) it contracts the Wigner-weighted v1 block against the v2 block over the
//surviving m pairs, applies the complex-to-real transform, accumulates the normalisation
//sum over every feature, and writes only the sparsified ones.
//
//Returns false if no device is present or the problem will not fit, and the caller keeps
//the CPU loop.

bool salted_gpu_available();

struct salted_gpu_problem {
	int natoms = 0;
	int nrad1 = 0;            //nspe1 * nrad1, the channel count of v1
	int nrad2 = 0;            //nspe2 * nrad2, the channel count of v2
	int llmax = 0;
	int lam = 0;
	int l21 = 0;              //2 * lam + 1
	int featsize = 0;         //nrad1 * nrad2 * llmax
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
	//sel[ifeat] is the output slot for that shell triple, or -1
	const int* sel = nullptr;

	double* p = nullptr;                 //[natoms * l21 * nfps], the caller's buffer
};

//Runs the whole lambda block. Returns false if the caller should fall back to the CPU.
bool salted_gpu_equicomb(const salted_gpu_problem& prob, int* empty_environments);
