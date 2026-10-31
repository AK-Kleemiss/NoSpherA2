#pragma once

#include "convenience.h"
#include "SALTED_utilities.h"
#include <complex>

//Per-atom normalisation of the sparse equicomb for every lambda = index into w3j, llvec
//(transposed, [2][llmax]) and c2r, as [lambda][atom]; 0 for an empty environment. One pass
//over the atoms serves all lambda, holding one atom's density matrices per thread at a time.
vec2 equicomb_norms(int natoms, int nrad1, int nrad2,
	const SALTEDDescriptors& v1, const SALTEDDescriptors& v2,
	const std::vector<const vec*>& w3j, const std::vector<ivec2>& llvec, const std::vector<cvec2>& c2r,
	bool v2_is_conj_of_v1 = false);

//Sparse implementation if sparsification is enabled
void equicomb(int natoms, int nrad1, int nrad2,
	const SALTEDDescriptors& v1,
	const SALTEDDescriptors& v2,
	const vec& w3j,
	const ivec2& llvec, const int& lam,
	const cvec2& c2r, const int& featsize,
	const int& nfps, const std::vector<int64_t>& vfps,
	vec& p,
	// v2 identical to conj(v1): read v1 and flip the sign instead of
	// holding a second copy of the same gigabytes
	bool v2_is_conj_of_v1 = false,
	// natoms values from equicomb_norms for this lambda; computed here when absent,
	// unused on the GPU, which computes its own
	const double* norms = nullptr);

//Normal implementation
void equicomb(int natoms, int nrad1, int nrad2,
	const SALTEDDescriptors& v1,
	const SALTEDDescriptors& v2,
	vec& w3j, int llmax,
	ivec2& llvec, int lam,
	cvec2& c2r, int featsize,
	vec& p,
	bool v2_is_conj_of_v1 = false);

//GPU path for the descriptor combination; -no_gpu_salted keeps it on the CPU
void equicomb_set_gpu(bool on);
bool equicomb_gpu_enabled();
