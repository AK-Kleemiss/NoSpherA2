#pragma once

#include "convenience.h"
#include "SALTED_utilities.h"
#include <complex>

//Per-atom, per-l density matrices A = sum_n x conj(x') and B = sum_n x x' of v1 and of the
//factor u the sparse equicomb multiplies by, for every l the descriptors hold. Its norm
//contracts them; they do not depend on lambda, so a caller looping over lambda builds them once.
struct equicomb_density {
	int natoms = 0, nrad1 = 0, nrad2 = 0;
	bool conj = false;
	int lmax1 = -1, lmax2 = -1;
	//[atom][are, aim, bre, bim][(2l+1)^2 blocks, l = 0..lmax]; m2 is empty when u's are v1's
	//(u = conj(v1) and nrad2 == nrad1), the norm then flips the sign of their imaginary parts
	vec m1, m2;
	bool matches(int na, int n1, int n2, bool c) const
	{
		return natoms == na && nrad1 == n1 && nrad2 == n2 && conj == c;
	}
};
equicomb_density equicomb_density_matrices(int natoms, int nrad1, int nrad2,
	const SALTEDDescriptors& v1, const SALTEDDescriptors& v2, bool v2_is_conj_of_v1 = false);

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
	// Built here when absent or for other dimensions; unused on the GPU, which keeps its own
	const equicomb_density* density = nullptr);

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
