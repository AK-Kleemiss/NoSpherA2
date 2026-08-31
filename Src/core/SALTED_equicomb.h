#pragma once

#include "convenience.h"
#include "SALTED_utilities.h"
#include <complex>

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
    bool v2_is_conj_of_v1 = false);

//Normal implementation
void equicomb(int natoms, int nrad1, int nrad2,
    const SALTEDDescriptors& v1,
    const SALTEDDescriptors& v2,
    vec& w3j, int llmax,
    ivec2& llvec, int lam,
    cvec2& c2r, int featsize,
    vec& p,
    bool v2_is_conj_of_v1 = false);

//Opt-in GPU path for the descriptor combination; off unless -gpu_salted is given
void equicomb_set_gpu(bool on);
bool equicomb_gpu_enabled();
