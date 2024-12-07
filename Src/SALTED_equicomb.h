#pragma once

#include "convenience.h"
#include <complex>
#include "SALTED_math.h"

//Sparse implementation if sparsification is enabled
void equicomb(int natoms, int nrad1, int nrad2,
    const cvec4& v1,
    const cvec4& v2,
    const vec& w3j,
    const ivec2& llvec, const int& lam,
    const cvec2& c2r, const int& featsize,
    const int& nfps, const std::vector<int64_t>& vfps,
    vec& p);

//Normal implementation
void equicomb(int natoms, int nrad1, int nrad2,
    cvec4& v1,
    cvec4& v2,
    vec& w3j, int llmax,
    ivec2& llvec, int lam,
    cvec2& c2r, int featsize,
    vec& p);