#pragma once

#include <array>
#include <complex>
#include <filesystem>
#include <string>
#include <vector>

#include "api.h"

#include "../core/pch.h"
#include "../core/isosurface.h"

#include "../../mdspan/include/mdspan/mdspan.hpp"

using RGB = std::array<int, 3>;
using cdouble = std::complex<double>;
using vec = std::vector<double>;
using vec2 = std::vector<vec>;
using vec3 = std::vector<vec2>;
using ivec = std::vector<int>;
using ivec2 = std::vector<ivec>;
using cvec = std::vector<cdouble>;
using cvec2 = std::vector<cvec>;
using cvec3 = std::vector<cvec2>;
using cvec4 = std::vector<std::vector<cvec2>>;
using bvec = std::vector<bool>;
using svec = std::vector<std::string>;
using pathvec = std::vector<std::filesystem::path>;

DLL_EXPORT std::vector<Triangle> __cdecl compute_Hirshfeld_suface_i(
    const std::filesystem::path& fn1,
    const std::filesystem::path& fn2,
    double resolution,
    double radius);
