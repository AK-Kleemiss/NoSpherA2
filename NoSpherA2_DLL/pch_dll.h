// pch_dll.h: This is a precompiled header file.
// Files listed below are compiled only once, improving build performance for future builds.
// This also affects IntelliSense performance, including code completion and many code browsing features.
// However, files listed here are ALL re-compiled if any one of them is updated between builds.
// Do not add files here that you will be updating frequently as this negates the performance advantage.

#pragma once

#include <math.h>
#include <complex>
#include <string>
#include <vector>
#include <array>
#include <filesystem>
#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Windows headers
// Windows Header Files
#define NOMINMAX
#include <windows.h>
#include "../Src/pch.h"
#include "../Src/isosurface.h"
#include "../mdspan/include/mdspan/mdspan.hpp"

#ifdef _WIN32
#ifdef BUILDING_DLL
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT __declspec(dllimport)
#endif
#else
#define DLL_EXPORT
#endif

// Suppose we store each face color as an (R,G,B) triple
typedef std::array<int, 3> RGB;
typedef std::complex<double> cdouble;
typedef std::vector<double> vec;
typedef std::vector<vec> vec2;
typedef std::vector<vec2> vec3;
typedef std::vector<int> ivec;
typedef std::vector<ivec> ivec2;
typedef std::vector<cdouble> cvec;
typedef std::vector<cvec> cvec2;
typedef std::vector<cvec2> cvec3;
typedef std::vector<std::vector<cvec2>> cvec4;
typedef std::vector<bool> bvec;
typedef std::vector<std::string> svec;
typedef std::vector<std::filesystem::path> pathvec;

// A simple triangle struct composed of three vertices
/*class Triangle {
public:
    Triangle(d3 v1, d3 v2, d3 v3, RGB colour, int colour_index);
    Triangle(d3 v1, d3 v2, d3 v3, RGB colour);
    Triangle(d3 v1, d3 v2, d3 v3);

    double calc_area() const;
    double calc_inner_volume() const;
    d3 calc_center() const;
    RGB get_colour() const;
    std::array<double, 3 > get_v(const int& nr) const;
    void set_colour(const RGB& given);
    void set_colour_index(const int& given);
    int get_colour_index() const;
};*/

//std::array<double, 3 > DLL_EXPORT Triangle::get_v(const int& nr) const;
//RGB DLL_EXPORT Triangle::get_colour() const;

std::vector<Triangle> DLL_EXPORT compute_Hirshfeld_suface_i(const std::filesystem::path& fn1, const std::filesystem::path& fn2, double resolution, double radius);
