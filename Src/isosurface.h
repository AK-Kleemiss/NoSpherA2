#pragma once
#include "convenience.h"

// Suppose we store each face color as an (R,G,B) triple
typedef std::array<int, 3> RGB;

// A simple triangle struct composed of three vertices
struct Triangle {
    std::array<double, 3> v1, v2, v3;
    RGB colour;
    int colour_index;
    double calc_area() const;
    double calc_inner_volume() const;
    std::array<double, 3> calc_center() const;
};
std::vector<Triangle> marchingCubes(const cube& volumeData, const double isoVal, const int subdivisionLevel = 2);
bool writeObj(const std::filesystem::path& filename, const std::vector<Triangle>& triangles);
bool writeColourObj(const std::filesystem::path& filename, std::vector<Triangle>& triangles);
bool writeMTL(const std::string& mtlFilename, std::vector<Triangle>& triangles);
RGB get_colour(const Triangle& t, const cube& volumeData, std::array<std::array<int, 3>, 3> Colourcode, double low_lim, double high_lim);
RGB get_colour(const Triangle& t, double(*func)(const std::array<double, 3>&, const WFN&), const WFN& wavy, std::array<std::array<int, 3>, 3> Colourcode, double low_lim, double high_lim);
double calc_d_i(const std::array<double, 3>& p_t, const WFN& wavy);