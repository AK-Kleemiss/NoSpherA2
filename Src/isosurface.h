#pragma once
#include "convenience.h"

// Suppose we store each face color as an (R,G,B) triple
typedef std::array<int, 3> RGB;

// A simple triangle struct composed of three vertices
struct Triangle {
    std::array<double, 3> v1, v2, v3;
    RGB colour;
    int colour_index;
    double calc_area() const {
        std::array<double, 3> a = { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
        std::array<double, 3> b = { v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2] };
        return 0.5 * array_length(cross(a, b));
    }
    double calc_inner_volume() const {
        std::array<double, 3> a = { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
        std::array<double, 3> b = { v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2] };
        return a_dot(cross(b,a), v1) / 6.;
    }
};
std::vector<Triangle> marchingCubes(const cube& volumeData, const double isoVal);
bool writeObj(const std::filesystem::path& filename, const std::vector<Triangle>& triangles);
bool writeColourObj(const std::filesystem::path& filename, std::vector<Triangle>& triangles);
bool writeMTL(const std::string& mtlFilename, std::vector<Triangle>& triangles);
RGB get_colour(const Triangle& t, const cube& volumeData, std::array<std::array<int, 3>, 3> Colourcode, double& low_lim, double& high_lim);
RGB get_colour(const Triangle& t, double(*func)(const double&, const double&, const double&, const WFN&), const WFN& wavy, std::array<std::array<int, 3>, 3> Colourcode, double& low_lim, double& high_lim);
double calc_d_i(const double& x, const double& y, const double& z, const WFN& wavy);