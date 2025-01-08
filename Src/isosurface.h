#pragma once
#include "convenience.h"

// A simple triangle struct composed of three vertices
struct Triangle {
    std::array<double, 3> v1, v2, v3;
};
// Suppose we store each face color as an (R,G,B) triple
struct RGB {
    float r, g, b;
};
std::vector<Triangle> marchingCubes(const cube& volumeData, const double isoVal);
bool writeObj(const std::filesystem::path& filename, const std::vector<Triangle>& triangles);
bool writeMTL(const std::string& mtlFilename, const std::vector<RGB>& faceColors);