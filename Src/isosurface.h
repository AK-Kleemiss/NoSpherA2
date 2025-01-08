#pragma once
#include "convenience.h"

// A simple triangle struct composed of three vertices
struct Triangle {
    std::array<double, 3> v1, v2, v3;
};
std::vector<Triangle> marchingCubes(const cube& volumeData);