// Force-included on MSVC ARM64 only, where Eigen 5 replaces 3.4 (see CMakeLists.txt). occ still
// relies on two things 3.4 had and 5.0 dropped: <cassert> pulled in transitively and Eigen::all.
#pragma once
#include <cassert>
#include <Eigen/Core>
namespace Eigen { using placeholders::all; }
