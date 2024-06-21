#pragma once
#include <vector>
#include <string>
#include <sstream>

#include "convenience.h"
#include "npy.h"
#include "properties.h"
//#include "equicomb.h"

#ifdef _WIN32
#include "rascaline.hpp"
#include "metatensor.h"
#endif

#include <complex>
#include <cmath>
#include <omp.h>

void predict();

