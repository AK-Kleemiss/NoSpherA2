#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <set>
#include <algorithm>

#include "convenience.h"
#include "wfn_class.h"
#include "npy.h"
#include "properties.h"
#include "JKFit.h"

#if defined(_WIN32) || defined(__RASCALINE__)
#include "H5Cpp.h"
#include "rascaline.hpp"
#include "metatensor.h"
#include <H5public.h>
#include <H5Dpublic.h>
#include <H5Opublic.h>
#endif

#include <complex>
#include <cmath>
#include <omp.h>

vec predict(const WFN& wavy, std::string model_folder);
vec gen_SALTED_densities(const WFN& wave, options opt, time_point& start, time_point& end_SALTED);

