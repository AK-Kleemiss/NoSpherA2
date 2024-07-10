#pragma once
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

#include "SALTED_math.h"

// Define types for simplicity
using cvec4 = std::vector<std::vector<cvec2>>;

vec predict(const WFN& wavy, std::string model_folder);
vec gen_SALTED_densities(const WFN& wave, options opt, time_point& start, time_point& end_SALTED);

