#pragma once
#include "convenience.h"
#include "wfn_class.h"
#include "properties.h"
#include "JKFit.h"

#include "SALTED_math.h"
#include "SALTED_io.h"
#include "SALTED_utilities.h"
#include "SALTED_equicomb.h"

vec predict(const WFN& wavy, std::string model_folder);
vec gen_SALTED_densities(const WFN& wave, options opt, time_point& start, time_point& end_SALTED);

