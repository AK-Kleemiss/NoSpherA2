#pragma once
#include "convenience.h"


int fixed_density_fit_test();

vec density_fit(const WFN& wavy, const WFN& wavy_aux, const double max_mem, const char metric);