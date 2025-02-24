#pragma once
#include "wfn_class.h"
#include "SALTED_predictor.h"
#include "convenience.h"

void solve_linear_problem(const vec2& A, vec& b);

int fixed_density_fit_test();

vec density_fit(const WFN& wavy, const std::string auxname, const double max_mem, const char metric);