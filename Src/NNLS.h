#pragma once

#include "convenience.h"
// A small struct to hold results
struct NNLSResult {
    std::vector<double> x;   // solution
    double rnorm;            // residual norm
    int status;              // 0 if success, -1 if iteration limit, or other codes
};


NNLSResult nnls(
    std::vector<double>& A, int m, int n,
    std::vector<double>& B,
    int maxiter = -1,
    double tol = -1);

