#pragma once

#include "convenience.h"
// A small struct to hold results
struct NNLSResult {
    std::vector<double> x;   // solution
    double rnorm;            // residual norm
    int status;              // 0 if success, -1 if iteration limit, or other codes
};

/*
 * Solves the non-negative least squares problem:
 *
 * min ||Ax - B||_2 subject to x >= 0
 *
 * Parameters:
 * - A: Matrix A in Col-major order (vector of doubles)
 * - m: Number of rows in matrix A
 * - n: Number of columns in matrix A
 * - B: Vector B (vector of doubles)
 * - maxiter: Maximum number of iterations (default is 3 * n)
 * - tol: Tolerance for convergence (default is 10 * max(m, n) * epsilon)
 *
 * Returns:
 * - NNLSResult: Struct containing the solution vector x, residual norm, and status code
 *
 * Note:
 * - Requires LAPACKE and CBLAS libraries for matrix operations
 * - If LAPACKE or CBLAS is not available, the function will terminate with an error message
 * - Bro, Rasmus and de Jong, Sijmen, "A Fast Non-Negativity-
       Constrained Least Squares Algorithm", Journal Of Chemometrics, 1997,
       :doi:`10.1002/(SICI)1099-128X(199709/10)11:5<393::AID-CEM483>3.0.CO;2-L`
 */
NNLSResult nnls(
    std::vector<double>& A, int m, int n,
    std::vector<double>& B,
    int maxiter = -1,
    double tol = -1);

