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
 */
NNLSResult nnls(
    std::vector<double>& A, int m, int n,
    std::vector<double>& B,
    int maxiter = -1,
    double tol = -1);

