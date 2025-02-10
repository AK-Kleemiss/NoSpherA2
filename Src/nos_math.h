#pragma once
#include "convenience.h"

struct Shape2D
{
    int rows;
    int cols;
    Shape2D() : rows(0), cols(0) {}
    Shape2D(int rows, int cols) : rows(rows), cols(cols) {
      err_checkf(rows >= 0, "cannot have negative size!", std::cout);
      err_checkf(cols >= 0, "cannot have negative cols!", std::cout);
    }
};

struct Shape3D
{
    int depth;
    int rows;
    int cols;
    Shape3D() : depth(0), rows(0), cols(0) {}
    Shape3D(int depth, int rows, int cols) : depth(depth), rows(rows), cols(cols) { 
      err_checkf(depth >= 0, "cannot have negative size!", std::cout);
      err_checkf(rows >= 0, "cannot have negative size!", std::cout);
      err_checkf(cols >= 0, "cannot have negative size!", std::cout);
    }
};

// To_2D
template <typename T>
std::vector<std::vector<T>> reshape(const std::vector<T>& flatVec, Shape2D sizes);
// To_3D
template <typename T>
std::vector<std::vector<std::vector<T>>> reshape(const std::vector<T>& flatVec, Shape3D sizes);

// TRANSPOSES
// 3D MATRIX
template <typename T>
std::vector<std::vector<std::vector<T>>> transpose(const std::vector<std::vector<std::vector<T>>> &originalVec);

// 2D MATRIX
template <class T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>> &mat);

// Flat 2D MATRIX
template <class T>
std::vector<T> transpose(const std::vector<T>& mat, const int rows, const int cols);

// vec -> 2D MATRIX
template <typename T>
std::vector<std::vector<T>> transpose(const std::vector<T>& vector);


// Reorder 3D std::vectors following a given order
template <typename T>
std::vector<std::vector<std::vector<T>>> reorder3D(const std::vector<std::vector<std::vector<T>>> &original);

// Function to collect rows from a vec2 based on a vector of indices
vec2 collectRows(const vec2 &matrix, const ivec &indices);

// Function to collect rows from a vec3 based on a vector of indices
vec3 collectRows(const vec3 &cube, const ivec &indices);

// FLATTEN Operation
//  Flatten vector 2D
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &vec2D);
// Flatten vector 3D
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<std::vector<T>>> &vec3D);

// SLICE Operation
std::vector<double> slice(const std::vector<double> &vec, size_t start, size_t length);

// Matrix multiplication



// 2D x 2D MATRIX MULTIPLICATION
//BLAS implementation of matrix multiplication
template <typename T>
std::vector<std::vector<T>> dot(const std::vector<std::vector<T>> &mat1,
                            const std::vector<std::vector<T>> &mat2,
                            bool transp1=false,
                            bool transp2=false);

//BLAS dot product from flattend matrices
template <typename T>
std::vector<std::vector<T>> dot(const std::vector<T>& faltMat1,
                                const std::vector<T>& flatMat2,
                                const int& mat1_d0,
                                const int& mat1_d1,
                                const int& mat2_d0,
                                const int& mat2_d1,
                                bool transp1=false,
                                bool transp2=false);

//BLAS implementation of matrix multiplication 2D x 2D
template <typename T>
std::vector<std::vector<T>> dot_BLAS(const std::vector<T>& flatMat1, const std::vector<T>& flatMat2, const int& m, const int& k1, const int& k2, const int& n, bool transp1 = false, bool transp2 = false);

//Own implementation of matrix multiplication
template <typename T>
std::vector<std::vector<T>> self_dot(const std::vector<std::vector<T>>& mat1, const std::vector<std::vector<T>>& mat2, bool transp1 = false, bool transp2 = false);




// 2D x 1D MATRIX MULTIPLICATION
// Wrapper for BLAS dot product
template <typename T>
std::vector<T> dot(const std::vector<std::vector<T>>& mat, 
                    const std::vector<T>& vec,
                    bool transp);

//BLAS implementation of matrix multiplication 2D x 1D
template <typename T>
std::vector<T> dot_BLAS(const std::vector<T>& flatMat, const std::vector<T>& vec, const int& m, const int& n, bool transp);

template <typename T>
std::vector<T> self_dot(const std::vector<std::vector<T>>& mat, const std::vector<T>& vec, bool transp1 = false);


// 1D x 1D MATRIX MULTIPLICATION
// Wrapper for BLAS dot product
template <typename T>
T dot(const std::vector<T>& vec1, const std::vector<T>& vec2, bool conjugate);

//BLAS implementation of matrix multiplication 1D x 1D
template <typename T>
T dot_BLAS(const std::vector<T>& vec1, const std::vector<T>& vec2, bool conjugate);

// Self written matrix multiplication with flat vectors
template <typename T>
T self_dot(const std::vector<T>& vec1, const std::vector<T>& vec2, bool conjugate);

template <typename T>
std::vector<std::vector<T>> diag_dot(const std::vector<std::vector<T>>& mat, const std::vector<T>& _vec, bool transp1 = false);

// Element-wise exponentiation of a matrix
vec2 elementWiseExponentiation(const vec2 &matrix, double exponent);

void _test_openblas();



//Small implementation of the non-negative least squares problem
// A small struct to hold results
struct NNLSResult {
    vec x;   // solution
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



void* math_load_BLAS(int num_threads);
void math_unload_BLAS(void* _hOpenBlas);
