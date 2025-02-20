#pragma once
#include "convenience.h"

struct Shape2D
{
    unsigned long long rows;
    unsigned long long cols;
    Shape2D() : rows(0), cols(0) {}
    Shape2D(unsigned long long rows, unsigned long long cols) : rows(rows), cols(cols) {
      err_checkf(rows >= 0, "cannot have negative size!", std::cout);
      err_checkf(cols >= 0, "cannot have negative cols!", std::cout);
    }
};

struct Shape3D
{
    unsigned long long depth;
    unsigned long long rows;
    unsigned long long cols;
    Shape3D() : depth(0), rows(0), cols(0) {}
    Shape3D(unsigned long long depth, unsigned long long rows, unsigned long long cols) : depth(depth), rows(rows), cols(cols) {
      err_checkf(depth >= 0, "cannot have negative size!", std::cout);
      err_checkf(rows >= 0, "cannot have negative size!", std::cout);
      err_checkf(cols >= 0, "cannot have negative size!", std::cout);
    }
};

struct Shape4D
{
    unsigned long long depth;
    unsigned long long rows;
    unsigned long long cols;
    unsigned long long time;
    Shape4D() : depth(0), rows(0), cols(0), time(0) {}
    Shape4D(unsigned long long depth, unsigned long long rows, unsigned long long cols, unsigned long long time) : depth(depth), rows(rows), cols(cols), time(time){
        err_checkf(depth >= 0, "cannot have negative size!", std::cout);
        err_checkf(rows >= 0, "cannot have negative size!", std::cout);
        err_checkf(cols >= 0, "cannot have negative size!", std::cout);
        err_checkf(time >= 0, "cannot have negative size!", std::cout);
    }
};

template <typename mat_t, typename vec_t, typename Shape_t>
mat_t reshape(vec_t& flatVec, const Shape_t size);

// TRANSPOSES
// 3D MATRIX
template <typename T>
std::vector<std::vector<std::vector<T>>> transpose(const std::vector<std::vector<std::vector<T>>> &originalVec);

// 2D MATRIX
template <class T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>> &mat);
template <typename T>
Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> transpose(const Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>>& mat);

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
std::vector<T> flatten(const std::vector<std::vector<T>>& vec2D);
// Flatten vector 3D
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<std::vector<T>>> &vec3D);

template <typename T1, typename T2>
T1 flatten(const T2& vecND);

// SLICE Operation
std::vector<double> slice(const std::vector<double> &vec, size_t start, size_t length);

// Matrix multiplication



// 2D x 2D MATRIX MULTIPLICATION
//BLAS implementation of matrix multiplication
template <typename T>
T dot(const T& mat1, const T& mat2, bool transp1 = false, bool transp2 = false);

//BLAS dot product from flattend matrices
template <typename T>
Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dot(const std::vector<T>& flatMat1,
                                                                                               const std::vector<T>& flatMat2,
                                                                                               const int& mat1_d0,
                                                                                               const int& mat1_d1,
                                                                                               const int& mat2_d0,
                                                                                               const int& mat2_d1,
                                                                                               bool transp1=false,
                                                                                               bool transp2=false);

//BLAS implementation of matrix multiplication 2D x 2D
template <typename T>
Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dot_BLAS(const std::vector<T>& flatMat1, const std::vector<T>& flatMat2, const int& m, const int& k1, const int& k2, const int& n, bool transp1 = false, bool transp2 = false);
template <typename T>
T dot_BLAS(const T& Mat1, const T& Mat2, const int& m, const int& k1, const int& k2, const int& n, bool transp1, bool transp2);

//Own implementation of matrix multiplication
template <typename T>
T self_dot(const T& mat1, const T& mat2, bool transp1 = false, bool transp2 = false);




// 2D x 1D MATRIX MULTIPLICATION
// Wrapper for BLAS dot product
template <typename T>
std::vector<T> dot(const std::vector<std::vector<T>>& mat, 
                    const std::vector<T>& vec,
                    bool transp = false);
template <typename T, typename T2>
T dot(const T2& mat,
    const T& vec,
    bool transp = false);

//BLAS implementation of matrix multiplication 2D x 1D
template <typename T>
std::vector<T> dot_BLAS(const std::vector<T>& flatMat, const std::vector<T>& vec, const int& m, const int& n, bool transp = false);
template <typename T, typename T2>
T dot_BLAS(const T2& Mat, const T& vec, bool transp = false);

template <typename T>
std::vector<T> self_dot(const std::vector<std::vector<T>>& mat, const std::vector<T>& vec, bool transp1 = false);


// 1D x 1D MATRIX MULTIPLICATION
// Wrapper for BLAS dot product
template <typename T>
T dot(const std::vector<T>& vec1, const std::vector<T>& vec2, bool conjugate = false);

//BLAS implementation of matrix multiplication 1D x 1D
template <typename T>
T dot_BLAS(const std::vector<T>& vec1, const std::vector<T>& vec2, bool conjugate = false);

// Self written matrix multiplication with flat vectors
template <typename T>
std::vector<T> dot(const std::vector<T> &mat, const std::vector<T> &vec, bool transp = false);

template <typename T>
Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> diag_dot(const Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>>& mat, const std::vector<T>& _vec, bool transp1 = false);

// Element-wise exponentiation of a matrix
vec2 elementWiseExponentiation(const vec2 &matrix, double exponent);
dMatrix2 elementWiseExponentiation(dMatrix2& matrix, double exponent);

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



void math_load_BLAS(int num_threads);
void math_unload_BLAS();
