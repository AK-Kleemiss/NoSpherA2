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
                            bool transp2=false,
	                        bool BLAS_enabled=false);

//BLAS dot product from flattend matrices
template <typename T>
std::vector<std::vector<T>> dot(const std::vector<T>& faltMat1,
                                const std::vector<T>& flatMat2,
                                const int& mat1_d0,
                                const int& mat1_d1,
                                const int& mat2_d0,
                                const int& mat2_d1,
                                bool transp1=false,
                                bool transp2=false,
                                bool BLAS_enabled=false);

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
                    bool transp,
                    bool BLAS_enabled = false);

//BLAS implementation of matrix multiplication 2D x 1D
template <typename T>
std::vector<T> dot_BLAS(const std::vector<T>& flatMat, const std::vector<T>& vec, const int& m, const int& n, bool transp);

template <typename T>
std::vector<T> self_dot(const std::vector<std::vector<T>>& mat, const std::vector<T>& vec);


// 1D x 1D MATRIX MULTIPLICATION
// Wrapper for BLAS dot product
template <typename T>
T dot(const std::vector<T>& vec1, const std::vector<T>& vec2, bool conjugate, bool BLAS_enabled = false);

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

void* math_load_BLAS(int num_threads);
void math_unload_BLAS(void* _hOpenBlas);
