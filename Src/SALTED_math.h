#pragma once
#include <vector>

#include "convenience.h"

struct Shape2D
{
    int rows;
    int cols;
};

struct Shape3D
{
    int depth;
    int rows;
    int cols;
};


// To_2D
template <typename T>
std::vector<std::vector<T>> reshape(std::vector<T> flatVec, Shape2D sizes);
// To_3D
template <typename T>
std::vector<std::vector<std::vector<T>>> reshape(std::vector<T> flatVec, Shape3D sizes);

// TRANSPOSES
// 3D MATRIX
template <typename T>
std::vector<std::vector<std::vector<T>>> transpose(const std::vector<std::vector<std::vector<T>>>& originalVec);

// 2D MATRIX
template <class T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>>& mat);

// Reorder 3D std::vectors following a given order
template <typename T>
std::vector<std::vector<std::vector<T>>> reorder3D(const std::vector<std::vector<std::vector<T>>>& original);

// Function to collect rows from a vec2 based on a vector of indices
vec2 collectRows(const vec2& matrix, const std::vector<int>& indices);

// Function to collect rows from a vec3 based on a vector of indices
vec3 collectRows(const vec3& cube, const std::vector<int>& indices);

//FLATTEN Operation
// Flatten vector 2D
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& vec2D);
// Flatten vector 3D
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<std::vector<T>>>& vec3D);

// SLICE Operation
std::vector<double> slice(const std::vector<double>& vec, size_t start, size_t length);

// Matrix multiplication
// 2D x 2D MATRIX MULTIPLICATION
template <typename T>
std::vector<std::vector<T>> dot(const std::vector<std::vector<T>>& mat1, const std::vector<std::vector<T>>& mat2);

// 2D x 1D MATRIX MULTIPLICATION
template <typename T>
std::vector<T> dot(const std::vector<std::vector<T>>& mat, const std::vector<T>& vec);

//Element-wise exponentiation of a matrix
vec2 elementWiseExponentiation(const vec2& matrix, double exponent);

void test_dot();
