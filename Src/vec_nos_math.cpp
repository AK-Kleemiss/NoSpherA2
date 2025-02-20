#include "pch.h"
#include "nos_math.h"
#include "lapacke.h" // for LAPACKE_xxx
#include "cblas.h"

// Flatten Vectors 2D
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &vec2D)
{
    std::vector<T> flatVec;
    size_t totalSize = 0;
    for (const auto &row : vec2D)
    {
        totalSize += row.size();
    }
    flatVec.reserve(totalSize);
    for (const std::vector<T> &row : vec2D)
    {
        // Use std::copy to copy the entire row at once
        flatVec.insert(flatVec.end(), row.begin(), row.end());
    }
    return flatVec;
}
template vec flatten(const vec2 &vec2D);
template cvec flatten(const cvec2 &vec2D);
template ivec flatten(const ivec2 &vec2D);

// Flatten Vectors 3D
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<std::vector<T>>> &vec3D)
{
    std::vector<T> flatVec;
    size_t totalSize = 0;
    for (const auto &row : vec3D)
    {
        for (const auto &innerRow : row)
        {
            totalSize += innerRow.size();
        }
    }
    flatVec.reserve(totalSize);
    for (const std::vector<std::vector<T>> &row : vec3D)
    {
        for (const std::vector<T> &innerRow : row)
        {
            // Use std::copy to copy the entire innerRow at once
            flatVec.insert(flatVec.end(), innerRow.begin(), innerRow.end());
        }
    }
    return flatVec;
}
template vec flatten(const vec3 &vec3D);
template cvec flatten(const cvec3 &vec3D);
template ivec flatten(const std::vector<ivec2> &vec3D);

// SLICE Operation
std::vector<double> slice(const std::vector<double> &vec, size_t start, size_t length)
{
    err_checkf(start + length < vec.size(), "Slice range is out of bounds.", std::cout);

    std::vector<double> result;
    for (size_t i = start; i < start + length; ++i)
    {
        result.push_back(vec[i]);
    }

    return result;
}

// Matrix multiplication
// 2D x 2D MATRIX MULTIPLICATION
// Slow direct implementation
template <typename T>
T self_dot(const T &mat1, const T &mat2, bool transp1, bool transp2)
{
    // if either of the matrices is empty, return a empty matrix
    if (mat1.empty() || mat2.empty())
    {
        return {};
    }

    T mat1Copy = mat1;
    T mat2Copy = mat2;

    if (transp1)
    {
        mat1Copy = transpose(mat1);
    }
    if (transp2)
    {
        mat2Copy = transpose(mat2);
    }

    int rows1 = (int)mat1Copy.size();
    int cols1 = (int)mat1Copy[0].size();
    int rows2 = (int)mat2Copy.size();
    int cols2 = (int)mat2Copy[0].size();

    // Check if matrix multiplication is possible
    err_checkf(cols1 == rows2, "Matrix dimensions do not match for multiplication", std::cout);

    using v_t = typename T::value_type;

    T result(rows1, v_t(cols2, 0.0));
    const long long int totalIterations = static_cast<long long int>(rows1 * cols2 * cols1);
    size_t total_size = rows1 * cols2;

#pragma omp parallel
    {
        v_t local_result(total_size, 0.0);
        int i, j, k, flatIndex;

#pragma omp for schedule(static) private(i, j, k, flatIndex) nowait
        for (long long int n = 0; n < totalIterations; ++n)
        {
            i = static_cast<int>(n / (cols2 * cols1));
            j = static_cast<int>((n / cols1) % cols2);
            k = static_cast<int>(n % cols1);
            flatIndex = i * cols2 + j;
            local_result[flatIndex] += mat1Copy[i][k] * mat2Copy[k][j];
        }

#pragma omp critical
        {
            for (i = 0; i < rows1; ++i)
            {
                for (j = 0; j < cols2; ++j)
                {
                    flatIndex = i * cols2 + j;
                    result[i][j] += local_result[flatIndex];
                }
            }
        }
    }

    return result;
}
template vec2 self_dot(const vec2 &mat1, const vec2 &mat2, bool transp1, bool transp2);
template cvec2 self_dot(const cvec2 &mat1, const cvec2 &mat2, bool transp1, bool transp2);

// 2D x 1D MATRIX MULTIPLICATION
template <typename T>
std::vector<T> self_dot(const std::vector<std::vector<T>> &mat, const std::vector<T> &_vec, bool transp1)
{
    std::vector<std::vector<T>> matCopy = mat;
    if (transp1)
    {
        matCopy = transpose(mat);
    }

    int mat_rows = static_cast<int>(matCopy.size());
    int mat_cols = static_cast<int>(matCopy[0].size());
    int vec_size = static_cast<int>(_vec.size());

    // Check if matrix multiplication is possible
    if (mat_cols != vec_size)
    {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication");
    }

    std::vector<T> result(mat_rows, 0.0);
    int totalIterations = mat_rows * mat_cols;
#pragma omp parallel
    {
        std::vector<T> local_result(mat_rows, 0.0);
#pragma omp for schedule(static)
        for (int n = 0; n < totalIterations; ++n)
        {
            int i = n / mat_cols;
            int j = n % mat_cols;
            local_result[i] += matCopy[i][j] * _vec[j];
        }
#pragma omp critical
        {
            for (int i = 0; i < mat_rows; ++i)
            {
                result[i] += local_result[i];
            }
        }
    }

    return result;
}
template vec self_dot(const vec2 &mat, const vec &_vec, bool transp1);
template cvec self_dot(const cvec2 &mat, const cvec &_vec, bool transp1);

// Base implementation of matrix-vector multiplication
template <typename T>
std::vector<T> dot(const std::vector<std::vector<T>> &mat, const std::vector<T> &vec, bool transp)
{
    int mat_rows = static_cast<int>(mat.size());
    int mat_cols = static_cast<int>(mat[0].size());
    int vec_size = static_cast<int>(vec.size());

    // Check if matrix multiplication is possible
    err_checkf(mat_cols == vec_size, "Matrix dimensions do not match for multiplication", std::cout);

    return dot_BLAS(flatten(mat), vec, mat_rows, mat_cols, transp);
}
template vec dot(const vec2 &mat, const vec &vec, bool transp);
template cvec dot(const cvec2 &mat, const cvec &vec, bool transp);

// mat x Vec
template <typename T>
std::vector<T> dot_BLAS(const std::vector<T> &flatMat, const std::vector<T> &vec, const int &m, const int &n, bool transp)
{
    std::vector<T> result(transp ? n : m, 0.0);
    if constexpr (std::is_same_v<T, double>)
    {
        // Call cblas_dgemv
        cblas_dgemv(CblasRowMajor,
                    transp ? CblasTrans : CblasNoTrans,
                    m, n,
                    1.0,
                    flatMat.data(), transp ? m : n,
                    vec.data(), 1,
                    0.0,
                    result.data(), 1);
    }
    else if constexpr (std::is_same_v<T, cdouble>)
    {
        cdouble one = cdouble(1.0, 0.0);
        cdouble zero = cdouble(0.0, 0.0);
        cblas_zgemv(CblasRowMajor,
                    transp ? CblasTrans : CblasNoTrans,
                    m, n,
                    &(one),
                    reinterpret_cast<const cdouble *>(flatMat.data()), transp ? m : n,
                    reinterpret_cast<const cdouble *>(vec.data()), 1,
                    &(zero),
                    reinterpret_cast<cdouble *>(result.data()), 1);
    }
    else
    {
        err_not_impl_f("Unsupported data type for matrix multiplication", std::cout);
        return {};
    }
    return result;
}
template vec dot_BLAS(const std::vector<double> &flatMat, const std::vector<double> &vec, const int &m, const int &n, bool transp);
template cvec dot_BLAS(const std::vector<cdouble> &flatMat, const std::vector<cdouble> &vec, const int &m, const int &n, bool transp);

// 1D x 1D Vector multiplication
template <typename T>
T self_dot(const std::vector<T> &vec1, const std::vector<T> &vec2, bool conjugate)
{
    T result{};
    if (conjugate)
    {
        for (size_t i = 0; i < vec1.size(); ++i)
        {
            result += conj(vec1[i]) * vec2[i];
        }
    }
    else
    {
        for (size_t i = 0; i < vec1.size(); ++i)
        {
            result += vec1[i] * vec2[i];
        }
    }
    return result;
}
template double self_dot(const std::vector<double> &vec1, const std::vector<double> &vec2, bool conjugate);
template cdouble self_dot(const std::vector<cdouble> &vec1, const std::vector<cdouble> &vec2, bool conjugate);

template <typename T>
T dot(const std::vector<T> &vec1, const std::vector<T> &vec2, bool conjugate)
{
    // if either of the vectors is empty, return 0
    if (vec1.empty() || vec2.empty())
    {
        return 0;
    }

    size_t size = vec1.size();
    // Check if vector dimensions match
    err_checkf(size == vec2.size(), "Vector dimensions do not match for multiplication", std::cout);

    return dot_BLAS(vec1, vec2, conjugate);
}
template double dot(const std::vector<double> &vec1, const std::vector<double> &vec2, bool conjugate);
template cdouble dot(const std::vector<cdouble> &vec1, const std::vector<cdouble> &vec2, bool conjugate);

template <typename T>
T dot_BLAS(const std::vector<T> &vec1, const std::vector<T> &vec2, bool conjugate)
{
    T result = 0.0;
    if constexpr (std::is_same_v<T, double>)
    {
        result = cblas_ddot((int)vec1.size(), vec1.data(), 1, vec2.data(), 1);
    }
    else if constexpr (std::is_same_v<T, cdouble>)
    {
        const openblas_complex_double t = conjugate ? cblas_zdotu((int)vec1.size(), reinterpret_cast<const cdouble *>(vec1.data()), 1, reinterpret_cast<const cdouble *>(vec2.data()), 1)
                                                    : cblas_zdotc((int)vec1.size(), reinterpret_cast<const cdouble *>(vec1.data()), 1, reinterpret_cast<const cdouble *>(vec2.data()), 1);
#if defined(_WIN32) || defined(__APPLE__)
        result = cdouble(t.real, t.imag);
#else
        result = t;
#endif
    }
    else
    {
        err_not_impl_f("Unsupported data type for vector multiplication", std::cout);
        return {};
    }
    return result;
}
template double dot_BLAS(const std::vector<double> &vec1, const std::vector<double> &vec2, bool conjugate);
template cdouble dot_BLAS(const std::vector<cdouble> &vec1, const std::vector<cdouble> &vec2, bool conjugate);


// TRANSPOSES
// 3D MATRIX
template <typename T>
std::vector<std::vector<std::vector<T>>> transpose(const std::vector<std::vector<std::vector<T>>> &originalVec)
{
    if (originalVec.empty() || originalVec[0].empty() || originalVec[0][0].empty())
    {
        return {}; // Return an empty vector if the original vector is empty or not 3D
    }

    size_t newDim1 = originalVec[0][0].size(); // New first dimension is the old third dimension
    size_t newDim2 = originalVec.size();       // New second dimension is the old first dimension
    size_t newDim3 = originalVec[0].size();    // New third dimension is the old second dimension

    std::vector<std::vector<std::vector<T>>> transposedVec(newDim1, std::vector<std::vector<T>>(newDim2, std::vector<T>(newDim3)));

#pragma omp parallel for
    for (int i = 0; i < int(originalVec.size()); ++i)
    {
        for (size_t j = 0; j < originalVec[i].size(); ++j)
        {
            for (size_t k = 0; k < originalVec[i][j].size(); ++k)
            {
                transposedVec[k][i][j] = originalVec[i][j][k];
            }
        }
    }

    return transposedVec;
}
template vec3 transpose(const vec3 &originalVec);
template cvec3 transpose(const cvec3 &originalVec);
template std::vector<ivec2> transpose(const std::vector<ivec2> &originalVec);

// 2D MATRIX
template <class T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>> &mat)
{
    int rows = static_cast<int>(mat.size());
    int cols = static_cast<int>(mat[0].size());
    std::vector<std::vector<T>> result(cols, std::vector<T>(rows));

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result[j][i] = mat[i][j];
        }
    }
    return result;
}
template vec2 transpose(const vec2 &mat);
template cvec2 transpose(const cvec2 &mat);
template ivec2 transpose(const ivec2 &mat);

// Flat 2D MATRIX
template <class T>
std::vector<T> transpose(const std::vector<T> &mat, const int rows, const int cols)
{
    std::vector<T> result(rows * cols);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result[j * rows + i] = mat[i * cols + j];
        }
    }

    return result;
}
template vec transpose(const vec &mat, const int rows, const int cols);
template cvec transpose(const cvec &mat, const int rows, const int cols);
template ivec transpose(const ivec &mat, const int rows, const int cols);

// vec -> 2D MATRIX
template <class T>
std::vector<std::vector<T>> transpose(const std::vector<T> &vector)
{
    int size = static_cast<int>(vector.size());
    std::vector<std::vector<T>> result(1, std::vector<T>(size));

    for (int i = 0; i < size; ++i)
    {
        result[0][i] = vector[i];
    }
    return result;
}
template vec2 transpose(const vec &vector);
template cvec2 transpose(const cvec &vector);
template ivec2 transpose(const ivec &vector);
