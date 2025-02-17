#include "pch.h"
#include "nos_math.h"

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include "lapacke.h" // for LAPACKE_xxx
#include "cblas.h"


template <typename mat_t, typename vec_t, typename Shape_t>
mat_t reshape(vec_t& flatVec, Shape_t size) {
    if constexpr (std::is_same_v<vec_t, std::vector<double>> || std::is_same_v<vec_t, std::vector<cdouble>> || std::is_same_v<vec_t, std::vector<int>>) {
        if constexpr (std::is_same_v<Shape_t, Shape2D>) {
            mat_t result(flatVec.data(), size.rows, size.cols);
            return result;
        }
        else if constexpr (std::is_same_v<Shape_t, Shape3D>) {
            mat_t result(flatVec.data(), size.depth, size.rows, size.cols);
            return result;
        }
        else if constexpr (std::is_same_v<Shape_t, Shape4D>) {
            mat_t result(flatVec.data(), size.depth, size.rows, size.cols, size.time);
            return result;
        }
        else {
            err_checkf(false, "Invalid Shape!", std::cout);
        }
    }
    else if constexpr (std::is_same_v < vec_t, dMatrix1> || std::is_same_v < vec_t, cMatrix1> || std::is_same_v < vec_t, iMatrix1>) {
        if constexpr (std::is_same_v<Shape_t, Shape2D>) {
            mat_t result(flatVec.data_handle(), size.rows, size.cols);
            return result;
        }
        else if constexpr (std::is_same_v<Shape_t, Shape3D>) {
            mat_t result(flatVec.data_handle(), size.depth, size.rows, size.cols);
            return result;
        }
        else if constexpr (std::is_same_v<Shape_t, Shape4D>) {
            mat_t result(flatVec.data_handle(), size.depth, size.rows, size.cols, size.time);
            return result;
        }
        else {
            err_checkf(false, "Invalid Shape!", std::cout);
        }
    }
    else {
        err_checkf(false, "Invalid Types!", std::cout);
    }
}
template Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> reshape(dMatrix1& fmat, Shape2D size);
template Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> reshape(cMatrix1& fmat, Shape2D size);
template Kokkos::mdspan<int, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> reshape(iMatrix1& fmat, Shape2D size);
template Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> reshape(vec& flatVec, Shape2D size);
template Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> reshape(cvec& flatVec, Shape2D size);
template Kokkos::mdspan<int, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> reshape(ivec& flatVec, Shape2D size);
template Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> reshape(vec& flatVec, Shape3D size);
template Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> reshape(cvec& flatVec, Shape3D size);
template Kokkos::mdspan<int, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> reshape(ivec& flatVec, Shape3D size);
template Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> reshape(vec& flatVec, Shape4D size);
template Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> reshape(cvec& flatVec, Shape4D size);
template Kokkos::mdspan<int, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> reshape(ivec& flatVec, Shape4D size);

// Flatten Vectors 2D
template <typename T>
T flatten(const std::vector<T> &vec2D)
{
    T flatVec;
    size_t totalSize = 0;
    for (const auto &row : vec2D)
    {
        totalSize += row.size();
    }
    flatVec.reserve(totalSize);
    for (const T &row : vec2D)
    {
        // Use std::copy to copy the entire row at once
        flatVec.insert(flatVec.end(), row.begin(), row.end());
    }
    return flatVec;
}
template vec flatten(const vec2 &vec2D);
template cvec flatten(const cvec2 &vec2D);
template ivec flatten(const ivec2 &vec2D);

// Flatten Matrix ND
template <typename T1, typename T2>
T1 flatten(const T2& vecND)
{
    auto DH = vecND.data_handle();
    auto size = vecND.size();
    T1 res (DH, size);
    return res;
}
template Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent>> flatten(const dMatrix2& vecND);
template Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent>> flatten(const cMatrix2& vecND);
template Kokkos::mdspan<int, Kokkos::extents<unsigned long long, std::dynamic_extent>> flatten(const iMatrix2& vecND);
template Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent>> flatten(const dMatrix3& vecND);
template Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent>> flatten(const cMatrix3& vecND);
template Kokkos::mdspan<int, Kokkos::extents<unsigned long long, std::dynamic_extent>> flatten(const iMatrix3& vecND);
template Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent>> flatten(const dMatrix4& vecND);
template Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent>> flatten(const cMatrix4& vecND);
template Kokkos::mdspan<int, Kokkos::extents<unsigned long long, std::dynamic_extent>> flatten(const iMatrix4& vecND);


// Flatten Vectors 3D
template <typename T>
T flatten(const std::vector<std::vector<T>> &vec3D)
{
    T flatVec;
    size_t totalSize = 0;
    for (const auto &row : vec3D)
    {
        for (const auto &innerRow : row)
        {
            totalSize += innerRow.size();
        }
    }
    flatVec.reserve(totalSize);
    for (const std::vector<T> &row : vec3D)
    {
        for (const T &innerRow : row)
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
template std::vector<std::vector<float>> self_dot(const std::vector<std::vector<float>> &mat1, const std::vector<std::vector<float>> &mat2, bool transp1, bool transp2);
template vec2 self_dot(const vec2 &mat1, const vec2 &mat2, bool transp1, bool transp2);
template cvec2 self_dot(const cvec2 &mat1, const cvec2 &mat2, bool transp1, bool transp2);

// typedef void (*ExampleFunctionType)(void);

// Fast 2Dx2D dot product using OpenBLAS
template <typename T>
T dot(const T &mat1, const T &mat2, bool transp1, bool transp2)
{
    // if either of the matrices is empty, return a empty matrix
    if (mat1.empty() || mat2.empty())
    {
        return {};
    }
    int m = transp1 ? (int)mat1.extent(1) : (int)mat1.extent(0);
    int k1 = transp1 ? (int)mat1.extent(0) : (int)mat1.extent(1);
    int k2 = transp2 ? (int)mat2.extent(1) : (int)mat2.extent(0);
    int n = transp2 ? (int)mat2.extent(0) : (int)mat2.extent(1);
    // The resulting matrix will have dimensions m x n

    // Check if matrix multiplication is possible
    err_checkf(k1 == k2, "Inner matrix dimensions must agree.", std::cout);

    if (has_BLAS)
    {
        // Flatten input matrices
        return dot_BLAS(mat1, mat2, m, k1, k2, n, transp1, transp2);
    }
    else
    {
        std::cout << "Something went wrong, using dot fallback." << std::endl;
        exit(-1);
        //return self_dot(mat1, mat2, transp1, transp2);
    }
    // return dot_BLAS(flatMat1, flatMat2, m, k1, k2, n, transp1, transp2);
}
//template std::vector<std::vector<float>> dot(const std::vector<std::vector<float>> &mat1, const std::vector<std::vector<float>> &mat2, bool transp1, bool transp2);
template dMatrix2 dot(const dMatrix2& mat1, const dMatrix2& mat2, bool transp1, bool transp2);
template cMatrix2 dot(const cMatrix2& mat1, const cMatrix2& mat2, bool transp1, bool transp2);

// When the matrices are given as flat vectors
template <typename T>
Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dot(const std::vector<T> &flatMat1, const std::vector<T> &flatMat2, const int &mat1_d0, const int &mat1_d1, const int &mat2_d0, const int &mat2_d1, bool transp1, bool transp2)
{
    // Check if flatMat1 and flatMat2 have the correct size
    err_checkf(flatMat1.size() == mat1_d0 * mat1_d1, "flat Matrix 1 has incorrect size", std::cout);
    err_checkf(flatMat2.size() == mat2_d0 * mat2_d1, "flat Matrix 2 has incorrect size", std::cout);

    // if either of the matrices is empty, return a empty matrix
    if (flatMat1.empty() || flatMat2.empty())
    {
        return {};
    }
    int m = transp1 ? mat1_d1 : mat1_d0;
    int k1 = transp1 ? mat1_d0 : mat1_d1;
    int k2 = transp2 ? mat2_d1 : mat2_d0;
    int n = transp2 ? mat2_d0 : mat2_d1;
    // The resulting matrix will have dimensions m x n

    // Check if matrix multiplication is possible
    err_checkf(k1 == k2, "Inner matrix dimensions must agree.", std::cout);

    if (has_BLAS)
    {
        return dot_BLAS(flatMat1, flatMat2, m, k1, k2, n, transp1, transp2);
    }
    else
    {
        std::cout << "Something went wrong, using dot fallback." << std::endl;
        exit(-1);
    }
}
//template std::vector<std::vector<float>> dot(const std::vector<float> &flatMat1, const std::vector<float> &flatMat2, const int &mat1_d0, const int &mat1_d1, const int &mat2_d0, const int &mat2_d1, bool transp1, bool transp2);
template Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dot(const vec &flatMat1, const vec &flatMat2, const int &mat1_d0, const int &mat1_d1, const int &mat2_d0, const int &mat2_d1, bool transp1, bool transp2);
template Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dot(const cvec &flatMat1, const cvec &flatMat2, const int &mat1_d0, const int &mat1_d1, const int &mat2_d0, const int &mat2_d1, bool transp1, bool transp2);

template <typename T>
Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dot_BLAS(const std::vector<T> &flatMat1, const std::vector<T> &flatMat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2)
{
    std::vector<T> result_flat(m * n, 0.0);
    if constexpr (std::is_same_v<T, double>)
    {
        // Call cblas_dgemm
        cblas_dgemm(CblasRowMajor,
                    transp1 ? CblasTrans : CblasNoTrans,
                    transp2 ? CblasTrans : CblasNoTrans,
                    m, n, k1,
                    1.0,
                    flatMat1.data(), transp1 ? m : k1,
                    flatMat2.data(), transp2 ? k2 : n,
                    0.0,
                    result_flat.data(), n);
    }
    else if constexpr (std::is_same_v<T, float>)
    {
        // Call cblas_sgemm
        cblas_sgemm(CblasRowMajor,
                    transp1 ? CblasTrans : CblasNoTrans,
                    transp2 ? CblasTrans : CblasNoTrans,
                    m, n, k1,
                    1.0f,
                    flatMat1.data(), transp1 ? m : k1,
                    flatMat2.data(), transp2 ? k2 : n,
                    0.0f,
                    result_flat.data(), n);
    }
    else if constexpr (std::is_same_v<T, cdouble>)
    {
        cdouble one = cdouble(1.0, 0.0);
        cdouble zero = cdouble(0.0, 0.0);
        cblas_zgemm(CblasRowMajor,
                    transp1 ? CblasTrans : CblasNoTrans,
                    transp2 ? CblasTrans : CblasNoTrans,
                    m, n, k1,
                    &(one),
                    reinterpret_cast<const cdouble *>(flatMat1.data()), transp1 ? m : k1,
                    reinterpret_cast<const cdouble *>(flatMat2.data()), transp2 ? k2 : n,
                    &(zero),
                    reinterpret_cast<cdouble *>(result_flat.data()), n);
    }
    else
    {
        err_not_impl_f("Unsupported data type for matrix multiplication", std::cout);
    }
    Shape2D result_shape({ (unsigned long long)m, (unsigned long long)n });
    return reshape< Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>>>(result_flat, result_shape);
}
//template std::vector<std::vector<float>> dot_BLAS(const std::vector<float> &mat1, const std::vector<float> &mat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2);
template Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dot_BLAS(const vec & flatMat1, const vec & flatMat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2);
template Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dot_BLAS(const cvec & flatMat1, const cvec & flatMat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2);

template <typename T>
T dot_BLAS(const T& Mat1, const T& Mat2, const int& m, const int& k1, const int& k2, const int& n, bool transp1, bool transp2)
{
    using v_t = T::value_type;
    std::vector<v_t> result_flat(m * n, 0.0);
    if constexpr (std::is_same_v<v_t, double>)
    {
        // Call cblas_dgemm
        cblas_dgemm(CblasRowMajor,
            transp1 ? CblasTrans : CblasNoTrans,
            transp2 ? CblasTrans : CblasNoTrans,
            m, n, k1,
            1.0,
            Mat1.data_handle(), transp1 ? m : k1,
            Mat2.data_handle(), transp2 ? k2 : n,
            0.0,
            result_flat.data(), n);
    }
    else if constexpr (std::is_same_v<v_t, float>)
    {
        // Call cblas_sgemm
        cblas_sgemm(CblasRowMajor,
            transp1 ? CblasTrans : CblasNoTrans,
            transp2 ? CblasTrans : CblasNoTrans,
            m, n, k1,
            1.0f,
            Mat1.data_handle(), transp1 ? m : k1,
            Mat2.data_handle(), transp2 ? k2 : n,
            0.0f,
            result_flat.data(), n);
    }
    else if constexpr (std::is_same_v<v_t, cdouble>)
    {
        cdouble one = cdouble(1.0, 0.0);
        cdouble zero = cdouble(0.0, 0.0);
        cblas_zgemm(CblasRowMajor,
            transp1 ? CblasTrans : CblasNoTrans,
            transp2 ? CblasTrans : CblasNoTrans,
            m, n, k1,
            &(one),
            reinterpret_cast<const cdouble*>(Mat1.data_handle()), transp1 ? m : k1,
            reinterpret_cast<const cdouble*>(Mat2.data_handle()), transp2 ? k2 : n,
            &(zero),
            reinterpret_cast<cdouble*>(result_flat.data()), n);
    }
    else
    {
        err_not_impl_f("Unsupported data type for matrix multiplication", std::cout);
    }
    Shape2D result_shape({ (unsigned long long)m, (unsigned long long)n });
    return reshape<T>(result_flat, result_shape);
}
//template std::vector<std::vector<float>> dot_BLAS(const std::vector<float> &mat1, const std::vector<float> &mat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2);
template dMatrix2 dot_BLAS(const dMatrix2& Mat1, const dMatrix2& Mat2, const int& m, const int& k1, const int& k2, const int& n, bool transp1, bool transp2);
template cMatrix2 dot_BLAS(const cMatrix2& Mat1, const cMatrix2& Mat2, const int& m, const int& k1, const int& k2, const int& n, bool transp1, bool transp2);

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
template std::vector<float> self_dot(const std::vector<std::vector<float>> &mat, const std::vector<float> &_vec, bool transp1);
template vec self_dot(const vec2 &mat, const vec &_vec, bool transp1);
template cvec self_dot(const cvec2 &mat, const cvec &_vec, bool transp1);

// 2D x 1D MATRIX MULTIPLICATION
template <typename T>
Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> diag_dot(const Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> &mat, const std::vector<T> &_vec, bool transp1)
{
    Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> matCopy = mat;
    if (transp1)
    {
        matCopy = transpose(mat);
    }

    int mat_rows = static_cast<int>(matCopy.extent(0));
    int mat_cols = static_cast<int>(matCopy.extent(1));
    int vec_size = static_cast<int>(_vec.size());

    // Check if matrix multiplication is possible
    err_checkf(mat_cols == vec_size || mat_rows == vec_size, "Matrix dimensions do not match for multiplication", std::cout);

    std::vector<T> result(mat_rows * mat_cols, 0.0);
    Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> result_m(result.data(), Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>(mat_rows, mat_cols));
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (int i = 0; i < mat_rows; i++)
        {
            for (int j = 0; j < mat_cols; j++)
            {
                result_m[std::array{ i,j }] = matCopy[std::array{ i,j }] * _vec[j];
            }
        }
    }

    return result_m;
}
//template std::vector<std::vector<float>> diag_dot(const std::vector<std::vector<float>> &mat, const std::vector<float> &_vec, bool transp1);
template Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> diag_dot(const Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>>&mat, const vec &_vec, bool transp1);
template Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> diag_dot(const Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>>&mat, const cvec &_vec, bool transp1);

// Base implementation of matrix-vector multiplication
template <typename T>
T dot(const std::vector<T> &mat, const T &vec, bool transp)
{
    int mat_rows = static_cast<int>(mat.size());
    int mat_cols = static_cast<int>(mat[0].size());
    int vec_size = static_cast<int>(vec.size());

    // Check if matrix multiplication is possible
    err_checkf(mat_cols == vec_size, "Matrix dimensions do not match for multiplication", std::cout);

    if (has_BLAS)
    {
        return dot_BLAS<T>(flatten(mat), vec, mat_rows, mat_cols, transp);
    }
    else
    {
        std::cout << "Something went wrong, using dot fallback." << std::endl;
        //return self_dot(mat, vec, transp);
    }
}
template std::vector<float> dot(const std::vector<std::vector<float>> &mat, const std::vector<float> &vec, bool transp);
template vec dot(const vec2 &mat, const vec &vec, bool transp);
template cvec dot(const cvec2 &mat, const cvec &vec, bool transp);

//mat x Vec
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
    else if constexpr (std::is_same_v<T, float>)
    {
        // Call cblas_sgemv
        cblas_sgemv(CblasRowMajor,
                    transp ? CblasTrans : CblasNoTrans,
                    m, n,
                    1.0f,
                    flatMat.data(), transp ? m : n,
                    vec.data(), 1,
                    0.0f,
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
template std::vector<float> dot_BLAS(const std::vector<float> &flatMat, const std::vector<float> &vec, const int &m, const int &n, bool transp);
template vec dot_BLAS(const std::vector<double> &flatMat, const std::vector<double> &vec, const int &m, const int &n, bool transp);
template cvec dot_BLAS(const std::vector<cdouble> &flatMat, const std::vector<cdouble> &vec, const int &m, const int &n, bool transp);

//mat x Vec
template <typename T, typename T2>
T dot_BLAS(const T2& Mat, const T& vec, bool transp)
{
    int n = Mat.extent(0);
    int m = Mat.extent(1);
    using DataType = typename T2::element_type;
    std::vector<DataType> result(transp ? n : m, 0.0);
    if constexpr (std::is_same_v<T, double>)
    {
        // Call cblas_dgemv
        cblas_dgemv(CblasRowMajor,
            transp ? CblasTrans : CblasNoTrans,
            m, n,
            1.0,
            Mat.data_handle(), transp ? m : n,
            vec.data_handle(), 1,
            0.0,
            result.data(), 1);
    }
    else if constexpr (std::is_same_v<T, float>)
    {
        // Call cblas_sgemv
        cblas_sgemv(CblasRowMajor,
            transp ? CblasTrans : CblasNoTrans,
            m, n,
            1.0f,
            Mat.data_handle(), transp ? m : n,
            vec.data_handle(), 1,
            0.0f,
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
            reinterpret_cast<const cdouble*>(Mat.data_handle()), transp ? m : n,
            reinterpret_cast<const cdouble*>(vec.data_handle()), 1,
            &(zero),
            reinterpret_cast<cdouble*>(result.data()), 1);
    }
    else
    {
        err_not_impl_f("Unsupported data type for matrix multiplication", std::cout);
        return {};
    }
    return T(result.data(), Kokkos::extents<unsigned long long, std::dynamic_extent>(result.size()));
}
template dMatrix1 dot_BLAS(const dMatrix2& Mat, const dMatrix1& vec, bool transp);
template cMatrix1 dot_BLAS(const cMatrix2& Mat, const cMatrix1& vec, bool transp);

template <typename T, typename T2>
T dot(const T2& mat, const T& vec, bool transp)
{
    unsigned long long mat_rows = mat.extent(0);
    unsigned long long mat_cols = mat.extent(1);
    unsigned long long vec_size = vec.extent(0);

    // Check if matrix multiplication is possible
    if (!transp)
        err_checkf(mat_cols == vec_size, "Matrix dimensions do not match for multiplication", std::cout);
    else
        err_checkf(mat_rows == vec_size, "Matrix dimensions do not match for multiplication", std::cout);

    if (has_BLAS)
    {
        return dot_BLAS(mat, vec, transp);
    }
    else
    {
        std::cout << "Something went wrong, using dot fallback." << std::endl;
        exit(-1);
        //return self_dot(mat, vec, transp);
    }
}
//template std::vector<float> dot(const std::vector<std::vector<float>>& mat, const std::vector<float>& vec, bool transp);
template dMatrix1 dot(const dMatrix2& mat, const dMatrix1& vec, bool transp);
template cMatrix1 dot(const cMatrix2& mat, const cMatrix1& vec, bool transp);

template <typename T>
T conj(const T &val)
{
    if constexpr (std::is_same_v<T, cdouble>)
    {
        return std::conj(val);
    }
    else
    {
        return val;
    }
}

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
template float self_dot(const std::vector<float> &vec1, const std::vector<float> &vec2, bool conjugate);
template double self_dot(const std::vector<double> &vec1, const std::vector<double> &vec2, bool conjugate);
template std::complex<double> self_dot(const std::vector<std::complex<double>> &vec1, const std::vector<std::complex<double>> &vec2, bool conjugate);

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

    if (has_BLAS)
    {
        return dot_BLAS(vec1, vec2, conjugate);
    }
    else
    {
        return self_dot(vec1, vec2, conjugate);
    }
}
template float dot(const std::vector<float> &vec1, const std::vector<float> &vec2, bool conjugate);
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
    else if constexpr (std::is_same_v<T, float>)
    {
        // Call cblas_sdot
        result = cblas_sdot((int)vec1.size(), vec1.data(), 1, vec2.data(), 1);
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
template float dot_BLAS(const std::vector<float> &vec1, const std::vector<float> &vec2, bool conjugate);
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

// 2D MATRIX
template <class T>
Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> transpose(const Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>>& mat)
{
    int rows = static_cast<int>(mat.extent(0));
    int cols = static_cast<int>(mat.extent(1));
    std::vector<T> result(cols*rows, 0.0);
    Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> result_m(result.data(), Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>(cols, rows));

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result_m[std::array{ j,i }] = mat[std::array{ i,j }];
        }
    }
    return result_m;
}
template dMatrix2 transpose(const dMatrix2& mat);
template cMatrix2 transpose(const cMatrix2& mat);
template iMatrix2 transpose(const iMatrix2& mat);

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

// Reorder 3D Vectors following a given order
template <typename T>
std::vector<std::vector<std::vector<T>>> reorder3D(const std::vector<std::vector<std::vector<T>>> &original)
{
    if (original.empty() || original[0].empty() || original[0][0].empty())
    {
        return {}; // Return an empty vector if the original vector is empty or not properly sized
    }

    size_t size1 = original.size();       // Original first dimension size
    size_t size2 = original[0].size();    // Original second dimension size
    size_t size3 = original[0][0].size(); // Original third dimension size

    // New vector with dimensions rearranged according to (2, 0, 1)
    std::vector<std::vector<std::vector<T>>> transposed(size3, std::vector<std::vector<T>>(size1, std::vector<T>(size2)));

    for (size_t i = 0; i < size1; ++i)
    {
        for (size_t j = 0; j < size2; ++j)
        {
            for (size_t k = 0; k < size3; ++k)
            {
                transposed[k][i][j] = original[i][j][k];
            }
        }
    }

    return transposed;
}
template vec3 reorder3D(const vec3 &original);

// Function to collect rows from a matrix based on a vector of indices
vec2 collectRows(const vec2 &matrix, const ivec &indices)
{
    // If no indices are provided, return empty matrix
    if (indices.empty())
    {
        return {};
    }
    vec2 result;
    for (int index : indices)
    {
        err_checkf(index < matrix.size(), "Index out of range", std::cout);
        result.push_back(matrix[index]);
    }
    return result;
}

// Function to collect rows from a Cube based on a vector of indices
vec3 collectRows(const vec3 &cube, const ivec &indices)
{
    // If no indices are provided, return empty matrix
    if (indices.empty())
    {
        return {};
    }
    vec3 result;
    for (int index : indices)
    {
        err_checkf(index < cube.size(), "Index out of range", std::cout);
        result.push_back(cube[index]);
    }
    return result;
}

// Element-wise exponentiation of a matrix
vec2 elementWiseExponentiation(const vec2 &matrix, double exponent)
{
    vec2 result = matrix; // Copy the original matrix to preserve its dimensions

    for (size_t i = 0; i < matrix.size(); ++i)
    { // Iterate over rows
        for (size_t j = 0; j < matrix[i].size(); ++j)
        {                                                    // Iterate over columns
            result[i][j] = std::pow(matrix[i][j], exponent); // Apply exponentiation
        }
    }

    return result;
}

// Element-wise exponentiation of a matrix
dMatrix2 elementWiseExponentiation(dMatrix2& matrix, double exponent)
{
    vec result(matrix.size(), 0.0);
    dMatrix2 result_m = reshape<dMatrix2>(result, Shape2D({matrix.extent(0), matrix.extent(1)}));

    for (size_t i = 0; i < matrix.extent(0); ++i)
    { // Iterate over rows
        for (size_t j = 0; j < matrix.extent(1); ++j)
        {                                                    // Iterate over columns
            result_m[std::array{ i,j }] = std::pow(matrix[std::array{ i,j }], exponent); // Apply exponentiation
        }
    }

    return result_m;
}

template <typename T>
void compare_matrices(const Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> &A, const std::vector<std::vector<T>> &B)
{
    std::cout << "Matrices have size " << A.extent(0) << "x" << A.extent(1) << std::endl;
    for (int i = 0; i < A.extent(0); i++)
    {
        for (int j = 0; j < A.extent(1); j++)
        {
            auto a = A[std::array{ i,j }];
            auto b = B[i][j];
            if (a != b)
            {
                std::cout << "Values not matching in comparison! " << i << "," << j << std::endl;
                std::cout << a << " != " << b << std::endl;
            }
            //err_checkf(a == b, "values not matching in comparison!", std::cout);
        }
    }
}

template <typename T>
void compare_matrices(const std::vector<std::vector<T>>& A, const Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>>& B)
{
    std::cout << "Matrices have size " << B.extent(0) << "x" << B.extent(1) << std::endl;
    for (int i = 0; i < A.extent(0); i++)
    {
        for (int j = 0; j < A.extent(1); j++)
        {
            auto a = B[std::array{ i,j }];
            auto b = A[i][j];
            if (a != b)
            {
                std::cout << "Values not matching in comparison! " << i << "," << j << std::endl;
                std::cout << a << " != " << b << std::endl;
            }
            //err_checkf(a == b, "values not matching in comparison!", std::cout);
        }
    }
}

template <typename T>
void compare_matrices(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B)
{
    for (int i = 0; i < A.extent(0); i++)
    {
        for (int j = 0; j < A.extent(1); j++)
        {
            assert(A[i][j] == B[i][j]);
        }
    }
}

void _test_openblas()
{

    vec2 A = { {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0} };
    vec2 B = { {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0} };
    // Init Mat A with some values as a 3x3 matrix
    vec fA = flatten<vec>(A);
    // Init Mat B with some values
    vec fB = flatten<vec>(B);
    Shape2D shape = { 3, 3 };
    dMatrix2 matA = reshape<dMatrix2>(fA, shape);
    dMatrix2 matB = reshape<dMatrix2>(fB, shape);
    math_load_BLAS(1);
    std::cout << "Testing matrices directly" << std::endl;
    compare_matrices(matA, A);
    compare_matrices(matB, B);

    std::cout << "Testing untransposed matrices" << std::endl;
    // First test regular dot-product
    compare_matrices(dot(matA, matB, false, false), self_dot(A, B));

    std::cout << "Testing transpose A" << std::endl;
    ////Second compare first transpose
    compare_matrices(dot(matA, matB, true, false), self_dot(transpose(A), B));

    std::cout << "Testing transpose B" << std::endl;
    ////Third comparte second transpose
    compare_matrices(dot(matA, matB, false, true), self_dot(A, transpose(B)));

    std::cout << "Testing transpose A and B" << std::endl;
    ////Fourth compare both transposed
    compare_matrices(dot(matA, matB, true, true), self_dot(transpose(A), transpose(B)));

    // Init Complex matrices
    cvec2 C = {{{1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}, {{4.0, 4.0}, {5.0, 5.0}, {6.0, 6.0}}, {{7.0, 7.0}, {8.0, 8.0}, {9.0, 9.0}}};
    cvec2 D = {{{1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}, {{4.0, 4.0}, {5.0, 5.0}, {6.0, 6.0}}, {{7.0, 7.0}, {8.0, 8.0}, {9.0, 9.0}}};
    cvec fC = flatten<cvec>(C);
    cvec fD = flatten<cvec>(D);
    shape = { 3, 3 };
    cMatrix2 matC = reshape<cMatrix2>(fC, shape);
    cMatrix2 matD = reshape<cMatrix2>(fD, shape);
    std::cout << "Testing C-matrices directly" << std::endl;
    compare_matrices(matC, C);
    compare_matrices(matD, D);

    std::cout << "Testing untransposed C-matrices" << std::endl;
    // First test regular dot-product
    compare_matrices(dot(matC, matD, false, false), self_dot(C, D));

    std::cout << "Testing transpose C" << std::endl;
    ////Second compare first transpose
    compare_matrices(dot(matC, matD, true, false), self_dot(transpose(C), D));

    std::cout << "Testing transpose D" << std::endl;
    ////Third comparte second transpose
    compare_matrices(dot(matC, matD, false, true), self_dot(C, transpose(D)));

    std::cout << "Testing transpose C and D" << std::endl;
    ////Fourth compare both transposed
    compare_matrices(dot(matC, matD, true, true), self_dot(transpose(C), transpose(D)));

    std::cout << "All BLAS tests passed!" << std::endl;
}

NNLSResult nnls(
    std::vector<double> &A, int m, int n,
    std::vector<double> &B,
    int maxiter,
    double tol)
{
    if (has_BLAS)
    {
        // Define output Variables
        vec X(n, 0);
        double RNORM = 0.0;
        int MODE = 0;

        // Check input dimensions
        if (A.size() != m * n)
        {
            std::cerr << "Error: Matrix A has incorrect dimensions in NNLS.\n";
            return NNLSResult{X, RNORM, 1};
        }

        // Define workspace variables
        std::vector<double> AtA(n * n, 0.0); // A^T * A
        std::vector<double> Atb(n, 0.0);     // A^T * b
        std::vector<double> W(n, 0.0);       // Dual vector
        std::vector<double> S(n, 0.0);       // Trial solution
        std::vector<bool> P(n, false);       // Active set (boolean)

        // Compute A^T * A (normal equations matrix)
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    n, n, m, 1.0, A.data(), m, A.data(), m,
                    0.0, AtA.data(), n);

        // Compute A^T * B (normal equations RHS)
        cblas_dgemv(CblasColMajor, CblasTrans, m, n,
                    1.0, A.data(), m, B.data(), 1,
                    0.0, Atb.data(), 1);
        // Set max iterations
        if (maxiter == -1)
            maxiter = 3 * n;
        if (tol == -1)
            tol = 10 * std::max(m, n) * std::numeric_limits<double>::epsilon();

        // Initialize W, S
        W = Atb; // Projected residual W = A^T * B

        int iter = 0;

        // Initialize workspace variables
        std::vector<double> AtA_active(n * n, 0.0);
        std::vector<double> Atb_active(n, 0);
        std::vector<int> ipiv(n);
        while (iter < maxiter)
        {
            // Step B: Find most active coefficient
            int k = -1;
            double maxW = -1e12;
            for (int i = 0; i < n; i++)
            {
                if (!P[i] && W[i] > tol && W[i] > maxW)
                {
                    k = i;
                    maxW = W[i];
                }
            }

            if (k == -1)
                break; // No positive residuals, terminate.

            // Step B.3: Move k to active set
            P[k] = true;

            // Solve least squares for active set (B.4)
            std::vector<int> activeIndices;
            for (int i = 0; i < n; i++)
            {
                if (P[i])
                    activeIndices.push_back(i);
            }

            int activeCount = activeIndices.size();

            // Extract submatrix AtA[P, P] and Atb[P]
            for (int i = 0; i < activeCount; i++)
            {
                int row = activeIndices[i];
                for (int j = 0; j < activeCount; j++)
                {
                    int col = activeIndices[j];
                    AtA_active[i * activeCount + j] = AtA[row * n + col];
                }
                Atb_active[i] = Atb[row];
            }

            // Solve AtA_active * S[P] = Atb_active
            int info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', activeCount, 1, AtA_active.data(), activeCount, Atb_active.data(), activeCount);

            if (info != 0)
            {
                std::cerr << "Warning: Ill-conditioned matrix detected in NNLS.\n";
                std::cerr << "Error Code: " << info << std::endl;
                MODE = 1;
                break;
            }

            // Assign solution to S
            for (int i = 0; i < activeCount; i++)
            {
                S[activeIndices[i]] = Atb_active[i];
            }

            // Step C: Check feasibility
            while (iter < maxiter)
            {
                iter++;
                double minS = 1e12;
                int minIdx = -1;

                for (int i = 0; i < activeCount; i++)
                {
                    if (S[activeIndices[i]] < 0 && S[activeIndices[i]] < minS)
                    {
                        minS = S[activeIndices[i]];
                        minIdx = activeIndices[i];
                    }
                }

                if (minS >= 0)
                    break; // All positive, proceed.

                // Compute alpha to move back in feasible space
                double alpha = 1.0;
                for (int i = 0; i < activeCount; i++)
                {
                    int idx = activeIndices[i];
                    if (S[idx] < 0)
                    {
                        alpha = std::min(alpha, X[idx] / (X[idx] - S[idx]));
                    }
                }

                // Adjust X and remove minIdx from active set
                for (int i = 0; i < n; i++)
                {
                    X[i] = X[i] + alpha * (S[i] - X[i]);
                }
                P[minIdx] = false; // Remove from active set
            }

            // Assign final solution
            for (int i = 0; i < n; i++)
            {
                X[i] = S[i];
            }
            // Compute residual W = Atb - AtA @ X
            cblas_dcopy(n, Atb.data(), 1, W.data(), 1);
            cblas_dgemv(CblasColMajor, CblasNoTrans, n, n,
                        -1.0, AtA.data(), n, X.data(), 1,
                        1.0, W.data(), 1);
        }

        // Compute residual norm ||A * X - B||
        std::vector<double> Ax(m, 0.0);
        cblas_dgemv(CblasColMajor, CblasNoTrans, m, n,
                    1.0, A.data(), m, X.data(), 1,
                    0.0, Ax.data(), 1);
        double sum_sq = 0.0;
        for (int i = 0; i < m; i++)
        {
            sum_sq += (Ax[i] - B[i]) * (Ax[i] - B[i]);
        }
        sum_sq = std::sqrt(sum_sq);
        NNLSResult resy({X, sum_sq, MODE});
        return resy;
    }
    else
    {
        std::cerr << "Error: NNLS requires LAPACKE and CBLAS.\n";
        exit(1);
    }
}

void math_load_BLAS(int num_threads)
{
    if (has_BLAS)
    {
        return;
    }
#ifdef _WIN32
    _putenv_s("OPENBLAS_NUM_THREADS", std::to_string(num_threads).c_str());
#else
    std::string nums = "OPENBLAS_NUM_THREADS=" + std::to_string(num_threads);
    char *env = strdup(nums.c_str());
    putenv(env);
#endif
    has_BLAS = true;
    return;
}

void math_unload_BLAS()
{
    has_BLAS = false;
}