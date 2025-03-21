#include "pch.h"
#include "nos_math.h"
#include "lapacke.h" // for LAPACKE_xxx
#include "cblas.h"

template <typename mat_t, typename vec_t, typename Shape_t>
mat_t reshape(vec_t& vec_in, const Shape_t size)
{
    using ext_t = typename mat_t::extents_type;
	using ele_t = typename mat_t::element_type;
    ext_t exts;
    if constexpr (std::is_same_v<Shape_t, Shape2D>)
    {
        exts = ext_t(size.rows, size.cols);
    }
    else if constexpr (std::is_same_v<Shape_t, Shape3D>)
    {
        exts = ext_t(size.depth, size.rows, size.cols);
    }
    else if constexpr (std::is_same_v<Shape_t, Shape4D>)
    {
        exts = ext_t(size.depth, size.rows, size.cols, size.time);
    }
    else
    {
        err_checkf(false, "Invalid Shape!", std::cout);
    }

	mat_t result(exts);
    if constexpr (std::is_same_v<vec_t, vec> || std::is_same_v<vec_t, cvec> || std::is_same_v<vec_t, ivec>)
    {
        std::copy(vec_in.begin(), vec_in.end(), result.data());
        return result;
    }
    else if constexpr (std::is_same_v<vec_t, dMatrix1> || std::is_same_v<vec_t, cMatrix1> || std::is_same_v<vec_t, iMatrix1>
        || std::is_same_v<vec_t, dMatrix2> || std::is_same_v<vec_t, cMatrix2> || std::is_same_v<vec_t, iMatrix2>
        || std::is_same_v<vec_t, dMatrix3> || std::is_same_v<vec_t, cMatrix3> || std::is_same_v<vec_t, iMatrix3>
        || std::is_same_v<vec_t, dMatrix4> || std::is_same_v<vec_t, cMatrix4> || std::is_same_v<vec_t, iMatrix4>)
    {
		std::copy(vec_in.data(), vec_in.data() + vec_in.size(), result.data());
        return result;
    }
    else
    {
        err_checkf(false, "Invalid Types!", std::cout);
    }
    mat_t res_0;
	return res_0;
}
template dMatrix2 reshape(dMatrix1 &vec_in, const Shape2D size);
template cMatrix2 reshape(cMatrix1 &vec_in, const Shape2D size);
template iMatrix2 reshape(iMatrix1 &vec_in, const Shape2D size);
template dMatrix2 reshape(dMatrix2& vec_in, const Shape2D size);
template cMatrix2 reshape(cMatrix2& vec_in, const Shape2D size);
template iMatrix2 reshape(iMatrix2& vec_in, const Shape2D size);
template dMatrix2 reshape(dMatrix3& vec_in, const Shape2D size);
template cMatrix2 reshape(cMatrix3& vec_in, const Shape2D size);
template iMatrix2 reshape(iMatrix3& vec_in, const Shape2D size);
template dMatrix2 reshape(dMatrix4& vec_in, const Shape2D size);
template cMatrix2 reshape(cMatrix4& vec_in, const Shape2D size);
template iMatrix2 reshape(iMatrix4& vec_in, const Shape2D size);
template dMatrix3 reshape(dMatrix1& vec_in, const Shape3D size);
template cMatrix3 reshape(cMatrix1& vec_in, const Shape3D size);
template iMatrix3 reshape(iMatrix1& vec_in, const Shape3D size);
template dMatrix3 reshape(dMatrix2& vec_in, const Shape3D size);
template cMatrix3 reshape(cMatrix2& vec_in, const Shape3D size);
template iMatrix3 reshape(iMatrix2& vec_in, const Shape3D size);
template dMatrix3 reshape(dMatrix3& vec_in, const Shape3D size);
template cMatrix3 reshape(cMatrix3& vec_in, const Shape3D size);
template iMatrix3 reshape(iMatrix3& vec_in, const Shape3D size);
template dMatrix3 reshape(dMatrix4& vec_in, const Shape3D size);
template cMatrix3 reshape(cMatrix4& vec_in, const Shape3D size);
template iMatrix3 reshape(iMatrix4& vec_in, const Shape3D size);
template dMatrix4 reshape(dMatrix1& vec_in, const Shape4D size);
template cMatrix4 reshape(cMatrix1& vec_in, const Shape4D size);
template iMatrix4 reshape(iMatrix1& vec_in, const Shape4D size);
template dMatrix4 reshape(dMatrix2& vec_in, const Shape4D size);
template cMatrix4 reshape(cMatrix2& vec_in, const Shape4D size);
template iMatrix4 reshape(iMatrix2& vec_in, const Shape4D size);
template dMatrix4 reshape(dMatrix3& vec_in, const Shape4D size);
template cMatrix4 reshape(cMatrix3& vec_in, const Shape4D size);
template iMatrix4 reshape(iMatrix3& vec_in, const Shape4D size);
template dMatrix4 reshape(dMatrix4& vec_in, const Shape4D size);
template cMatrix4 reshape(cMatrix4& vec_in, const Shape4D size);
template iMatrix4 reshape(iMatrix4& vec_in, const Shape4D size);
template dMatrix2 reshape(vec&      vec_in, const Shape2D size);
template cMatrix2 reshape(cvec&     vec_in, const Shape2D size);
template iMatrix2 reshape(ivec &    vec_in, const Shape2D size);
template dMatrix3 reshape(vec &     vec_in, const Shape3D size);
template cMatrix3 reshape(cvec &    vec_in, const Shape3D size);
template iMatrix3 reshape(ivec &    vec_in, const Shape3D size);
template dMatrix4 reshape(vec &     vec_in, const Shape4D size);
template cMatrix4 reshape(cvec &    vec_in, const Shape4D size);
template iMatrix4 reshape(ivec &    vec_in, const Shape4D size);

// Flatten Matrix ND
template <typename T1, typename T2>
T1 flatten(const T2 &vecND)
{
    auto size = vecND.size();
    T1 res(size);
	std::copy(vecND.data(), vecND.data() + size, res.data());
    return res;
}
template dMatrix1 flatten(const dMatrix2 &vecND);
template cMatrix1 flatten(const cMatrix2 &vecND);
template iMatrix1 flatten(const iMatrix2 &vecND);
template dMatrix1 flatten(const dMatrix3 &vecND);
template cMatrix1 flatten(const cMatrix3 &vecND);
template iMatrix1 flatten(const iMatrix3 &vecND);
template dMatrix1 flatten(const dMatrix4 &vecND);
template cMatrix1 flatten(const cMatrix4 &vecND);
template iMatrix1 flatten(const iMatrix4 &vecND);

// Fast 2Dx2D dot product using OpenBLAS
template <typename T>
T dot(const T &mat1, const T &mat2, bool transp1, bool transp2)
{
    // if either of the matrices is empty, return a empty matrix
    if ((mat1.size() == 0) || (mat2.size() == 0))
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

    // Flatten input matrices
    return dot_BLAS(mat1, mat2, m, k1, k2, n, transp1, transp2);
}
template dMatrix2 dot(const dMatrix2 &mat1, const dMatrix2 &mat2, bool transp1, bool transp2);
template cMatrix2 dot(const cMatrix2 &mat1, const cMatrix2 &mat2, bool transp1, bool transp2);


template <typename T>
Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dot_BLAS(const std::vector<T> &flatMat1, const std::vector<T> &flatMat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2)
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
    Shape2D result_shape({(unsigned long long)m, (unsigned long long)n});
    return reshape<Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>>>(result_flat, result_shape);
}
template dMatrix2 dot_BLAS(const vec &flatMat1, const vec &flatMat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2);
template cMatrix2 dot_BLAS(const cvec &flatMat1, const cvec &flatMat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2);

template <typename T>
T dot_BLAS(const T &Mat1, const T &Mat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2)
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
                    Mat1.data(), transp1 ? m : k1,
                    Mat2.data(), transp2 ? k2 : n,
                    0.0,
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
                    reinterpret_cast<const cdouble *>(Mat1.data()), transp1 ? m : k1,
                    reinterpret_cast<const cdouble *>(Mat2.data()), transp2 ? k2 : n,
                    &(zero),
                    reinterpret_cast<cdouble *>(result_flat.data()), n);
    }
    else
    {
        err_not_impl_f("Unsupported data type for matrix multiplication", std::cout);
    }
    Shape2D result_shape({(unsigned long long)m, (unsigned long long)n});
    return reshape<T>(result_flat, result_shape);
}
template dMatrix2 dot_BLAS(const dMatrix2 &Mat1, const dMatrix2 &Mat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2);
template cMatrix2 dot_BLAS(const cMatrix2 &Mat1, const cMatrix2 &Mat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2);

// 2D x 1D MATRIX MULTIPLICATION
template <typename T>
Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> diag_dot(const Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> &mat, const std::vector<T> &_vec, bool transp1)
{
    Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> matCopy = mat;
    if (transp1)
    {
        matCopy = transpose(mat);
    }

    int mat_rows = static_cast<int>(matCopy.extent(0));
    int mat_cols = static_cast<int>(matCopy.extent(1));
    int vec_size = static_cast<int>(_vec.size());

    // Check if matrix multiplication is possible
    err_checkf(mat_cols == vec_size || mat_rows == vec_size, "Matrix dimensions do not match for multiplication", std::cout);

    Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> result(mat_rows, mat_cols);
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (int i = 0; i < mat_rows; i++)
        {
            for (int j = 0; j < mat_cols; j++)
            {
                result(i, j) = matCopy(i, j) * _vec[j];
            }
        }
    }

    return result;
}
template dMatrix2 diag_dot(const dMatrix2 &mat, const vec &_vec, bool transp1);
template cMatrix2 diag_dot(const cMatrix2 &mat, const cvec &_vec, bool transp1);

// 2D MATRIX
template <class T>
Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> transpose(const Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> &mat)
{
    int rows = static_cast<int>(mat.extent(0));
    int cols = static_cast<int>(mat.extent(1));
    Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> result(cols, rows);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result(j, i) = mat(i, j);
        }
    }
    return result;
}
template dMatrix2 transpose(const dMatrix2 &mat);
template cMatrix2 transpose(const cMatrix2 &mat);
template iMatrix2 transpose(const iMatrix2 &mat);

// mat x Vec
template <typename T, typename T2>
T dot_BLAS(const T2 &Mat, const T &vec, bool transp)
{
    int m = (int)Mat.extent(0);
    int n = (int)Mat.extent(1);
    using DataType = typename T2::element_type;
    std::vector<DataType> result(transp ? n : m, 0.0);
    if constexpr (std::is_same_v<DataType, double>)
    {
        // Call cblas_dgemv
        cblas_dgemv(CblasRowMajor,
                    transp ? CblasTrans : CblasNoTrans,
                    m, n,
                    1.0,
                    Mat.data(), transp ? m : n,
                    vec.data(), 1,
                    0.0,
                    result.data(), 1);
    }
    else if constexpr (std::is_same_v<DataType, cdouble>)
    {
        cdouble one = cdouble(1.0, 0.0);
        cdouble zero = cdouble(0.0, 0.0);
        cblas_zgemv(CblasRowMajor,
                    transp ? CblasTrans : CblasNoTrans,
                    m, n,
                    &(one),
                    reinterpret_cast<const cdouble *>(Mat.data()), transp ? m : n,
                    reinterpret_cast<const cdouble *>(vec.data()), 1,
                    &(zero),
                    reinterpret_cast<cdouble *>(result.data()), 1);
    }
    else
    {
        err_not_impl_f("Unsupported data type for matrix multiplication", std::cout);
        return {};
    }
    T result_m(result.size());
	std::move(result.begin(), result.end(), result_m.data());
    return result_m;
}
template dMatrix1 dot_BLAS(const dMatrix2 &Mat, const dMatrix1 &vec, bool transp);
template cMatrix1 dot_BLAS(const cMatrix2 &Mat, const cMatrix1 &vec, bool transp);

template <typename T, typename T2>
T dot(const T2 &mat, const T &vec, bool transp)
{
    unsigned long long mat_rows = mat.extent(0);
    unsigned long long mat_cols = mat.extent(1);
    unsigned long long vec_size = vec.extent(0);

    // Check if matrix multiplication is possible
    if (!transp)
        err_checkf(mat_cols == vec_size, "Matrix dimensions do not match for multiplication", std::cout);
    else
        err_checkf(mat_rows == vec_size, "Matrix dimensions do not match for multiplication", std::cout);

    return dot_BLAS(mat, vec, transp);
}
template dMatrix1 dot(const dMatrix2 &mat, const dMatrix1 &vec, bool transp);
template cMatrix1 dot(const cMatrix2 &mat, const cMatrix1 &vec, bool transp);