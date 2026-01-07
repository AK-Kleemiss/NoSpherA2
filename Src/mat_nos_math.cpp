#include "pch.h"
#include "nos_math.h"


#if defined(__APPLE__)
// On macOS we are using Accelerate for BLAS/LAPACK
#include <Accelerate/Accelerate.h>
#else
// Linux/Windows with oneMKL
#include <mkl.h>
#endif

template <typename mat_t, typename vec_t, typename Shape_t>
mat_t reshape(const vec_t& vec_in, const Shape_t size)
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
template dMatrix2 reshape(const dMatrix1& vec_in, const Shape2D size);
template cMatrix2 reshape(const cMatrix1& vec_in, const Shape2D size);
template iMatrix2 reshape(const iMatrix1& vec_in, const Shape2D size);
template dMatrix2 reshape(const dMatrix2& vec_in, const Shape2D size);
template cMatrix2 reshape(const cMatrix2& vec_in, const Shape2D size);
template iMatrix2 reshape(const iMatrix2& vec_in, const Shape2D size);
template dMatrix2 reshape(const dMatrix3& vec_in, const Shape2D size);
template cMatrix2 reshape(const cMatrix3& vec_in, const Shape2D size);
template iMatrix2 reshape(const iMatrix3& vec_in, const Shape2D size);
template dMatrix2 reshape(const dMatrix4& vec_in, const Shape2D size);
template cMatrix2 reshape(const cMatrix4& vec_in, const Shape2D size);
template iMatrix2 reshape(const iMatrix4& vec_in, const Shape2D size);
template dMatrix3 reshape(const dMatrix1& vec_in, const Shape3D size);
template cMatrix3 reshape(const cMatrix1& vec_in, const Shape3D size);
template iMatrix3 reshape(const iMatrix1& vec_in, const Shape3D size);
template dMatrix3 reshape(const dMatrix2& vec_in, const Shape3D size);
template cMatrix3 reshape(const cMatrix2& vec_in, const Shape3D size);
template iMatrix3 reshape(const iMatrix2& vec_in, const Shape3D size);
template dMatrix3 reshape(const dMatrix3& vec_in, const Shape3D size);
template cMatrix3 reshape(const cMatrix3& vec_in, const Shape3D size);
template iMatrix3 reshape(const iMatrix3& vec_in, const Shape3D size);
template dMatrix3 reshape(const dMatrix4& vec_in, const Shape3D size);
template cMatrix3 reshape(const cMatrix4& vec_in, const Shape3D size);
template iMatrix3 reshape(const iMatrix4& vec_in, const Shape3D size);
template dMatrix4 reshape(const dMatrix1& vec_in, const Shape4D size);
template cMatrix4 reshape(const cMatrix1& vec_in, const Shape4D size);
template iMatrix4 reshape(const iMatrix1& vec_in, const Shape4D size);
template dMatrix4 reshape(const dMatrix2& vec_in, const Shape4D size);
template cMatrix4 reshape(const cMatrix2& vec_in, const Shape4D size);
template iMatrix4 reshape(const iMatrix2& vec_in, const Shape4D size);
template dMatrix4 reshape(const dMatrix3& vec_in, const Shape4D size);
template cMatrix4 reshape(const cMatrix3& vec_in, const Shape4D size);
template iMatrix4 reshape(const iMatrix3& vec_in, const Shape4D size);
template dMatrix4 reshape(const dMatrix4& vec_in, const Shape4D size);
template cMatrix4 reshape(const cMatrix4& vec_in, const Shape4D size);
template iMatrix4 reshape(const iMatrix4& vec_in, const Shape4D size);
template dMatrix2 reshape(const vec& vec_in, const Shape2D size);
template cMatrix2 reshape(const cvec& vec_in, const Shape2D size);
template iMatrix2 reshape(const ivec& vec_in, const Shape2D size);
template dMatrix3 reshape(const vec& vec_in, const Shape3D size);
template cMatrix3 reshape(const cvec& vec_in, const Shape3D size);
template iMatrix3 reshape(const ivec& vec_in, const Shape3D size);
template dMatrix4 reshape(const vec& vec_in, const Shape4D size);
template cMatrix4 reshape(const cvec& vec_in, const Shape4D size);
template iMatrix4 reshape(const ivec& vec_in, const Shape4D size);

// Flatten Matrix ND
template <typename T1, typename T2>
T1 flatten(const T2& vecND)
{
    auto size = vecND.size();
    T1 res(size);
    std::copy(vecND.data(), vecND.data() + size, res.data());
    return res;
}
template dMatrix1 flatten(const dMatrix2& vecND);
template cMatrix1 flatten(const cMatrix2& vecND);
template iMatrix1 flatten(const iMatrix2& vecND);
template dMatrix1 flatten(const dMatrix3& vecND);
template cMatrix1 flatten(const cMatrix3& vecND);
template iMatrix1 flatten(const iMatrix3& vecND);
template dMatrix1 flatten(const dMatrix4& vecND);
template cMatrix1 flatten(const cMatrix4& vecND);
template iMatrix1 flatten(const iMatrix4& vecND);

// Fast 2Dx2D dot product using OpenBLAS
template <typename T>
T dot(const T& mat1, const T& mat2, bool transp1, bool transp2)
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
template dMatrix2 dot(const dMatrix2& mat1, const dMatrix2& mat2, bool transp1, bool transp2);
template cMatrix2 dot(const cMatrix2& mat1, const cMatrix2& mat2, bool transp1, bool transp2);


template <typename T>
Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dot_BLAS(const std::vector<T>& flatMat1, const std::vector<T>& flatMat2, const int& m, const int& k1, const int& k2, const int& n, bool transp1, bool transp2)
{
    std::vector<T> result_flat((size_t)m * (size_t)n, 0.0);
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
            reinterpret_cast<const cdouble*>(flatMat1.data()), transp1 ? m : k1,
            reinterpret_cast<const cdouble*>(flatMat2.data()), transp2 ? k2 : n,
            &(zero),
            reinterpret_cast<cdouble*>(result_flat.data()), n);
    }
    else
    {
        err_not_impl_f("Unsupported data type for matrix multiplication", std::cout);
    }
    Shape2D result_shape({ (unsigned long long)m, (unsigned long long)n });
    return reshape<Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>>>(result_flat, result_shape);
}
template dMatrix2 dot_BLAS(const vec& flatMat1, const vec& flatMat2, const int& m, const int& k1, const int& k2, const int& n, bool transp1, bool transp2);
template cMatrix2 dot_BLAS(const cvec& flatMat1, const cvec& flatMat2, const int& m, const int& k1, const int& k2, const int& n, bool transp1, bool transp2);

template <typename T>
T dot_BLAS(const T& Mat1, const T& Mat2, const int& m, const int& k1, const int& k2, const int& n, bool transp1, bool transp2)
{
    using v_t = T::value_type;
    std::vector<v_t> result_flat((size_t)m * (size_t)n, 0.0);
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
            reinterpret_cast<const cdouble*>(Mat1.data()), transp1 ? m : k1,
            reinterpret_cast<const cdouble*>(Mat2.data()), transp2 ? k2 : n,
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
template dMatrix2 dot_BLAS(const dMatrix2& Mat1, const dMatrix2& Mat2, const int& m, const int& k1, const int& k2, const int& n, bool transp1, bool transp2);
template cMatrix2 dot_BLAS(const cMatrix2& Mat1, const cMatrix2& Mat2, const int& m, const int& k1, const int& k2, const int& n, bool transp1, bool transp2);

// 2D x 1D MATRIX MULTIPLICATION
template <typename T>
Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> diag_dot(const Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>>& mat, const std::vector<T>& _vec, bool transp1)
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
template dMatrix2 diag_dot(const dMatrix2& mat, const vec& _vec, bool transp1);
template cMatrix2 diag_dot(const cMatrix2& mat, const cvec& _vec, bool transp1);

// 2D MATRIX
template <class T>
Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> transpose(const Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>>& mat)
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
template dMatrix2 transpose(const dMatrix2& mat);
template cMatrix2 transpose(const cMatrix2& mat);
template iMatrix2 transpose(const iMatrix2& mat);

// mat x Vec
template <typename T, typename T2>
T dot_BLAS(const T2& Mat, const T& vec, bool transp)
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
            reinterpret_cast<const cdouble*>(Mat.data()), transp ? m : n,
            reinterpret_cast<const cdouble*>(vec.data()), 1,
            &(zero),
            reinterpret_cast<cdouble*>(result.data()), 1);
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

    return dot_BLAS(mat, vec, transp);
}
template dMatrix1 dot(const dMatrix2& mat, const dMatrix1& vec, bool transp);
template cMatrix1 dot(const cMatrix2& mat, const cMatrix1& vec, bool transp);

template <typename T, typename T2>
T trace_product(const T2& a, const T2& b, const bool transp) {
    const int rows_a = a.extent(0);
    const int cols_a = a.extent(1);
    const int rows_b = b.extent(0);
    const int cols_b = b.extent(1);
    err_checkf(cols_a == rows_b, "Incompatible sizes in trace product", std::cout);
    err_checkf(rows_a == cols_b, "Incompatible sizes in trace product", std::cout);
    err_checkf(cols_a == cols_b, "Incompatible matrix sizes in trace product!", std::cout);
    using Datatype = typename T2::value_type;
    if constexpr (std::is_same_v<Datatype, double>)
    {
        double res = 0;
        double* b_data = (double*)b.data();
        double* a_data = (double*)a.data();
        if (!transp) {
#pragma omp parallel for reduction(+:res) schedule(static)
            for (int i = 0; i < cols_a; i++) {
                double* ina = a_data + (i * rows_a);
                for (int j = 0; j < rows_a; j++) {
                    res += *(ina + j) * b_data[j * cols_a + i];
                }
            }
        }
        else {
#pragma omp parallel for reduction(+:res) schedule(static)
            for (int i = 0; i < rows_a; i++) {
                for (int j = 0; j < cols_a; j++) {
                    res += a(i, j) * b(i, j);
                }
            }
        }
        return res;
    }
    else if constexpr (std::is_same_v<Datatype, cdouble>)
    {
        cdouble res = cdouble(0.0, 0.0);
        if (!transp) {
            for (int i = 0; i < cols_a; i++) {
                cdouble* ina = (cdouble*)a.data() + (i * rows_a);
                for (int j = 0; j < rows_a; j++) {
                    res += *(ina + j) * b.data()[j * cols_a + i];
                }
            }
        }
        else {
            for (int i = 0; i < rows_a; i++) {
                for (int j = 0; j < cols_a; j++) {
                    res += a(i, j) * b(i, j);
                }
            }
        }
        return res;
    }
    else if constexpr (std::is_same_v<Datatype, int>)
    {
        int res = 0;
        int* b_data = (int*)b.data();
        int* a_data = (int*)a.data();
        if (!transp) {
#pragma omp parallel for reduction(+:res) schedule(static)
            for (int i = 0; i < cols_a; i++) {
                int* ina = a_data + (i * rows_a);
                for (int j = 0; j < rows_a; j++) {
                    res += *(ina + j) * b_data[j * cols_a + i];
                }
            }
        }
        else {
#pragma omp parallel for reduction(+:res) schedule(static)
            for (int i = 0; i < rows_a; i++) {
                for (int j = 0; j < cols_a; j++) {
                    res += a(i, j) * b(i, j);
                }
            }
        }
        return res;
    }
}
template int trace_product(const iMatrix2& a, const iMatrix2& b, const bool transp);
template double trace_product(const dMatrix2& a, const dMatrix2& b, const bool transp);
template cdouble trace_product(const cMatrix2& a, const cMatrix2& b, const bool transp);

template <typename T>
T get_rectangle(const T& a, const ivec& rows) {
    //IMPORTANT: rows will be filled up as rows dictated, hence a shuffle is possible
    //const int rows_a = a.extent(0);
    const int cols_a = a.extent(1);
    size_t runny = 0;
    T res(rows.size(), cols_a);
    for (auto row : rows) {
        for (int j = 0; j < cols_a; j++) {
            res(runny, j) = a(row, j);
        }
        runny++;
    }
    return res;
}
template iMatrix2 get_rectangle(const iMatrix2& a, const ivec& rows);
template dMatrix2 get_rectangle(const dMatrix2& a, const ivec& rows);
template cMatrix2 get_rectangle(const cMatrix2& a, const ivec& rows);

template <typename T, typename T2>
void get_submatrix(const T2& full,
    T& sub,
    const ivec& indices) {
    err_checkf(indices.size() > 0, "Atom indices list is empty.", std::cout);
    const int n = static_cast<int>(indices.size());
    err_checkf(full.extent(0) == full.extent(1), "Matrix must be square.", std::cout);
    err_checkf(sub.size() == n * n, "Submatrix has incorrect size.", std::cout);

    for (int i = 0; i < n; ++i) {
        const int global_i = indices[i];
        const int in = i * n;
        for (int j = 0; j < n; ++j) {
            const int global_j = indices[j];
            sub[in + j] = full(global_i, global_j);
        }
    }
}
template void get_submatrix(const iMatrix2& full, ivec& sub, const ivec& indices);
template void get_submatrix(const dMatrix2& full, vec& sub, const ivec& indices);
template void get_submatrix(const cMatrix2& full, cvec& sub, const ivec& indices);

template <typename T, typename T2>
void get_submatrix(const T2& full,
    T& sub,
    const ivec& val_indices,
    const ivec& vec_indices) {
    const int n1 = static_cast<int>(val_indices.size());
    const int n2 = static_cast<int>(vec_indices.size());
    err_checkf(n2 > 0, "Val indices list is empty.", std::cout);
    err_checkf(n2 > 0, "Vec indices list is empty.", std::cout);
    err_checkf(sub.size() == n1 * n2, "Submatrix has incorrect size.", std::cout);

    for (int i = 0; i < n1; ++i) {
        const int global_i = val_indices[i];
        const int in = i * n2;
        for (int j = 0; j < n2; ++j) {
            const int global_j = vec_indices[j];
            sub[in + j] = full(global_i, global_j);
        }
    }
}
template void get_submatrix(const iMatrix2& full, ivec& sub, const ivec& val_indices, const ivec& vec_indices);
template void get_submatrix(const dMatrix2& full, vec& sub, const ivec& val_indices, const ivec& vec_indices);
template void get_submatrix(const cMatrix2& full, cvec& sub, const ivec& val_indices, const ivec& vec_indices);

template<typename T>
bool isSymmetricViaEigenvalues(const T& A, int n, double tol) {
    T A_copy = A;  // d/zgeev destroys input
    T wr(n), wi(n);
    using Datatype = typename T::value_type;
    if constexpr (std::is_same_v<Datatype, double>)
    {
#ifdef __APPLE__
        // Use Accelerate framework on macOS for eigenvalue computation
        __CLPK_integer info = 0;
        LAPACK_DoubleComplex* a_data = reinterpret_cast<LAPACK_DoubleComplex*>(A_copy.data());
        double* wr_data = reinterpret_cast<double*>(wr.data());
        double* wi_data = reinterpret_cast<double*>(wi.data());
        __CLPK_integer n_clpk = static_cast<__CLPK_integer>(n);
        __CLPK_integer lda = static_cast<__CLPK_integer>(n);
        __CLPK_integer ldvl = static_cast<__CLPK_integer>(n);
        __CLPK_integer ldvr = static_cast<__CLPK_integer>(n);
        // Since we don't need left/right eigenvectors, we can pass nullptr
        dgeev_((char*)"N", (char*)"N", &n_clpk, a_data, &lda,
            wr_data, wi_data, nullptr, &ldvl, nullptr, &ldvr,
            nullptr, &info);
#else
        LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'N',
            n, A_copy.data(),
            n, wr.data(),
            wi.data(), nullptr, n,
            nullptr, n);
#endif
    }
    else if constexpr (std::is_same_v<Datatype, cdouble>)
    {
        // int matrix_layout, char jobvl, char jobvr,
        //     lapack_int n, lapack_complex_double* a,
        //     lapack_int lda, lapack_complex_double* w,
        //     lapack_complex_double* vl, lapack_int ldvl,
        //     lapack_complex_double* vr, lapack_int ldvr
        // //ISSUE WITH MKL COMPLEX TYPES I AM CURRENTL Y NOT INTERESTED IN SOLVING!
           //LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'N', 
           //    n, reinterpret_cast<const cdouble*>(A_copy.data()), 
           //    n, reinterpret_cast<const cdouble*>(wr.data()), 
           //    reinterpret_cast<const cdouble*>(wi.data()), n, nullptr, n);
        return false;
    }

    for (int i = 0; i < n; ++i) {
        if (std::abs(wi[i]) > tol) {
            return false;
        }
    }
    return true;
}
template bool isSymmetricViaEigenvalues(const vec& A, int n, double tol);

template <typename T, typename T2>
void get_submatrices(const T2& D_full,
    const T2& S_full,
    T& D_sub,
    T& S_sub,
    const ivec& indices) {
    err_checkf(indices.size() > 0, "Atom indices list is empty.", std::cout);
    const int n = static_cast<int>(indices.size());
    err_checkf(D_full.extent(0) == D_full.extent(1), "Density matrix D must be square.", std::cout);
    err_checkf(S_full.extent(0) == S_full.extent(1), "Overlap matrix S must be square.", std::cout);
    err_checkf(D_full.extent(0) == S_full.extent(0), "Density and Overlap matrices must be of the same size.", std::cout);
    err_checkf(D_sub.size() == n * n, "Density submatrix has incorrect size.", std::cout);
    err_checkf(S_sub.size() == n * n, "Overlap submatrix has incorrect size.", std::cout);

    for (int i = 0; i < n; ++i) {
        const int global_i = indices[i];
        const int in = i * n;
        for (int j = 0; j < n; ++j) {
            const int global_j = indices[j];
            D_sub[in + j] = D_full(global_i, global_j);
            S_sub[in + j] = S_full(global_i, global_j);
        }
    }
}
template void get_submatrices(const dMatrix2& D_full, const dMatrix2& S_full, vec& D_sub, vec& S_sub, const ivec& indices);
//template void get_submatrices(const cMatrix2& D_full, const cMatrix2& S_full, cvec& D_sub, cvec& S_sub, const ivec& indices);

//calculates the Moore-Penrose pseudo-inverse of a matrix A using SVD
dMatrix2 LAPACKE_invert(const dMatrix2& A, const double cutoff) {
    const int m = static_cast<int>(A.extent(0)); // rows
    const int n = static_cast<int>(A.extent(1)); // cols
    const int k = std::min(m, n);

    // 1. Allocate memory for SVD results
    vec S(k);                 // Singular values
    vec U(m * k);             // Left singular vectors (m x k)
    vec Vt(k * n);            // Right singular vectors transposed (k x n)
    vec superb(k - 1);        // Workspace for dgesvd

    // Make a copy of A because dgesvd destroys the input matrix
    vec A_copy = A.container();
    // 2. Compute SVD: A = U * S * Vt
    // use 'S' for jobu/jobvt to compute the "thin" SVD (only the first k columns/rows)
#ifdef __APPLE__
    // Use Accelerate framework on macOS for SVD computation
    __CLPK_integer info = 0;
    __CLPK_integer lwork = -1;
    double work_query = 0.0;

    // Query optimal work size
    dgesvd_((char*)"S", (char*)"S", (__CLPK_integer*)&m, (__CLPK_integer*)&n,
        A_copy.data(), (__CLPK_integer*)&n,
        S.data(),
        U.data(), (__CLPK_integer*)&k,
        Vt.data(), (__CLPK_integer*)&n,
        &work_query, (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);

    lwork = static_cast<__CLPK_integer>(work_query);
    vec work(lwork);

    // Perform SVD: A = U * S * Vt
    dgesvd_((char*)"S", (char*)"S", (__CLPK_integer*)&m, (__CLPK_integer*)&n,
        A_copy.data(), (__CLPK_integer*)&n,
        S.data(),
        U.data(), (__CLPK_integer*)&k,
        Vt.data(), (__CLPK_integer*)&n,
        work.data(), (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
#else
    err_checkf(LAPACKE_dgesvd(
        LAPACK_ROW_MAJOR,
        'S', 'S',
        m, n,
        A_copy.data(), n,
        S.data(),
        U.data(), k,
        Vt.data(), n,
        superb.data()) == 0,
        "LAPACKE_dgesvd failed with info unequal 0!", std::cout);
#endif

    // 3. Invert Singular Values (Sigma^+)
    // Filter out small singular values
    for (int i = 0; i < k; ++i)
        S[i] = S[i] < cutoff ? 0.0 : 1.0 / S[i];

    // 4. Compute A^+ = V * S^+ * U^T
    A_copy = vec(k * m, 0.0); // Reuse A_copy as W to save memory
    for (int i = 0; i < k; ++i) {
        const int im = i * m;
        if (S[i] == 0.0)
            continue;
#pragma omp simd
        for (int j = 0; j < m; ++j) {
            A_copy[im + j] = S[i] * U[j * k + i]; // U is row-major (m x k)
        }
    }

    // Perform Matrix Multiplication: C = alpha * A * B + beta * C
    return dot<dMatrix2>(reshape<dMatrix2>(Vt, Shape2D(k, n)), reshape<dMatrix2>(A_copy, Shape2D(k, m)), true, false);
}

void make_Eigenvalues(vec& A, vec& W) {
    const int n = static_cast<int>(W.size());
#ifdef __APPLE__
    // Use Accelerate framework on macOS for eigenvalue computation
    __CLPK_integer info = 0;
    __CLPK_integer lwork = -1;
    double work_query = 0.0;
    // Query optimal work size
    dsyev_((char*)"V", (char*)"U", (__CLPK_integer*)&n,
        A.data(), (__CLPK_integer*)&n,
        W.data(),
        &work_query, (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
    lwork = static_cast<__CLPK_integer>(work_query);
    vec work(lwork);
    // Perform eigenvalue decomposition
    dsyev_((char*)"V", (char*)"U", (__CLPK_integer*)&n,
        A.data(), (__CLPK_integer*)&n,
        W.data(),
        work.data(), (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
    err_checkf(info == 0, "The algorithm failed to compute eigenvalues.", std::cout);
#else
    err_checkf(LAPACKE_dsyev(
        LAPACK_ROW_MAJOR,
        'V', 'U',
        n,
        A.data(),
        n,
        W.data()) == 0, "The algorithm failed to compute eigenvalues.", std::cout);
#endif
}

vec mat_sqrt(vec& A, vec& W, const double cutoff) {
    const int n = static_cast<int>(W.size());
    vec Temp(n * n, 0.0);

    make_Eigenvalues(A, W);

    for (int i = 0; i < n; ++i)
        W[i] = abs(W[i]) < cutoff ? 0.0 : std::sqrt(abs(W[i]));

    double* T;
    int in, jn;
    for (int i = 0; i < n; i++) {
        in = i * n;
        for (int j = 0; j < n; j++) {
            jn = j * n;
            T = &Temp[in + j];
#pragma omp simd
            for (int k = 0; k < n; k++)
                *T += A[in + k] * W[k] * A[jn + k];
        }
    }

    return Temp;
}

//Swippedy swappdy column/rows in a symmteric matrix
template <typename T>
void swap_rows_cols_symm(T& mat, const int i, const int j) {
    err_checkf(mat.extent(0) == mat.extent(1), "Matrix must be square to swap rows/cols symmetrically!", std::cout);
    const int n = static_cast<int>(mat.extent(0));
    using Datatype = typename T::value_type;
    if constexpr (std::is_same_v<Datatype, double>)
    {
        // swap rows i and j
        cblas_dswap(n,               // number of elements in a row
            mat.data() + i * n, 1,       // row i, contiguous
            mat.data() + j * n, 1);      // row j, contiguous

        // swap columns i and j
        cblas_dswap(n,               // number of elements in a column
            mat.data() + i, n,      // column i, stride = lda
            mat.data() + j, n);     // column j, stride = lda
    }
    else if constexpr (std::is_same_v<Datatype, cdouble>)
    {
        // swap rows i and j
        cblas_zswap(n,                                       // number of elements in a row
            reinterpret_cast<cdouble*>(mat.data() + i * n), 1,       // row i, contiguous
            reinterpret_cast<cdouble*>(mat.data() + j * n), 1);      // row j, contiguous
        // swap columns i and j
        cblas_zswap(n,                                       // number of elements in a column
            reinterpret_cast<cdouble*>(mat.data() + i), n,      // column i, stride = lda
            reinterpret_cast<cdouble*>(mat.data() + j), n);     // column j, stride = lda
    }
    else
    {
        err_not_impl_f("Unsupported data type for swapping rows/cols symmetrically", std::cout);
    }
}
template void swap_rows_cols_symm(dMatrix2& mat, const int row1, const int row2);
template void swap_rows_cols_symm(cMatrix2& mat, const int row1, const int row2);
