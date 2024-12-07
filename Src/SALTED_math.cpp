#include "SALTED_math.h"
#if has_RAS
#include "cblas.h"
#include "openblas_config.h"
#endif

template <typename T>
std::vector<std::vector<T>> reshape(const std::vector<T> &flatVec, Shape2D sizes)
{
    std::vector<std::vector<T>> reshapedVec(sizes.rows, std::vector<T>(sizes.cols));
    for (int i = 0; i < sizes.rows; ++i)
    {
        for (int j = 0; j < sizes.cols; ++j)
        {
            reshapedVec[i][j] = flatVec[i * sizes.cols + j];
        }
    }
    return reshapedVec;
}
template vec2 reshape(const vec &flatVec, Shape2D sizes);
template cvec2 reshape(const cvec &flatVec, Shape2D sizes);
template ivec2 reshape(const ivec &flatVec, Shape2D sizes);

// To_3D
template <typename T>
std::vector<std::vector<std::vector<T>>> reshape(const std::vector<T> &flatVec, Shape3D sizes)
{
    std::vector<std::vector<std::vector<T>>> reshapedVec(sizes.depth, std::vector<std::vector<T>>(sizes.rows, std::vector<T>(sizes.cols)));
    for (int i = 0; i < sizes.depth; ++i)
    {
        for (int j = 0; j < sizes.rows; ++j)
        {
            for (int k = 0; k < sizes.cols; ++k)
            {
                reshapedVec[i][j][k] = flatVec[i * sizes.rows * sizes.cols + j * sizes.cols + k];
            }
        }
    }
    return reshapedVec;
}
template vec3 reshape(const vec &flatVec, Shape3D sizes);
template cvec3 reshape(const cvec &flatVec, Shape3D sizes);
template std::vector<ivec2> reshape(const ivec &flatVec, Shape3D sizes);

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
std::vector<std::vector<T>> self_dot(const std::vector<std::vector<T>> &mat1, const std::vector<std::vector<T>> &mat2, bool transp1, bool transp2)
{
    // if either of the matrices is empty, return a empty matrix
    if (mat1.empty() || mat2.empty())
    {
        return {};
    }

    std::vector<std::vector<T>> mat1Copy = mat1;
    std::vector<std::vector<T>> mat2Copy = mat2;

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

    std::vector<std::vector<T>> result(rows1, std::vector<T>(cols2, 0.0));
    const long long int totalIterations = static_cast<long long int>(rows1 * cols2 * cols1);
    size_t total_size = rows1 * cols2;

#pragma omp parallel
    {
        std::vector<T> local_result(total_size, 0.0);
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

// Fast dot product using OpenBLAS
template <typename T>
std::vector<std::vector<T>> dot(const std::vector<std::vector<T>> &mat1, const std::vector<std::vector<T>> &mat2, bool transp1, bool transp2, bool BLAS_enabled)
{
    // if either of the matrices is empty, return a empty matrix
    if (mat1.empty() || mat2.empty())
    {
        return {};
    }
    int m = transp1 ? (int)mat1[0].size() : (int)mat1.size();
    int k1 = transp1 ? (int)mat1.size() : (int)mat1[0].size();
    int k2 = transp2 ? (int)mat2[0].size() : (int)mat2.size();
    int n = transp2 ? (int)mat2.size() : (int)mat2[0].size();
    // The resulting matrix will have dimensions m x n

    // Check if matrix multiplication is possible
    err_checkf(k1 == k2, "Inner matrix dimensions must agree.", std::cout);

    if (BLAS_enabled)
    {
        // Flatten input matrices
        std::vector<T> flatMat1 = flatten(mat1);
        std::vector<T> flatMat2 = flatten(mat2);
        return dot_BLAS(flatMat1, flatMat2, m, k1, k2, n, transp1, transp2);
    }
    else
    {
        std::cout << "Something went wrong, using dot fallback." << std::endl;
        return self_dot(mat1, mat2, transp1, transp2);
    }
    // return dot_BLAS(flatMat1, flatMat2, m, k1, k2, n, transp1, transp2);
}
template std::vector<std::vector<float>> dot(const std::vector<std::vector<float>> &mat1, const std::vector<std::vector<float>> &mat2, bool transp1, bool transp2, bool BLAS_enabled);
template vec2 dot(const vec2 &mat1, const vec2 &mat2, bool transp1, bool transp2, bool BLAS_enabled);
template cvec2 dot(const cvec2 &mat1, const cvec2 &mat2, bool transp1, bool transp2, bool BLAS_enabled);

// When the matrices are given as flat vectors
template <typename T>
std::vector<std::vector<T>> dot(const std::vector<T> &flatMat1, const std::vector<T> &flatMat2, const int &mat1_d0, const int &mat1_d1, const int &mat2_d0, const int &mat2_d1, bool transp1, bool transp2, bool BLAS_enabled)
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

    if (BLAS_enabled)
    {
        return dot_BLAS(flatMat1, flatMat2, m, k1, k2, n, transp1, transp2);
    }
    else
    {
        std::cout << "Something went wrong, using dot fallback." << std::endl;
        std::vector<T> result_flat(m * n, 0.0);
        std::vector<std::vector<T>> mat1_2D = reshape(flatMat1, {m, k1});
        std::vector<std::vector<T>> mat2_2D = reshape(flatMat2, {n, k2});
        return self_dot(mat1_2D, mat2_2D, transp1, transp2);
    }
    // return dot_BLAS(flatMat1, flatMat2, m, k1, k2, n, transp1, transp2);
}
template std::vector<std::vector<float>> dot(const std::vector<float> &flatMat1, const std::vector<float> &flatMat2, const int &mat1_d0, const int &mat1_d1, const int &mat2_d0, const int &mat2_d1, bool transp1, bool transp2, bool BLAS_enabled);
template vec2 dot(const vec &flatMat1, const vec &flatMat2, const int &mat1_d0, const int &mat1_d1, const int &mat2_d0, const int &mat2_d1, bool transp1, bool transp2, bool BLAS_enabled);
template cvec2 dot(const cvec &flatMat1, const cvec &flatMat2, const int &mat1_d0, const int &mat1_d1, const int &mat2_d0, const int &mat2_d1, bool transp1, bool transp2, bool BLAS_enabled);

template <typename T>
std::vector<std::vector<T>> dot_BLAS(const std::vector<T> &flatMat1, const std::vector<T> &flatMat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2)
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

    Shape2D sizes = {m, n};
    return reshape(result_flat, sizes);
}
template std::vector<std::vector<float>> dot_BLAS(const std::vector<float> &mat1, const std::vector<float> &mat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2);
template vec2 dot_BLAS(const std::vector<double> &mat1, const std::vector<double> &mat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2);
template cvec2 dot_BLAS(const std::vector<cdouble> &mat1, const std::vector<cdouble> &mat2, const int &m, const int &k1, const int &k2, const int &n, bool transp1, bool transp2);

// 2D x 1D MATRIX MULTIPLICATION
template <typename T>
std::vector<T> self_dot(const std::vector<std::vector<T>> &mat, const std::vector<T> &_vec, bool transp1 = false)
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

// Base implementation of matrix-vector multiplication
template <typename T>
std::vector<T> dot(const std::vector<std::vector<T>> &mat, const std::vector<T> &vec, bool transp, bool BLAS_enabled)
{
    int mat_rows = static_cast<int>(mat.size());
    int mat_cols = static_cast<int>(mat[0].size());
    int vec_size = static_cast<int>(vec.size());

    // Check if matrix multiplication is possible
    err_checkf(mat_cols == vec_size, "Matrix dimensions do not match for multiplication", std::cout);

    if (BLAS_enabled)
    {
        return dot_BLAS(flatten(mat), vec, mat_rows, mat_cols, transp);
    }
    else
    {
        std::cout << "Something went wrong, using dot fallback." << std::endl;
        return self_dot(mat, vec, transp);
    }
}
template std::vector<float> dot(const std::vector<std::vector<float>> &mat, const std::vector<float> &vec, bool transp, bool BLAS_enabled);
template vec dot(const vec2 &mat, const vec &vec, bool transp, bool BLAS_enabled);
template cvec dot(const cvec2 &mat, const cvec &vec, bool transp, bool BLAS_enabled);

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

// 1D x 1D Vector multiplication
// Here only for non-complex values
template <typename T>
T self_dot(const std::vector<T> &vec1, const std::vector<T> &vec2, bool conjugate)
{
    T result{};
    for (size_t i = 0; i < vec1.size(); ++i)
    {
        result += vec1[i] * vec2[i];
    }
    return result;
}
template float self_dot(const std::vector<float> &vec1, const std::vector<float> &vec2, bool conjugate);
template double self_dot(const std::vector<double> &vec1, const std::vector<double> &vec2, bool conjugate);

// Special implementation for complex values
cdouble self_dot(const std::vector<cdouble> &vec1, const std::vector<cdouble> &vec2, bool conjugate)
{
    cdouble result{};
    if (!conjugate)
    {
        for (size_t i = 0; i < vec1.size(); ++i)
        {
            result += vec1[i] * vec2[i];
        }
    }
    else
    {
        for (size_t i = 0; i < vec1.size(); ++i)
        {
            result += std::conj(vec1[i]) * vec2[i];
        }
    }
    return result;
}

template <typename T>
T dot(const std::vector<T> &vec1, const std::vector<T> &vec2, bool conjugate, bool BLAS_enabled)
{
    // if either of the vectors is empty, return 0
    if (vec1.empty() || vec2.empty())
    {
        return 0;
    }

    size_t size = vec1.size();
    // Check if vector dimensions match
    err_checkf(size == vec2.size(), "Vector dimensions do not match for multiplication", std::cout);

    if (BLAS_enabled)
    {
        return dot_BLAS(vec1, vec2, conjugate);
    }
    else
    {
        return self_dot(vec1, vec2, conjugate);
    }
}
template float dot(const std::vector<float> &vec1, const std::vector<float> &vec2, bool conjugate, bool BLAS_enabled);
template double dot(const std::vector<double> &vec1, const std::vector<double> &vec2, bool conjugate, bool BLAS_enabled);
template cdouble dot(const std::vector<cdouble> &vec1, const std::vector<cdouble> &vec2, bool conjugate, bool BLAS_enabled);

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
#ifdef _WIN32
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
// #if has_RAS
// #ifdef _WIN32
////Windows specific implementation
//{
//	T result = 0.0;
//	HMODULE hOpenBlas = LoadLibrary(TEXT("libopenblas.dll"));
//	if (hOpenBlas != NULL)
//	{
//		ExampleFunctionType eF = (ExampleFunctionType)GetProcAddress(hOpenBlas, "cblas_sdot");
//		if (eF != NULL)
//		{
//			if constexpr (std::is_same_v<T, double>)
//			{
//                result = cblas_ddot((int)vec1.size(), vec1.data(), 1, vec2.data(), 1);
//			}
//			else if constexpr (std::is_same_v<T, float>)
//			{
//				// Call cblas_sdot
//				result = cblas_sdot((int)vec1.size(), vec1.data(), 1, vec2.data(), 1);
//			}
//			else if constexpr (std::is_same_v<T, cdouble>)
//			{
//                if (!conjugate) {
//                    result = reinterpret_cast<cdouble&>(cblas_zdotu((int)vec1.size(), reinterpret_cast<const cdouble*>(vec1.data()), 1, reinterpret_cast<const cdouble*>(vec2.data()), 1));
//                }
//				else {
//					result = reinterpret_cast<cdouble&>(cblas_zdotc((int)vec1.size(), reinterpret_cast<const cdouble*>(vec1.data()), 1, reinterpret_cast<const cdouble*>(vec2.data()), 1));
//				}
//			}
//			else
//			{
//				err_not_impl_f("Unsupported data type for vector multiplication", std::cout);
//                return {};
//			}
//			return result;
//		}
//		else
//		{
//			// DLL found but function not found
//			err_not_impl_f("OpenBLAS DLL found but function not found!", std::cout);
//            return {};
//		}
//	}
//	else
//	{
//		std::cout << "OpenBLAS DLL not found, using fallback." << std::endl;
//        if (!conjugate) {
//            for (size_t i = 0; i < vec1.size(); ++i)
//            {
//                result += vec1[i] * vec2[i];
//            }
//		}
//		else {
//            err_not_impl_f("Conjugate dot product not implemented for this data type", std::cout);
//            return {};
//        }
//		return result;
//	}
//}
// #else
////Linux specific implementation
//{
//	T result = 0.0;
//	if constexpr (std::is_same_v<T, double>)
//	{
//		// Call cblas_ddot
//		result = cblas_ddot(vec1.size(), vec1.data(), 1, vec2.data(), 1);
//	}
//	else if constexpr (std::is_same_v<T, float>)
//	{
//		// Call cblas_sdot
//		result = cblas_sdot(vec1.size(), vec1.data(), 1, vec2.data(), 1);
//	}
//	else if constexpr (std::is_same_v<T, cdouble>)
//	{
//		if (!conjugate) {
//			result = cblas_zdotu(vec1.size(), reinterpret_cast<const cdouble*>(vec1.data()), 1, reinterpret_cast<const cdouble*>(vec2.data()), 1);
//		}
//		else {
//			result = cblas_zdotc(vec1.size(), reinterpret_cast<const cdouble*>(vec1.data()), 1, reinterpret_cast<const cdouble*>(vec2.data()), 1);
//		}
//	}
//	else
//	{
//		err_not_impl_f("Unsupported data type for vector multiplication", std::cout);
//	}
//	return result;
//}
// #endif
// #else
//// Fallback implementation
//{
//	T result = 0.0;
//	if (!conjugate) {
//		for (size_t i = 0; i < vec1.size(); ++i)
//		{
//			result += vec1[i] * vec2[i];
//		}
//	}
//	else {
//		for (size_t i = 0; i < vec1.size(); ++i)
//		{
//			result += conj(vec1[i]) * vec2[i];
//		}
//	}
//	return result;
//}
// #endif
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

template <typename T>
void compare_matrices(const std::vector<std::vector<T>> &A, const std::vector<std::vector<T>> &B)
{
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A[0].size(); j++)
        {
            assert(A[i][j] == B[i][j]);
        }
    }
}

void _test_openblas()
{
    // Init Mat A with some values as a 3x3 matrix
    vec2 A = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    // Init Mat B with some values
    vec2 B = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};

  //First test regular dot-product
  compare_matrices(dot(A, B, false, false, true), self_dot(A, B));

  ////Second compare first transpose
  compare_matrices(dot(A, B, true, false, true), self_dot(transpose(A), B));

  ////Third comparte second transpose
  compare_matrices(dot(A, B, false, true, true), self_dot(A, transpose(B)));

  ////Fourth compare both transposed
  compare_matrices(dot(A, B, true, true, true), self_dot(transpose(A), transpose(B)));

    // Init Complex matrices
    cvec2 C = {{{1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}, {{4.0, 4.0}, {5.0, 5.0}, {6.0, 6.0}}, {{7.0, 7.0}, {8.0, 8.0}, {9.0, 9.0}}};
    cvec2 D = {{{1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}, {{4.0, 4.0}, {5.0, 5.0}, {6.0, 6.0}}, {{7.0, 7.0}, {8.0, 8.0}, {9.0, 9.0}}};

  //First test regular dot-product
  compare_matrices(dot(C, D, false, false, true), self_dot(C, D));

  ////Second compare first transpose
  compare_matrices(dot(C, D, true, false, true), self_dot(transpose(C), D));

  ////Third comparte second transpose
  compare_matrices(dot(C, D, false, true, true), self_dot(C, transpose(D)));

  ////Fourth compare both transposed
  compare_matrices(dot(C, D, true, true, true), self_dot(transpose(C), transpose(D)));

  std::cout << "All BLAS tests passed!" << std::endl;
}
#ifdef _WIN32
#include "DLL_Helper.h"
#endif

void* math_load_BLAS(int num_threads)
{
#if has_RAS
#ifdef _WIN32
    _putenv_s("OPENBLAS_NUM_THREADS", std::to_string(num_threads).c_str());
    typedef void (*ExampleFunctionType)(void);
    void* _hOpenBlas = static_cast<void*>(LoadLibrary(TEXT("libopenblas.dll")));
    if (_hOpenBlas != NULL)
    {
        ExampleFunctionType eF = (ExampleFunctionType)GetProcAddress((HMODULE)_hOpenBlas, "cblas_sgemm");
        if (eF == NULL)
            return NULL;
    }
    return _hOpenBlas;
#else
    std::string nums = "OPENBLAS_NUM_THREADS=" + std::to_string(_opt.threads);
    char* env = strdup(nums.c_str());
    putenv(env);
    _blas_enabled = true;
#endif
#endif
}

void math_unload_BLAS(void* _hOpenBlas)
{
#ifdef _WIN32
    if (_hOpenBlas != NULL)
    {
        int ret;
        int max_iterations = 150;
        while (max_iterations > 0)
        {
            ret = FreeLibrary((HMODULE)_hOpenBlas);
            if (ret == 0)
            {
                break;
            }
            max_iterations--;
        }
        if (max_iterations == 0)
        {
            std::cout << "Could not free the OpenBLAS library" << std::endl;
        }
        else
        {
            _hOpenBlas = NULL;
        }
    }
#endif
}