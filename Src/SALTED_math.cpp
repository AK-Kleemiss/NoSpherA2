#include "SALTED_math.h"
#ifdef _WIN32
#include "DLL_Helper.h"
#endif
#if has_RAS
#include "cblas.h"
#endif

using namespace std;

template <typename T>
vector<vector<T>> reshape(vector<T> flatVec, Shape2D sizes)
{
    vector<vector<T>> reshapedVec(sizes.rows, vector<T>(sizes.cols));
    for (int i = 0; i < sizes.rows; ++i)
    {
        for (int j = 0; j < sizes.cols; ++j)
        {
            reshapedVec[i][j] = flatVec[i * sizes.cols + j];
        }
    }
    return reshapedVec;
}
template vector<vector<double>> reshape(vector<double> flatVec, Shape2D sizes);
template vector<vector<cdouble>> reshape(vector<cdouble> flatVec, Shape2D sizes);
template vector<vector<int>> reshape(vector<int> flatVec, Shape2D sizes);

// To_3D
template <typename T>
vector<vector<vector<T>>> reshape(vector<T> flatVec, Shape3D sizes)
{
    vector<vector<vector<T>>> reshapedVec(sizes.depth, vector<vector<T>>(sizes.rows, vector<T>(sizes.cols)));
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
template vector<vector<vector<double>>> reshape(vector<double> flatVec, Shape3D sizes);
template vector<vector<vector<cdouble>>> reshape(vector<cdouble> flatVec, Shape3D sizes);
template vector<vector<vector<int>>> reshape(vector<int> flatVec, Shape3D sizes);


// Flatten Vectors 2D
template <typename T>
vector<T> flatten(const vector<vector<T>> &vec2D)
{
    vector<T> flatVec;
    size_t totalSize = 0;
    for (const auto &row : vec2D)
    {
        totalSize += row.size();
    }
    flatVec.reserve(totalSize);
    for (const vector<T> &row : vec2D)
    {
        // Use std::copy to copy the entire row at once
        flatVec.insert(flatVec.end(), row.begin(), row.end());
    }
    return flatVec;
}
template vector<double> flatten(const vector<vector<double>> &vec2D);
template vector<cdouble> flatten(const vector<vector<cdouble>> &vec2D);
template vector<int> flatten(const vector<vector<int>> &vec2D);

// Flatten Vectors 3D
template <typename T>
vector<T> flatten(const vector<vector<vector<T>>> &vec3D)
{
    vector<T> flatVec;
    size_t totalSize = 0;
    for (const auto &row : vec3D)
    {
        for (const auto &innerRow : row)
        {
            totalSize += innerRow.size();
        }
    }
    flatVec.reserve(totalSize);
    for (const vector<vector<T>> &row : vec3D)
    {
        for (const vector<T> &innerRow : row)
        {
            // Use std::copy to copy the entire innerRow at once
            flatVec.insert(flatVec.end(), innerRow.begin(), innerRow.end());
        }
    }
    return flatVec;
}
template vector<double> flatten(const vector<vector<vector<double>>> &vec3D);
template vector<cdouble> flatten(const vector<vector<vector<cdouble>>> &vec3D);
template vector<int> flatten(const vector<vector<vector<int>>> &vec3D);

// SLICE Operation
vector<double> slice(const vector<double> &vec, size_t start, size_t length)
{
    if (start + length > vec.size())
    {
        throw std::out_of_range("Slice range is out of bounds.");
    }

    std::vector<double> result;
    for (size_t i = start; i < start + length; ++i)
    {
        result.push_back(vec[i]);
    }

    return result;
}

// Matrix multiplication
// 2D x 2D MATRIX MULTIPLICATION
//Slow direct implementation
template <typename T>
vector<vector<T>> self_dot(const vector<vector<T>> &mat1, const vector<vector<T>> &mat2)
{
    // if either of the matrices is empty, return a empty matrix
    if (mat1.empty() || mat2.empty())
    {
        return {};
    }

    size_t rows1 = mat1.size();
    size_t cols1 = mat1[0].size();
    size_t rows2 = mat2.size();
    size_t cols2 = mat2[0].size();

    // Check if matrix multiplication is possible
    if (cols1 != rows2)
    {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication");
    }

    vector<vector<T>> result(rows1, vector<T>(cols2, 0.0));
    const long long int totalIterations = static_cast<long long int>(rows1 * cols2 * cols1);
    size_t total_size = rows1 * cols2;
#pragma omp parallel
    {
        vector<T> local_result(total_size, 0.0);
        int i, j, k, flatIndex;

#pragma omp for schedule(static) private(i, j, k, flatIndex) nowait
        for (long long int n = 0; n < totalIterations; ++n)
        {
            i = static_cast<int>(n / (cols2 * cols1));
            j = static_cast<int>((n / cols1) % cols2);
            k = static_cast<int>(n % cols1);
            flatIndex = i * cols2 + j;
            local_result[flatIndex] += mat1[i][k] * mat2[k][j];
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
template vector<vector<float>> self_dot(const vector<vector<float>> &mat1, const vector<vector<float>> &mat2);
template vector<vector<double>> self_dot(const vector<vector<double>> &mat1, const vector<vector<double>> &mat2);
template vector<vector<cdouble>> self_dot(const vector<vector<cdouble>> &mat1, const vector<vector<cdouble>> &mat2);

typedef void (*ExampleFunctionType)(void);

//Fast dot product using OpenBLAS
template <typename T>
vector<vector<T>> dot(const vector<vector<T>>& mat1, const vector<vector<T>>& mat2, bool transp1, bool transp2) {
    // if either of the matrices is empty, return a empty matrix
    if (mat1.empty() || mat2.empty())
    {
        return {};
    }
    size_t m = transp1 ? mat1[0].size() : mat1.size();
    size_t k1 = transp1 ? mat1.size() : mat1[0].size();
    size_t k2 = transp2 ? mat2[0].size() : mat2.size();
    size_t n = transp2 ? mat2.size() : mat2[0].size();
    // The resulting matrix will have dimensions m x n

    // Check if matrix multiplication is possible
    err_checkf(k1 == k2, "Inner matrix dimensions must agree.", std::cout);

    // Flatten input matrices
    vector<T> flatMat1 = flatten(mat1);
    vector<T> flatMat2 = flatten(mat2);
    
    return dot_BLAS(flatMat1, flatMat2, m, k1, k2, n, transp1, transp2);
}
template vector<vector<float>> dot(const vector<vector<float>>& mat1, const vector<vector<float>>& mat2, bool transp1, bool transp2);
template vector<vector<double>> dot(const vector<vector<double>> &mat1, const vector<vector<double>> &mat2, bool transp1, bool transp2);
template vector<vector<cdouble>> dot(const vector<vector<cdouble>> &mat1, const vector<vector<cdouble>> &mat2, bool transp1, bool transp2);

//When the matrices are given as flat vectors
template <typename T>
vector<vector<T>> dot(const vector<T>& flatMat1, const vector<T>& flatMat2, size_t mat1_d0, size_t mat1_d1, size_t mat2_d0, size_t mat2_d1, bool transp1, bool transp2) {
    //Check if flatMat1 and flatMat2 have the correct size
    err_checkf(flatMat1.size() == mat1_d0 * mat1_d1, "flat Matrix 1 has incorrect size", std::cout);
    err_checkf(flatMat2.size() == mat2_d0 * mat2_d1, "flat Matrix 2 has incorrect size", std::cout);

    // if either of the matrices is empty, return a empty matrix
	if (flatMat1.empty() || flatMat2.empty())
	{
		return {};
	}
	size_t m = transp1 ? mat1_d1 : mat1_d0;
    size_t k1 = transp1 ? mat1_d0 : mat1_d1;
    size_t k2 = transp2 ? mat2_d1 : mat2_d0;
    size_t n = transp2 ? mat2_d0 : mat2_d1;
	// The resulting matrix will have dimensions m x n

	// Check if matrix multiplication is possible
	err_checkf(k1 == k2, "Inner matrix dimensions must agree.", std::cout);

	return dot_BLAS(flatMat1, flatMat2, m, k1, k2, n, transp1, transp2);
}
template vector<vector<float>> dot(const vector<float> & flatMat1, const vector<float> & flatMat2, size_t mat1_d0, size_t mat1_d1, size_t mat2_d0, size_t mat2_d1, bool transp1, bool transp2);
template vector<vector<double>> dot(const vector<double> & flatMat1, const vector<double> & flatMat2, size_t mat1_d0, size_t mat1_d1, size_t mat2_d0, size_t mat2_d1, bool transp1, bool transp2);
template vector<vector<cdouble>> dot(const vector<cdouble> & flatMat1, const vector<cdouble> & flatMat2, size_t mat1_d0, size_t mat1_d1, size_t mat2_d0, size_t mat2_d1, bool transp1, bool transp2);

template <typename T>
std::vector<std::vector<T>> dot_BLAS(const std::vector<T>& flatMat1, const std::vector<T>& flatMat2, const size_t m, const size_t k1, const size_t k2, const size_t n, bool transp1, bool transp2)
#if has_RAS
#ifdef _WIN32
//Windows specific implementation
{
    vector<T> result_flat(m * n, 0.0);
    HMODULE hOpenBlas = LoadLibrary(TEXT("libopenblas.dll"));
    if (hOpenBlas != NULL)
    {
        ExampleFunctionType eF = (ExampleFunctionType)GetProcAddress(hOpenBlas, "cblas_sgemm");
        if (eF != NULL)
        {
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
                    reinterpret_cast<const cdouble*>(flatMat1.data()), transp1 ? m : k1,
                    reinterpret_cast<const cdouble*>(flatMat2.data()), transp2 ? k2 : n,
                    &(zero),
                    reinterpret_cast<cdouble*>(result_flat.data()), n);
            }
            else
            {
                err_not_impl_f("Unsupported data type for matrix multiplication", std::cout);
            }

            Shape2D sizes = { m, n };
            return reshape(result_flat, sizes);
        }
        else
        {
            // DLL found but function not found
            err_not_impl_f("OpenBLAS DLL found but function not found!", std::cout);
        }
    }
    else
    {
        // DLL not found, fallback
        std::cout << "OpenBLAS DLL not found, using fallback." << std::endl;
        std::vector<std::vector<T>> mat1_2D = reshape(flatMat1, { m, k1 });
        std::vector<std::vector<T>> mat2_2D = reshape(flatMat2, { k2, n });
        if (transp1 && !transp2)
            return self_dot(transpose(mat1_2D), mat2_2D);
        else if (transp1 && transp2)
            return self_dot(transpose(mat1_2D), transpose(mat2_2D));
        else if (!transp1 && transp2)
            return self_dot(mat1_2D, transpose(mat2_2D));
        else
            return self_dot(mat1_2D, mat2_2D);
    }
}
#else
//Linux specific implementation
{
    vector<T> result_flat(m * n, 0.0);
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
            reinterpret_cast<const cdouble*>(flatMat1.data()), transp1 ? m : k1,
            reinterpret_cast<const cdouble*>(flatMat2.data()), transp2 ? k2 : n,
            &(zero),
            reinterpret_cast<cdouble*>(result_flat.data()), n);
    }
    else
    {
        err_not_impl_f("Unsupported data type for matrix multiplication", std::cout);
    }

    Shape2D sizes = { m, n };
    return reshape(result_flat, sizes);
}
#endif
#else
// Fallback implementation
{
    std::cout << "Something went wrong, using dot fallback." << std::endl;
    vector<T> result_flat(m * n, 0.0);
    std::vector<std::vector<T>> mat1_2D = reshape(flatMat1, { m, k1 });
    std::vector<std::vector<T>> mat2_2D = reshape(flatMat2, { k2, n });
    if (transp1 && !transp2)
        return self_dot(transpose(mat1_2D), mat2_2D);
    else if (transp1 && transp2)
        return self_dot(transpose(mat1_2D), transpose(mat2_2D));
    else if (!transp1 && transp2)
        return self_dot(mat1_2D, transpose(mat2_2D));
    else
        return self_dot(mat1_2D, mat2_2D);
}
#endif
template vector<vector<float>> dot_BLAS(const std::vector<float>& mat1, const std::vector<float>& mat2, const size_t m, const size_t k1, const size_t k2, const size_t n, bool transp1, bool transp2);
template vector<vector<double>> dot_BLAS(const std::vector<double>& mat1, const std::vector<double>& mat2, const size_t m, const size_t k1, const size_t k2, const size_t n, bool transp1, bool transp2);
template vector<vector<cdouble>> dot_BLAS(const std::vector<cdouble>& mat1, const std::vector<cdouble>& mat2, const size_t m, const size_t k1, const size_t k2, const size_t n, bool transp1, bool transp2);

// 2D x 1D MATRIX MULTIPLICATION
template <typename T>
vector<T> self_dot(const vector<vector<T>>& mat, const vector<T>& vec)
{
    int mat_rows = static_cast<int>(mat.size());
    int mat_cols = static_cast<int>(mat[0].size());
    int vec_size = static_cast<int>(vec.size());

    // Check if matrix multiplication is possible
    if (mat_cols != vec_size)
    {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication");
    }

    vector<T> result(mat_rows, 0.0);
    int totalIterations = mat_rows * mat_cols;
#pragma omp parallel
    {
        vector<T> local_result(mat_rows, 0.0);
#pragma omp for schedule(static)
        for (int n = 0; n < totalIterations; ++n)
        {
            int i = n / mat_cols;
            int j = n % mat_cols;
            local_result[i] += mat[i][j] * vec[j];
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
template vector<float> self_dot(const vector<vector<float>>& mat, const vector<float>& vec);
template vector<double> self_dot(const vector<vector<double>>& mat, const vector<double>& vec);
template vector<cdouble> self_dot(const vector<vector<cdouble>>& mat, const vector<cdouble>& vec);

// Base implementation of matrix-vector multiplication
template <typename T>
vector<T> dot(const vector<vector<T>>& mat, const vector<T>& vec, bool transp)
{
	int mat_rows = static_cast<int>(mat.size());
	int mat_cols = static_cast<int>(mat[0].size());
	int vec_size = static_cast<int>(vec.size());

	// Check if matrix multiplication is possible
	if (mat_cols != vec_size)
	{
		throw std::invalid_argument("Matrix dimensions do not match for multiplication");
	}
    
    return dot_BLAS(flatten(mat), vec, mat_rows, mat_cols, transp);
}
template vector<float> dot(const vector<vector<float>>& mat, const vector<float>& vec, bool transp);
template vector<double> dot(const vector<vector<double>>& mat, const vector<double>& vec, bool transp);
template vector<cdouble> dot(const vector<vector<cdouble>>& mat, const vector<cdouble>& vec, bool transp);

template <typename T>
std::vector<T> dot_BLAS(const std::vector<T>& flatMat, const std::vector<T>& vec, const size_t m, const size_t n,  bool transp)
#if has_RAS
#ifdef _WIN32
//Windows specific implementation
{
    vector<T> result(transp ? n : m, 0.0);
    HMODULE hOpenBlas = LoadLibrary(TEXT("libopenblas.dll"));
    if (hOpenBlas != NULL)
    {
        ExampleFunctionType eF = (ExampleFunctionType)GetProcAddress(hOpenBlas, "cblas_sgemm");
        if (eF != NULL)
        {
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
					reinterpret_cast<const cdouble*>(flatMat.data()), transp ? m : n,
					reinterpret_cast<const cdouble*>(vec.data()), 1,
					&(zero),
					reinterpret_cast<cdouble*>(result.data()), 1);
			}
			else
			{
				err_not_impl_f("Unsupported data type for matrix multiplication", std::cout);
			}
            return result;
 
        }
        else
        {
            // DLL found but function not found
            err_not_impl_f("OpenBLAS DLL found but function not found!", std::cout);
        }
    }
    else
    {
        std::cout << "OpenBLAS DLL not found, using fallback." << std::endl;
        std::vector<std::vector<T>> mat1_2D = reshape(flatMat, { m, n });
        if (transp) {
            return self_dot(transpose(mat1_2D), vec);
        }
        else {
            return self_dot(mat1_2D, vec);
        }
    }
}
#else
//Linux specific implementation
{
    vector<T> result_flat(m * n, 0.0);
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
            result_flat.data(), 1);
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
			result_flat.data(), 1);
	}
	else if constexpr (std::is_same_v<T, cdouble>)
        {
		cdouble one = cdouble(1.0, 0.0);
		cdouble zero = cdouble(0.0, 0.0);
		cblas_zgemv(CblasRowMajor,
			transp ? CblasTrans : CblasNoTrans,
			m, n,
			&(one),
			reinterpret_cast<const cdouble*>(flatMat.data()), transp ? m : n,
			reinterpret_cast<const cdouble*>(vec.data()), 1,
			&(zero),
			reinterpret_cast<cdouble*>(result_flat.data()), 1);
	}
	else
	{
		err_not_impl_f("Unsupported data type for matrix multiplication", std::cout);
	}
	return result_flat;
}
#endif
#else
// Fallback implementation
{
    std::vector<std::vector<T>> mat1_2D = reshape(flatMat, { m, n });
    if transp
        return self_dot(transpose(mat1_2D), vec);
    else
        return self_dot(mat1_2D, vec);
}
#endif
template vector<float> dot_BLAS(const std::vector<float>& flatMat, const std::vector<float>& vec, const size_t m, const size_t n, bool transp);
template vector<double> dot_BLAS(const std::vector<double>& flatMat, const std::vector<double>& vec, const size_t m, const size_t n, bool transp);
template vector<cdouble> dot_BLAS(const std::vector<cdouble>& flatMat, const std::vector<cdouble>& vec, const size_t m, const size_t n, bool transp);

//1D x 1D Vector multiplication
template <typename T>
T dot(const std::vector<T>& vec1, const std::vector<T>& vec2, bool conjugate)
{
	// if either of the vectors is empty, return 0
	if (vec1.empty() || vec2.empty())
	{
		return 0;
	}

	size_t size = vec1.size();
	// Check if vector dimensions match
	if (size != vec2.size())
	{
		throw std::invalid_argument("Vector dimensions do not match for multiplication");
	}

    return dot_BLAS(vec1, vec2, conjugate);
}
template float dot(const std::vector<float>& vec1, const std::vector<float>& vec2, bool conjugate);
template double dot(const std::vector<double>& vec1, const std::vector<double>& vec2, bool conjugate);
template cdouble dot(const std::vector<cdouble>& vec1, const std::vector<cdouble>& vec2, bool conjugate);


template <typename T>
T dot_BLAS(const std::vector<T>& vec1, const std::vector<T>& vec2, bool conjugate)
#if has_RAS
#ifdef _WIN32
//Windows specific implementation
{
	T result = 0.0;
	HMODULE hOpenBlas = LoadLibrary(TEXT("libopenblas.dll"));
	if (hOpenBlas != NULL)
	{
		ExampleFunctionType eF = (ExampleFunctionType)GetProcAddress(hOpenBlas, "cblas_sdot");
		if (eF != NULL)
		{
			if constexpr (std::is_same_v<T, double>)
			{
				// Call cblas_ddot
                if (!conjugate) {
                    result = cblas_ddot(vec1.size(), vec1.data(), 1, vec2.data(), 1);
                }
                else {
                    result = cblas_ddot(vec1.size(), vec1.data(), 1, vec2.data(), 1);
                }
				
			}
			else if constexpr (std::is_same_v<T, float>)
			{
				// Call cblas_sdot
				result = cblas_sdot(vec1.size(), vec1.data(), 1, vec2.data(), 1);
			}
			else if constexpr (std::is_same_v<T, cdouble>)
			{
                if (!conjugate) {
                    result = reinterpret_cast<cdouble&>(cblas_zdotu(vec1.size(), reinterpret_cast<const cdouble*>(vec1.data()), 1, reinterpret_cast<const cdouble*>(vec2.data()), 1));
                }
				else {
					result = reinterpret_cast<cdouble&>(cblas_zdotc(vec1.size(), reinterpret_cast<const cdouble*>(vec1.data()), 1, reinterpret_cast<const cdouble*>(vec2.data()), 1));
				}
			}
			else
			{
				err_not_impl_f("Unsupported data type for vector multiplication", std::cout);
			}
			return result;
		}
		else
		{
			// DLL found but function not found
			err_not_impl_f("OpenBLAS DLL found but function not found!", std::cout);
		}
	}
	else
	{
		std::cout << "OpenBLAS DLL not found, using fallback." << std::endl;
        if (!conjugate) {
            for (size_t i = 0; i < vec1.size(); ++i)
            {
                result += vec1[i] * vec2[i];
            }
		}
		else {
            err_not_impl_f("Conjugate dot product not implemented for this data type", std::cout);
        }
		return result;
	}
}
#else
//Linux specific implementation
{
	T result = 0.0;
	if constexpr (std::is_same_v<T, double>)
	{
		// Call cblas_ddot
		result = cblas_ddot(vec1.size(), vec1.data(), 1, vec2.data(), 1);
	}
	else if constexpr (std::is_same_v<T, float>)
	{
		// Call cblas_sdot
		result = cblas_sdot(vec1.size(), vec1.data(), 1, vec2.data(), 1);
	}
	else if constexpr (std::is_same_v<T, cdouble>)
	{
		if (!conjugate) {
			result = cblas_zdotu(vec1.size(), reinterpret_cast<const cdouble*>(vec1.data()), 1, reinterpret_cast<const cdouble*>(vec2.data()), 1);
		}
		else {
			result = cblas_zdotc(vec1.size(), reinterpret_cast<const cdouble*>(vec1.data()), 1, reinterpret_cast<const cdouble*>(vec2.data()), 1);
		}
	}
	else
	{
		err_not_impl_f("Unsupported data type for vector multiplication", std::cout);
	}
	return result;
}
#endif
#else
// Fallback implementation
{
	T result = 0.0;
	if (!conjugate) {
		for (size_t i = 0; i < vec1.size(); ++i)
		{
			result += vec1[i] * vec2[i];
		}
	}
	else {
		for (size_t i = 0; i < vec1.size(); ++i)
		{
			result += conj(vec1[i]) * vec2[i];
		}
	}
	return result;
}
#endif
template float dot_BLAS(const std::vector<float>& vec1, const std::vector<float>& vec2, bool conjugate);    
template double dot_BLAS(const std::vector<double>& vec1, const std::vector<double>& vec2, bool conjugate);
template cdouble dot_BLAS(const std::vector<cdouble>& vec1, const std::vector<cdouble>& vec2, bool conjugate);

// TRANSPOSES
// 3D MATRIX
template <typename T>
vector<vector<vector<T>>> transpose(const vector<vector<vector<T>>> &originalVec)
{
    if (originalVec.empty() || originalVec[0].empty() || originalVec[0][0].empty())
    {
        return {}; // Return an empty vector if the original vector is empty or not 3D
    }

    size_t newDim1 = originalVec[0][0].size(); // New first dimension is the old third dimension
    size_t newDim2 = originalVec.size();       // New second dimension is the old first dimension
    size_t newDim3 = originalVec[0].size();    // New third dimension is the old second dimension

    vector<vector<vector<T>>> transposedVec(newDim1, vector<vector<T>>(newDim2, std::vector<T>(newDim3)));

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
template vector<vector<vector<double>>> transpose(const vector<vector<vector<double>>>& originalVec);
template vector<vector<vector<cdouble>>> transpose(const vector<vector<vector<cdouble>>>& originalVec);
template vector<vector<vector<int>>> transpose(const vector<vector<vector<int>>>& originalVec);

// 2D MATRIX
template <class T>
vector<vector<T>> transpose(const vector<vector<T>> &mat)
{
    int rows = static_cast<int>(mat.size());
    int cols = static_cast<int>(mat[0].size());
    vector<vector<T>> result(cols, vector<T>(rows));

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result[j][i] = mat[i][j];
        }
    }
    return result;
}
template vector<vector<double>> transpose(const vector<vector<double>> &mat);
template vector<vector<cdouble>> transpose(const vector<vector<cdouble>> &mat);
template vector<vector<int>> transpose(const vector<vector<int>> &mat);

// Reorder 3D Vectors following a given order
template <typename T>
vector<vector<vector<T>>> reorder3D(const vector<vector<vector<T>>> &original)
{
    if (original.empty() || original[0].empty() || original[0][0].empty())
    {
        return {}; // Return an empty vector if the original vector is empty or not properly sized
    }

    size_t size1 = original.size();       // Original first dimension size
    size_t size2 = original[0].size();    // Original second dimension size
    size_t size3 = original[0][0].size(); // Original third dimension size

    // New vector with dimensions rearranged according to (2, 0, 1)
    vector<vector<vector<T>>> transposed(size3, vector<vector<T>>(size1, vector<T>(size2)));

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
template vector<vector<vector<double>>> reorder3D(const vector<vector<vector<double>>> &original);

// Function to collect rows from a matrix based on a vector of indices
vec2 collectRows(const vec2 &matrix, const vector<int> &indices)
{
    // If no indices are provided, return empty matrix
    if (indices.empty())
    {
        return {};
    }
    vec2 result;
    for (int index : indices)
    {
        if (index < matrix.size())
        {
            result.push_back(matrix[index]);
        }
        else
        {
            // Handle the case where index is out of bounds
            throw std::out_of_range("Index out of range");
        }
    }
    return result;
}

// Function to collect rows from a Cube based on a vector of indices
vec3 collectRows(const vec3 &cube, const vector<int> &indices)
{
    // If no indices are provided, return empty matrix
    if (indices.empty())
    {
        return {};
    }
    vec3 result;
    for (int index : indices)
    {
        if (index < cube.size())
        {
            result.push_back(cube[index]);
        }
        else
        {
            // Handle the case where index is out of bounds
            throw std::out_of_range("Index out of range");
        }
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
void compare_matrices(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B)
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
  //Init Mat A with some values as a 3x3 matrix
  vec2 A = { {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0} };
  //Init Mat B with some values
  vec2 B = { {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0} };

  //First test regular dot-product
  compare_matrices(dot(A, B, false, false), self_dot(A, B));

  ////Second compare first transpose
  compare_matrices(dot(A, B, true, false), self_dot(transpose(A), B));

  ////Third comparte second transpose
  compare_matrices(dot(A, B, false, true), self_dot(A, transpose(B)));

  ////Fourth compare both transposed
  compare_matrices(dot(A, B, true, true), self_dot(transpose(A), transpose(B)));

  //Init Complex matrices
  cvec2 C = { {{1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}, {{4.0, 4.0}, {5.0, 5.0}, {6.0, 6.0}}, {{7.0, 7.0}, {8.0, 8.0}, {9.0, 9.0}} };
  cvec2 D = { {{1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}, {{4.0, 4.0}, {5.0, 5.0}, {6.0, 6.0}}, {{7.0, 7.0}, {8.0, 8.0}, {9.0, 9.0}} };

  //First test regular dot-product
  compare_matrices(dot(C, D, false, false), self_dot(C, D));

  ////Second compare first transpose
  compare_matrices(dot(C, D, true, false), self_dot(transpose(C), D));

  ////Third comparte second transpose
  compare_matrices(dot(C, D, false, true), self_dot(C, transpose(D)));

  ////Fourth compare both transposed
  compare_matrices(dot(C, D, true, true), self_dot(transpose(C), transpose(D)));

  std::cout << "All BLAS tests passed!" << std::endl;

}