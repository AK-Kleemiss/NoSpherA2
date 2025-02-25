#include "pch.h"
#include "nos_math.h"

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include "lapacke.h" // for LAPACKE_xxx
#include "cblas.h"
#include <execution>

template <typename T>
T conj(const T &val)
{
    if constexpr (std::is_same_v<T, cdouble> || std::is_same_v<T, std::complex<int>>)
    {
        return std::conj(val);
    }
    else
    {
        return val;
    }
}

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
dMatrix2 elementWiseExponentiation(dMatrix2 &matrix, double exponent)
{
    vec result(matrix.size(), 0.0);
    dMatrix2 result_m = reshape<dMatrix2>(result, Shape2D({matrix.extent(0), matrix.extent(1)}));

	std::transform(std::execution::par, matrix.container().begin(), matrix.container().end(), result_m.data(), [exponent](double val) { return std::pow(val, exponent); });

    return result_m;
}

template <typename T>
void compare_matrices(const Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> &A, const std::vector<std::vector<T>> &B)
{
    std::cout << "Matrices have size " << A.extent(0) << "x" << A.extent(1) << std::endl;
    for (int i = 0; i < A.extent(0); i++)
    {
        for (int j = 0; j < A.extent(1); j++)
        {
            auto a = A(i, j);
            auto b = B[i][j];
            if (a != b)
            {
                std::cout << "Values not matching in comparison! " << i << "," << j << std::endl;
                std::cout << a << " != " << b << std::endl;
            }
            // err_checkf(a == b, "values not matching in comparison!", std::cout);
        }
    }
}

template <typename T>
void compare_matrices(const Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long,  std::dynamic_extent>>& A, const std::vector<T>& B)
{
    std::cout << "Matrices have size " << A.extent(0) << std::endl;
    for (int i = 0; i < A.extent(0); i++)
    {
        auto a = A(i);
        auto b = B[i];
        if (a != b)
        {
            std::cout << "Values not matching in comparison! " << i << std::endl;
            std::cout << a << " != " << b << std::endl;
        }
        // err_checkf(a == b, "values not matching in comparison!", std::cout);
    }
}

template <typename T>
void compare_matrices(const std::vector<std::vector<T>> &A, const Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> &B)
{
    std::cout << "Matrices have size " << B.extent(0) << "x" << B.extent(1) << std::endl;
    for (int i = 0; i < A.extent(0); i++)
    {
        for (int j = 0; j < A.extent(1); j++)
        {
            auto a = B(i, j);
            auto b = A[i][j];
            if (a != b)
            {
                std::cout << "Values not matching in comparison! " << i << "," << j << std::endl;
                std::cout << a << " != " << b << std::endl;
            }
            // err_checkf(a == b, "values not matching in comparison!", std::cout);
        }
    }
}

template <typename T>
void compare_matrices(const std::vector<std::vector<T>> &A, const std::vector<std::vector<T>> &B)
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
    ivec dims = { 10, 10 };
    // Init Mat A with some values as a 3x3 matrix
    vec2 A(dims[0], vec(dims[1]));
    vec2 B(dims[0], vec(dims[1]));
    //Init A and B with random values between -100 and 100
	for (int i = 0; i < dims[0]; i++)
	{
		for (int j = 0; j < dims[1]; j++)
		{
			A[i][j] = rand() % 200 - 100;
			B[i][j] = rand() % 200 - 100;
		}
	}

	vec fA = flatten<double>(A);
	vec fB = flatten<double>(B);

    dMatrix2 matA(dims[0], dims[1]);
	std::copy(fA.data(), fA.data() + fA.size(), matA.data());
    dMatrix2 matB(dims[0], dims[1]);
	std::copy(fB.data(), fB.data() + fB.size(), matB.data());
     //Init Mat A and Mat B as 3x3 matrices


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
	cvec2 C(dims[0], cvec(dims[1])), D(dims[0], cvec(dims[1]));
	for (int i = 0; i < dims[0]; i++)
	{
		for (int j = 0; j < dims[1]; j++)
		{
			C[i][j] = cdouble(rand() % 200 - 100, rand() % 200 - 100);
			D[i][j] = cdouble(rand() % 200 - 100, rand() % 200 - 100);
		}
	}

    cvec fC = flatten<cdouble>(C);
    cvec fD = flatten<cdouble>(D);
    cMatrix2 matC(dims[0], dims[1]);
    std::copy(fC.data(), fC.data() + fC.size(), matC.data());
    cMatrix2 matD(dims[0], dims[1]);
    std::copy(fD.data(), fD.data() + fD.size(), matD.data());


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

	// Test 2D x 1D matrix multiplication
    dims[0] = 12;
	vec E(dims[1]);
	vec2 F(dims[0], vec(dims[1]));
    for (int i = 0; i < dims[1]; i++) {
		E[i] = rand() % 200 - 100;
		for (int j = 0; j < dims[0]; j++) {
			F[j][i] = rand() % 200 - 100;
		}
    }
	vec fE = flatten<double>(F);
	dMatrix1 matE(dims[1]);
	dMatrix2 matF(dims[0], dims[1]);
	std::copy(E.data(), E.data() + E.size(), matE.data());
	std::copy(fE.data(), fE.data() + fE.size(), matF.data());

    //For matrix just reuse matA
	std::cout << "Testing 2D x 1D matrix multiplication" << std::endl;
	compare_matrices(dot(matF, matE, false), self_dot(F, E));

    std::cout << "All BLAS tests passed!" << std::endl;
}


NNLSResult nnls(
    dMatrix2& A, dMatrix1& B,
    int maxiter,
    double tol) {
	int m = A.extent(0);
	int n = A.extent(1);
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
            return NNLSResult{ X, RNORM, 1 };
        }

        // Define workspace variables
        std::vector<double> AtA(n * n, 0.0); // A^T * A
        std::vector<double> Atb(n, 0.0);     // A^T * b
        std::vector<double> W(n, 0.0);       // Dual vector
        std::vector<double> S(n, 0.0);       // Trial solution
        std::vector<bool> P(n, false);       // Active set (boolean)

        // Compute A^T * A (normal equations matrix)
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
            n, n, m, 1.0, A.data(), n, A.data(), n,
            0.0, AtA.data(), n);

        // Compute A^T * B (normal equations RHS)
        cblas_dgemv(CblasRowMajor, CblasTrans, m, n,
            1.0, A.data(), n, B.data(), 1,
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
                int col = activeIndices[i];
                for (int j = 0; j < activeCount; j++)
                {
                    int row = activeIndices[j];
                    AtA_active[i * activeCount + j] = AtA[col * n + row];
                }
                Atb_active[i] = Atb[col]; 
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
        cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n,
            1.0, A.data(), n, X.data(), 1,
            0.0, Ax.data(), 1);
        double sum_sq = 0.0;
        for (int i = 0; i < m; i++)
        {
            sum_sq += (Ax[i] - B(i)) * (Ax[i] - B(i));
        }
        sum_sq = std::sqrt(sum_sq);
        NNLSResult resy({ X, sum_sq, MODE });
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


//Remeaining functions which were separated for readibility
#include "mat_nos_math.cpp"
#include "vec_nos_math.cpp"