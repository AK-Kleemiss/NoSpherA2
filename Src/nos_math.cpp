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

#ifdef __APPLE__
    std::transform(matrix.container().begin(), matrix.container().end(), result_m.data(), [exponent](double val)
                   { return std::pow(val, exponent); });
#else
    std::transform(std::execution::par, matrix.container().begin(), matrix.container().end(), result_m.data(), [exponent](double val)
                   { return std::pow(val, exponent); });
#endif

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
void compare_matrices(const Kokkos::Experimental::mdarray<T, Kokkos::extents<unsigned long long, std::dynamic_extent>> &A, const std::vector<T> &B)
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
    ivec dims = {10, 10};
    // Init Mat A with some values as a 3x3 matrix
    vec2 A(dims[0], vec(dims[1]));
    vec2 B(dims[0], vec(dims[1]));
    // Init A and B with random values between -100 and 100
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
    // Init Mat A and Mat B as 3x3 matrices

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
    for (int i = 0; i < dims[1]; i++)
    {
        E[i] = rand() % 200 - 100;
        for (int j = 0; j < dims[0]; j++)
        {
            F[j][i] = rand() % 200 - 100;
        }
    }
    vec fE = flatten<double>(F);
    dMatrix1 matE(dims[1]);
    dMatrix2 matF(dims[0], dims[1]);
    std::copy(E.data(), E.data() + E.size(), matE.data());
    std::copy(fE.data(), fE.data() + fE.size(), matF.data());

    // For matrix just reuse matA
    std::cout << "Testing 2D x 1D matrix multiplication" << std::endl;
    compare_matrices(dot(matF, matE, false), self_dot(F, E));

    std::cout << "All BLAS tests passed!" << std::endl;
}

//#include <linalg.hpp>
//void _test_lahva() {
//    lahva::cpu::Vector<double> p(5, 2.0);
//
//    using namespace lahva::cpu;
//    // construct a 5 by 5 matrix, using the Shape struct and initializing the values to 1.0
//    Matrix<float> s(lahva::Shape(5, 5), 1.0);
//
//
//    //lahva::gpu::Vector<double> p(5, 2.0);
//
//    //using namespace lahva::gpu;
//    //// similar to the CPU Matrix, we have a quadratic 5 x 5 matrix
//    //// here we explicitly give the template parameters for the Allocators instead of relying on default values. 
//    //Matrix<float, CudaHostAllocator<float>, CudaDeviceAsyncAllocator<float>> s(5, 1.0);
//
//}


NNLSResult nnls(dMatrix2& A,
    dMatrix1& b,
    int maxiter) {

    int m = A.extent(0), n = A.extent(1);

    if (maxiter == -1) maxiter = 3 * n;

    ivec inds(n);
    vec w(n), x(n), work(m), zz(m);

    for (int i = 0; i < n; ++i) inds[i] = i;

    int iteration = 0, iz1 = 0, nrow = 0, nsetp = 0, jj = 0;
    double tau = 0.0, unorm = 0.0, alpha, beta, cc, ss, wmax, T, tmp;
    bool skip = false;

    while (iz1 < n && nsetp < m) {
        // simulating a goto from col independence check
        if (skip) {
            skip = false;
        }
        else {
            std::fill(w.begin() + iz1, w.end(), 0.0);
            for (int i = iz1; i < n; ++i) {
                for (int j = nrow; j < m; ++j) {
                    w[i] += b(j) * A(j, inds[i]);
                }
            }
        }

        //Find the largest w[j] and its index.
        vec::iterator max_it = std::max_element(w.begin() + iz1, w.end());
        wmax = *max_it;
        int izmax = std::distance(w.begin(), max_it);

        int iz = izmax;
        int j = inds[iz];

        // If wmax <= 0.0, terminate since this is a KKT certificate.
        if (wmax <= 0.0) break;

        //# The sign of wmax is OK for j to be moved to set p.Begin the transformation
        for (int i = nrow; i < m; ++i) {
            work[i] = A(i, j);
        }
        int tmpint = m - nrow;

        LAPACKE_dlarfg(tmpint, &work[nrow], &work[nrow + 1], 1, &tau);
        beta = work[nrow];
        work[nrow] = 1.0;
        unorm = 0.0;
        if (nsetp > 0) {
            for (int i = 0; i < nsetp; ++i) {
                unorm += A(i, j) * A(i, j);
            }
            unorm = std::sqrt(unorm);
        }
		if (unorm + std::abs(beta) * 0.01 - unorm > 0.0) {
			// Column j is sufficiently independent.Copy b into zz and solve for
			// ztest which is the new prospective value for x[j].
			std::copy(b.data(), b.data() + m, zz.begin());

			LAPACKE_dlarfx(LAPACK_COL_MAJOR, 'L', tmpint, 1.0, &work[nrow], tau, &zz[nrow], tmpint, &tmp);
			if (zz[nrow] / beta <= 0.0) {
				// reject column j as a candidate to be moved from set z to set p.
				// Set w[j] to 0.0 and move to the next greatest entry in w.
				w[j] = 0.0;
				continue;
			}
		}
        else {
			// Column j is not numerically independent, reject column j
            w[j] = 0.0;
            continue;
        }
        // column j accepted
        A(nrow, j) = beta;
		std::copy(zz.begin(),zz.end(), b.data());
        inds[iz] = inds[iz1];
        inds[iz1] = j;
        iz1 += 1;
        nsetp += 1;

		if (iz1 < n) {
			for (int i = iz1; i < n; ++i) {
				int col = inds[i];
				for (int j = nrow; j < m; ++j) {
					zz[j] = A(j, col);
				}
				LAPACKE_dlarfx(LAPACK_COL_MAJOR, 'L', tmpint, 1.0, &work[nrow], tau, &zz[nrow], tmpint, &tmp);
				for (int j = nrow; j < m; ++j) {
					A(j, col) = zz[j];
				}
			}
		}
		nrow += 1;

        if (nsetp < m - 1) {
			for (int i = nrow; i < m; ++i) {
				A(i, j) = 0.0;
			}
        }
        w[j] = 0.0;

		std::copy(b.container().begin(), b.container().end(), zz.begin());
		for (int k = 0; k < nsetp; ++k) {
			int ip = nsetp - k - 1;
			if (k != 0) {
				for (int ii = 0; ii < ip + 1; ++ii) {
					zz[ii] -= A(ii, jj) * zz[ip + 1];
				}
			}
			jj = inds[ip];
			zz[ip] /= A(ip, jj);
		}

        while (true) {
            iteration++;
			if (iteration > maxiter) {
				std::cerr << "NNLS did not converge after " << maxiter << " iterations.\n";
				return NNLSResult{ x, 0.0, 1 };
			}

            alpha = 2.0;
			for (int ip = 0; ip < nsetp; ++ip) {
                int k = inds[ip];
                if (zz[ip] <= 0.0) {
					T = -x[k] / (zz[ip] - x[k]);
                    if (alpha > T) {
                        alpha = T;
                        jj = ip;
                    }
                }
			}
			if (alpha == 2.0) break;

            for (int i = 0; i < nsetp; ++i) {
                x[inds[i]] = (1 - alpha) * x[inds[i]] + alpha * zz[i];
            }

            // Modify A, B, and the indices to move coefficient
            // i from set p to set z.While loop simulates a goto
            int i = inds[jj];
            while (true)
            {
                x[i] = 0.0;
                if (jj != nsetp) {
                    jj += 1;
                    for (int j = jj; j < nsetp; ++j) {
                        int ii = inds[j];
                        inds[j - 1] = ii;
                        LAPACKE_dlartgp(A(j - 1, ii), A(j, ii), &cc, &ss, &A(j - 1, ii));
                        A(j, ii) = 0.0;
                        for (int col = 0; col < n; ++col) {
                            if (col != ii) {
                                tmp = A(j - 1, col);
                                A(j - 1, col) = cc * tmp + ss * A(j, col);
                                A(j, col) = -ss * tmp + cc * A(j, col);
                            }
                        }
                        tmp = b(j - 1);
                        b(j - 1) = cc * tmp + ss * b(j);
                        b(j) = -ss * tmp + cc * b(j);
                    }
                }
				nrow -= 1;
				nsetp -= 1;
				iz1 -= 1;
				inds[iz1] = i;
				bool loop_broken = false;
                for (int jj = 0; jj < nsetp; ++jj) {
                    i = inds[jj];
                    if (x[i] <= 0.0) {
						loop_broken = true;
                        break;
                    }
                }
                if (!loop_broken) break;
            }
			std::copy(b.container().begin(), b.container().end(), zz.begin());
            for (int k = 0; k < nsetp; ++k) {
				int ip = nsetp - k - 1;
                if (k != 0) {
                    for (int ii = 0; ii < ip + 1; ++ii) {
                        zz[ii] -= A(ii, jj) * zz[ip + 1];
                    }
                }
				jj = inds[ip];
				zz[ip] /= A(ip, jj);
            }
        }
		for (int i = 0; i < nsetp; ++i) {
			x[inds[i]] = zz[i];
		}

    }
    //Calculate the residual np.linalg.norm(b[nrow:])
    double res = cblas_dnrm2(m - nrow, &b(nrow), 1);
	return NNLSResult{ x, res, 0 };
}

//NNLSResult nnls(dMatrix2& A,
//    dMatrix1& b,
//    int maxiter, double tol) {
//    int m = A.extent(0), n = A.extent(1);
//
//    if (maxiter == -1) maxiter = 3 * n;
//
//    std::vector<int> inds(n);
//    std::vector<double> w(n, 0.0), x(n, 0.0), work(m, 0.0), zz(m, 0.0);
//
//    for (int i = 0; i < n; ++i) inds[i] = i;
//
//    int i = 0, ii = 0, ip = 0, iteration = 0, iz = 0, iz1 = 0, izmax = 0, j = 0, jj = 0, k = 0, col = 0, nrow = 0, nsetp = 0, one = 1, tmpint = 0;
//    double tau = 0.0, unorm = 0.0, ztest, tmp, alpha, beta, cc, ss, wmax, T;
//    bool skip = false;
//
//    while (iz1 < n && nsetp < m) {
//        // simulating a goto from col independence check
//        if (skip) {
//            skip = false;
//        }
//        else {
//            std::fill(w.begin() + iz1, w.end(), 0.0);
//            for (int i = iz1; i < n; ++i) {
//                for (int j = nrow; j < m; ++j) {
//                    w[i] += b(j) * A(j, inds[i]);
//                }
//            }
//        }
//
//        //Find the largest w[j] and its index.
//        wmax = 0.0;
//        for (int i = iz1; i < n; ++i) {
//            if (w[i] > wmax) {
//                wmax = w[i];
//                izmax = i;
//            }
//        }
//        iz = izmax;
//        j = inds[iz];
//
//        // If wmax <= 0.0, terminate since this is a KKT certificate.
//        if (wmax <= 0.0) break;
//
//        //# The sign of wmax is OK for j to be moved to set p.Begin the transformation
//        for (int i = nrow; i < m; ++i) {
//            work[i] = A(i, j);
//        }
//        int tmpint = m - nrow;
//
//        //DLARFGP(N, ALPHA, X, INCX, TAU)
//        LAPACKE_dlarfg(tmpint, &work[nrow], &work[nrow + 1], 1, &tau);
//        beta = work[nrow];
//        work[nrow] = 1.0;
//        unorm = 0.0;
//        if (nsetp > 0) {
//            for (int i = 0; i < nsetp; ++i) {
//                unorm += A(i, j) * A(i, j);
//            }
//            unorm = std::sqrt(unorm);
//        }
//        if (unorm + std::abs(beta) * 0.01 - unorm > 0.0) {
//            // Column j is sufficiently independent.Copy b into zz and solve for
//            // ztest which is the new prospective value for x[j].
//            for (int i = 0; i < m; ++i) {
//                zz[i] = b(i);
//            }
//            LAPACKE_dlarfx(LAPACK_COL_MAJOR, 'L', tmpint, one, &work[nrow], tau, &zz[nrow], tmpint, &tmp);
//            ztest = zz[nrow] / beta;
//            if (ztest <= 0.0) {
//                // reject column j as a candidate to be moved from set z to set p.
//                // Set w[j] to 0.0 and move to the next greatest entry in w.
//                w[j] = 0.0;
//                skip = true;
//                continue;
//            }
//        }
//        else {
//            // Column j is not numerically independent, reject column j
//            w[j] = 0.0;
//            skip = true;
//            continue;
//        }
//        // column j accepted
//        A(nrow, j) = beta;
//        std::copy(zz.begin(), zz.end(), b.data());
//        inds[iz] = inds[iz1];
//        inds[iz1] = j;
//        iz1 += 1;
//        nsetp += 1;
//
//        if (iz1 < n) {
//            for (int i = iz1; i < n; ++i) {
//                col = inds[i];
//                for (int j = nrow; j < m; ++j) {
//                    zz[j] = A(j, col);
//                }
//                LAPACKE_dlarfx(LAPACK_COL_MAJOR, 'L', tmpint, one, &work[nrow], tau, &zz[nrow], tmpint, &tmp);
//                for (int j = nrow; j < m; ++j) {
//                    A(j, col) = zz[j];
//                }
//            }
//        }
//        nrow += 1;
//
//        if (nsetp < m - 1) {
//            for (int i = nrow; i < m; ++i) {
//                A(i, j) = 0.0;
//            }
//        }
//        w[j] = 0.0;
//
//        std::copy(b.container().begin(), b.container().end(), zz.begin());
//        for (int k = 0; k < nsetp; ++k) {
//            ip = nsetp - k - 1;
//            if (k != 0) {
//                for (int ii = 0; ii < ip + 1; ++ii) {
//                    zz[ii] -= A(ii, jj) * zz[ip + 1];
//                }
//            }
//            jj = inds[ip];
//            zz[ip] /= A(ip, jj);
//        }
//        while (true) {
//            iteration++;
//            if (iteration > maxiter) {
//                std::cerr << "NNLS did not converge after " << maxiter << " iterations.\n";
//                return NNLSResult{ x, 0.0, 1 };
//            }
//
//            alpha = 2.0;
//            for (int ip = 0; ip < nsetp; ++ip) {
//                k = inds[ip];
//                if (zz[ip] <= 0.0) {
//                    T = -x[k] / (zz[ip] - x[k]);
//                    if (alpha > T) {
//                        alpha = T;
//                        jj = ip;
//                    }
//                }
//            }
//            if (alpha == 2.0) break;
//
//            for (int i = 0; i < nsetp; ++i) {
//                x[inds[i]] = (1 - alpha) * x[inds[i]] + alpha * zz[i];
//            }
//
//            // Modify A, B, and the indices to move coefficient
//            // i from set p to set z.While loop simulates a goto
//            i = inds[jj];
//            while (true)
//            {
//                x[i] = 0.0;
//                if (jj != nsetp) {
//                    jj += 1;
//                    for (int j = jj; j < nsetp; ++j) {
//                        ii = inds[j];
//                        inds[j - 1] = ii;
//                        LAPACKE_dlartgp(A(j - 1, ii), A(j, ii), &cc, &ss, &A(j - 1, ii));
//                        A(j, ii) = 0.0;
//                        for (int col = 0; col < n; ++col) {
//                            if (col != ii) {
//                                tmp = A(j - 1, col);
//                                A(j - 1, col) = cc * tmp + ss * A(j, col);
//                                A(j, col) = -ss * tmp + cc * A(j, col);
//                            }
//                        }
//                        tmp = b(j - 1);
//                        b(j - 1) = cc * tmp + ss * b(j);
//                        b(j) = -ss * tmp + cc * b(j);
//                    }
//                }
//                nrow -= 1;
//                nsetp -= 1;
//                iz1 -= 1;
//                inds[iz1] = i;
//                bool loop_broken = false;
//                for (int jj = 0; jj < nsetp; ++jj) {
//                    i = inds[jj];
//                    if (x[i] <= 0.0) {
//                        loop_broken = true;
//                        break;
//                    }
//                }
//                if (!loop_broken) break;
//            }
//            std::copy(b.container().begin(), b.container().end(), zz.begin());
//            for (int k = 0; k < nsetp; ++k) {
//                ip = nsetp - k - 1;
//                if (k != 0) {
//                    for (int ii = 0; ii < ip + 1; ++ii) {
//                        zz[ii] -= A(ii, jj) * zz[ip + 1];
//                    }
//                }
//                jj = inds[ip];
//                zz[ip] /= A(ip, jj);
//            }
//        }
//        for (int i = 0; i < nsetp; ++i) {
//            x[inds[i]] = zz[i];
//        }
//
//    }
//    //Calculate the residual np.linalg.norm(b[nrow:])
//    double res = 0.0;
//    for (int i = nrow; i < m; ++i) {
//        res += b(i) * b(i);
//    }
//    res = std::sqrt(res);
//    return NNLSResult{ x, res, 0 };
//}


void set_BLAS_threads(int num_threads)
{
#ifdef _WIN32
    _putenv_s("OPENBLAS_NUM_THREADS", std::to_string(num_threads).c_str());
#else
    std::string nums = "OPENBLAS_NUM_THREADS=" + std::to_string(num_threads);
    char* env = strdup(nums.c_str());
    putenv(env);
#endif
}

// Remeaining functions which were separated for readibility
#include "mat_nos_math.cpp"
#include "vec_nos_math.cpp"