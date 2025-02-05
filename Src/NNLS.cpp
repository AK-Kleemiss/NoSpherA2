
#include "NNLS.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

// You may need to adjust these to match your LAPACKE/CBLAS headers: 
#if has_RAS == 1
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include "lapacke.h" // for LAPACKE_xxx
#include "cblas.h" // for cblas_dnrm2, etc.
#endif



NNLSResult nnls(
    std::vector<double>& A, int m, int n,
    std::vector<double>& B,
    int maxiter,
    double tol)
{
#if has_RAS != 1
	std::cerr << "Error: NNLS requires LAPACKE and CBLAS.\n";
	exit(1);
#endif
    // Define output Variables
    vec X(n,0);
	double RNORM = 0.0;
	int MODE = 0;

	// Check input dimensions
	if (A.size() != m * n) {
		std::cerr << "Error: Matrix A has incorrect dimensions in NNLS.\n";
		return NNLSResult{ X, RNORM, 1 };
	}


    // Define workspace variables
    std::vector<double> AtA(n * n, 0.0); // A^T * A
    std::vector<double> Atb(n, 0.0);     // A^T * b
    std::vector<double> W(n, 0.0);       // Dual vector
    std::vector<double> S(n, 0.0);       // Trial solution
    std::vector<bool> P(n, false);       // Active set (boolean)

#if has_RAS == 1
    // Compute A^T * A (normal equations matrix)
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
        n, n, m, 1.0, A.data(), m, A.data(), m,
        0.0, AtA.data(), n);

    // Compute A^T * B (normal equations RHS)
    cblas_dgemv(CblasColMajor, CblasTrans, m, n,
        1.0, A.data(), m, B.data(), 1,
        0.0, Atb.data(), 1);
#endif
    // Set max iterations
    if (maxiter == -1) maxiter = 3 * n;
    if (tol == -1) tol = 10 * std::max(m, n) * std::numeric_limits<double>::epsilon();

    // Initialize W, S
    W = Atb; // Projected residual W = A^T * B

    int iter = 0;

    while (iter < maxiter) {
        // Step B: Find most active coefficient
        int k = -1;
        double maxW = -1e12;
        for (int i = 0; i < n; i++) {
            if (!P[i] && W[i] > tol && W[i] > maxW) {
                k = i;
                maxW = W[i];
            }
        }

        if (k == -1) break; // No positive residuals, terminate.

        // Step B.3: Move k to active set
        P[k] = true;

        // Solve least squares for active set (B.4)
        std::vector<int> activeIndices;
        for (int i = 0; i < n; i++) {
            if (P[i]) activeIndices.push_back(i);
        }

        int activeCount = activeIndices.size();
        std::vector<double> AtA_active(activeCount * activeCount, 0.0);
        std::vector<double> Atb_active(activeCount, 0.0);
        std::vector<int> ipiv(activeCount);

        // Extract submatrix AtA[P, P] and Atb[P]
        for (int i = 0; i < activeCount; i++) {
            int row = activeIndices[i];
            for (int j = 0; j < activeCount; j++) {
                int col = activeIndices[j];
                AtA_active[i * activeCount + j] = AtA[row * n + col];
            }
            Atb_active[i] = Atb[row];
        }

        // Solve AtA_active * S[P] = Atb_active
#if has_RAS == 1
        int info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', activeCount, 1,
            AtA_active.data(), activeCount, Atb_active.data(), activeCount);
#endif

        if (info != 0) {
            std::cerr << "Warning: Ill-conditioned matrix detected in NNLS.\n";
            MODE = 1;
            break;
        }

        // Assign solution to S
        for (int i = 0; i < activeCount; i++) {
            S[activeIndices[i]] = Atb_active[i];
        }

        // Step C: Check feasibility
        while (iter < maxiter) {
            iter++;
            double minS = 1e12;
            int minIdx = -1;

            for (int i = 0; i < activeCount; i++) {
                if (S[activeIndices[i]] < 0 && S[activeIndices[i]] < minS) {
                    minS = S[activeIndices[i]];
                    minIdx = activeIndices[i];
                }
            }

            if (minS >= 0) break; // All positive, proceed.

            // Compute alpha to move back in feasible space
            double alpha = 1.0;
            for (int i = 0; i < activeCount; i++) {
                int idx = activeIndices[i];
                if (S[idx] < 0) {
                    alpha = std::min(alpha, X[idx] / (X[idx] - S[idx]));
                }
            }

            // Adjust X and remove minIdx from active set
            for (int i = 0; i < n; i++) {
                X[i] = X[i] + alpha * (S[i] - X[i]);
            }
            P[minIdx] = false; // Remove from active set
        }

        // Assign final solution
        for (int i = 0; i < n; i++) {
            X[i] = S[i];
        }
#if has_RAS == 1
        // Compute residual W = Atb - AtA @ X
        cblas_dcopy(n, Atb.data(), 1, W.data(), 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, n, n,
            -1.0, AtA.data(), n, X.data(), 1,
            1.0, W.data(), 1);
#endif
    }

    // Compute residual norm ||A * X - B||
    std::vector<double> Ax(m, 0.0);
#if has_RAS == 1
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, n,
        1.0, A.data(), m, X.data(), 1,
        0.0, Ax.data(), 1);
#endif
    double sum_sq = 0.0;
    for (int i = 0; i < m; i++) {
        sum_sq += (Ax[i] - B[i]) * (Ax[i] - B[i]);
    }

	return NNLSResult(X, std::sqrt(sum_sq), MODE);
}