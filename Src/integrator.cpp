#include "integrator.h"
#include "integration_params.h"
#include "JKFit.h"
#include "int2e.h"
#include "math.h"

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <memory>
#include <cstddef>
#if has_RAS == 1
#include "lapacke.h"
#include "cblas.h"
#endif



vec einsum_ijk_ij_p(const vec3& v1, const vec2& v2)
{
    const int I = v1.size();
    const int J = v1[0].size();
    const int P = v1[0][0].size();
    // Initialize the result vector
    vec rho(P, 0.0);

    // Perform the summation
    for (int p = 0; p < P; ++p)
    {
        for (int i = 0; i < I; ++i)
        {
            for (int j = 0; j < J; ++j)
            {
                rho[p] += v1[i][j][p] * v2[i][j];
            }
        }
    }
    return rho;
}

vec einsum_ijp_ij_p(const vec& v1, const vec& v2, const int I, const int J, const int P) {
    // Initialize the result vector
    vec rho(P, 0.0);

    //v1 of size I * J * P
    //v2 of size I * J

    // Perform the summation
    for (int p = 0; p < P; ++p)
    {
        for (int i = 0; i < I; ++i)
        {
            for (int j = 0; j < J; ++j)
            {
                int inda = p * I * J + i * J + j;
                rho[p] += v1[p * I * J + i * J + j] * v2[i * J + j];
            }
        }
    }

    return rho;
}



void solve_linear_system(const vec2& A, vec& b)
{
    err_checkf(A.size() == b.size(), "Inconsitent size of arrays in linear_solve", std::cout);
    // LAPACK variables
    const int n = A.size(); // The order of the matrix eri2c
    const int nrhs = 1;     // Number of right-hand sides (columns of rho and )
    const int lda = n;      // Leading dimension of eri2c
    const int ldb = n;      // Leading dimension of rho
    ivec ipiv(n, 0);        // Pivot indices
    vec temp = flatten(transpose(A));

#if has_RAS == 1
    // Call LAPACK function to solve the system
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, temp.data(), lda, ipiv.data(), b.data(), ldb);

    if (info != 0)
    {
        std::cout << "Error: LAPACKE_dgesv returned " << info << std::endl;
    }
#endif
}

void solve_linear_system(const vec& A, const int& size_A, vec& b)
{
    err_checkf(size_A == b.size(), "Inconsitent size of arrays in linear_solve", std::cout);
    // LAPACK variables
    const int n = size_A; // The order of the matrix eri2c
    const int nrhs = 1;     // Number of right-hand sides (columns of rho and )
    const int lda = n;      // Leading dimension of eri2c
    const int ldb = n;      // Leading dimension of rho
    ivec ipiv(n, 0);        // Pivot indices

    //Manually transpose the flattened matrix in A of size(size_A, size_A)
    vec temp = transpose(A, n, n);

#if has_RAS == 1
    // Call LAPACK function to solve the system
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, temp.data(), lda, ipiv.data(), b.data(), ldb);

    if (info != 0)
    {
        std::cout << "Error: LAPACKE_dgesv returned " << info << std::endl;
    }
#endif
}

int density_fit(WFN& wavy, const std::string auxname)
{
    WFN wavy_aux;
    load_basis_into_WFN(wavy_aux, BasisSetLibrary().get_basis_set(auxname));
    std::vector<double> eri2c;
    std::vector<double> eri3c;

    // Initialize basis functions (qmBasis and auxBasis)

    wavy.get_DensityMatrix();

    Int_Params normal_basis(wavy);
    Int_Params aux_basis(wavy);

    // Compute integrals
    computeEri2c(aux_basis, eri2c);
    computeEri3c(normal_basis, aux_basis, eri3c);

    // Convert eri3c to matrix form and perform contractions using BLAS

    std::cout << "Done!" << std::endl;
    return 0;
}

int fixed_density_fit_test()
{
#ifdef _WIN32
    void* blas = math_load_BLAS(4);
    if (blas == NULL)
    {
        ExtractDLL("libopenblas.dll");
        blas = math_load_BLAS(4);
    }
    err_checkf(blas != NULL, "BLAS NOT LOADED CORRECTLY!", std::cout);
#endif // __WIN32

    WFN wavy_gbw("H2.gbw");
    Int_Params normal_basis(wavy_gbw);

    WFN wavy_aux(0);
    wavy_aux.atoms = wavy_gbw.atoms;
    wavy_aux.set_ncen(wavy_gbw.get_ncen());
    wavy_aux.delete_basis_set();
    //Int_Params aux_basis(wavy_aux, "cc-pvqz-jkfit");
    Int_Params aux_basis(wavy_aux, "6-31G**_rifit");

    vec eri2c_test_test;
    vec eri3c_test_test;
    computeEri2c(aux_basis, eri2c_test_test);


    std::ofstream file("eri2c_nospherA2.txt");
    if (file.is_open())
    {
    	for (int i = 0; i < eri2c_test_test.size(); i++)
    	{
    		file << std::fixed << std::showpoint << std::setprecision(12) <<eri2c_test_test[i] << std::endl;
    	}
    	file.close();
    }

    computeEri3c(normal_basis, aux_basis, eri3c_test_test);
    std::ofstream file2("eri3c_nospherA2.txt");
    if (file2.is_open())
    {
        for (int i = 0; i < eri3c_test_test.size(); i++)
        {
            file2 << std::fixed << std::showpoint << std::setprecision(12) << eri3c_test_test[i] << std::endl;
        }
        file2.close();
    }

    //   vec dm_test = flatten(wavy_gbw.get_dm());
       //vec rho_test = einsum_ijp_ij_p(eri3c_test_test, dm_test, normal_basis.get_nao(), normal_basis.get_nao(), aux_basis.get_nao());
    //   vec2 intermediate = reshape(eri2c_test_test, { aux_basis.get_nao(), aux_basis.get_nao() });
    //   solve_linear_system(intermediate, rho_test);

#ifdef _WIN32
    math_unload_BLAS(blas);
#endif
    return 0;
}