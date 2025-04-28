#include "pch.h"
#include "integrator.h"
#include "integration_params.h"
#include "JKFit.h"
#include "libCintMain.h"
#include "nos_math.h"

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#include "lapacke.h"
#include "cblas.h"

vec einsum_ijk_ij_p(const dMatrix3 &v1, const dMatrix2 &v2)
{
    const int I = (int)v1.extent(0);
    const int J = (int)v1.extent(1);
    const int P = (int)v1.extent(2);
    // Initialize the result vector
    vec rho(P, 0.0);

    // Perform the summation
    for (int p = 0; p < P; ++p)
    {
        for (int i = 0; i < I; ++i)
        {
            for (int j = 0; j < J; ++j)
            {
                rho[p] += v1(i, j, p) * v2(i, j);
            }
        }
    }
    return rho;
}

void solve_linear_system(const vec2 &A, vec &b)
{
    err_checkf(A.size() == b.size(), "Inconsitent size of arrays in linear_solve", std::cout);
    // LAPACK variables
    const int n = (int) A.size(); // The order of the matrix eri2c
    const int nrhs = 1;     // Number of right-hand sides (columns of rho and )
    const int lda = n;      // Leading dimension of eri2c
    const int ldb = n;      // Leading dimension of rho
    ivec ipiv(n, 0);        // Pivot indices
    vec temp = flatten<double>(transpose(A));
    // Call LAPACK function to solve the system
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, temp.data(), lda, ipiv.data(), b.data(), ldb);

    if (info != 0)
    {
        std::cout << "Error: LAPACKE_dgesv returned " << info << std::endl;
    }
}

void solve_linear_system(const vec &A, const int &size_A, vec &b)
{
    err_checkf(size_A == b.size(), "Inconsitent size of arrays in linear_solve", std::cout);
    // LAPACK variables
    const int n = size_A; // The order of the matrix eri2c
    const int nrhs = 1;   // Number of right-hand sides (columns of rho and )
    const int lda = n;    // Leading dimension of eri2c
    const int ldb = n;    // Leading dimension of rho
    ivec ipiv(n, 0);      // Pivot indices

    // Manually transpose the flattened matrix in A of size(size_A, size_A)
    vec temp = transpose(A, n, n);

    // Call LAPACK function to solve the system
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, temp.data(), lda, ipiv.data(), b.data(), ldb);

    if (info != 0)
    {
        std::cout << "Error: LAPACKE_dgesv returned " << info << std::endl;
    }
}

// Reorder p-Orbitals to the SALTED convention of:
//  L = 1 components following the - 1, 0, +1 convention
// Meaning for every p-Orbital we swap the first and last component
vec reorder_p(vec coefs_in, WFN aux_basis)
{
    vec coefs_out = coefs_in;
    int coef_idx = 0;
    for (int atm_idx = 0; atm_idx < aux_basis.get_ncen(); atm_idx++)
    {
        for (int shell = 0; shell < aux_basis.get_atom_shell_count(atm_idx); shell++)
        {
            int type = aux_basis.get_shell_type(atm_idx, shell); // Guessing only NoSpherA2 basis sets are used! l starts at 0!!!!!
            if (type != 1)
            {
                coef_idx += 2 * type + 1;
                continue;
            }
            coefs_out[coef_idx] = coefs_in[coef_idx + 1];
            coefs_out[coef_idx + 1] = coefs_in[coef_idx + 2];
            coefs_out[coef_idx + 2] = coefs_in[coef_idx];
            coef_idx += 3;
        }
    }

    return coefs_out;
}

vec density_fit(const WFN &wavy, const WFN& wavy_aux, const double max_mem, const char metric)
{
    vec eri2c;
    vec eri3c;
    vec rho;

    // Initialize basis functions (qmBasis and auxBasis)

    Int_Params normal_basis(wavy);

    Int_Params aux_basis(wavy_aux);

    dMatrix2 dm = wavy.get_dm();

    // Compute integrals
    if (metric == 'C')
    {
        computeEri2c(aux_basis, eri2c);
        computeRho_Coulomb(normal_basis, aux_basis, dm, rho, max_mem);
    }
    else if (metric == 'O')
    {
        compute2c_Overlap(aux_basis, eri2c);
        computeRho_Overlap(normal_basis, aux_basis, dm, rho, max_mem);
    }

    solve_linear_system(eri2c, aux_basis.get_nao(), rho);

    rho = reorder_p(rho, wavy_aux);

    return rho;
}

int fixed_density_fit_test()
{
    WFN wavy_gbw("TESTMOL.gbw");
    Int_Params normal_basis(wavy_gbw);

    WFN wavy_aux(0);
    wavy_aux.set_atoms(wavy_gbw.get_atoms());
    wavy_aux.set_ncen(wavy_gbw.get_ncen());
    wavy_aux.delete_basis_set();
    load_basis_into_WFN(wavy_aux, BasisSetLibrary().get_basis_set("TESTING"));

    Int_Params aux_basis(wavy_aux);

    //   normal_basis.print_data("normal_basis");
    // aux_basis.print_data("aux_basis");

    vec eri2c;
    vec eri3c;
    computeEri2c(aux_basis, eri2c);
    computeEri3c(normal_basis, aux_basis, eri3c);

    // std::ofstream file("eri2c_nospherA2.txt");
    // if (file.is_open())
    //{
    //	for (int i = 0; i < eri2c.size(); i++)
    //	{
    //		file << std::fixed << std::showpoint << std::setprecision(12) <<eri2c[i] << std::endl;
    //	}
    //	file.close();
    // }
    // std::ofstream file2("eri3c_nospherA2.txt");
    // if (file2.is_open())
    //{
    //     for (int i = 0; i < eri3c.size(); i++)
    //     {
    //         file2 << std::fixed << std::showpoint << std::setprecision(12) << eri3c[i] << std::endl;
    //     }
    //     file2.close();
    // }

    dMatrix2 dm = wavy_gbw.get_dm();
    dMatrix3 eri3c_3d = reshape<dMatrix3>(eri3c, Shape3D({ normal_basis.get_nao(), normal_basis.get_nao(), aux_basis.get_nao() }));
    // Perform contractions using BLAS
    vec rho_test = einsum_ijk_ij_p(eri3c_3d, dm);

    solve_linear_system(eri2c, aux_basis.get_nao(), rho_test);

    vec coefs = reorder_p(rho_test, wavy_aux);
    //
    //    std::ofstream file2("coeffs_nospherA2.txt");
    //    if (file2.is_open())
    //{
    //    for (int i = 0; i < coefs.size(); i++)
    //    {
    //        file2 << std::fixed << std::showpoint << std::setprecision(12) << coefs[i] << std::endl;
    //    }
    //    file2.close();
    //}

    for (int i = 0; i < rho_test.size(); i++)
    {
        std::cout << std::fixed << std::showpoint << std::setprecision(12) << coefs[i] << "\t";
        // New line every 10 vals
        if ((i + 1) % 10 == 0)
            std::cout << std::endl;
    }
    return 0;
}