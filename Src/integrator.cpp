#include "integrator.h"
#include "integration_params.h"
#include "JKFit.h"
#include "libCintMain.h"
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

vec einsum_ijk_ij_p(const vec& v1, const vec& v2, const int I, const int J, const int P) {
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

//Reorder p-Orbitals to the SALTED convention of:
// L = 1 components following the - 1, 0, +1 convention
//Meaning for every p-Orbital we swap the first and last component
vec reorder_p(vec coefs_in, WFN aux_basis) {
    vec coefs_out = coefs_in;
    int coef_idx = 0;
    for (int atm_idx = 0; atm_idx < aux_basis.get_ncen(); atm_idx++) {
        for (int shell = 0; shell < aux_basis.get_atom_shell_count(atm_idx); shell++) {
            int type = aux_basis.get_shell_type(atm_idx, shell); //Guessing only NoSpherA2 basis sets are used! l starts at 0!!!!!
            if (type != 1) {
                coef_idx += 2 * type + 1;
                continue;
            }
            coefs_out[coef_idx] = coefs_in[coef_idx +1];
            coefs_out[coef_idx + 1] = coefs_in[coef_idx + 2];
            coefs_out[coef_idx + 2] = coefs_in[coef_idx];
            coef_idx += 3;   
        }
    }


	return coefs_out;
}


vec density_fit(const WFN& wavy, const std::string auxname, const double max_mem)
{
    vec eri2c;
    vec eri3c;
    vec rho;

    // Initialize basis functions (qmBasis and auxBasis)

    Int_Params normal_basis(wavy);

    WFN wavy_aux(0);
    wavy_aux.set_atoms(wavy.get_atoms());
    wavy_aux.set_ncen(wavy.get_ncen());
    wavy_aux.delete_basis_set();
    Int_Params aux_basis(wavy_aux, auxname);

    // Compute integrals
    computeEri2c(aux_basis, eri2c);
    //computeEri3c(normal_basis, aux_basis, eri3c);


    vec2 dm = wavy.get_dm();
    //vec3 eri3c_3d = reshape(eri3c, { normal_basis.get_nao(), normal_basis.get_nao(), aux_basis.get_nao() });
	computeRho(normal_basis, aux_basis, dm, rho, max_mem);
    // Perform contractions using BLAS
    //vec rho = einsum_ijk_ij_p(eri3c_3d, dm);
    solve_linear_system(eri2c, aux_basis.get_nao(), rho);

    vec coefs = reorder_p(rho, wavy_aux);
    return coefs;
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

    WFN wavy_gbw("TESTMOL.gbw");
    Int_Params normal_basis(wavy_gbw);
	


    WFN wavy_aux(0);
	wavy_aux.set_atoms(wavy_gbw.get_atoms());
    wavy_aux.set_ncen(wavy_gbw.get_ncen());
    wavy_aux.delete_basis_set();

    Int_Params aux_basis(wavy_aux, "TESTING");


 //   normal_basis.print_data("normal_basis");
	//aux_basis.print_data("aux_basis");


    vec eri2c;
    vec eri3c;
    computeEri2c(aux_basis, eri2c);
    computeEri3c(normal_basis, aux_basis, eri3c);

    //std::ofstream file("eri2c_nospherA2.txt");
    //if (file.is_open())
    //{
    //	for (int i = 0; i < eri2c.size(); i++)
    //	{
    //		file << std::fixed << std::showpoint << std::setprecision(12) <<eri2c[i] << std::endl;
    //	}
    //	file.close();
    //}
    //std::ofstream file2("eri3c_nospherA2.txt");
    //if (file2.is_open())
    //{
    //    for (int i = 0; i < eri3c.size(); i++)
    //    {
    //        file2 << std::fixed << std::showpoint << std::setprecision(12) << eri3c[i] << std::endl;
    //    }
    //    file2.close();
    //}

    vec2 dm = wavy_gbw.get_dm();
    vec3 eri3c_3d = reshape(eri3c, { normal_basis.get_nao(), normal_basis.get_nao(), aux_basis.get_nao() });
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
        //New line every 10 vals
		if ((i + 1) % 10 == 0) std::cout << std::endl;
    }

#ifdef _WIN32
    math_unload_BLAS(blas);
#endif
    return 0;
}