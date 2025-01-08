#include "integrator.h"
#include "convenience.h"
#include "constants.h"
#include "JKFit.h"

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <memory>
#include <cstddef>
#if has_RAS == 1
#include "lapacke.h"
#include "cblas.h"
#endif
#include "math.h"

// This is an implementation of libcint from PySCF in C++ for use during density fitting calculations
// Libcint is published under Apache License 2.0
//  The original code can be found at
//  https://github.com/sunqm/libcint

#define bas(SLOT, I) bas[8 * (I) + (SLOT)] // Basis set data for atom I
#define atm(SLOT, I) atm[6 * (I) + (SLOT)] // Atom data for atom I


Int_Params::Int_Params() {
    wfn_origin = 0;
    ncen = 0;
}

Int_Params::Int_Params(const WFN& wavy) {
    atoms = wavy.atoms;
	wfn_origin = wavy.get_origin();
	ncen = wavy.get_ncen();
	calc_integration_parameters();
}

Int_Params::Int_Params(WFN& wavy, const std::string auxname) {
    atoms = wavy.atoms;
    wfn_origin = wavy.get_origin();
    ncen = wavy.get_ncen();
    calc_integration_parameters();
    load_basis_into_WFN(wavy, BasisSetLibrary().get_basis_set(auxname));
}

vec Int_Params::normalize_gto(vec coef, vec exp, int l) {
    //GTO norm Ref: H. B. Schlegel and M. J. Frisch, Int. J. Quant.  Chem., 54(1995), 83-87.
    for (int i = 0; i < coef.size(); i++) {
        coef[i] *= 1.0 / std::sqrt(gaussian_int(l * 2 + 2, 2 * exp[i]));
    }

    //Normalize contracted GTO
    //#ee = numpy.empty((nprim, nprim))
    //    #for i in range(nprim) :
    //    #    for j in range(i + 1) :
    //    #        ee[i, j] = ee[j, i] = gaussian_int(angl * 2 + 2, es[i] + es[j])
    //    #s1 = 1 / numpy.sqrt(numpy.einsum('pi,pq,qi->i', cs, ee, cs))
    //    return numpy.einsum('pi,i->pi', cs, s1)
	vec2 ee(coef.size(), vec(coef.size(), 0.0));
	for (int i = 0; i < coef.size(); i++) {
		for (int j = 0; j < i + 1; j++) {
            ee[i][j] = ee[j][i] = gaussian_int(l * 2 + 2, exp[i] + exp[j]);
		}
	}
	//Do the einsum by hand
	vec s1(coef.size(), 0.0);
	for (int i = 0; i < coef.size(); i++) {
		for (int j = 0; j < coef.size(); j++) {
			for (int k = 0; k < coef.size(); k++) {
				s1[i] += coef[k] * ee[k][j] * coef[j];
			}
		}
	}
	for (int i = 0; i < s1.size(); i++) {
		s1[i] = 1.0/ std::sqrt(s1[i]);
	}
	for (int i = 0; i < coef.size(); i++) {
		coef[i] *= s1[i];
	}
	return coef;
}

void Int_Params::collect_basis_data() {
    for (int atom_idx = 0; atom_idx < ncen; atom_idx++) {
        nbas += atoms[atom_idx].shellcount.size();
        if (basis_sets.find(atoms[atom_idx].charge) != basis_sets.end()) {
            continue;
        }
        std::vector<basis_set_entry> basis = atoms[atom_idx].basis_set;
        vec coefficients, exponents;
        for (int shell = 0; shell < basis.size(); shell++){
            coefficients.push_back(basis[shell].coefficient);
			exponents.push_back(basis[shell].exponent);
		}
        if (wfn_origin == 9) {
            for (int i = 0; i < coefficients.size(); i++) {
				coefficients[i] *= constants::sqr_pi * 2;  //Conversion factor from GBW to libcint
            }
        }
        if (wfn_origin == 0) {
            int coef_idx = 0;
            for (int shell = 0; shell < atoms[atom_idx].shellcount.size(); shell++) {
				vec shell_coefs(coefficients.begin() + coef_idx, coefficients.begin() + atoms[atom_idx].shellcount[shell] + coef_idx);
				vec shell_exp(exponents.begin() + coef_idx, exponents.begin() + atoms[atom_idx].shellcount[shell] + coef_idx);
				

                shell_coefs = normalize_gto(shell_coefs, shell_exp, basis[shell].type);
                

				//Place the new coefs at the correct place in the coefficients vector
				std::copy(shell_coefs.begin(), shell_coefs.end(), coefficients.begin() + coef_idx);
                coef_idx += atoms[atom_idx].shellcount[shell];
            }
        }


        basis_sets.insert({ atoms[atom_idx].charge, {coefficients, exponents, {0.0}} });
    }
}

//void Int_Params::collect_basis_data_from_gbw() {
//    for (int atom_idx = 0; atom_idx < ncen; atom_idx++) {
//        nbas += atoms[atom_idx].shellcount.size();
//        if (basis_sets.find(atoms[atom_idx].charge) != basis_sets.end()) {
//            continue;
//        }
//        std::vector<basis_set_entry> basis = atoms[atom_idx].basis_set;
//        vec coefficients, exponents;
//        for (int shell = 0; shell < basis.size(); shell++){
//			double coef = basis[shell].coefficient * constants::sqr_pi * 2;  //Conversion factor from GBW to libcint
//            coefficients.push_back(coef);
//    		exponents.push_back(basis[shell].exponent);
//    	}
//
//        basis_sets.insert({ atoms[atom_idx].charge, {coefficients, exponents, {0.0}} });
//    }
//}
//
//void Int_Params::collect_basis_data_internal() {
//    for (int atom_idx = 0; atom_idx < ncen; atom_idx++) {
//        nbas += atoms[atom_idx].shellcount.size();
//        if (basis_sets.find(atoms[atom_idx].charge) != basis_sets.end()) {
//            continue;
//        }
//        std::vector<basis_set_entry> basis = atoms[atom_idx].basis_set;
//        vec coefficients, exponents;
//        for (int shell = 0; shell < basis.size(); shell++) {
//            double coef = basis[shell].coefficient * constants::sqr_pi * 2;  //Conversion factor from GBW to libcint
//            coefficients.push_back(coef);
//            exponents.push_back(basis[shell].exponent);
//        }
//
//        basis_sets.insert({ atoms[atom_idx].charge, {coefficients, exponents, {0.0}} });
//    }
//}

void Int_Params::populate_atm() {
    //atom: (Z, ptr, nuclear_model=1, ptr+3, 0,0)  nuclear_model seems to be 1 for all atoms  USURE!!! Nested arrays per atom
	_atm.resize(ncen * 6, 0);
	int ptr = 20;
	for (int atom_idx = 0; atom_idx < ncen; atom_idx++)
	{
		//Assign_atoms
		_atm[atom_idx * 6 + 0] = atoms[atom_idx].charge; // Z
		_atm[atom_idx * 6 + 1] = ptr; // ptr
		_atm[atom_idx * 6 + 2] = 1; // nuclear_model
		ptr += 3;
		_atm[atom_idx * 6 + 3] = ptr; // ptr+3
		ptr += 1;
	}
}

void Int_Params::populate_env() {
    //env: (leading 20*0, (x,y,z,zeta=0)*ncen, (coefficients_l=0, exponentsl=0, ...)*ncen )  flat array for everything
    _env.resize(20 + ncen * 4, 0);
	int env_offset = 20;  //Do not know if this is needed
    for (int atom_idx = 0; atom_idx < ncen; atom_idx++) {
        _env[env_offset + atom_idx * 4] = atoms[atom_idx].x;
        _env[env_offset + atom_idx * 4 + 1] = atoms[atom_idx].y;
        _env[env_offset + atom_idx * 4 + 2] = atoms[atom_idx].z;
		_env[env_offset + atom_idx * 4 + 3] = 0;  //zeta
    }
	ivec seen_elements;
    int bas_ptr = _env.size();
	for (int atom_idx = 0; atom_idx < ncen; atom_idx++) {
		if (std::find(seen_elements.begin(), seen_elements.end(), atoms[atom_idx].charge) != seen_elements.end()) {continue;}
		seen_elements.push_back(atoms[atom_idx].charge);

        basis_sets[atoms[atom_idx].charge][2][0] = (double)bas_ptr;

		vec2 basis_data = basis_sets[atoms[atom_idx].charge];
		vec coefficients = basis_data[0];
		vec exponents = basis_data[1];

		int func_count = 0;
        for (int shell = 0; shell < atoms[atom_idx].shellcount.size(); shell++) {
            int n_funcs = atoms[atom_idx].shellcount[shell];
            for (int func = 0; func < n_funcs; func++){
				_env.push_back(exponents[func + func_count]);
			}
            for (int func = 0; func < n_funcs; func++) {
				_env.push_back(coefficients[func + func_count]);
            }
            func_count += n_funcs;
        }
        bas_ptr += func_count;

	}

}

void Int_Params::populate_bas() {
    //basis: (atom_id, l, nprim, ncentr, kappa=0, ptr, ptr+nprim, 0)  atom_id has to be consistent with order in other arrays  |  Nested arrays for each contracted basis function and then stacked for all atoms
    _bas.resize(8 * nbas, 0);
    int index = 0;
    for (int atom_idx = 0; atom_idx < ncen; atom_idx++){
		int bas_ptr = basis_sets[atoms[atom_idx].charge][2][0];
        for (int shell = 0; shell < atoms[atom_idx].shellcount.size(); shell += 1) {
			_bas[index * 8 + 0] = atom_idx; // atom_id
            if (wfn_origin == 0) {
                _bas[index * 8 + 1] = atoms[atom_idx].basis_set[shell].type; // l s=0, p=1, d=2 ...
			}
			else {
				_bas[index * 8 + 1] = atoms[atom_idx].basis_set[shell].type - 1; // l s=0, p=1, d=2 ...
			}		
			_bas[index * 8 + 2] = atoms[atom_idx].shellcount[shell];  // nprim
			_bas[index * 8 + 3] = 1; // ncentr    Not sure
			_bas[index * 8 + 5] = bas_ptr;  //Pointer to the end of the _env vector
			bas_ptr += _bas[index * 8 + 2];
			_bas[index * 8 + 6] = bas_ptr;
			//_bas 8 is set to 0 
			bas_ptr += _bas[index * 8 + 2] * _bas[atom_idx * 8 + 3];
            index++;
        }
    }

}

void Int_Params::calc_integration_parameters() {
    collect_basis_data();
    populate_atm();
    populate_env();
    populate_bas();
}

Int_Params Int_Params::operator+(const Int_Params& other) {
    //off = len(env1)
    //natm_off = len(atm1)
    //atm2 = numpy.copy(atm2)
    //bas2 = numpy.copy(bas2)
    //atm2[:, PTR_COORD] += off
    //atm2[:, PTR_ZETA] += off
    //bas2[:, ATOM_OF] += natm_off
    //bas2[:, PTR_EXP] += off
    //bas2[:, PTR_COEFF] += off
	//atm = numpy.vstack((atm1, atm2))
	//bas = numpy.vstack((bas1, bas2))
	//env = numpy.hstack((env1, env2))
    
	Int_Params combined;
	combined._atm.resize(_atm.size() + other._atm.size(), 0);
	combined._bas.resize(_bas.size() + other._bas.size(), 0);
	combined._env.resize(_env.size() + other._env.size(), 0);
    combined.ncen = ncen;
	combined.nbas = nbas + other.nbas;

	std::copy(_atm.begin(), _atm.end(), combined._atm.begin());
	std::copy(_bas.begin(), _bas.end(), combined._bas.begin());
	std::copy(_env.begin(), _env.end(), combined._env.begin());
	std::copy(other._env.begin(), other._env.end(), combined._env.begin() + _env.size());

	int off = _env.size();
	int natm_off = _atm.size() / 6.0;
	ivec atm2 = other._atm;
	ivec bas2 = other._bas;
    for (int a = 0; a < natm_off; a++) {
		atm2[a * 6 + 1] += off;
		atm2[a * 6 + 3] += off;
    }
    for (int b = 0; b < other.nbas; b++) {
		bas2[b * 8 + 0] += natm_off;
		bas2[b * 8 + 5] += off;
		bas2[b * 8 + 6] += off;
    }
	std::copy(atm2.begin(), atm2.end(), combined._atm.begin() + _atm.size());
	std::copy(bas2.begin(), bas2.end(), combined._bas.begin() + _bas.size());
	return combined;
}


int* make_loc(int* bas, int nbas) {
    constexpr int ANG_OF = 1; // Column index for angular momentum
    constexpr int NCTR_OF = 3; // Column index for number of centers

    ivec dims(nbas, 0);
    // Calculate (2*l + 1) * nctr for spherical harmonics
    for (size_t i = 0; i < nbas; i++) {
        dims[i] = (2 * bas(ANG_OF, i) + 1) * bas(NCTR_OF, i);
    }

    // Create the ao_loc array
    int* ao_loc = new int[nbas + 1];
    ao_loc[0] = 0;

    // Compute the cumulative sum
    std::partial_sum(dims.begin(), dims.end(), ao_loc+1);
    
    return ao_loc;
}

// Function to compute two-center two-electron integrals (eri2c)
void computeEri2c(Int_Params &params, std::vector<double> &eri2c)
{
    int* bas = params.get_ptr_bas();
	int* atm = params.get_ptr_atm();
	double* env = params.get_ptr_env();

    int nbas = params.get_nbas();

	int shl_slice[] = { 0, nbas, 0, nbas };
	int nat = 1;  //No idea what this has to be
	int *aoloc = make_loc(bas, nbas);

	int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
	int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];

	eri2c.resize(naoi * naoj, 0.0);
	Opt opty = int2c2e_optimizer(atm, nat, bas, nbas, env);
	// Compute integrals
	GTOint2c(int2c2e_sph, eri2c.data(), 1, 0, shl_slice, aoloc, &opty, atm, nat, bas, nbas, env);
	free(aoloc);
}

// Function to compute three-center two-electron integrals (eri3c)
void computeEri3c(Int_Params &param1,
                  Int_Params &param2,
                  std::vector<double> &flat_eri3c)
{   
    int nQM = param1.get_nbas();
    int nAux = param2.get_nbas();

    Int_Params combined = param1 + param2;

    int* bas = combined.get_ptr_bas();
    int* atm = combined.get_ptr_atm();
    double* env = combined.get_ptr_env();

    int shl_slice[] = {
        0,
        nQM,
        0,
        nQM,
        nQM,
        nQM + nAux,
    };
    int nat = 1;
    int nbas = combined.get_nbas();
    int *aoloc = make_loc(bas, nbas);
    int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
    int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];
    int naok = aoloc[shl_slice[5]] - aoloc[shl_slice[4]];

    flat_eri3c.resize(naoi * naoj * naok, 0.0);
    Opt opty = int3c2e_optimizer(atm, nat, bas, nbas, env);
    // Compute integrals
    GTOnr3c_drv(int3c2e_sph, GTOnr3c_fill_s1, flat_eri3c.data(), 1, shl_slice, aoloc, &opty, atm, nat, bas, nbas, env);
    free(aoloc);
}

vec einsum_ijk_ij_p(const vec3 &v1, const vec2 &v2)
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


void solve_linear_system(const vec2 &A, vec &b)
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

int density_fit( WFN &wavy, const std::string auxname)
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
    void *blas = math_load_BLAS(4);
    if (blas == NULL)
    {
        ExtractDLL("libopenblas.dll");
        blas = math_load_BLAS(4);
    }
    err_checkf(blas != NULL, "BLAS NOT LOADED CORRECTLY!", std::cout);
#endif // __WIN32
    vec2 dm = {
        {
            0.60245569,
            0.60245569,
        },
        {
            0.60245569,
            0.60245569,
        },
    };

    vec2 eri2c_test{{1., 2.}, {3., 4.}};
    vec3 eri3c_test{{{1., 2.}, {3., 4.}}, {{5., 6.}, {7., 8.}}};

    // vec rho_t(2, 0.0);
    vec rho = einsum_ijk_ij_p(eri3c_test, dm);
    solve_linear_system(eri2c_test, rho);
    std::cout << "Test done!" << std::endl;

    WFN wavy_gbw("H2.gbw");
	Int_Params normal_basis(wavy_gbw);

    WFN wavy_aux(0);
	wavy_aux.atoms = wavy_gbw.atoms;
    wavy_aux.set_ncen(wavy_gbw.get_ncen());
    wavy_aux.delete_basis_set();
	Int_Params aux_basis(wavy_aux, "cc-pvqz-jkfit");

    vec eri2c_test_test;
    vec eri3c_test_test;
	computeEri2c(aux_basis, eri2c_test_test);
    computeEri3c(normal_basis, aux_basis, eri3c_test_test);

#ifdef _WIN32
    math_unload_BLAS(blas);
#endif
	return 0;

    WFN wavy("H2.molden");

    vec2 eri2c{
        {5.19870815, 9.45570563, 12.35839308, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 3.01064326, 7.22843174, 10.70109626, 0., 0., -1.81824542, 0., 0., -3.16522589, 0., 0., 0.90587703, 0., 0.},
        {9.45570563, 21.11552664, 30.43158895, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 7.22843174, 17.64603115, 27.15849891, 0., 0., -2.71979468, 0., 0., -5.41825866, 0., 0., 0.7165021, 0., 0.},
        {12.35839308, 30.43158895, 47.51332266, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 10.70109626, 27.15849891, 43.71681781, 0., 0., -2.40501078, 0., 0., -5.45538584, 0., 0., 0.34262858, 0., 0.},
        {0., 0., 0., 2.98672814, 0., 0., 4.01618213, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.40570429, 0., -0., 2.52095233, 0., -0., 0., 0., 0., -1.02543666, 0.},
        {0., 0., 0., 0., 2.98672814, 0., 0., 4.01618213, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.40570429, -0., 0., 2.52095233, -0., 0., -1.02543666, 0., 0., 0.},
        {0., 0., 0., 0., 0., 2.98672814, 0., 0., 4.01618213, 0., 0., 0., 0., 0., 1.81824542, 2.71979468, 2.40501078, 0., 0., -0.51875067, 0., 0., 0.29259347, 0., 0., -0.07438893, 0., 0.},
        {0., 0., 0., 4.01618213, 0., 0., 7.0373349, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.52095233, 0., -0., 5.02200529, 0., -0., 0., 0., 0., -1.15041631, 0.},
        {0., 0., 0., 0., 4.01618213, 0., 0., 7.0373349, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.52095233, -0., 0., 5.02200529, -0., 0., -1.15041631, 0., 0., 0.},
        {0., 0., 0., 0., 0., 4.01618213, 0., 0., 7.0373349, 0., 0., 0., 0., 0., 3.16522589, 5.41825866, 5.45538584, 0., 0., 0.29259347, 0., 0., 1.75312933, 0., 0., -0.56790607, 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 1.57696834, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.54855118, 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.57696834, 0., 0., 0., 0., 0., 0., 0., 1.02543666, 0., 0., 1.15041631, 0., 0., -0.54294546, 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.57696834, 0., 0., 0.90587703, 0.7165021, 0.34262858, 0., 0., 0.07438893, 0., 0., 0.56790607, 0., 0., -0.08810792, 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.57696834, 0., 0., 0., 0., 1.02543666, 0., 0., 1.15041631, 0., 0., 0., 0., 0., -0.54294546, 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.57696834, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.54855118},
        {3.01064326, 7.22843174, 10.70109626, 0., 0., 1.81824542, 0., 0., 3.16522589, 0., 0., 0.90587703, 0., 0., 5.19870815, 9.45570563, 12.35839308, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {7.22843174, 17.64603115, 27.15849891, 0., 0., 2.71979468, 0., 0., 5.41825866, 0., 0., 0.7165021, 0., 0., 9.45570563, 21.11552664, 30.43158895, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {10.70109626, 27.15849891, 43.71681781, 0., 0., 2.40501078, 0., 0., 5.45538584, 0., 0., 0.34262858, 0., 0., 12.35839308, 30.43158895, 47.51332266, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 1.40570429, 0., 0., 2.52095233, 0., 0., 0., 0., 0., 1.02543666, 0., 0., 0., 0., 2.98672814, 0., 0., 4.01618213, 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 1.40570429, 0., 0., 2.52095233, 0., 0., 1.02543666, 0., 0., 0., 0., 0., 0., 0., 2.98672814, 0., 0., 4.01618213, 0., 0., 0., 0., 0., 0.},
        {-1.81824542, -2.71979468, -2.40501078, -0., -0., -0.51875067, -0., -0., 0.29259347, 0., 0., 0.07438893, 0., 0., 0., 0., 0., 0., 0., 2.98672814, 0., 0., 4.01618213, 0., 0., 0., 0., 0.},
        {0., 0., 0., 2.52095233, 0., 0., 5.02200529, 0., 0., 0., 0., 0., 1.15041631, 0., 0., 0., 0., 4.01618213, 0., 0., 7.0373349, 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 2.52095233, 0., 0., 5.02200529, 0., 0., 1.15041631, 0., 0., 0., 0., 0., 0., 0., 4.01618213, 0., 0., 7.0373349, 0., 0., 0., 0., 0., 0.},
        {-3.16522589, -5.41825866, -5.45538584, -0., -0., 0.29259347, -0., -0., 1.75312933, 0., 0., 0.56790607, 0., 0., 0., 0., 0., 0., 0., 4.01618213, 0., 0., 7.0373349, 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.54855118, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.57696834, 0., 0., 0., 0.},
        {0., 0., 0., 0., -1.02543666, 0., 0., -1.15041631, 0., 0., -0.54294546, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.57696834, 0., 0., 0.},
        {0.90587703, 0.7165021, 0.34262858, 0., 0., -0.07438893, 0., 0. - 0.56790607, 0., 0., -0.08810792, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.57696834, 0., 0.},
        {0., 0., 0., -1.02543666, 0., 0., -1.15041631, 0., 0., 0., 0., 0., -0.54294546, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.57696834, 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.54855118, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.57696834}};
    vec3 eri3c{
        {{1.95647032, 3.94550015, 5.44434985, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.3058732, 3.16766661, 4.79834086, 0., 0., -0.6241081, 0., 0., -1.15224588, 0., 0., 0.24996002, 0., 0.},
         {1.07410572, 2.36148332, 3.3968143, 0., 0., 0.27120783, 0., 0., 0.44664828, 0., 0., 0.06726545, 0., 0., 1.07410572, 2.36148332, 3.3968143, 0., 0., -0.27120783, 0., 0., -0.44664828, 0., 0., 0.06726545, 0., 0.}},
        {{1.07410572, 2.36148332, 3.3968143, 0., 0., 0.27120783, 0., 0., 0.44664828, 0., 0., 0.06726545, 0., 0., 1.07410572, 2.36148332, 3.3968143, 0., 0., -0.27120783, 0., 0., -0.44664828, 0., 0., 0.06726545, 0., 0.},
         {1.3058732, 3.16766661, 4.79834086, 0., 0., 0.6241081, 0., 0., 1.15224588, 0., 0., 0.24996002, 0., 0., 1.95647032, 3.94550015, 5.44434985, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}}};
    
    vec result = einsum_ijk_ij_p(eri3c, dm);
    solve_linear_system(eri2c, result);

    vec result_gbw = einsum_ijk_ij_p(eri3c, wavy_gbw.get_dm());
    solve_linear_system(eri2c, result_gbw);

    vec result_molden = einsum_ijk_ij_p(eri3c, wavy.get_dm());
    solve_linear_system(eri2c, result_molden);
    vec diff_gbw(result.size());
    vec diff_molden(result.size());
    err_checkf(result.size() == result_gbw.size(), "Error, Size mismatch", std::cout);
    err_checkf(result.size() == result_molden.size(), "Error, Size mismatch", std::cout);
    for (int i = 0; i < result.size(); i++) {
        diff_gbw[i] = result[i] - result_gbw[i];
        diff_molden[i] = result[i] - result_molden[i];
    }
    std::cout << "Min/Max Difference (gbw): " << *std::min_element(diff_gbw.begin(), diff_gbw.end()) << "/" << *std::min_element(diff_gbw.begin(), diff_gbw.end()) << std::endl;
    std::cout << "Min/Max Difference (molden): " << *std::min_element(diff_molden.begin(), diff_molden.end()) << "/" << *std::min_element(diff_molden.begin(), diff_molden.end()) << std::endl;
    std::cout << "Done!" << std::endl;
#ifdef _WIN32
    math_unload_BLAS(blas);
#endif
    return 0;
}

constexpr double common_fac_sp(int l)
{
    switch (l)
    {
    case 0:
        return 0.282094791773878143;
    case 1:
        return 0.488602511902919921;
    default:
        return 1;
    }
}

inline void _g0_2d4d_0000(double *g)
{
    g[0] = 1;
    g[1] = 1;
}

inline void _g0_2d4d_0001(double *g, Rys2eT *bc)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    g[0] = 1;
    g[1] = cpx[0];
    g[2] = 1;
    g[3] = cpy[0];
    g[5] = cpz[0] * g[4];
}

inline void _g0_2d4d_0002(double *g, Rys2eT *bc)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[6] = 1;
    g[7] = 1;
    g[8] = cpy[0];
    g[9] = cpy[1];
    g[10] = cpy[0] * cpy[0] + b01[0];
    g[11] = cpy[1] * cpy[1] + b01[1];
    g[14] = cpz[0] * g[12];
    g[15] = cpz[1] * g[13];
    g[16] = cpz[0] * g[14] + b01[0] * g[12];
    g[17] = cpz[1] * g[15] + b01[1] * g[13];
}

inline void _g0_2d4d_0003(double *g, Rys2eT *bc)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[6] = cpx[0] * (g[4] + 2 * b01[0]);
    g[7] = cpx[1] * (g[5] + 2 * b01[1]);
    g[8] = 1;
    g[9] = 1;
    g[10] = cpy[0];
    g[11] = cpy[1];
    g[12] = cpy[0] * cpy[0] + b01[0];
    g[13] = cpy[1] * cpy[1] + b01[1];
    g[14] = cpy[0] * (g[12] + 2 * b01[0]);
    g[15] = cpy[1] * (g[13] + 2 * b01[1]);
    g[18] = cpz[0] * g[16];
    g[19] = cpz[1] * g[17];
    g[20] = cpz[0] * g[18] + b01[0] * g[16];
    g[21] = cpz[1] * g[19] + b01[1] * g[17];
    g[22] = cpz[0] * g[20] + 2 * b01[0] * g[18];
    g[23] = cpz[1] * g[21] + 2 * b01[1] * g[19];
}

inline void _g0_2d4d_0010(double *g, Rys2eT *bc)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    g[0] = 1;
    g[1] = cpx[0];
    g[2] = 1;
    g[3] = cpy[0];
    g[5] = cpz[0] * g[4];
}

inline void _g0_2d4d_0011(double *g, Rys2eT *bc, Env *envs)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[7] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[2] = xkxl + cpx[0];
    g[3] = xkxl + cpx[1];
    g[12] = 1;
    g[13] = 1;
    g[16] = cpy[0];
    g[17] = cpy[1];
    g[18] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[19] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[14] = ykyl + cpy[0];
    g[15] = ykyl + cpy[1];
    g[28] = cpz[0] * g[24];
    g[29] = cpz[1] * g[25];
    g[30] = g[28] * (zkzl + cpz[0]) + b01[0] * g[24];
    g[31] = g[29] * (zkzl + cpz[1]) + b01[1] * g[25];
    g[26] = g[24] * (zkzl + cpz[0]);
    g[27] = g[25] * (zkzl + cpz[1]);
}

inline void _g0_2d4d_0012(double *g, Rys2eT *bc, Env *envs)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = g[8] * (xkxl + cpx[0]) + cpx[0] * 2 * b01[0];
    g[11] = g[9] * (xkxl + cpx[1]) + cpx[1] * 2 * b01[1];
    g[6] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[7] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[2] = xkxl + cpx[0];
    g[3] = xkxl + cpx[1];
    g[16] = 1;
    g[17] = 1;
    g[20] = cpy[0];
    g[21] = cpy[1];
    g[24] = cpy[0] * cpy[0] + b01[0];
    g[25] = cpy[1] * cpy[1] + b01[1];
    g[26] = g[24] * (ykyl + cpy[0]) + cpy[0] * 2 * b01[0];
    g[27] = g[25] * (ykyl + cpy[1]) + cpy[1] * 2 * b01[1];
    g[22] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[23] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[18] = ykyl + cpy[0];
    g[19] = ykyl + cpy[1];
    g[36] = cpz[0] * g[32];
    g[37] = cpz[1] * g[33];
    g[40] = cpz[0] * g[36] + b01[0] * g[32];
    g[41] = cpz[1] * g[37] + b01[1] * g[33];
    g[42] = g[40] * (zkzl + cpz[0]) + 2 * b01[0] * g[36];
    g[43] = g[41] * (zkzl + cpz[1]) + 2 * b01[1] * g[37];
    g[38] = g[36] * (zkzl + cpz[0]) + b01[0] * g[32];
    g[39] = g[37] * (zkzl + cpz[1]) + b01[1] * g[33];
    g[34] = g[32] * (zkzl + cpz[0]);
    g[35] = g[33] * (zkzl + cpz[1]);
}

inline void _g0_2d4d_0020(double *g, Rys2eT *bc)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[6] = 1;
    g[7] = 1;
    g[8] = cpy[0];
    g[9] = cpy[1];
    g[10] = cpy[0] * cpy[0] + b01[0];
    g[11] = cpy[1] * cpy[1] + b01[1];
    g[14] = cpz[0] * g[12];
    g[15] = cpz[1] * g[13];
    g[16] = cpz[0] * g[14] + b01[0] * g[12];
    g[17] = cpz[1] * g[15] + b01[1] * g[13];
}

inline void _g0_2d4d_0021(double *g, Rys2eT *bc, Env *envs)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[8] = xkxl + cpx[0];
    g[9] = xkxl + cpx[1];
    g[10] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[11] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[12] = g[4] * (xkxl + cpx[0]) + cpx[0] * 2 * b01[0];
    g[13] = g[5] * (xkxl + cpx[1]) + cpx[1] * 2 * b01[1];
    g[16] = 1;
    g[17] = 1;
    g[18] = cpy[0];
    g[19] = cpy[1];
    g[20] = cpy[0] * cpy[0] + b01[0];
    g[21] = cpy[1] * cpy[1] + b01[1];
    g[24] = ykyl + cpy[0];
    g[25] = ykyl + cpy[1];
    g[26] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[27] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[28] = g[20] * (ykyl + cpy[0]) + cpy[0] * 2 * b01[0];
    g[29] = g[21] * (ykyl + cpy[1]) + cpy[1] * 2 * b01[1];
    g[34] = cpz[0] * g[32];
    g[35] = cpz[1] * g[33];
    g[36] = cpz[0] * g[34] + b01[0] * g[32];
    g[37] = cpz[1] * g[35] + b01[1] * g[33];
    g[40] = g[32] * (zkzl + cpz[0]);
    g[41] = g[33] * (zkzl + cpz[1]);
    g[42] = g[34] * (zkzl + cpz[0]) + b01[0] * g[32];
    g[43] = g[35] * (zkzl + cpz[1]) + b01[1] * g[33];
    g[44] = g[36] * (zkzl + cpz[0]) + 2 * b01[0] * g[34];
    g[45] = g[37] * (zkzl + cpz[1]) + 2 * b01[1] * g[35];
}

inline void _g0_2d4d_0030(double *g, Rys2eT *bc)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[6] = cpx[0] * (g[4] + 2 * b01[0]);
    g[7] = cpx[1] * (g[5] + 2 * b01[1]);
    g[8] = 1;
    g[9] = 1;
    g[10] = cpy[0];
    g[11] = cpy[1];
    g[12] = cpy[0] * cpy[0] + b01[0];
    g[13] = cpy[1] * cpy[1] + b01[1];
    g[14] = cpy[0] * (g[12] + 2 * b01[0]);
    g[15] = cpy[1] * (g[13] + 2 * b01[1]);
    g[18] = cpz[0] * g[16];
    g[19] = cpz[1] * g[17];
    g[20] = cpz[0] * g[18] + b01[0] * g[16];
    g[21] = cpz[1] * g[19] + b01[1] * g[17];
    g[22] = cpz[0] * g[20] + 2 * b01[0] * g[18];
    g[23] = cpz[1] * g[21] + 2 * b01[1] * g[19];
}

inline void _g0_2d4d_0100(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    g[0] = 1;
    g[1] = c0x[0];
    g[2] = 1;
    g[3] = c0y[0];
    g[5] = c0z[0] * g[4];
}

inline void _g0_2d4d_0101(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = 1;
    g[9] = 1;
    g[10] = cpy[0];
    g[11] = cpy[1];
    g[12] = c0y[0];
    g[13] = c0y[1];
    g[14] = cpy[0] * c0y[0] + b00[0];
    g[15] = cpy[1] * c0y[1] + b00[1];
    g[18] = cpz[0] * g[16];
    g[19] = cpz[1] * g[17];
    g[20] = c0z[0] * g[16];
    g[21] = c0z[1] * g[17];
    g[22] = cpz[0] * g[20] + b00[0] * g[16];
    g[23] = cpz[1] * g[21] + b00[1] * g[17];
}

inline void _g0_2d4d_0102(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[6] = c0x[0];
    g[7] = c0x[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[8] = cpx[0] * c0x[0] + b00[0];
    g[9] = cpx[1] * c0x[1] + b00[1];
    g[10] = cpx[0] * (g[8] + b00[0]) + b01[0] * c0x[0];
    g[11] = cpx[1] * (g[9] + b00[1]) + b01[1] * c0x[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = cpy[0];
    g[15] = cpy[1];
    g[18] = c0y[0];
    g[19] = c0y[1];
    g[16] = cpy[0] * cpy[0] + b01[0];
    g[17] = cpy[1] * cpy[1] + b01[1];
    g[20] = cpy[0] * c0y[0] + b00[0];
    g[21] = cpy[1] * c0y[1] + b00[1];
    g[22] = cpy[0] * (g[20] + b00[0]) + b01[0] * c0y[0];
    g[23] = cpy[1] * (g[21] + b00[1]) + b01[1] * c0y[1];
    g[26] = cpz[0] * g[24];
    g[27] = cpz[1] * g[25];
    g[30] = c0z[0] * g[24];
    g[31] = c0z[1] * g[25];
    g[28] = cpz[0] * g[26] + b01[0] * g[24];
    g[29] = cpz[1] * g[27] + b01[1] * g[25];
    g[32] = cpz[0] * g[30] + b00[0] * g[24];
    g[33] = cpz[1] * g[31] + b00[1] * g[25];
    g[34] = cpz[0] * g[32] + b01[0] * g[30] + b00[0] * g[26];
    g[35] = cpz[1] * g[33] + b01[1] * g[31] + b00[1] * g[27];
}

inline void _g0_2d4d_0110(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = 1;
    g[9] = 1;
    g[10] = cpy[0];
    g[11] = cpy[1];
    g[12] = c0y[0];
    g[13] = c0y[1];
    g[14] = cpy[0] * c0y[0] + b00[0];
    g[15] = cpy[1] * c0y[1] + b00[1];
    g[18] = cpz[0] * g[16];
    g[19] = cpz[1] * g[17];
    g[20] = c0z[0] * g[16];
    g[21] = c0z[1] * g[17];
    g[22] = cpz[0] * g[20] + b00[0] * g[16];
    g[23] = cpz[1] * g[21] + b00[1] * g[17];
}

inline void _g0_2d4d_0111(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[12] = c0x[0];
    g[13] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[16] = cpx[0] * c0x[0] + b00[0];
    g[17] = cpx[1] * c0x[1] + b00[1];
    g[6] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[7] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[18] = g[16] * (xkxl + cpx[0]) + cpx[0] * b00[0] + b01[0] * c0x[0];
    g[19] = g[17] * (xkxl + cpx[1]) + cpx[1] * b00[1] + b01[1] * c0x[1];
    g[2] = xkxl + cpx[0];
    g[3] = xkxl + cpx[1];
    g[14] = c0x[0] * (xkxl + cpx[0]) + b00[0];
    g[15] = c0x[1] * (xkxl + cpx[1]) + b00[1];
    g[24] = 1;
    g[25] = 1;
    g[36] = c0y[0];
    g[37] = c0y[1];
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[40] = cpy[0] * c0y[0] + b00[0];
    g[41] = cpy[1] * c0y[1] + b00[1];
    g[30] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[31] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[42] = g[40] * (ykyl + cpy[0]) + cpy[0] * b00[0] + b01[0] * c0y[0];
    g[43] = g[41] * (ykyl + cpy[1]) + cpy[1] * b00[1] + b01[1] * c0y[1];
    g[26] = ykyl + cpy[0];
    g[27] = ykyl + cpy[1];
    g[38] = c0y[0] * (ykyl + cpy[0]) + b00[0];
    g[39] = c0y[1] * (ykyl + cpy[1]) + b00[1];
    g[60] = c0z[0] * g[48];
    g[61] = c0z[1] * g[49];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[64] = cpz[0] * g[60] + b00[0] * g[48];
    g[65] = cpz[1] * g[61] + b00[1] * g[49];
    g[54] = g[52] * (zkzl + cpz[0]) + b01[0] * g[48];
    g[55] = g[53] * (zkzl + cpz[1]) + b01[1] * g[49];
    g[66] = g[64] * (zkzl + cpz[0]) + b01[0] * g[60] + b00[0] * g[52];
    g[67] = g[65] * (zkzl + cpz[1]) + b01[1] * g[61] + b00[1] * g[53];
    g[50] = g[48] * (zkzl + cpz[0]);
    g[51] = g[49] * (zkzl + cpz[1]);
    g[62] = g[60] * (zkzl + cpz[0]) + b00[0] * g[48];
    g[63] = g[61] * (zkzl + cpz[1]) + b00[1] * g[49];
}

inline void _g0_2d4d_0120(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[6] = c0x[0];
    g[7] = c0x[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[8] = cpx[0] * c0x[0] + b00[0];
    g[9] = cpx[1] * c0x[1] + b00[1];
    g[10] = cpx[0] * (g[8] + b00[0]) + b01[0] * c0x[0];
    g[11] = cpx[1] * (g[9] + b00[1]) + b01[1] * c0x[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = cpy[0];
    g[15] = cpy[1];
    g[18] = c0y[0];
    g[19] = c0y[1];
    g[16] = cpy[0] * cpy[0] + b01[0];
    g[17] = cpy[1] * cpy[1] + b01[1];
    g[20] = cpy[0] * c0y[0] + b00[0];
    g[21] = cpy[1] * c0y[1] + b00[1];
    g[22] = cpy[0] * (g[20] + b00[0]) + b01[0] * c0y[0];
    g[23] = cpy[1] * (g[21] + b00[1]) + b01[1] * c0y[1];
    g[26] = cpz[0] * g[24];
    g[27] = cpz[1] * g[25];
    g[30] = c0z[0] * g[24];
    g[31] = c0z[1] * g[25];
    g[28] = cpz[0] * g[26] + b01[0] * g[24];
    g[29] = cpz[1] * g[27] + b01[1] * g[25];
    g[32] = cpz[0] * g[30] + b00[0] * g[24];
    g[33] = cpz[1] * g[31] + b00[1] * g[25];
    g[34] = cpz[0] * g[32] + b01[0] * g[30] + b00[0] * g[26];
    g[35] = cpz[1] * g[33] + b01[1] * g[31] + b00[1] * g[27];
}

inline void _g0_2d4d_0200(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[6] = 1;
    g[7] = 1;
    g[8] = c0y[0];
    g[9] = c0y[1];
    g[10] = c0y[0] * c0y[0] + b10[0];
    g[11] = c0y[1] * c0y[1] + b10[1];
    g[14] = c0z[0] * g[12];
    g[15] = c0z[1] * g[13];
    g[16] = c0z[0] * g[14] + b10[0] * g[12];
    g[17] = c0z[1] * g[15] + b10[1] * g[13];
}

inline void _g0_2d4d_0201(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[10] = c0x[0] * (g[6] + b00[0]) + b10[0] * cpx[0];
    g[11] = c0x[1] * (g[7] + b00[1]) + b10[1] * cpx[1];
    g[12] = 1;
    g[13] = 1;
    g[16] = c0y[0];
    g[17] = c0y[1];
    g[20] = c0y[0] * c0y[0] + b10[0];
    g[21] = c0y[1] * c0y[1] + b10[1];
    g[14] = cpy[0];
    g[15] = cpy[1];
    g[18] = cpy[0] * c0y[0] + b00[0];
    g[19] = cpy[1] * c0y[1] + b00[1];
    g[22] = c0y[0] * (g[18] + b00[0]) + b10[0] * cpy[0];
    g[23] = c0y[1] * (g[19] + b00[1]) + b10[1] * cpy[1];
    g[28] = c0z[0] * g[24];
    g[29] = c0z[1] * g[25];
    g[32] = c0z[0] * g[28] + b10[0] * g[24];
    g[33] = c0z[1] * g[29] + b10[1] * g[25];
    g[26] = cpz[0] * g[24];
    g[27] = cpz[1] * g[25];
    g[30] = cpz[0] * g[28] + b00[0] * g[24];
    g[31] = cpz[1] * g[29] + b00[1] * g[25];
    g[34] = c0z[0] * g[30] + b10[0] * g[26] + b00[0] * g[28];
    g[35] = c0z[1] * g[31] + b10[1] * g[27] + b00[1] * g[29];
}

inline void _g0_2d4d_0210(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[0] * (g[6] + b00[0]) + b10[0] * cpx[0];
    g[11] = c0x[1] * (g[7] + b00[1]) + b10[1] * cpx[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = cpy[0];
    g[15] = cpy[1];
    g[16] = c0y[0];
    g[17] = c0y[1];
    g[18] = cpy[0] * c0y[0] + b00[0];
    g[19] = cpy[1] * c0y[1] + b00[1];
    g[20] = c0y[0] * c0y[0] + b10[0];
    g[21] = c0y[1] * c0y[1] + b10[1];
    g[22] = c0y[0] * (g[18] + b00[0]) + b10[0] * cpy[0];
    g[23] = c0y[1] * (g[19] + b00[1]) + b10[1] * cpy[1];
    g[26] = cpz[0] * g[24];
    g[27] = cpz[1] * g[25];
    g[28] = c0z[0] * g[24];
    g[29] = c0z[1] * g[25];
    g[30] = cpz[0] * g[28] + b00[0] * g[24];
    g[31] = cpz[1] * g[29] + b00[1] * g[25];
    g[32] = c0z[0] * g[28] + b10[0] * g[24];
    g[33] = c0z[1] * g[29] + b10[1] * g[25];
    g[34] = c0z[0] * g[30] + b10[0] * g[26] + b00[0] * g[28];
    g[35] = c0z[1] * g[31] + b10[1] * g[27] + b00[1] * g[29];
}

inline void _g0_2d4d_0300(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[6] = c0x[0] * (g[4] + 2 * b10[0]);
    g[7] = c0x[1] * (g[5] + 2 * b10[1]);
    g[8] = 1;
    g[9] = 1;
    g[10] = c0y[0];
    g[11] = c0y[1];
    g[12] = c0y[0] * c0y[0] + b10[0];
    g[13] = c0y[1] * c0y[1] + b10[1];
    g[14] = c0y[0] * (g[12] + 2 * b10[0]);
    g[15] = c0y[1] * (g[13] + 2 * b10[1]);
    g[18] = c0z[0] * g[16];
    g[19] = c0z[1] * g[17];
    g[20] = c0z[0] * g[18] + b10[0] * g[16];
    g[21] = c0z[1] * g[19] + b10[1] * g[17];
    g[22] = c0z[0] * g[20] + 2 * b10[0] * g[18];
    g[23] = c0z[1] * g[21] + 2 * b10[1] * g[19];
}

inline void _g0_2d4d_1000(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    g[0] = 1;
    g[1] = c0x[0];
    g[2] = 1;
    g[3] = c0y[0];
    g[5] = c0z[0] * g[4];
}

inline void _g0_2d4d_1001(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = 1;
    g[9] = 1;
    g[10] = c0y[0];
    g[11] = c0y[1];
    g[12] = cpy[0];
    g[13] = cpy[1];
    g[14] = cpy[0] * c0y[0] + b00[0];
    g[15] = cpy[1] * c0y[1] + b00[1];
    g[18] = c0z[0] * g[16];
    g[19] = c0z[1] * g[17];
    g[20] = cpz[0] * g[16];
    g[21] = cpz[1] * g[17];
    g[22] = cpz[0] * g[18] + b00[0] * g[16];
    g[23] = cpz[1] * g[19] + b00[1] * g[17];
}

inline void _g0_2d4d_1002(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[0] * (g[6] + b00[0]) + b01[0] * c0x[0];
    g[11] = cpx[1] * (g[7] + b00[1]) + b01[1] * c0x[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = c0y[0];
    g[15] = c0y[1];
    g[16] = cpy[0];
    g[17] = cpy[1];
    g[18] = cpy[0] * c0y[0] + b00[0];
    g[19] = cpy[1] * c0y[1] + b00[1];
    g[20] = cpy[0] * cpy[0] + b01[0];
    g[21] = cpy[1] * cpy[1] + b01[1];
    g[22] = cpy[0] * (g[18] + b00[0]) + b01[0] * c0y[0];
    g[23] = cpy[1] * (g[19] + b00[1]) + b01[1] * c0y[1];
    g[26] = c0z[0] * g[24];
    g[27] = c0z[1] * g[25];
    g[28] = cpz[0] * g[24];
    g[29] = cpz[1] * g[25];
    g[30] = cpz[0] * g[26] + b00[0] * g[24];
    g[31] = cpz[1] * g[27] + b00[1] * g[25];
    g[32] = cpz[0] * g[28] + b01[0] * g[24];
    g[33] = cpz[1] * g[29] + b01[1] * g[25];
    g[34] = cpz[0] * g[30] + b01[0] * g[26] + b00[0] * g[28];
    g[35] = cpz[1] * g[31] + b01[1] * g[27] + b00[1] * g[29];
}

inline void _g0_2d4d_1010(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = 1;
    g[9] = 1;
    g[10] = c0y[0];
    g[11] = c0y[1];
    g[12] = cpy[0];
    g[13] = cpy[1];
    g[14] = cpy[0] * c0y[0] + b00[0];
    g[15] = cpy[1] * c0y[1] + b00[1];
    g[18] = c0z[0] * g[16];
    g[19] = c0z[1] * g[17];
    g[20] = cpz[0] * g[16];
    g[21] = cpz[1] * g[17];
    g[22] = cpz[0] * g[18] + b00[0] * g[16];
    g[23] = cpz[1] * g[19] + b00[1] * g[17];
}

inline void _g0_2d4d_1011(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[0] * c0x[0] + b00[0];
    g[11] = cpx[1] * c0x[1] + b00[1];
    g[12] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[13] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[14] = g[10] * (xkxl + cpx[0]) + cpx[0] * b00[0] + b01[0] * c0x[0];
    g[15] = g[11] * (xkxl + cpx[1]) + cpx[1] * b00[1] + b01[1] * c0x[1];
    g[4] = xkxl + cpx[0];
    g[5] = xkxl + cpx[1];
    g[6] = c0x[0] * (xkxl + cpx[0]) + b00[0];
    g[7] = c0x[1] * (xkxl + cpx[1]) + b00[1];
    g[24] = 1;
    g[25] = 1;
    g[26] = c0y[0];
    g[27] = c0y[1];
    g[32] = cpy[0];
    g[33] = cpy[1];
    g[34] = cpy[0] * c0y[0] + b00[0];
    g[35] = cpy[1] * c0y[1] + b00[1];
    g[36] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[37] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[38] = g[34] * (ykyl + cpy[0]) + cpy[0] * b00[0] + b01[0] * c0y[0];
    g[39] = g[35] * (ykyl + cpy[1]) + cpy[1] * b00[1] + b01[1] * c0y[1];
    g[28] = ykyl + cpy[0];
    g[29] = ykyl + cpy[1];
    g[30] = c0y[0] * (ykyl + cpy[0]) + b00[0];
    g[31] = c0y[1] * (ykyl + cpy[1]) + b00[1];
    g[50] = c0z[0] * g[48];
    g[51] = c0z[1] * g[49];
    g[56] = cpz[0] * g[48];
    g[57] = cpz[1] * g[49];
    g[58] = cpz[0] * g[50] + b00[0] * g[48];
    g[59] = cpz[1] * g[51] + b00[1] * g[49];
    g[60] = g[56] * (zkzl + cpz[0]) + b01[0] * g[48];
    g[61] = g[57] * (zkzl + cpz[1]) + b01[1] * g[49];
    g[62] = g[58] * (zkzl + cpz[0]) + b01[0] * g[50] + b00[0] * g[56];
    g[63] = g[59] * (zkzl + cpz[1]) + b01[1] * g[51] + b00[1] * g[57];
    g[52] = g[48] * (zkzl + cpz[0]);
    g[53] = g[49] * (zkzl + cpz[1]);
    g[54] = g[50] * (zkzl + cpz[0]) + b00[0] * g[48];
    g[55] = g[51] * (zkzl + cpz[1]) + b00[1] * g[49];
}

inline void _g0_2d4d_1020(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[0] * (g[6] + b00[0]) + b01[0] * c0x[0];
    g[11] = cpx[1] * (g[7] + b00[1]) + b01[1] * c0x[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = c0y[0];
    g[15] = c0y[1];
    g[16] = cpy[0];
    g[17] = cpy[1];
    g[18] = cpy[0] * c0y[0] + b00[0];
    g[19] = cpy[1] * c0y[1] + b00[1];
    g[20] = cpy[0] * cpy[0] + b01[0];
    g[21] = cpy[1] * cpy[1] + b01[1];
    g[22] = cpy[0] * (g[18] + b00[0]) + b01[0] * c0y[0];
    g[23] = cpy[1] * (g[19] + b00[1]) + b01[1] * c0y[1];
    // g[24] = w[0];
    // g[25] = w[1];
    g[26] = c0z[0] * g[24];
    g[27] = c0z[1] * g[25];
    g[28] = cpz[0] * g[24];
    g[29] = cpz[1] * g[25];
    g[30] = cpz[0] * g[26] + b00[0] * g[24];
    g[31] = cpz[1] * g[27] + b00[1] * g[25];
    g[32] = cpz[0] * g[28] + b01[0] * g[24];
    g[33] = cpz[1] * g[29] + b01[1] * g[25];
    g[34] = cpz[0] * g[30] + b01[0] * g[26] + b00[0] * g[28];
    g[35] = cpz[1] * g[31] + b01[1] * g[27] + b00[1] * g[29];
}

inline void _g0_2d4d_1100(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[7] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[2] = xixj + c0x[0];
    g[3] = xixj + c0x[1];
    g[12] = 1;
    g[13] = 1;
    g[16] = c0y[0];
    g[17] = c0y[1];
    g[18] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[19] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[14] = yiyj + c0y[0];
    g[15] = yiyj + c0y[1];
    g[28] = c0z[0] * g[24];
    g[29] = c0z[1] * g[25];
    g[30] = g[28] * (zizj + c0z[0]) + b10[0] * g[24];
    g[31] = g[29] * (zizj + c0z[1]) + b10[1] * g[25];
    g[26] = g[24] * (zizj + c0z[0]);
    g[27] = g[25] * (zizj + c0z[1]);
}

inline void _g0_2d4d_1101(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[10] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[11] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[2] = xixj + c0x[0];
    g[3] = xixj + c0x[1];
    g[14] = g[12] * (xixj + c0x[0]) + c0x[0] * b00[0] + b10[0] * cpx[0];
    g[15] = g[13] * (xixj + c0x[1]) + c0x[1] * b00[1] + b10[1] * cpx[1];
    g[6] = cpx[0] * (xixj + c0x[0]) + b00[0];
    g[7] = cpx[1] * (xixj + c0x[1]) + b00[1];
    g[24] = 1;
    g[25] = 1;
    g[32] = c0y[0];
    g[33] = c0y[1];
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[36] = cpy[0] * c0y[0] + b00[0];
    g[37] = cpy[1] * c0y[1] + b00[1];
    g[34] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[35] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[26] = yiyj + c0y[0];
    g[27] = yiyj + c0y[1];
    g[38] = g[36] * (yiyj + c0y[0]) + c0y[0] * b00[0] + b10[0] * cpy[0];
    g[39] = g[37] * (yiyj + c0y[1]) + c0y[1] * b00[1] + b10[1] * cpy[1];
    g[30] = cpy[0] * (yiyj + c0y[0]) + b00[0];
    g[31] = cpy[1] * (yiyj + c0y[1]) + b00[1];
    g[56] = c0z[0] * g[48];
    g[57] = c0z[1] * g[49];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[60] = cpz[0] * g[56] + b00[0] * g[48];
    g[61] = cpz[1] * g[57] + b00[1] * g[49];
    g[58] = g[56] * (zizj + c0z[0]) + b10[0] * g[48];
    g[59] = g[57] * (zizj + c0z[1]) + b10[1] * g[49];
    g[50] = g[48] * (zizj + c0z[0]);
    g[51] = g[49] * (zizj + c0z[1]);
    g[62] = g[60] * (zizj + c0z[0]) + b10[0] * g[52] + b00[0] * g[56];
    g[63] = g[61] * (zizj + c0z[1]) + b10[1] * g[53] + b00[1] * g[57];
    g[54] = zizj * g[52] + cpz[0] * g[56] + b00[0] * g[48];
    g[55] = zizj * g[53] + cpz[1] * g[57] + b00[1] * g[49];
}

inline void _g0_2d4d_1110(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[10] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[11] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[2] = xixj + c0x[0];
    g[3] = xixj + c0x[1];
    g[14] = g[12] * (xixj + c0x[0]) + c0x[0] * b00[0] + b10[0] * cpx[0];
    g[15] = g[13] * (xixj + c0x[1]) + c0x[1] * b00[1] + b10[1] * cpx[1];
    g[6] = cpx[0] * (xixj + c0x[0]) + b00[0];
    g[7] = cpx[1] * (xixj + c0x[1]) + b00[1];
    g[24] = 1;
    g[25] = 1;
    g[32] = c0y[0];
    g[33] = c0y[1];
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[36] = cpy[0] * c0y[0] + b00[0];
    g[37] = cpy[1] * c0y[1] + b00[1];
    g[34] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[35] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[26] = yiyj + c0y[0];
    g[27] = yiyj + c0y[1];
    g[38] = g[36] * (yiyj + c0y[0]) + c0y[0] * b00[0] + b10[0] * cpy[0];
    g[39] = g[37] * (yiyj + c0y[1]) + c0y[1] * b00[1] + b10[1] * cpy[1];
    g[30] = cpy[0] * (yiyj + c0y[0]) + b00[0];
    g[31] = cpy[1] * (yiyj + c0y[1]) + b00[1];
    g[56] = c0z[0] * g[48];
    g[57] = c0z[1] * g[49];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[60] = cpz[0] * g[56] + b00[0] * g[48];
    g[61] = cpz[1] * g[57] + b00[1] * g[49];
    g[58] = g[56] * (zizj + c0z[0]) + b10[0] * g[48];
    g[59] = g[57] * (zizj + c0z[1]) + b10[1] * g[49];
    g[50] = g[48] * (zizj + c0z[0]);
    g[51] = g[49] * (zizj + c0z[1]);
    g[62] = g[60] * (zizj + c0z[0]) + b10[0] * g[52] + b00[0] * g[56];
    g[63] = g[61] * (zizj + c0z[1]) + b10[1] * g[53] + b00[1] * g[57];
    g[54] = zizj * g[52] + cpz[0] * g[56] + b00[0] * g[48];
    g[55] = zizj * g[53] + cpz[1] * g[57] + b00[1] * g[49];
}

inline void _g0_2d4d_1200(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = g[8] * (xixj + c0x[0]) + c0x[0] * 2 * b10[0];
    g[11] = g[9] * (xixj + c0x[1]) + c0x[1] * 2 * b10[1];
    g[6] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[7] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[2] = xixj + c0x[0];
    g[3] = xixj + c0x[1];
    g[16] = 1;
    g[17] = 1;
    g[20] = c0y[0];
    g[21] = c0y[1];
    g[24] = c0y[0] * c0y[0] + b10[0];
    g[25] = c0y[1] * c0y[1] + b10[1];
    g[26] = g[24] * (yiyj + c0y[0]) + c0y[0] * 2 * b10[0];
    g[27] = g[25] * (yiyj + c0y[1]) + c0y[1] * 2 * b10[1];
    g[22] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[23] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[18] = yiyj + c0y[0];
    g[19] = yiyj + c0y[1];
    g[36] = c0z[0] * g[32];
    g[37] = c0z[1] * g[33];
    g[40] = c0z[0] * g[36] + b10[0] * g[32];
    g[41] = c0z[1] * g[37] + b10[1] * g[33];
    g[42] = g[40] * (zizj + c0z[0]) + 2 * b10[0] * g[36];
    g[43] = g[41] * (zizj + c0z[1]) + 2 * b10[1] * g[37];
    g[38] = g[36] * (zizj + c0z[0]) + b10[0] * g[32];
    g[39] = g[37] * (zizj + c0z[1]) + b10[1] * g[33];
    g[34] = g[32] * (zizj + c0z[0]);
    g[35] = g[33] * (zizj + c0z[1]);
}

inline void _g0_2d4d_2000(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[6] = 1;
    g[7] = 1;
    g[8] = c0y[0];
    g[9] = c0y[1];
    g[10] = c0y[0] * c0y[0] + b10[0];
    g[11] = c0y[1] * c0y[1] + b10[1];
    g[14] = c0z[0] * g[12];
    g[15] = c0z[1] * g[13];
    g[16] = c0z[0] * g[14] + b10[0] * g[12];
    g[17] = c0z[1] * g[15] + b10[1] * g[13];
}

inline void _g0_2d4d_2001(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[6] = cpx[0];
    g[7] = cpx[1];
    g[8] = cpx[0] * c0x[0] + b00[0];
    g[9] = cpx[1] * c0x[1] + b00[1];
    g[10] = c0x[0] * (g[8] + b00[0]) + b10[0] * cpx[0];
    g[11] = c0x[1] * (g[9] + b00[1]) + b10[1] * cpx[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = c0y[0];
    g[15] = c0y[1];
    g[16] = c0y[0] * c0y[0] + b10[0];
    g[17] = c0y[1] * c0y[1] + b10[1];
    g[18] = cpy[0];
    g[19] = cpy[1];
    g[20] = cpy[0] * c0y[0] + b00[0];
    g[21] = cpy[1] * c0y[1] + b00[1];
    g[22] = c0y[0] * (g[20] + b00[0]) + b10[0] * cpy[0];
    g[23] = c0y[1] * (g[21] + b00[1]) + b10[1] * cpy[1];
    g[26] = c0z[0] * g[24];
    g[27] = c0z[1] * g[25];
    g[28] = c0z[0] * g[26] + b10[0] * g[24];
    g[29] = c0z[1] * g[27] + b10[1] * g[25];
    g[30] = cpz[0] * g[24];
    g[31] = cpz[1] * g[25];
    g[32] = cpz[0] * g[26] + b00[0] * g[24];
    g[33] = cpz[1] * g[27] + b00[1] * g[25];
    g[34] = c0z[0] * g[32] + b10[0] * g[30] + b00[0] * g[26];
    g[35] = c0z[1] * g[33] + b10[1] * g[31] + b00[1] * g[27];
}

inline void _g0_2d4d_2010(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[6] = cpx[0];
    g[7] = cpx[1];
    g[8] = cpx[0] * c0x[0] + b00[0];
    g[9] = cpx[1] * c0x[1] + b00[1];
    g[10] = c0x[0] * (g[8] + b00[0]) + b10[0] * cpx[0];
    g[11] = c0x[1] * (g[9] + b00[1]) + b10[1] * cpx[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = c0y[0];
    g[15] = c0y[1];
    g[16] = c0y[0] * c0y[0] + b10[0];
    g[17] = c0y[1] * c0y[1] + b10[1];
    g[18] = cpy[0];
    g[19] = cpy[1];
    g[20] = cpy[0] * c0y[0] + b00[0];
    g[21] = cpy[1] * c0y[1] + b00[1];
    g[22] = c0y[0] * (g[20] + b00[0]) + b10[0] * cpy[0];
    g[23] = c0y[1] * (g[21] + b00[1]) + b10[1] * cpy[1];
    g[26] = c0z[0] * g[24];
    g[27] = c0z[1] * g[25];
    g[28] = c0z[0] * g[26] + b10[0] * g[24];
    g[29] = c0z[1] * g[27] + b10[1] * g[25];
    g[30] = cpz[0] * g[24];
    g[31] = cpz[1] * g[25];
    g[32] = cpz[0] * g[26] + b00[0] * g[24];
    g[33] = cpz[1] * g[27] + b00[1] * g[25];
    g[34] = c0z[0] * g[32] + b10[0] * g[30] + b00[0] * g[26];
    g[35] = c0z[1] * g[33] + b10[1] * g[31] + b00[1] * g[27];
}

inline void _g0_2d4d_2100(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[12] = g[4] * (xixj + c0x[0]) + c0x[0] * 2 * b10[0];
    g[13] = g[5] * (xixj + c0x[1]) + c0x[1] * 2 * b10[1];
    g[10] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[11] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[8] = xixj + c0x[0];
    g[9] = xixj + c0x[1];
    g[16] = 1;
    g[17] = 1;
    g[18] = c0y[0];
    g[19] = c0y[1];
    g[20] = c0y[0] * c0y[0] + b10[0];
    g[21] = c0y[1] * c0y[1] + b10[1];
    g[28] = g[20] * (yiyj + c0y[0]) + c0y[0] * 2 * b10[0];
    g[29] = g[21] * (yiyj + c0y[1]) + c0y[1] * 2 * b10[1];
    g[26] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[27] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[24] = yiyj + c0y[0];
    g[25] = yiyj + c0y[1];
    g[34] = c0z[0] * g[32];
    g[35] = c0z[1] * g[33];
    g[36] = c0z[0] * g[34] + b10[0] * g[32];
    g[37] = c0z[1] * g[35] + b10[1] * g[33];
    g[44] = g[36] * (zizj + c0z[0]) + 2 * b10[0] * g[34];
    g[45] = g[37] * (zizj + c0z[1]) + 2 * b10[1] * g[35];
    g[42] = g[34] * (zizj + c0z[0]) + b10[0] * g[32];
    g[43] = g[35] * (zizj + c0z[1]) + b10[1] * g[33];
    g[40] = g[32] * (zizj + c0z[0]);
    g[41] = g[33] * (zizj + c0z[1]);
}

inline void _g0_2d4d_3000(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[6] = c0x[0] * (g[4] + 2 * b10[0]);
    g[7] = c0x[1] * (g[5] + 2 * b10[1]);
    g[8] = 1;
    g[9] = 1;
    g[10] = c0y[0];
    g[11] = c0y[1];
    g[12] = c0y[0] * c0y[0] + b10[0];
    g[13] = c0y[1] * c0y[1] + b10[1];
    g[14] = c0y[0] * (g[12] + 2 * b10[0]);
    g[15] = c0y[1] * (g[13] + 2 * b10[1]);
    g[18] = c0z[0] * g[16];
    g[19] = c0z[1] * g[17];
    g[20] = c0z[0] * g[18] + b10[0] * g[16];
    g[21] = c0z[1] * g[19] + b10[1] * g[17];
    g[22] = c0z[0] * g[20] + 2 * b10[0] * g[18];
    g[23] = c0z[1] * g[21] + 2 * b10[1] * g[19];
}

void g0_2e_2d4d_unrolled(double *g, Rys2eT *bc, Env *envs)
{
    int type_ijkl = ((envs->li_ceil << 6) | (envs->lj_ceil << 4) |
                     (envs->lk_ceil << 2) | (envs->ll_ceil));
    switch (type_ijkl)
    {
    case 0b00000000:
        _g0_2d4d_0000(g);
        return;
    case 0b00000001:
        _g0_2d4d_0001(g, bc);
        return;
    case 0b00000010:
        _g0_2d4d_0002(g, bc);
        return;
    case 0b00000011:
        _g0_2d4d_0003(g, bc);
        return;
    case 0b00000100:
        _g0_2d4d_0010(g, bc);
        return;
    case 0b00000101:
        _g0_2d4d_0011(g, bc, envs);
        return;
    case 0b00000110:
        _g0_2d4d_0012(g, bc, envs);
        return;
    case 0b00001000:
        _g0_2d4d_0020(g, bc);
        return;
    case 0b00001001:
        _g0_2d4d_0021(g, bc, envs);
        return;
    case 0b00001100:
        _g0_2d4d_0030(g, bc);
        return;
    case 0b00010000:
        _g0_2d4d_0100(g, bc);
        return;
    case 0b00010001:
        _g0_2d4d_0101(g, bc);
        return;
    case 0b00010010:
        _g0_2d4d_0102(g, bc);
        return;
    case 0b00010100:
        _g0_2d4d_0110(g, bc);
        return;
    case 0b00010101:
        _g0_2d4d_0111(g, bc, envs);
        return;
    case 0b00011000:
        _g0_2d4d_0120(g, bc);
        return;
    case 0b00100000:
        _g0_2d4d_0200(g, bc);
        return;
    case 0b00100001:
        _g0_2d4d_0201(g, bc);
        return;
    case 0b00100100:
        _g0_2d4d_0210(g, bc);
        return;
    case 0b00110000:
        _g0_2d4d_0300(g, bc);
        return;
    case 0b01000000:
        _g0_2d4d_1000(g, bc);
        return;
    case 0b01000001:
        _g0_2d4d_1001(g, bc);
        return;
    case 0b01000010:
        _g0_2d4d_1002(g, bc);
        return;
    case 0b01000100:
        _g0_2d4d_1010(g, bc);
        return;
    case 0b01000101:
        _g0_2d4d_1011(g, bc, envs);
        return;
    case 0b01001000:
        _g0_2d4d_1020(g, bc);
        return;
    case 0b01010000:
        _g0_2d4d_1100(g, bc, envs);
        return;
    case 0b01010001:
        _g0_2d4d_1101(g, bc, envs);
        return;
    case 0b01010100:
        _g0_2d4d_1110(g, bc, envs);
        return;
    case 0b01100000:
        _g0_2d4d_1200(g, bc, envs);
        return;
    case 0b10000000:
        _g0_2d4d_2000(g, bc);
        return;
    case 0b10000001:
        _g0_2d4d_2001(g, bc);
        return;
    case 0b10000100:
        _g0_2d4d_2010(g, bc);
        return;
    case 0b10010000:
        _g0_2d4d_2100(g, bc, envs);
        return;
    case 0b11000000:
        _g0_2d4d_3000(g, bc);
        return;
    }
    fprintf(stderr, "Dimension error for CINTg0_2e_lj2d4d: iklj = %d %d %d %d",
            (int)envs->li_ceil, (int)envs->lk_ceil,
            (int)envs->ll_ceil, (int)envs->lj_ceil);
}

inline void _srg0_2d4d_0000(double *g)
{
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
}

inline void _srg0_2d4d_0001(double *g, Rys2eT *bc)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = 1;
    g[5] = 1;
    g[6] = cpy[0];
    g[7] = cpy[1];
    g[10] = cpz[0] * g[8];
    g[11] = cpz[1] * g[9];
}

inline void _srg0_2d4d_0002(double *g, Rys2eT *bc)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[12] = 1;
    g[13] = 1;
    g[14] = 1;
    g[15] = 1;
    g[16] = cpy[0];
    g[17] = cpy[1];
    g[18] = cpy[2];
    g[19] = cpy[3];
    g[20] = cpy[0] * cpy[0] + b01[0];
    g[21] = cpy[1] * cpy[1] + b01[1];
    g[22] = cpy[2] * cpy[2] + b01[2];
    g[23] = cpy[3] * cpy[3] + b01[3];
    g[28] = cpz[0] * g[24];
    g[29] = cpz[1] * g[25];
    g[30] = cpz[2] * g[26];
    g[31] = cpz[3] * g[27];
    g[32] = cpz[0] * g[28] + b01[0] * g[24];
    g[33] = cpz[1] * g[29] + b01[1] * g[25];
    g[34] = cpz[2] * g[30] + b01[2] * g[26];
    g[35] = cpz[3] * g[31] + b01[3] * g[27];
}

inline void _srg0_2d4d_0003(double *g, Rys2eT *bc)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[12] = cpx[0] * (g[8] + 2 * b01[0]);
    g[13] = cpx[1] * (g[9] + 2 * b01[1]);
    g[14] = cpx[2] * (g[10] + 2 * b01[2]);
    g[15] = cpx[3] * (g[11] + 2 * b01[3]);
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = cpy[0];
    g[21] = cpy[1];
    g[22] = cpy[2];
    g[23] = cpy[3];
    g[24] = cpy[0] * cpy[0] + b01[0];
    g[25] = cpy[1] * cpy[1] + b01[1];
    g[26] = cpy[2] * cpy[2] + b01[2];
    g[27] = cpy[3] * cpy[3] + b01[3];
    g[28] = cpy[0] * (g[24] + 2 * b01[0]);
    g[29] = cpy[1] * (g[25] + 2 * b01[1]);
    g[30] = cpy[2] * (g[26] + 2 * b01[2]);
    g[31] = cpy[3] * (g[27] + 2 * b01[3]);
    g[36] = cpz[0] * g[32];
    g[37] = cpz[1] * g[33];
    g[38] = cpz[2] * g[34];
    g[39] = cpz[3] * g[35];
    g[40] = cpz[0] * g[36] + b01[0] * g[32];
    g[41] = cpz[1] * g[37] + b01[1] * g[33];
    g[42] = cpz[2] * g[38] + b01[2] * g[34];
    g[43] = cpz[3] * g[39] + b01[3] * g[35];
    g[44] = cpz[0] * g[40] + 2 * b01[0] * g[36];
    g[45] = cpz[1] * g[41] + 2 * b01[1] * g[37];
    g[46] = cpz[2] * g[42] + 2 * b01[2] * g[38];
    g[47] = cpz[3] * g[43] + 2 * b01[3] * g[39];
}

inline void _srg0_2d4d_0010(double *g, Rys2eT *bc)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = 1;
    g[5] = 1;
    g[6] = cpy[0];
    g[7] = cpy[1];
    g[10] = cpz[0] * g[8];
    g[11] = cpz[1] * g[9];
}

inline void _srg0_2d4d_0011(double *g, Rys2eT *bc, Env *envs)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[12] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[13] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[14] = cpx[2] * (xkxl + cpx[2]) + b01[2];
    g[15] = cpx[3] * (xkxl + cpx[3]) + b01[3];
    g[4] = xkxl + cpx[0];
    g[5] = xkxl + cpx[1];
    g[6] = xkxl + cpx[2];
    g[7] = xkxl + cpx[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[32] = cpy[0];
    g[33] = cpy[1];
    g[34] = cpy[2];
    g[35] = cpy[3];
    g[36] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[37] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[38] = cpy[2] * (ykyl + cpy[2]) + b01[2];
    g[39] = cpy[3] * (ykyl + cpy[3]) + b01[3];
    g[28] = ykyl + cpy[0];
    g[29] = ykyl + cpy[1];
    g[30] = ykyl + cpy[2];
    g[31] = ykyl + cpy[3];
    g[56] = cpz[0] * g[48];
    g[57] = cpz[1] * g[49];
    g[58] = cpz[2] * g[50];
    g[59] = cpz[3] * g[51];
    g[60] = g[56] * (zkzl + cpz[0]) + b01[0] * g[48];
    g[61] = g[57] * (zkzl + cpz[1]) + b01[1] * g[49];
    g[62] = g[58] * (zkzl + cpz[2]) + b01[2] * g[50];
    g[63] = g[59] * (zkzl + cpz[3]) + b01[3] * g[51];
    g[52] = g[48] * (zkzl + cpz[0]);
    g[53] = g[49] * (zkzl + cpz[1]);
    g[54] = g[50] * (zkzl + cpz[2]);
    g[55] = g[51] * (zkzl + cpz[3]);
}

inline void _srg0_2d4d_0012(double *g, Rys2eT *bc, Env *envs)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[16] = cpx[0] * cpx[0] + b01[0];
    g[17] = cpx[1] * cpx[1] + b01[1];
    g[18] = cpx[2] * cpx[2] + b01[2];
    g[19] = cpx[3] * cpx[3] + b01[3];
    g[20] = g[16] * (xkxl + cpx[0]) + cpx[0] * 2 * b01[0];
    g[21] = g[17] * (xkxl + cpx[1]) + cpx[1] * 2 * b01[1];
    g[22] = g[18] * (xkxl + cpx[2]) + cpx[2] * 2 * b01[2];
    g[23] = g[19] * (xkxl + cpx[3]) + cpx[3] * 2 * b01[3];
    g[12] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[13] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[14] = cpx[2] * (xkxl + cpx[2]) + b01[2];
    g[15] = cpx[3] * (xkxl + cpx[3]) + b01[3];
    g[4] = xkxl + cpx[0];
    g[5] = xkxl + cpx[1];
    g[6] = xkxl + cpx[2];
    g[7] = xkxl + cpx[3];
    g[32] = 1;
    g[33] = 1;
    g[34] = 1;
    g[35] = 1;
    g[40] = cpy[0];
    g[41] = cpy[1];
    g[42] = cpy[2];
    g[43] = cpy[3];
    g[48] = cpy[0] * cpy[0] + b01[0];
    g[49] = cpy[1] * cpy[1] + b01[1];
    g[50] = cpy[2] * cpy[2] + b01[2];
    g[51] = cpy[3] * cpy[3] + b01[3];
    g[52] = g[48] * (ykyl + cpy[0]) + cpy[0] * 2 * b01[0];
    g[53] = g[49] * (ykyl + cpy[1]) + cpy[1] * 2 * b01[1];
    g[54] = g[50] * (ykyl + cpy[2]) + cpy[2] * 2 * b01[2];
    g[55] = g[51] * (ykyl + cpy[3]) + cpy[3] * 2 * b01[3];
    g[44] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[45] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[46] = cpy[2] * (ykyl + cpy[2]) + b01[2];
    g[47] = cpy[3] * (ykyl + cpy[3]) + b01[3];
    g[36] = ykyl + cpy[0];
    g[37] = ykyl + cpy[1];
    g[38] = ykyl + cpy[2];
    g[39] = ykyl + cpy[3];
    g[72] = cpz[0] * g[64];
    g[73] = cpz[1] * g[65];
    g[74] = cpz[2] * g[66];
    g[75] = cpz[3] * g[67];
    g[80] = cpz[0] * g[72] + b01[0] * g[64];
    g[81] = cpz[1] * g[73] + b01[1] * g[65];
    g[82] = cpz[2] * g[74] + b01[2] * g[66];
    g[83] = cpz[3] * g[75] + b01[3] * g[67];
    g[84] = g[80] * (zkzl + cpz[0]) + 2 * b01[0] * g[72];
    g[85] = g[81] * (zkzl + cpz[1]) + 2 * b01[1] * g[73];
    g[86] = g[82] * (zkzl + cpz[2]) + 2 * b01[2] * g[74];
    g[87] = g[83] * (zkzl + cpz[3]) + 2 * b01[3] * g[75];
    g[76] = g[72] * (zkzl + cpz[0]) + b01[0] * g[64];
    g[77] = g[73] * (zkzl + cpz[1]) + b01[1] * g[65];
    g[78] = g[74] * (zkzl + cpz[2]) + b01[2] * g[66];
    g[79] = g[75] * (zkzl + cpz[3]) + b01[3] * g[67];
    g[68] = g[64] * (zkzl + cpz[0]);
    g[69] = g[65] * (zkzl + cpz[1]);
    g[70] = g[66] * (zkzl + cpz[2]);
    g[71] = g[67] * (zkzl + cpz[3]);
}

inline void _srg0_2d4d_0020(double *g, Rys2eT *bc)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[12] = 1;
    g[13] = 1;
    g[14] = 1;
    g[15] = 1;
    g[16] = cpy[0];
    g[17] = cpy[1];
    g[18] = cpy[2];
    g[19] = cpy[3];
    g[20] = cpy[0] * cpy[0] + b01[0];
    g[21] = cpy[1] * cpy[1] + b01[1];
    g[22] = cpy[2] * cpy[2] + b01[2];
    g[23] = cpy[3] * cpy[3] + b01[3];
    g[28] = cpz[0] * g[24];
    g[29] = cpz[1] * g[25];
    g[30] = cpz[2] * g[26];
    g[31] = cpz[3] * g[27];
    g[32] = cpz[0] * g[28] + b01[0] * g[24];
    g[33] = cpz[1] * g[29] + b01[1] * g[25];
    g[34] = cpz[2] * g[30] + b01[2] * g[26];
    g[35] = cpz[3] * g[31] + b01[3] * g[27];
}

inline void _srg0_2d4d_0021(double *g, Rys2eT *bc, Env *envs)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[16] = xkxl + cpx[0];
    g[17] = xkxl + cpx[1];
    g[18] = xkxl + cpx[2];
    g[19] = xkxl + cpx[3];
    g[20] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[21] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[22] = cpx[2] * (xkxl + cpx[2]) + b01[2];
    g[23] = cpx[3] * (xkxl + cpx[3]) + b01[3];
    g[24] = g[8] * (xkxl + cpx[0]) + cpx[0] * 2 * b01[0];
    g[25] = g[9] * (xkxl + cpx[1]) + cpx[1] * 2 * b01[1];
    g[26] = g[10] * (xkxl + cpx[2]) + cpx[2] * 2 * b01[2];
    g[27] = g[11] * (xkxl + cpx[3]) + cpx[3] * 2 * b01[3];
    g[32] = 1;
    g[33] = 1;
    g[34] = 1;
    g[35] = 1;
    g[36] = cpy[0];
    g[37] = cpy[1];
    g[38] = cpy[2];
    g[39] = cpy[3];
    g[40] = cpy[0] * cpy[0] + b01[0];
    g[41] = cpy[1] * cpy[1] + b01[1];
    g[42] = cpy[2] * cpy[2] + b01[2];
    g[43] = cpy[3] * cpy[3] + b01[3];
    g[48] = ykyl + cpy[0];
    g[49] = ykyl + cpy[1];
    g[50] = ykyl + cpy[2];
    g[51] = ykyl + cpy[3];
    g[52] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[53] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[54] = cpy[2] * (ykyl + cpy[2]) + b01[2];
    g[55] = cpy[3] * (ykyl + cpy[3]) + b01[3];
    g[56] = g[40] * (ykyl + cpy[0]) + cpy[0] * 2 * b01[0];
    g[57] = g[41] * (ykyl + cpy[1]) + cpy[1] * 2 * b01[1];
    g[58] = g[42] * (ykyl + cpy[2]) + cpy[2] * 2 * b01[2];
    g[59] = g[43] * (ykyl + cpy[3]) + cpy[3] * 2 * b01[3];
    g[68] = cpz[0] * g[64];
    g[69] = cpz[1] * g[65];
    g[70] = cpz[2] * g[66];
    g[71] = cpz[3] * g[67];
    g[72] = cpz[0] * g[68] + b01[0] * g[64];
    g[73] = cpz[1] * g[69] + b01[1] * g[65];
    g[74] = cpz[2] * g[70] + b01[2] * g[66];
    g[75] = cpz[3] * g[71] + b01[3] * g[67];
    g[80] = g[64] * (zkzl + cpz[0]);
    g[81] = g[65] * (zkzl + cpz[1]);
    g[82] = g[66] * (zkzl + cpz[2]);
    g[83] = g[67] * (zkzl + cpz[3]);
    g[84] = g[68] * (zkzl + cpz[0]) + b01[0] * g[64];
    g[85] = g[69] * (zkzl + cpz[1]) + b01[1] * g[65];
    g[86] = g[70] * (zkzl + cpz[2]) + b01[2] * g[66];
    g[87] = g[71] * (zkzl + cpz[3]) + b01[3] * g[67];
    g[88] = g[72] * (zkzl + cpz[0]) + 2 * b01[0] * g[68];
    g[89] = g[73] * (zkzl + cpz[1]) + 2 * b01[1] * g[69];
    g[90] = g[74] * (zkzl + cpz[2]) + 2 * b01[2] * g[70];
    g[91] = g[75] * (zkzl + cpz[3]) + 2 * b01[3] * g[71];
}

inline void _srg0_2d4d_0030(double *g, Rys2eT *bc)
{
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[12] = cpx[0] * (g[8] + 2 * b01[0]);
    g[13] = cpx[1] * (g[9] + 2 * b01[1]);
    g[14] = cpx[2] * (g[10] + 2 * b01[2]);
    g[15] = cpx[3] * (g[11] + 2 * b01[3]);
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = cpy[0];
    g[21] = cpy[1];
    g[22] = cpy[2];
    g[23] = cpy[3];
    g[24] = cpy[0] * cpy[0] + b01[0];
    g[25] = cpy[1] * cpy[1] + b01[1];
    g[26] = cpy[2] * cpy[2] + b01[2];
    g[27] = cpy[3] * cpy[3] + b01[3];
    g[28] = cpy[0] * (g[24] + 2 * b01[0]);
    g[29] = cpy[1] * (g[25] + 2 * b01[1]);
    g[30] = cpy[2] * (g[26] + 2 * b01[2]);
    g[31] = cpy[3] * (g[27] + 2 * b01[3]);
    g[36] = cpz[0] * g[32];
    g[37] = cpz[1] * g[33];
    g[38] = cpz[2] * g[34];
    g[39] = cpz[3] * g[35];
    g[40] = cpz[0] * g[36] + b01[0] * g[32];
    g[41] = cpz[1] * g[37] + b01[1] * g[33];
    g[42] = cpz[2] * g[38] + b01[2] * g[34];
    g[43] = cpz[3] * g[39] + b01[3] * g[35];
    g[44] = cpz[0] * g[40] + 2 * b01[0] * g[36];
    g[45] = cpz[1] * g[41] + 2 * b01[1] * g[37];
    g[46] = cpz[2] * g[42] + 2 * b01[2] * g[38];
    g[47] = cpz[3] * g[43] + 2 * b01[3] * g[39];
}

inline void _srg0_2d4d_0100(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = 1;
    g[5] = 1;
    g[6] = c0y[0];
    g[7] = c0y[1];
    g[10] = c0z[0] * g[8];
    g[11] = c0z[1] * g[9];
}

inline void _srg0_2d4d_0101(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[10] = c0x[2];
    g[11] = c0x[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = cpy[0];
    g[21] = cpy[1];
    g[22] = cpy[2];
    g[23] = cpy[3];
    g[24] = c0y[0];
    g[25] = c0y[1];
    g[26] = c0y[2];
    g[27] = c0y[3];
    g[28] = cpy[0] * c0y[0] + b00[0];
    g[29] = cpy[1] * c0y[1] + b00[1];
    g[30] = cpy[2] * c0y[2] + b00[2];
    g[31] = cpy[3] * c0y[3] + b00[3];
    g[36] = cpz[0] * g[32];
    g[37] = cpz[1] * g[33];
    g[38] = cpz[2] * g[34];
    g[39] = cpz[3] * g[35];
    g[40] = c0z[0] * g[32];
    g[41] = c0z[1] * g[33];
    g[42] = c0z[2] * g[34];
    g[43] = c0z[3] * g[35];
    g[44] = cpz[0] * g[40] + b00[0] * g[32];
    g[45] = cpz[1] * g[41] + b00[1] * g[33];
    g[46] = cpz[2] * g[42] + b00[2] * g[34];
    g[47] = cpz[3] * g[43] + b00[3] * g[35];
}

inline void _srg0_2d4d_0102(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[12] = c0x[0];
    g[13] = c0x[1];
    g[14] = c0x[2];
    g[15] = c0x[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[16] = cpx[0] * c0x[0] + b00[0];
    g[17] = cpx[1] * c0x[1] + b00[1];
    g[18] = cpx[2] * c0x[2] + b00[2];
    g[19] = cpx[3] * c0x[3] + b00[3];
    g[20] = cpx[0] * (g[16] + b00[0]) + b01[0] * c0x[0];
    g[21] = cpx[1] * (g[17] + b00[1]) + b01[1] * c0x[1];
    g[22] = cpx[2] * (g[18] + b00[2]) + b01[2] * c0x[2];
    g[23] = cpx[3] * (g[19] + b00[3]) + b01[3] * c0x[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[30] = cpy[2];
    g[31] = cpy[3];
    g[36] = c0y[0];
    g[37] = c0y[1];
    g[38] = c0y[2];
    g[39] = c0y[3];
    g[32] = cpy[0] * cpy[0] + b01[0];
    g[33] = cpy[1] * cpy[1] + b01[1];
    g[34] = cpy[2] * cpy[2] + b01[2];
    g[35] = cpy[3] * cpy[3] + b01[3];
    g[40] = cpy[0] * c0y[0] + b00[0];
    g[41] = cpy[1] * c0y[1] + b00[1];
    g[42] = cpy[2] * c0y[2] + b00[2];
    g[43] = cpy[3] * c0y[3] + b00[3];
    g[44] = cpy[0] * (g[40] + b00[0]) + b01[0] * c0y[0];
    g[45] = cpy[1] * (g[41] + b00[1]) + b01[1] * c0y[1];
    g[46] = cpy[2] * (g[42] + b00[2]) + b01[2] * c0y[2];
    g[47] = cpy[3] * (g[43] + b00[3]) + b01[3] * c0y[3];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[54] = cpz[2] * g[50];
    g[55] = cpz[3] * g[51];
    g[60] = c0z[0] * g[48];
    g[61] = c0z[1] * g[49];
    g[62] = c0z[2] * g[50];
    g[63] = c0z[3] * g[51];
    g[56] = cpz[0] * g[52] + b01[0] * g[48];
    g[57] = cpz[1] * g[53] + b01[1] * g[49];
    g[58] = cpz[2] * g[54] + b01[2] * g[50];
    g[59] = cpz[3] * g[55] + b01[3] * g[51];
    g[64] = cpz[0] * g[60] + b00[0] * g[48];
    g[65] = cpz[1] * g[61] + b00[1] * g[49];
    g[66] = cpz[2] * g[62] + b00[2] * g[50];
    g[67] = cpz[3] * g[63] + b00[3] * g[51];
    g[68] = cpz[0] * g[64] + b01[0] * g[60] + b00[0] * g[52];
    g[69] = cpz[1] * g[65] + b01[1] * g[61] + b00[1] * g[53];
    g[70] = cpz[2] * g[66] + b01[2] * g[62] + b00[2] * g[54];
    g[71] = cpz[3] * g[67] + b01[3] * g[63] + b00[3] * g[55];
}

inline void _srg0_2d4d_0110(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[10] = c0x[2];
    g[11] = c0x[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = cpy[0];
    g[21] = cpy[1];
    g[22] = cpy[2];
    g[23] = cpy[3];
    g[24] = c0y[0];
    g[25] = c0y[1];
    g[26] = c0y[2];
    g[27] = c0y[3];
    g[28] = cpy[0] * c0y[0] + b00[0];
    g[29] = cpy[1] * c0y[1] + b00[1];
    g[30] = cpy[2] * c0y[2] + b00[2];
    g[31] = cpy[3] * c0y[3] + b00[3];
    g[36] = cpz[0] * g[32];
    g[37] = cpz[1] * g[33];
    g[38] = cpz[2] * g[34];
    g[39] = cpz[3] * g[35];
    g[40] = c0z[0] * g[32];
    g[41] = c0z[1] * g[33];
    g[42] = c0z[2] * g[34];
    g[43] = c0z[3] * g[35];
    g[44] = cpz[0] * g[40] + b00[0] * g[32];
    g[45] = cpz[1] * g[41] + b00[1] * g[33];
    g[46] = cpz[2] * g[42] + b00[2] * g[34];
    g[47] = cpz[3] * g[43] + b00[3] * g[35];
}

inline void _srg0_2d4d_0111(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[24] = c0x[0];
    g[25] = c0x[1];
    g[26] = c0x[2];
    g[27] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[32] = cpx[0] * c0x[0] + b00[0];
    g[33] = cpx[1] * c0x[1] + b00[1];
    g[34] = cpx[2] * c0x[2] + b00[2];
    g[35] = cpx[3] * c0x[3] + b00[3];
    g[12] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[13] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[14] = cpx[2] * (xkxl + cpx[2]) + b01[2];
    g[15] = cpx[3] * (xkxl + cpx[3]) + b01[3];
    g[36] = g[32] * (xkxl + cpx[0]) + cpx[0] * b00[0] + b01[0] * c0x[0];
    g[37] = g[33] * (xkxl + cpx[1]) + cpx[1] * b00[1] + b01[1] * c0x[1];
    g[38] = g[34] * (xkxl + cpx[2]) + cpx[2] * b00[2] + b01[2] * c0x[2];
    g[39] = g[35] * (xkxl + cpx[3]) + cpx[3] * b00[3] + b01[3] * c0x[3];
    g[4] = xkxl + cpx[0];
    g[5] = xkxl + cpx[1];
    g[6] = xkxl + cpx[2];
    g[7] = xkxl + cpx[3];
    g[28] = c0x[0] * (xkxl + cpx[0]) + b00[0];
    g[29] = c0x[1] * (xkxl + cpx[1]) + b00[1];
    g[30] = c0x[2] * (xkxl + cpx[2]) + b00[2];
    g[31] = c0x[3] * (xkxl + cpx[3]) + b00[3];
    g[48] = 1;
    g[49] = 1;
    g[50] = 1;
    g[51] = 1;
    g[72] = c0y[0];
    g[73] = c0y[1];
    g[74] = c0y[2];
    g[75] = c0y[3];
    g[56] = cpy[0];
    g[57] = cpy[1];
    g[58] = cpy[2];
    g[59] = cpy[3];
    g[80] = cpy[0] * c0y[0] + b00[0];
    g[81] = cpy[1] * c0y[1] + b00[1];
    g[82] = cpy[2] * c0y[2] + b00[2];
    g[83] = cpy[3] * c0y[3] + b00[3];
    g[60] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[61] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[62] = cpy[2] * (ykyl + cpy[2]) + b01[2];
    g[63] = cpy[3] * (ykyl + cpy[3]) + b01[3];
    g[84] = g[80] * (ykyl + cpy[0]) + cpy[0] * b00[0] + b01[0] * c0y[0];
    g[85] = g[81] * (ykyl + cpy[1]) + cpy[1] * b00[1] + b01[1] * c0y[1];
    g[86] = g[82] * (ykyl + cpy[2]) + cpy[2] * b00[2] + b01[2] * c0y[2];
    g[87] = g[83] * (ykyl + cpy[3]) + cpy[3] * b00[3] + b01[3] * c0y[3];
    g[52] = ykyl + cpy[0];
    g[53] = ykyl + cpy[1];
    g[54] = ykyl + cpy[2];
    g[55] = ykyl + cpy[3];
    g[76] = c0y[0] * (ykyl + cpy[0]) + b00[0];
    g[77] = c0y[1] * (ykyl + cpy[1]) + b00[1];
    g[78] = c0y[2] * (ykyl + cpy[2]) + b00[2];
    g[79] = c0y[3] * (ykyl + cpy[3]) + b00[3];
    g[120] = c0z[0] * g[96];
    g[121] = c0z[1] * g[97];
    g[122] = c0z[2] * g[98];
    g[123] = c0z[3] * g[99];
    g[104] = cpz[0] * g[96];
    g[105] = cpz[1] * g[97];
    g[106] = cpz[2] * g[98];
    g[107] = cpz[3] * g[99];
    g[128] = cpz[0] * g[120] + b00[0] * g[96];
    g[129] = cpz[1] * g[121] + b00[1] * g[97];
    g[130] = cpz[2] * g[122] + b00[2] * g[98];
    g[131] = cpz[3] * g[123] + b00[3] * g[99];
    g[108] = g[104] * (zkzl + cpz[0]) + b01[0] * g[96];
    g[109] = g[105] * (zkzl + cpz[1]) + b01[1] * g[97];
    g[110] = g[106] * (zkzl + cpz[2]) + b01[2] * g[98];
    g[111] = g[107] * (zkzl + cpz[3]) + b01[3] * g[99];
    g[132] = g[128] * (zkzl + cpz[0]) + b01[0] * g[120] + b00[0] * g[104];
    g[133] = g[129] * (zkzl + cpz[1]) + b01[1] * g[121] + b00[1] * g[105];
    g[134] = g[130] * (zkzl + cpz[2]) + b01[2] * g[122] + b00[2] * g[106];
    g[135] = g[131] * (zkzl + cpz[3]) + b01[3] * g[123] + b00[3] * g[107];
    g[100] = g[96] * (zkzl + cpz[0]);
    g[101] = g[97] * (zkzl + cpz[1]);
    g[102] = g[98] * (zkzl + cpz[2]);
    g[103] = g[99] * (zkzl + cpz[3]);
    g[124] = g[120] * (zkzl + cpz[0]) + b00[0] * g[96];
    g[125] = g[121] * (zkzl + cpz[1]) + b00[1] * g[97];
    g[126] = g[122] * (zkzl + cpz[2]) + b00[2] * g[98];
    g[127] = g[123] * (zkzl + cpz[3]) + b00[3] * g[99];
}

inline void _srg0_2d4d_0120(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[12] = c0x[0];
    g[13] = c0x[1];
    g[14] = c0x[2];
    g[15] = c0x[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[16] = cpx[0] * c0x[0] + b00[0];
    g[17] = cpx[1] * c0x[1] + b00[1];
    g[18] = cpx[2] * c0x[2] + b00[2];
    g[19] = cpx[3] * c0x[3] + b00[3];
    g[20] = cpx[0] * (g[16] + b00[0]) + b01[0] * c0x[0];
    g[21] = cpx[1] * (g[17] + b00[1]) + b01[1] * c0x[1];
    g[22] = cpx[2] * (g[18] + b00[2]) + b01[2] * c0x[2];
    g[23] = cpx[3] * (g[19] + b00[3]) + b01[3] * c0x[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[30] = cpy[2];
    g[31] = cpy[3];
    g[36] = c0y[0];
    g[37] = c0y[1];
    g[38] = c0y[2];
    g[39] = c0y[3];
    g[32] = cpy[0] * cpy[0] + b01[0];
    g[33] = cpy[1] * cpy[1] + b01[1];
    g[34] = cpy[2] * cpy[2] + b01[2];
    g[35] = cpy[3] * cpy[3] + b01[3];
    g[40] = cpy[0] * c0y[0] + b00[0];
    g[41] = cpy[1] * c0y[1] + b00[1];
    g[42] = cpy[2] * c0y[2] + b00[2];
    g[43] = cpy[3] * c0y[3] + b00[3];
    g[44] = cpy[0] * (g[40] + b00[0]) + b01[0] * c0y[0];
    g[45] = cpy[1] * (g[41] + b00[1]) + b01[1] * c0y[1];
    g[46] = cpy[2] * (g[42] + b00[2]) + b01[2] * c0y[2];
    g[47] = cpy[3] * (g[43] + b00[3]) + b01[3] * c0y[3];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[54] = cpz[2] * g[50];
    g[55] = cpz[3] * g[51];
    g[60] = c0z[0] * g[48];
    g[61] = c0z[1] * g[49];
    g[62] = c0z[2] * g[50];
    g[63] = c0z[3] * g[51];
    g[56] = cpz[0] * g[52] + b01[0] * g[48];
    g[57] = cpz[1] * g[53] + b01[1] * g[49];
    g[58] = cpz[2] * g[54] + b01[2] * g[50];
    g[59] = cpz[3] * g[55] + b01[3] * g[51];
    g[64] = cpz[0] * g[60] + b00[0] * g[48];
    g[65] = cpz[1] * g[61] + b00[1] * g[49];
    g[66] = cpz[2] * g[62] + b00[2] * g[50];
    g[67] = cpz[3] * g[63] + b00[3] * g[51];
    g[68] = cpz[0] * g[64] + b01[0] * g[60] + b00[0] * g[52];
    g[69] = cpz[1] * g[65] + b01[1] * g[61] + b00[1] * g[53];
    g[70] = cpz[2] * g[66] + b01[2] * g[62] + b00[2] * g[54];
    g[71] = cpz[3] * g[67] + b01[3] * g[63] + b00[3] * g[55];
}

inline void _srg0_2d4d_0200(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[12] = 1;
    g[13] = 1;
    g[14] = 1;
    g[15] = 1;
    g[16] = c0y[0];
    g[17] = c0y[1];
    g[18] = c0y[2];
    g[19] = c0y[3];
    g[20] = c0y[0] * c0y[0] + b10[0];
    g[21] = c0y[1] * c0y[1] + b10[1];
    g[22] = c0y[2] * c0y[2] + b10[2];
    g[23] = c0y[3] * c0y[3] + b10[3];
    g[28] = c0z[0] * g[24];
    g[29] = c0z[1] * g[25];
    g[30] = c0z[2] * g[26];
    g[31] = c0z[3] * g[27];
    g[32] = c0z[0] * g[28] + b10[0] * g[24];
    g[33] = c0z[1] * g[29] + b10[1] * g[25];
    g[34] = c0z[2] * g[30] + b10[2] * g[26];
    g[35] = c0z[3] * g[31] + b10[3] * g[27];
}

inline void _srg0_2d4d_0201(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[10] = c0x[2];
    g[11] = c0x[3];
    g[16] = c0x[0] * c0x[0] + b10[0];
    g[17] = c0x[1] * c0x[1] + b10[1];
    g[18] = c0x[2] * c0x[2] + b10[2];
    g[19] = c0x[3] * c0x[3] + b10[3];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[20] = c0x[0] * (g[12] + b00[0]) + b10[0] * cpx[0];
    g[21] = c0x[1] * (g[13] + b00[1]) + b10[1] * cpx[1];
    g[22] = c0x[2] * (g[14] + b00[2]) + b10[2] * cpx[2];
    g[23] = c0x[3] * (g[15] + b00[3]) + b10[3] * cpx[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[32] = c0y[0];
    g[33] = c0y[1];
    g[34] = c0y[2];
    g[35] = c0y[3];
    g[40] = c0y[0] * c0y[0] + b10[0];
    g[41] = c0y[1] * c0y[1] + b10[1];
    g[42] = c0y[2] * c0y[2] + b10[2];
    g[43] = c0y[3] * c0y[3] + b10[3];
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[30] = cpy[2];
    g[31] = cpy[3];
    g[36] = cpy[0] * c0y[0] + b00[0];
    g[37] = cpy[1] * c0y[1] + b00[1];
    g[38] = cpy[2] * c0y[2] + b00[2];
    g[39] = cpy[3] * c0y[3] + b00[3];
    g[44] = c0y[0] * (g[36] + b00[0]) + b10[0] * cpy[0];
    g[45] = c0y[1] * (g[37] + b00[1]) + b10[1] * cpy[1];
    g[46] = c0y[2] * (g[38] + b00[2]) + b10[2] * cpy[2];
    g[47] = c0y[3] * (g[39] + b00[3]) + b10[3] * cpy[3];
    g[56] = c0z[0] * g[48];
    g[57] = c0z[1] * g[49];
    g[58] = c0z[2] * g[50];
    g[59] = c0z[3] * g[51];
    g[64] = c0z[0] * g[56] + b10[0] * g[48];
    g[65] = c0z[1] * g[57] + b10[1] * g[49];
    g[66] = c0z[2] * g[58] + b10[2] * g[50];
    g[67] = c0z[3] * g[59] + b10[3] * g[51];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[54] = cpz[2] * g[50];
    g[55] = cpz[3] * g[51];
    g[60] = cpz[0] * g[56] + b00[0] * g[48];
    g[61] = cpz[1] * g[57] + b00[1] * g[49];
    g[62] = cpz[2] * g[58] + b00[2] * g[50];
    g[63] = cpz[3] * g[59] + b00[3] * g[51];
    g[68] = c0z[0] * g[60] + b10[0] * g[52] + b00[0] * g[56];
    g[69] = c0z[1] * g[61] + b10[1] * g[53] + b00[1] * g[57];
    g[70] = c0z[2] * g[62] + b10[2] * g[54] + b00[2] * g[58];
    g[71] = c0z[3] * g[63] + b10[3] * g[55] + b00[3] * g[59];
}

inline void _srg0_2d4d_0210(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[10] = c0x[2];
    g[11] = c0x[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = c0x[0] * c0x[0] + b10[0];
    g[17] = c0x[1] * c0x[1] + b10[1];
    g[18] = c0x[2] * c0x[2] + b10[2];
    g[19] = c0x[3] * c0x[3] + b10[3];
    g[20] = c0x[0] * (g[12] + b00[0]) + b10[0] * cpx[0];
    g[21] = c0x[1] * (g[13] + b00[1]) + b10[1] * cpx[1];
    g[22] = c0x[2] * (g[14] + b00[2]) + b10[2] * cpx[2];
    g[23] = c0x[3] * (g[15] + b00[3]) + b10[3] * cpx[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[30] = cpy[2];
    g[31] = cpy[3];
    g[32] = c0y[0];
    g[33] = c0y[1];
    g[34] = c0y[2];
    g[35] = c0y[3];
    g[36] = cpy[0] * c0y[0] + b00[0];
    g[37] = cpy[1] * c0y[1] + b00[1];
    g[38] = cpy[2] * c0y[2] + b00[2];
    g[39] = cpy[3] * c0y[3] + b00[3];
    g[40] = c0y[0] * c0y[0] + b10[0];
    g[41] = c0y[1] * c0y[1] + b10[1];
    g[42] = c0y[2] * c0y[2] + b10[2];
    g[43] = c0y[3] * c0y[3] + b10[3];
    g[44] = c0y[0] * (g[36] + b00[0]) + b10[0] * cpy[0];
    g[45] = c0y[1] * (g[37] + b00[1]) + b10[1] * cpy[1];
    g[46] = c0y[2] * (g[38] + b00[2]) + b10[2] * cpy[2];
    g[47] = c0y[3] * (g[39] + b00[3]) + b10[3] * cpy[3];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[54] = cpz[2] * g[50];
    g[55] = cpz[3] * g[51];
    g[56] = c0z[0] * g[48];
    g[57] = c0z[1] * g[49];
    g[58] = c0z[2] * g[50];
    g[59] = c0z[3] * g[51];
    g[60] = cpz[0] * g[56] + b00[0] * g[48];
    g[61] = cpz[1] * g[57] + b00[1] * g[49];
    g[62] = cpz[2] * g[58] + b00[2] * g[50];
    g[63] = cpz[3] * g[59] + b00[3] * g[51];
    g[64] = c0z[0] * g[56] + b10[0] * g[48];
    g[65] = c0z[1] * g[57] + b10[1] * g[49];
    g[66] = c0z[2] * g[58] + b10[2] * g[50];
    g[67] = c0z[3] * g[59] + b10[3] * g[51];
    g[68] = c0z[0] * g[60] + b10[0] * g[52] + b00[0] * g[56];
    g[69] = c0z[1] * g[61] + b10[1] * g[53] + b00[1] * g[57];
    g[70] = c0z[2] * g[62] + b10[2] * g[54] + b00[2] * g[58];
    g[71] = c0z[3] * g[63] + b10[3] * g[55] + b00[3] * g[59];
}

inline void _srg0_2d4d_0300(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[12] = c0x[0] * (g[8] + 2 * b10[0]);
    g[13] = c0x[1] * (g[9] + 2 * b10[1]);
    g[14] = c0x[2] * (g[10] + 2 * b10[2]);
    g[15] = c0x[3] * (g[11] + 2 * b10[3]);
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = c0y[0];
    g[21] = c0y[1];
    g[22] = c0y[2];
    g[23] = c0y[3];
    g[24] = c0y[0] * c0y[0] + b10[0];
    g[25] = c0y[1] * c0y[1] + b10[1];
    g[26] = c0y[2] * c0y[2] + b10[2];
    g[27] = c0y[3] * c0y[3] + b10[3];
    g[28] = c0y[0] * (g[24] + 2 * b10[0]);
    g[29] = c0y[1] * (g[25] + 2 * b10[1]);
    g[30] = c0y[2] * (g[26] + 2 * b10[2]);
    g[31] = c0y[3] * (g[27] + 2 * b10[3]);
    g[36] = c0z[0] * g[32];
    g[37] = c0z[1] * g[33];
    g[38] = c0z[2] * g[34];
    g[39] = c0z[3] * g[35];
    g[40] = c0z[0] * g[36] + b10[0] * g[32];
    g[41] = c0z[1] * g[37] + b10[1] * g[33];
    g[42] = c0z[2] * g[38] + b10[2] * g[34];
    g[43] = c0z[3] * g[39] + b10[3] * g[35];
    g[44] = c0z[0] * g[40] + 2 * b10[0] * g[36];
    g[45] = c0z[1] * g[41] + 2 * b10[1] * g[37];
    g[46] = c0z[2] * g[42] + 2 * b10[2] * g[38];
    g[47] = c0z[3] * g[43] + 2 * b10[3] * g[39];
}

inline void _srg0_2d4d_1000(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = 1;
    g[5] = 1;
    g[6] = c0y[0];
    g[7] = c0y[1];
    g[10] = c0z[0] * g[8];
    g[11] = c0z[1] * g[9];
}

inline void _srg0_2d4d_1001(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = c0y[0];
    g[21] = c0y[1];
    g[22] = c0y[2];
    g[23] = c0y[3];
    g[24] = cpy[0];
    g[25] = cpy[1];
    g[26] = cpy[2];
    g[27] = cpy[3];
    g[28] = cpy[0] * c0y[0] + b00[0];
    g[29] = cpy[1] * c0y[1] + b00[1];
    g[30] = cpy[2] * c0y[2] + b00[2];
    g[31] = cpy[3] * c0y[3] + b00[3];
    g[36] = c0z[0] * g[32];
    g[37] = c0z[1] * g[33];
    g[38] = c0z[2] * g[34];
    g[39] = c0z[3] * g[35];
    g[40] = cpz[0] * g[32];
    g[41] = cpz[1] * g[33];
    g[42] = cpz[2] * g[34];
    g[43] = cpz[3] * g[35];
    g[44] = cpz[0] * g[36] + b00[0] * g[32];
    g[45] = cpz[1] * g[37] + b00[1] * g[33];
    g[46] = cpz[2] * g[38] + b00[2] * g[34];
    g[47] = cpz[3] * g[39] + b00[3] * g[35];
}

inline void _srg0_2d4d_1002(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = cpx[0] * cpx[0] + b01[0];
    g[17] = cpx[1] * cpx[1] + b01[1];
    g[18] = cpx[2] * cpx[2] + b01[2];
    g[19] = cpx[3] * cpx[3] + b01[3];
    g[20] = cpx[0] * (g[12] + b00[0]) + b01[0] * c0x[0];
    g[21] = cpx[1] * (g[13] + b00[1]) + b01[1] * c0x[1];
    g[22] = cpx[2] * (g[14] + b00[2]) + b01[2] * c0x[2];
    g[23] = cpx[3] * (g[15] + b00[3]) + b01[3] * c0x[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = c0y[0];
    g[29] = c0y[1];
    g[30] = c0y[2];
    g[31] = c0y[3];
    g[32] = cpy[0];
    g[33] = cpy[1];
    g[34] = cpy[2];
    g[35] = cpy[3];
    g[36] = cpy[0] * c0y[0] + b00[0];
    g[37] = cpy[1] * c0y[1] + b00[1];
    g[38] = cpy[2] * c0y[2] + b00[2];
    g[39] = cpy[3] * c0y[3] + b00[3];
    g[40] = cpy[0] * cpy[0] + b01[0];
    g[41] = cpy[1] * cpy[1] + b01[1];
    g[42] = cpy[2] * cpy[2] + b01[2];
    g[43] = cpy[3] * cpy[3] + b01[3];
    g[44] = cpy[0] * (g[36] + b00[0]) + b01[0] * c0y[0];
    g[45] = cpy[1] * (g[37] + b00[1]) + b01[1] * c0y[1];
    g[46] = cpy[2] * (g[38] + b00[2]) + b01[2] * c0y[2];
    g[47] = cpy[3] * (g[39] + b00[3]) + b01[3] * c0y[3];
    g[52] = c0z[0] * g[48];
    g[53] = c0z[1] * g[49];
    g[54] = c0z[2] * g[50];
    g[55] = c0z[3] * g[51];
    g[56] = cpz[0] * g[48];
    g[57] = cpz[1] * g[49];
    g[58] = cpz[2] * g[50];
    g[59] = cpz[3] * g[51];
    g[60] = cpz[0] * g[52] + b00[0] * g[48];
    g[61] = cpz[1] * g[53] + b00[1] * g[49];
    g[62] = cpz[2] * g[54] + b00[2] * g[50];
    g[63] = cpz[3] * g[55] + b00[3] * g[51];
    g[64] = cpz[0] * g[56] + b01[0] * g[48];
    g[65] = cpz[1] * g[57] + b01[1] * g[49];
    g[66] = cpz[2] * g[58] + b01[2] * g[50];
    g[67] = cpz[3] * g[59] + b01[3] * g[51];
    g[68] = cpz[0] * g[60] + b01[0] * g[52] + b00[0] * g[56];
    g[69] = cpz[1] * g[61] + b01[1] * g[53] + b00[1] * g[57];
    g[70] = cpz[2] * g[62] + b01[2] * g[54] + b00[2] * g[58];
    g[71] = cpz[3] * g[63] + b01[3] * g[55] + b00[3] * g[59];
}

inline void _srg0_2d4d_1010(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = c0y[0];
    g[21] = c0y[1];
    g[22] = c0y[2];
    g[23] = c0y[3];
    g[24] = cpy[0];
    g[25] = cpy[1];
    g[26] = cpy[2];
    g[27] = cpy[3];
    g[28] = cpy[0] * c0y[0] + b00[0];
    g[29] = cpy[1] * c0y[1] + b00[1];
    g[30] = cpy[2] * c0y[2] + b00[2];
    g[31] = cpy[3] * c0y[3] + b00[3];
    g[36] = c0z[0] * g[32];
    g[37] = c0z[1] * g[33];
    g[38] = c0z[2] * g[34];
    g[39] = c0z[3] * g[35];
    g[40] = cpz[0] * g[32];
    g[41] = cpz[1] * g[33];
    g[42] = cpz[2] * g[34];
    g[43] = cpz[3] * g[35];
    g[44] = cpz[0] * g[36] + b00[0] * g[32];
    g[45] = cpz[1] * g[37] + b00[1] * g[33];
    g[46] = cpz[2] * g[38] + b00[2] * g[34];
    g[47] = cpz[3] * g[39] + b00[3] * g[35];
}

inline void _srg0_2d4d_1011(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[16] = cpx[0];
    g[17] = cpx[1];
    g[18] = cpx[2];
    g[19] = cpx[3];
    g[20] = cpx[0] * c0x[0] + b00[0];
    g[21] = cpx[1] * c0x[1] + b00[1];
    g[22] = cpx[2] * c0x[2] + b00[2];
    g[23] = cpx[3] * c0x[3] + b00[3];
    g[24] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[25] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[26] = cpx[2] * (xkxl + cpx[2]) + b01[2];
    g[27] = cpx[3] * (xkxl + cpx[3]) + b01[3];
    g[28] = g[20] * (xkxl + cpx[0]) + cpx[0] * b00[0] + b01[0] * c0x[0];
    g[29] = g[21] * (xkxl + cpx[1]) + cpx[1] * b00[1] + b01[1] * c0x[1];
    g[30] = g[22] * (xkxl + cpx[2]) + cpx[2] * b00[2] + b01[2] * c0x[2];
    g[31] = g[23] * (xkxl + cpx[3]) + cpx[3] * b00[3] + b01[3] * c0x[3];
    g[8] = xkxl + cpx[0];
    g[9] = xkxl + cpx[1];
    g[10] = xkxl + cpx[2];
    g[11] = xkxl + cpx[3];
    g[12] = c0x[0] * (xkxl + cpx[0]) + b00[0];
    g[13] = c0x[1] * (xkxl + cpx[1]) + b00[1];
    g[14] = c0x[2] * (xkxl + cpx[2]) + b00[2];
    g[15] = c0x[3] * (xkxl + cpx[3]) + b00[3];
    g[48] = 1;
    g[49] = 1;
    g[50] = 1;
    g[51] = 1;
    g[52] = c0y[0];
    g[53] = c0y[1];
    g[54] = c0y[2];
    g[55] = c0y[3];
    g[64] = cpy[0];
    g[65] = cpy[1];
    g[66] = cpy[2];
    g[67] = cpy[3];
    g[68] = cpy[0] * c0y[0] + b00[0];
    g[69] = cpy[1] * c0y[1] + b00[1];
    g[70] = cpy[2] * c0y[2] + b00[2];
    g[71] = cpy[3] * c0y[3] + b00[3];
    g[72] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[73] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[74] = cpy[2] * (ykyl + cpy[2]) + b01[2];
    g[75] = cpy[3] * (ykyl + cpy[3]) + b01[3];
    g[76] = g[68] * (ykyl + cpy[0]) + cpy[0] * b00[0] + b01[0] * c0y[0];
    g[77] = g[69] * (ykyl + cpy[1]) + cpy[1] * b00[1] + b01[1] * c0y[1];
    g[78] = g[70] * (ykyl + cpy[2]) + cpy[2] * b00[2] + b01[2] * c0y[2];
    g[79] = g[71] * (ykyl + cpy[3]) + cpy[3] * b00[3] + b01[3] * c0y[3];
    g[56] = ykyl + cpy[0];
    g[57] = ykyl + cpy[1];
    g[58] = ykyl + cpy[2];
    g[59] = ykyl + cpy[3];
    g[60] = c0y[0] * (ykyl + cpy[0]) + b00[0];
    g[61] = c0y[1] * (ykyl + cpy[1]) + b00[1];
    g[62] = c0y[2] * (ykyl + cpy[2]) + b00[2];
    g[63] = c0y[3] * (ykyl + cpy[3]) + b00[3];
    g[100] = c0z[0] * g[96];
    g[101] = c0z[1] * g[97];
    g[102] = c0z[2] * g[98];
    g[103] = c0z[3] * g[99];
    g[112] = cpz[0] * g[96];
    g[113] = cpz[1] * g[97];
    g[114] = cpz[2] * g[98];
    g[115] = cpz[3] * g[99];
    g[116] = cpz[0] * g[100] + b00[0] * g[96];
    g[117] = cpz[1] * g[101] + b00[1] * g[97];
    g[118] = cpz[2] * g[102] + b00[2] * g[98];
    g[119] = cpz[3] * g[103] + b00[3] * g[99];
    g[120] = g[112] * (zkzl + cpz[0]) + b01[0] * g[96];
    g[121] = g[113] * (zkzl + cpz[1]) + b01[1] * g[97];
    g[122] = g[114] * (zkzl + cpz[2]) + b01[2] * g[98];
    g[123] = g[115] * (zkzl + cpz[3]) + b01[3] * g[99];
    g[124] = g[116] * (zkzl + cpz[0]) + b01[0] * g[100] + b00[0] * g[112];
    g[125] = g[117] * (zkzl + cpz[1]) + b01[1] * g[101] + b00[1] * g[113];
    g[126] = g[118] * (zkzl + cpz[2]) + b01[2] * g[102] + b00[2] * g[114];
    g[127] = g[119] * (zkzl + cpz[3]) + b01[3] * g[103] + b00[3] * g[115];
    g[104] = g[96] * (zkzl + cpz[0]);
    g[105] = g[97] * (zkzl + cpz[1]);
    g[106] = g[98] * (zkzl + cpz[2]);
    g[107] = g[99] * (zkzl + cpz[3]);
    g[108] = g[100] * (zkzl + cpz[0]) + b00[0] * g[96];
    g[109] = g[101] * (zkzl + cpz[1]) + b00[1] * g[97];
    g[110] = g[102] * (zkzl + cpz[2]) + b00[2] * g[98];
    g[111] = g[103] * (zkzl + cpz[3]) + b00[3] * g[99];
}

inline void _srg0_2d4d_1020(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = cpx[0] * cpx[0] + b01[0];
    g[17] = cpx[1] * cpx[1] + b01[1];
    g[18] = cpx[2] * cpx[2] + b01[2];
    g[19] = cpx[3] * cpx[3] + b01[3];
    g[20] = cpx[0] * (g[12] + b00[0]) + b01[0] * c0x[0];
    g[21] = cpx[1] * (g[13] + b00[1]) + b01[1] * c0x[1];
    g[22] = cpx[2] * (g[14] + b00[2]) + b01[2] * c0x[2];
    g[23] = cpx[3] * (g[15] + b00[3]) + b01[3] * c0x[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = c0y[0];
    g[29] = c0y[1];
    g[30] = c0y[2];
    g[31] = c0y[3];
    g[32] = cpy[0];
    g[33] = cpy[1];
    g[34] = cpy[2];
    g[35] = cpy[3];
    g[36] = cpy[0] * c0y[0] + b00[0];
    g[37] = cpy[1] * c0y[1] + b00[1];
    g[38] = cpy[2] * c0y[2] + b00[2];
    g[39] = cpy[3] * c0y[3] + b00[3];
    g[40] = cpy[0] * cpy[0] + b01[0];
    g[41] = cpy[1] * cpy[1] + b01[1];
    g[42] = cpy[2] * cpy[2] + b01[2];
    g[43] = cpy[3] * cpy[3] + b01[3];
    g[44] = cpy[0] * (g[36] + b00[0]) + b01[0] * c0y[0];
    g[45] = cpy[1] * (g[37] + b00[1]) + b01[1] * c0y[1];
    g[46] = cpy[2] * (g[38] + b00[2]) + b01[2] * c0y[2];
    g[47] = cpy[3] * (g[39] + b00[3]) + b01[3] * c0y[3];
    g[52] = c0z[0] * g[48];
    g[53] = c0z[1] * g[49];
    g[54] = c0z[2] * g[50];
    g[55] = c0z[3] * g[51];
    g[56] = cpz[0] * g[48];
    g[57] = cpz[1] * g[49];
    g[58] = cpz[2] * g[50];
    g[59] = cpz[3] * g[51];
    g[60] = cpz[0] * g[52] + b00[0] * g[48];
    g[61] = cpz[1] * g[53] + b00[1] * g[49];
    g[62] = cpz[2] * g[54] + b00[2] * g[50];
    g[63] = cpz[3] * g[55] + b00[3] * g[51];
    g[64] = cpz[0] * g[56] + b01[0] * g[48];
    g[65] = cpz[1] * g[57] + b01[1] * g[49];
    g[66] = cpz[2] * g[58] + b01[2] * g[50];
    g[67] = cpz[3] * g[59] + b01[3] * g[51];
    g[68] = cpz[0] * g[60] + b01[0] * g[52] + b00[0] * g[56];
    g[69] = cpz[1] * g[61] + b01[1] * g[53] + b00[1] * g[57];
    g[70] = cpz[2] * g[62] + b01[2] * g[54] + b00[2] * g[58];
    g[71] = cpz[3] * g[63] + b01[3] * g[55] + b00[3] * g[59];
}

inline void _srg0_2d4d_1100(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[10] = c0x[2];
    g[11] = c0x[3];
    g[12] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[13] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[14] = c0x[2] * (xixj + c0x[2]) + b10[2];
    g[15] = c0x[3] * (xixj + c0x[3]) + b10[3];
    g[4] = xixj + c0x[0];
    g[5] = xixj + c0x[1];
    g[6] = xixj + c0x[2];
    g[7] = xixj + c0x[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[32] = c0y[0];
    g[33] = c0y[1];
    g[34] = c0y[2];
    g[35] = c0y[3];
    g[36] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[37] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[38] = c0y[2] * (yiyj + c0y[2]) + b10[2];
    g[39] = c0y[3] * (yiyj + c0y[3]) + b10[3];
    g[28] = yiyj + c0y[0];
    g[29] = yiyj + c0y[1];
    g[30] = yiyj + c0y[2];
    g[31] = yiyj + c0y[3];
    g[56] = c0z[0] * g[48];
    g[57] = c0z[1] * g[49];
    g[58] = c0z[2] * g[50];
    g[59] = c0z[3] * g[51];
    g[60] = g[56] * (zizj + c0z[0]) + b10[0] * g[48];
    g[61] = g[57] * (zizj + c0z[1]) + b10[1] * g[49];
    g[62] = g[58] * (zizj + c0z[2]) + b10[2] * g[50];
    g[63] = g[59] * (zizj + c0z[3]) + b10[3] * g[51];
    g[52] = g[48] * (zizj + c0z[0]);
    g[53] = g[49] * (zizj + c0z[1]);
    g[54] = g[50] * (zizj + c0z[2]);
    g[55] = g[51] * (zizj + c0z[3]);
}

inline void _srg0_2d4d_1101(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[16] = c0x[0];
    g[17] = c0x[1];
    g[18] = c0x[2];
    g[19] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[24] = cpx[0] * c0x[0] + b00[0];
    g[25] = cpx[1] * c0x[1] + b00[1];
    g[26] = cpx[2] * c0x[2] + b00[2];
    g[27] = cpx[3] * c0x[3] + b00[3];
    g[20] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[21] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[22] = c0x[2] * (xixj + c0x[2]) + b10[2];
    g[23] = c0x[3] * (xixj + c0x[3]) + b10[3];
    g[4] = xixj + c0x[0];
    g[5] = xixj + c0x[1];
    g[6] = xixj + c0x[2];
    g[7] = xixj + c0x[3];
    g[28] = g[24] * (xixj + c0x[0]) + c0x[0] * b00[0] + b10[0] * cpx[0];
    g[29] = g[25] * (xixj + c0x[1]) + c0x[1] * b00[1] + b10[1] * cpx[1];
    g[30] = g[26] * (xixj + c0x[2]) + c0x[2] * b00[2] + b10[2] * cpx[2];
    g[31] = g[27] * (xixj + c0x[3]) + c0x[3] * b00[3] + b10[3] * cpx[3];
    g[12] = cpx[0] * (xixj + c0x[0]) + b00[0];
    g[13] = cpx[1] * (xixj + c0x[1]) + b00[1];
    g[14] = cpx[2] * (xixj + c0x[2]) + b00[2];
    g[15] = cpx[3] * (xixj + c0x[3]) + b00[3];
    g[48] = 1;
    g[49] = 1;
    g[50] = 1;
    g[51] = 1;
    g[64] = c0y[0];
    g[65] = c0y[1];
    g[66] = c0y[2];
    g[67] = c0y[3];
    g[56] = cpy[0];
    g[57] = cpy[1];
    g[58] = cpy[2];
    g[59] = cpy[3];
    g[72] = cpy[0] * c0y[0] + b00[0];
    g[73] = cpy[1] * c0y[1] + b00[1];
    g[74] = cpy[2] * c0y[2] + b00[2];
    g[75] = cpy[3] * c0y[3] + b00[3];
    g[68] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[69] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[70] = c0y[2] * (yiyj + c0y[2]) + b10[2];
    g[71] = c0y[3] * (yiyj + c0y[3]) + b10[3];
    g[52] = yiyj + c0y[0];
    g[53] = yiyj + c0y[1];
    g[54] = yiyj + c0y[2];
    g[55] = yiyj + c0y[3];
    g[76] = g[72] * (yiyj + c0y[0]) + c0y[0] * b00[0] + b10[0] * cpy[0];
    g[77] = g[73] * (yiyj + c0y[1]) + c0y[1] * b00[1] + b10[1] * cpy[1];
    g[78] = g[74] * (yiyj + c0y[2]) + c0y[2] * b00[2] + b10[2] * cpy[2];
    g[79] = g[75] * (yiyj + c0y[3]) + c0y[3] * b00[3] + b10[3] * cpy[3];
    g[60] = cpy[0] * (yiyj + c0y[0]) + b00[0];
    g[61] = cpy[1] * (yiyj + c0y[1]) + b00[1];
    g[62] = cpy[2] * (yiyj + c0y[2]) + b00[2];
    g[63] = cpy[3] * (yiyj + c0y[3]) + b00[3];
    g[112] = c0z[0] * g[96];
    g[113] = c0z[1] * g[97];
    g[114] = c0z[2] * g[98];
    g[115] = c0z[3] * g[99];
    g[104] = cpz[0] * g[96];
    g[105] = cpz[1] * g[97];
    g[106] = cpz[2] * g[98];
    g[107] = cpz[3] * g[99];
    g[120] = cpz[0] * g[112] + b00[0] * g[96];
    g[121] = cpz[1] * g[113] + b00[1] * g[97];
    g[122] = cpz[2] * g[114] + b00[2] * g[98];
    g[123] = cpz[3] * g[115] + b00[3] * g[99];
    g[116] = g[112] * (zizj + c0z[0]) + b10[0] * g[96];
    g[117] = g[113] * (zizj + c0z[1]) + b10[1] * g[97];
    g[118] = g[114] * (zizj + c0z[2]) + b10[2] * g[98];
    g[119] = g[115] * (zizj + c0z[3]) + b10[3] * g[99];
    g[100] = g[96] * (zizj + c0z[0]);
    g[101] = g[97] * (zizj + c0z[1]);
    g[102] = g[98] * (zizj + c0z[2]);
    g[103] = g[99] * (zizj + c0z[3]);
    g[124] = g[120] * (zizj + c0z[0]) + b10[0] * g[104] + b00[0] * g[112];
    g[125] = g[121] * (zizj + c0z[1]) + b10[1] * g[105] + b00[1] * g[113];
    g[126] = g[122] * (zizj + c0z[2]) + b10[2] * g[106] + b00[2] * g[114];
    g[127] = g[123] * (zizj + c0z[3]) + b10[3] * g[107] + b00[3] * g[115];
    g[108] = zizj * g[104] + cpz[0] * g[112] + b00[0] * g[96];
    g[109] = zizj * g[105] + cpz[1] * g[113] + b00[1] * g[97];
    g[110] = zizj * g[106] + cpz[2] * g[114] + b00[2] * g[98];
    g[111] = zizj * g[107] + cpz[3] * g[115] + b00[3] * g[99];
}

inline void _srg0_2d4d_1110(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[16] = c0x[0];
    g[17] = c0x[1];
    g[18] = c0x[2];
    g[19] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[24] = cpx[0] * c0x[0] + b00[0];
    g[25] = cpx[1] * c0x[1] + b00[1];
    g[26] = cpx[2] * c0x[2] + b00[2];
    g[27] = cpx[3] * c0x[3] + b00[3];
    g[20] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[21] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[22] = c0x[2] * (xixj + c0x[2]) + b10[2];
    g[23] = c0x[3] * (xixj + c0x[3]) + b10[3];
    g[4] = xixj + c0x[0];
    g[5] = xixj + c0x[1];
    g[6] = xixj + c0x[2];
    g[7] = xixj + c0x[3];
    g[28] = g[24] * (xixj + c0x[0]) + c0x[0] * b00[0] + b10[0] * cpx[0];
    g[29] = g[25] * (xixj + c0x[1]) + c0x[1] * b00[1] + b10[1] * cpx[1];
    g[30] = g[26] * (xixj + c0x[2]) + c0x[2] * b00[2] + b10[2] * cpx[2];
    g[31] = g[27] * (xixj + c0x[3]) + c0x[3] * b00[3] + b10[3] * cpx[3];
    g[12] = cpx[0] * (xixj + c0x[0]) + b00[0];
    g[13] = cpx[1] * (xixj + c0x[1]) + b00[1];
    g[14] = cpx[2] * (xixj + c0x[2]) + b00[2];
    g[15] = cpx[3] * (xixj + c0x[3]) + b00[3];
    g[48] = 1;
    g[49] = 1;
    g[50] = 1;
    g[51] = 1;
    g[64] = c0y[0];
    g[65] = c0y[1];
    g[66] = c0y[2];
    g[67] = c0y[3];
    g[56] = cpy[0];
    g[57] = cpy[1];
    g[58] = cpy[2];
    g[59] = cpy[3];
    g[72] = cpy[0] * c0y[0] + b00[0];
    g[73] = cpy[1] * c0y[1] + b00[1];
    g[74] = cpy[2] * c0y[2] + b00[2];
    g[75] = cpy[3] * c0y[3] + b00[3];
    g[68] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[69] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[70] = c0y[2] * (yiyj + c0y[2]) + b10[2];
    g[71] = c0y[3] * (yiyj + c0y[3]) + b10[3];
    g[52] = yiyj + c0y[0];
    g[53] = yiyj + c0y[1];
    g[54] = yiyj + c0y[2];
    g[55] = yiyj + c0y[3];
    g[76] = g[72] * (yiyj + c0y[0]) + c0y[0] * b00[0] + b10[0] * cpy[0];
    g[77] = g[73] * (yiyj + c0y[1]) + c0y[1] * b00[1] + b10[1] * cpy[1];
    g[78] = g[74] * (yiyj + c0y[2]) + c0y[2] * b00[2] + b10[2] * cpy[2];
    g[79] = g[75] * (yiyj + c0y[3]) + c0y[3] * b00[3] + b10[3] * cpy[3];
    g[60] = cpy[0] * (yiyj + c0y[0]) + b00[0];
    g[61] = cpy[1] * (yiyj + c0y[1]) + b00[1];
    g[62] = cpy[2] * (yiyj + c0y[2]) + b00[2];
    g[63] = cpy[3] * (yiyj + c0y[3]) + b00[3];
    g[112] = c0z[0] * g[96];
    g[113] = c0z[1] * g[97];
    g[114] = c0z[2] * g[98];
    g[115] = c0z[3] * g[99];
    g[104] = cpz[0] * g[96];
    g[105] = cpz[1] * g[97];
    g[106] = cpz[2] * g[98];
    g[107] = cpz[3] * g[99];
    g[120] = cpz[0] * g[112] + b00[0] * g[96];
    g[121] = cpz[1] * g[113] + b00[1] * g[97];
    g[122] = cpz[2] * g[114] + b00[2] * g[98];
    g[123] = cpz[3] * g[115] + b00[3] * g[99];
    g[116] = g[112] * (zizj + c0z[0]) + b10[0] * g[96];
    g[117] = g[113] * (zizj + c0z[1]) + b10[1] * g[97];
    g[118] = g[114] * (zizj + c0z[2]) + b10[2] * g[98];
    g[119] = g[115] * (zizj + c0z[3]) + b10[3] * g[99];
    g[100] = g[96] * (zizj + c0z[0]);
    g[101] = g[97] * (zizj + c0z[1]);
    g[102] = g[98] * (zizj + c0z[2]);
    g[103] = g[99] * (zizj + c0z[3]);
    g[124] = g[120] * (zizj + c0z[0]) + b10[0] * g[104] + b00[0] * g[112];
    g[125] = g[121] * (zizj + c0z[1]) + b10[1] * g[105] + b00[1] * g[113];
    g[126] = g[122] * (zizj + c0z[2]) + b10[2] * g[106] + b00[2] * g[114];
    g[127] = g[123] * (zizj + c0z[3]) + b10[3] * g[107] + b00[3] * g[115];
    g[108] = zizj * g[104] + cpz[0] * g[112] + b00[0] * g[96];
    g[109] = zizj * g[105] + cpz[1] * g[113] + b00[1] * g[97];
    g[110] = zizj * g[106] + cpz[2] * g[114] + b00[2] * g[98];
    g[111] = zizj * g[107] + cpz[3] * g[115] + b00[3] * g[99];
}

inline void _srg0_2d4d_1200(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[10] = c0x[2];
    g[11] = c0x[3];
    g[16] = c0x[0] * c0x[0] + b10[0];
    g[17] = c0x[1] * c0x[1] + b10[1];
    g[18] = c0x[2] * c0x[2] + b10[2];
    g[19] = c0x[3] * c0x[3] + b10[3];
    g[20] = g[16] * (xixj + c0x[0]) + c0x[0] * 2 * b10[0];
    g[21] = g[17] * (xixj + c0x[1]) + c0x[1] * 2 * b10[1];
    g[22] = g[18] * (xixj + c0x[2]) + c0x[2] * 2 * b10[2];
    g[23] = g[19] * (xixj + c0x[3]) + c0x[3] * 2 * b10[3];
    g[12] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[13] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[14] = c0x[2] * (xixj + c0x[2]) + b10[2];
    g[15] = c0x[3] * (xixj + c0x[3]) + b10[3];
    g[4] = xixj + c0x[0];
    g[5] = xixj + c0x[1];
    g[6] = xixj + c0x[2];
    g[7] = xixj + c0x[3];
    g[32] = 1;
    g[33] = 1;
    g[34] = 1;
    g[35] = 1;
    g[40] = c0y[0];
    g[41] = c0y[1];
    g[42] = c0y[2];
    g[43] = c0y[3];
    g[48] = c0y[0] * c0y[0] + b10[0];
    g[49] = c0y[1] * c0y[1] + b10[1];
    g[50] = c0y[2] * c0y[2] + b10[2];
    g[51] = c0y[3] * c0y[3] + b10[3];
    g[52] = g[48] * (yiyj + c0y[0]) + c0y[0] * 2 * b10[0];
    g[53] = g[49] * (yiyj + c0y[1]) + c0y[1] * 2 * b10[1];
    g[54] = g[50] * (yiyj + c0y[2]) + c0y[2] * 2 * b10[2];
    g[55] = g[51] * (yiyj + c0y[3]) + c0y[3] * 2 * b10[3];
    g[44] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[45] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[46] = c0y[2] * (yiyj + c0y[2]) + b10[2];
    g[47] = c0y[3] * (yiyj + c0y[3]) + b10[3];
    g[36] = yiyj + c0y[0];
    g[37] = yiyj + c0y[1];
    g[38] = yiyj + c0y[2];
    g[39] = yiyj + c0y[3];
    g[72] = c0z[0] * g[64];
    g[73] = c0z[1] * g[65];
    g[74] = c0z[2] * g[66];
    g[75] = c0z[3] * g[67];
    g[80] = c0z[0] * g[72] + b10[0] * g[64];
    g[81] = c0z[1] * g[73] + b10[1] * g[65];
    g[82] = c0z[2] * g[74] + b10[2] * g[66];
    g[83] = c0z[3] * g[75] + b10[3] * g[67];
    g[84] = g[80] * (zizj + c0z[0]) + 2 * b10[0] * g[72];
    g[85] = g[81] * (zizj + c0z[1]) + 2 * b10[1] * g[73];
    g[86] = g[82] * (zizj + c0z[2]) + 2 * b10[2] * g[74];
    g[87] = g[83] * (zizj + c0z[3]) + 2 * b10[3] * g[75];
    g[76] = g[72] * (zizj + c0z[0]) + b10[0] * g[64];
    g[77] = g[73] * (zizj + c0z[1]) + b10[1] * g[65];
    g[78] = g[74] * (zizj + c0z[2]) + b10[2] * g[66];
    g[79] = g[75] * (zizj + c0z[3]) + b10[3] * g[67];
    g[68] = g[64] * (zizj + c0z[0]);
    g[69] = g[65] * (zizj + c0z[1]);
    g[70] = g[66] * (zizj + c0z[2]);
    g[71] = g[67] * (zizj + c0z[3]);
}

inline void _srg0_2d4d_2000(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[12] = 1;
    g[13] = 1;
    g[14] = 1;
    g[15] = 1;
    g[16] = c0y[0];
    g[17] = c0y[1];
    g[18] = c0y[2];
    g[19] = c0y[3];
    g[20] = c0y[0] * c0y[0] + b10[0];
    g[21] = c0y[1] * c0y[1] + b10[1];
    g[22] = c0y[2] * c0y[2] + b10[2];
    g[23] = c0y[3] * c0y[3] + b10[3];
    g[28] = c0z[0] * g[24];
    g[29] = c0z[1] * g[25];
    g[30] = c0z[2] * g[26];
    g[31] = c0z[3] * g[27];
    g[32] = c0z[0] * g[28] + b10[0] * g[24];
    g[33] = c0z[1] * g[29] + b10[1] * g[25];
    g[34] = c0z[2] * g[30] + b10[2] * g[26];
    g[35] = c0z[3] * g[31] + b10[3] * g[27];
}

inline void _srg0_2d4d_2001(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[12] = cpx[0];
    g[13] = cpx[1];
    g[14] = cpx[2];
    g[15] = cpx[3];
    g[16] = cpx[0] * c0x[0] + b00[0];
    g[17] = cpx[1] * c0x[1] + b00[1];
    g[18] = cpx[2] * c0x[2] + b00[2];
    g[19] = cpx[3] * c0x[3] + b00[3];
    g[20] = c0x[0] * (g[16] + b00[0]) + b10[0] * cpx[0];
    g[21] = c0x[1] * (g[17] + b00[1]) + b10[1] * cpx[1];
    g[22] = c0x[2] * (g[18] + b00[2]) + b10[2] * cpx[2];
    g[23] = c0x[3] * (g[19] + b00[3]) + b10[3] * cpx[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = c0y[0];
    g[29] = c0y[1];
    g[30] = c0y[2];
    g[31] = c0y[3];
    g[32] = c0y[0] * c0y[0] + b10[0];
    g[33] = c0y[1] * c0y[1] + b10[1];
    g[34] = c0y[2] * c0y[2] + b10[2];
    g[35] = c0y[3] * c0y[3] + b10[3];
    g[36] = cpy[0];
    g[37] = cpy[1];
    g[38] = cpy[2];
    g[39] = cpy[3];
    g[40] = cpy[0] * c0y[0] + b00[0];
    g[41] = cpy[1] * c0y[1] + b00[1];
    g[42] = cpy[2] * c0y[2] + b00[2];
    g[43] = cpy[3] * c0y[3] + b00[3];
    g[44] = c0y[0] * (g[40] + b00[0]) + b10[0] * cpy[0];
    g[45] = c0y[1] * (g[41] + b00[1]) + b10[1] * cpy[1];
    g[46] = c0y[2] * (g[42] + b00[2]) + b10[2] * cpy[2];
    g[47] = c0y[3] * (g[43] + b00[3]) + b10[3] * cpy[3];
    g[52] = c0z[0] * g[48];
    g[53] = c0z[1] * g[49];
    g[54] = c0z[2] * g[50];
    g[55] = c0z[3] * g[51];
    g[56] = c0z[0] * g[52] + b10[0] * g[48];
    g[57] = c0z[1] * g[53] + b10[1] * g[49];
    g[58] = c0z[2] * g[54] + b10[2] * g[50];
    g[59] = c0z[3] * g[55] + b10[3] * g[51];
    g[60] = cpz[0] * g[48];
    g[61] = cpz[1] * g[49];
    g[62] = cpz[2] * g[50];
    g[63] = cpz[3] * g[51];
    g[64] = cpz[0] * g[52] + b00[0] * g[48];
    g[65] = cpz[1] * g[53] + b00[1] * g[49];
    g[66] = cpz[2] * g[54] + b00[2] * g[50];
    g[67] = cpz[3] * g[55] + b00[3] * g[51];
    g[68] = c0z[0] * g[64] + b10[0] * g[60] + b00[0] * g[52];
    g[69] = c0z[1] * g[65] + b10[1] * g[61] + b00[1] * g[53];
    g[70] = c0z[2] * g[66] + b10[2] * g[62] + b00[2] * g[54];
    g[71] = c0z[3] * g[67] + b10[3] * g[63] + b00[3] * g[55];
}

inline void _srg0_2d4d_2010(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *cpx = bc->c0px;
    double *cpy = bc->c0py;
    double *cpz = bc->c0pz;
    double *b00 = bc->b00;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[12] = cpx[0];
    g[13] = cpx[1];
    g[14] = cpx[2];
    g[15] = cpx[3];
    g[16] = cpx[0] * c0x[0] + b00[0];
    g[17] = cpx[1] * c0x[1] + b00[1];
    g[18] = cpx[2] * c0x[2] + b00[2];
    g[19] = cpx[3] * c0x[3] + b00[3];
    g[20] = c0x[0] * (g[16] + b00[0]) + b10[0] * cpx[0];
    g[21] = c0x[1] * (g[17] + b00[1]) + b10[1] * cpx[1];
    g[22] = c0x[2] * (g[18] + b00[2]) + b10[2] * cpx[2];
    g[23] = c0x[3] * (g[19] + b00[3]) + b10[3] * cpx[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = c0y[0];
    g[29] = c0y[1];
    g[30] = c0y[2];
    g[31] = c0y[3];
    g[32] = c0y[0] * c0y[0] + b10[0];
    g[33] = c0y[1] * c0y[1] + b10[1];
    g[34] = c0y[2] * c0y[2] + b10[2];
    g[35] = c0y[3] * c0y[3] + b10[3];
    g[36] = cpy[0];
    g[37] = cpy[1];
    g[38] = cpy[2];
    g[39] = cpy[3];
    g[40] = cpy[0] * c0y[0] + b00[0];
    g[41] = cpy[1] * c0y[1] + b00[1];
    g[42] = cpy[2] * c0y[2] + b00[2];
    g[43] = cpy[3] * c0y[3] + b00[3];
    g[44] = c0y[0] * (g[40] + b00[0]) + b10[0] * cpy[0];
    g[45] = c0y[1] * (g[41] + b00[1]) + b10[1] * cpy[1];
    g[46] = c0y[2] * (g[42] + b00[2]) + b10[2] * cpy[2];
    g[47] = c0y[3] * (g[43] + b00[3]) + b10[3] * cpy[3];
    g[52] = c0z[0] * g[48];
    g[53] = c0z[1] * g[49];
    g[54] = c0z[2] * g[50];
    g[55] = c0z[3] * g[51];
    g[56] = c0z[0] * g[52] + b10[0] * g[48];
    g[57] = c0z[1] * g[53] + b10[1] * g[49];
    g[58] = c0z[2] * g[54] + b10[2] * g[50];
    g[59] = c0z[3] * g[55] + b10[3] * g[51];
    g[60] = cpz[0] * g[48];
    g[61] = cpz[1] * g[49];
    g[62] = cpz[2] * g[50];
    g[63] = cpz[3] * g[51];
    g[64] = cpz[0] * g[52] + b00[0] * g[48];
    g[65] = cpz[1] * g[53] + b00[1] * g[49];
    g[66] = cpz[2] * g[54] + b00[2] * g[50];
    g[67] = cpz[3] * g[55] + b00[3] * g[51];
    g[68] = c0z[0] * g[64] + b10[0] * g[60] + b00[0] * g[52];
    g[69] = c0z[1] * g[65] + b10[1] * g[61] + b00[1] * g[53];
    g[70] = c0z[2] * g[66] + b10[2] * g[62] + b00[2] * g[54];
    g[71] = c0z[3] * g[67] + b10[3] * g[63] + b00[3] * g[55];
}

inline void _srg0_2d4d_2100(double *g, Rys2eT *bc, Env *envs)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[24] = g[8] * (xixj + c0x[0]) + c0x[0] * 2 * b10[0];
    g[25] = g[9] * (xixj + c0x[1]) + c0x[1] * 2 * b10[1];
    g[26] = g[10] * (xixj + c0x[2]) + c0x[2] * 2 * b10[2];
    g[27] = g[11] * (xixj + c0x[3]) + c0x[3] * 2 * b10[3];
    g[20] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[21] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[22] = c0x[2] * (xixj + c0x[2]) + b10[2];
    g[23] = c0x[3] * (xixj + c0x[3]) + b10[3];
    g[16] = xixj + c0x[0];
    g[17] = xixj + c0x[1];
    g[18] = xixj + c0x[2];
    g[19] = xixj + c0x[3];
    g[32] = 1;
    g[33] = 1;
    g[34] = 1;
    g[35] = 1;
    g[36] = c0y[0];
    g[37] = c0y[1];
    g[38] = c0y[2];
    g[39] = c0y[3];
    g[40] = c0y[0] * c0y[0] + b10[0];
    g[41] = c0y[1] * c0y[1] + b10[1];
    g[42] = c0y[2] * c0y[2] + b10[2];
    g[43] = c0y[3] * c0y[3] + b10[3];
    g[56] = g[40] * (yiyj + c0y[0]) + c0y[0] * 2 * b10[0];
    g[57] = g[41] * (yiyj + c0y[1]) + c0y[1] * 2 * b10[1];
    g[58] = g[42] * (yiyj + c0y[2]) + c0y[2] * 2 * b10[2];
    g[59] = g[43] * (yiyj + c0y[3]) + c0y[3] * 2 * b10[3];
    g[52] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[53] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[54] = c0y[2] * (yiyj + c0y[2]) + b10[2];
    g[55] = c0y[3] * (yiyj + c0y[3]) + b10[3];
    g[48] = yiyj + c0y[0];
    g[49] = yiyj + c0y[1];
    g[50] = yiyj + c0y[2];
    g[51] = yiyj + c0y[3];
    g[68] = c0z[0] * g[64];
    g[69] = c0z[1] * g[65];
    g[70] = c0z[2] * g[66];
    g[71] = c0z[3] * g[67];
    g[72] = c0z[0] * g[68] + b10[0] * g[64];
    g[73] = c0z[1] * g[69] + b10[1] * g[65];
    g[74] = c0z[2] * g[70] + b10[2] * g[66];
    g[75] = c0z[3] * g[71] + b10[3] * g[67];
    g[88] = g[72] * (zizj + c0z[0]) + 2 * b10[0] * g[68];
    g[89] = g[73] * (zizj + c0z[1]) + 2 * b10[1] * g[69];
    g[90] = g[74] * (zizj + c0z[2]) + 2 * b10[2] * g[70];
    g[91] = g[75] * (zizj + c0z[3]) + 2 * b10[3] * g[71];
    g[84] = g[68] * (zizj + c0z[0]) + b10[0] * g[64];
    g[85] = g[69] * (zizj + c0z[1]) + b10[1] * g[65];
    g[86] = g[70] * (zizj + c0z[2]) + b10[2] * g[66];
    g[87] = g[71] * (zizj + c0z[3]) + b10[3] * g[67];
    g[80] = g[64] * (zizj + c0z[0]);
    g[81] = g[65] * (zizj + c0z[1]);
    g[82] = g[66] * (zizj + c0z[2]);
    g[83] = g[67] * (zizj + c0z[3]);
}

inline void _srg0_2d4d_3000(double *g, Rys2eT *bc)
{
    double *c0x = bc->c00x;
    double *c0y = bc->c00y;
    double *c0z = bc->c00z;
    double *b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[12] = c0x[0] * (g[8] + 2 * b10[0]);
    g[13] = c0x[1] * (g[9] + 2 * b10[1]);
    g[14] = c0x[2] * (g[10] + 2 * b10[2]);
    g[15] = c0x[3] * (g[11] + 2 * b10[3]);
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = c0y[0];
    g[21] = c0y[1];
    g[22] = c0y[2];
    g[23] = c0y[3];
    g[24] = c0y[0] * c0y[0] + b10[0];
    g[25] = c0y[1] * c0y[1] + b10[1];
    g[26] = c0y[2] * c0y[2] + b10[2];
    g[27] = c0y[3] * c0y[3] + b10[3];
    g[28] = c0y[0] * (g[24] + 2 * b10[0]);
    g[29] = c0y[1] * (g[25] + 2 * b10[1]);
    g[30] = c0y[2] * (g[26] + 2 * b10[2]);
    g[31] = c0y[3] * (g[27] + 2 * b10[3]);
    g[36] = c0z[0] * g[32];
    g[37] = c0z[1] * g[33];
    g[38] = c0z[2] * g[34];
    g[39] = c0z[3] * g[35];
    g[40] = c0z[0] * g[36] + b10[0] * g[32];
    g[41] = c0z[1] * g[37] + b10[1] * g[33];
    g[42] = c0z[2] * g[38] + b10[2] * g[34];
    g[43] = c0z[3] * g[39] + b10[3] * g[35];
    g[44] = c0z[0] * g[40] + 2 * b10[0] * g[36];
    g[45] = c0z[1] * g[41] + 2 * b10[1] * g[37];
    g[46] = c0z[2] * g[42] + 2 * b10[2] * g[38];
    g[47] = c0z[3] * g[43] + 2 * b10[3] * g[39];
}

void srg0_2e_2d4d_unrolled(double *g, Rys2eT *bc, Env *envs)
{
    int type_ijkl = ((envs->li_ceil << 6) | (envs->lj_ceil << 4) |
                     (envs->lk_ceil << 2) | (envs->ll_ceil));
    switch (type_ijkl)
    {
    case 0b00000000:
        _srg0_2d4d_0000(g);
        return;
    case 0b00000001:
        _srg0_2d4d_0001(g, bc);
        return;
    case 0b00000010:
        _srg0_2d4d_0002(g, bc);
        return;
    case 0b00000011:
        _srg0_2d4d_0003(g, bc);
        return;
    case 0b00000100:
        _srg0_2d4d_0010(g, bc);
        return;
    case 0b00000101:
        _srg0_2d4d_0011(g, bc, envs);
        return;
    case 0b00000110:
        _srg0_2d4d_0012(g, bc, envs);
        return;
    case 0b00001000:
        _srg0_2d4d_0020(g, bc);
        return;
    case 0b00001001:
        _srg0_2d4d_0021(g, bc, envs);
        return;
    case 0b00001100:
        _srg0_2d4d_0030(g, bc);
        return;
    case 0b00010000:
        _srg0_2d4d_0100(g, bc);
        return;
    case 0b00010001:
        _srg0_2d4d_0101(g, bc);
        return;
    case 0b00010010:
        _srg0_2d4d_0102(g, bc);
        return;
    case 0b00010100:
        _srg0_2d4d_0110(g, bc);
        return;
    case 0b00010101:
        _srg0_2d4d_0111(g, bc, envs);
        return;
    case 0b00011000:
        _srg0_2d4d_0120(g, bc);
        return;
    case 0b00100000:
        _srg0_2d4d_0200(g, bc);
        return;
    case 0b00100001:
        _srg0_2d4d_0201(g, bc);
        return;
    case 0b00100100:
        _srg0_2d4d_0210(g, bc);
        return;
    case 0b00110000:
        _srg0_2d4d_0300(g, bc);
        return;
    case 0b01000000:
        _srg0_2d4d_1000(g, bc);
        return;
    case 0b01000001:
        _srg0_2d4d_1001(g, bc);
        return;
    case 0b01000010:
        _srg0_2d4d_1002(g, bc);
        return;
    case 0b01000100:
        _srg0_2d4d_1010(g, bc);
        return;
    case 0b01000101:
        _srg0_2d4d_1011(g, bc, envs);
        return;
    case 0b01001000:
        _srg0_2d4d_1020(g, bc);
        return;
    case 0b01010000:
        _srg0_2d4d_1100(g, bc, envs);
        return;
    case 0b01010001:
        _srg0_2d4d_1101(g, bc, envs);
        return;
    case 0b01010100:
        _srg0_2d4d_1110(g, bc, envs);
        return;
    case 0b01100000:
        _srg0_2d4d_1200(g, bc, envs);
        return;
    case 0b10000000:
        _srg0_2d4d_2000(g, bc);
        return;
    case 0b10000001:
        _srg0_2d4d_2001(g, bc);
        return;
    case 0b10000100:
        _srg0_2d4d_2010(g, bc);
        return;
    case 0b10010000:
        _srg0_2d4d_2100(g, bc, envs);
        return;
    case 0b11000000:
        _srg0_2d4d_3000(g, bc);
        return;
    }
    fprintf(stderr, "Dimension error for CINTg0_2e_lj2d4d: iklj = %d %d %d %d",
            (int)envs->li_ceil, (int)envs->lk_ceil,
            (int)envs->ll_ceil, (int)envs->lj_ceil);
}

void g0_2e_2d(double *g, Rys2eT *bc, Env *envs)
{
    const int nroots = envs->nrys_roots;
    const int nmax = envs->li_ceil + envs->lj_ceil;
    const int mmax = envs->lk_ceil + envs->ll_ceil;
    const int dm = envs->g2d_klmax;
    const int dn = envs->g2d_ijmax;
    int i, j, m, n, off;
    double *gx = g;
    double *gy = gx + envs->g_size;
    double *gz = gy + 2 * envs->g_size;

    for (i = 0; i < nroots; i++)
    {
        gx[i] = 1;
        gy[i] = 1;
        // gz[i] = w[i];
    }

    double s0x, s1x, s2x;
    double s0y, s1y, s2y;
    double s0z, s1z, s2z;
    double c00x, c00y, c00z, c0px, c0py, c0pz, b10, b01, b00;
    for (i = 0; i < nroots; i++)
    {
        c00x = bc->c00x[i];
        c00y = bc->c00y[i];
        c00z = bc->c00z[i];
        c0px = bc->c0px[i];
        c0py = bc->c0py[i];
        c0pz = bc->c0pz[i];
        b10 = bc->b10[i];
        b01 = bc->b01[i];
        b00 = bc->b00[i];
        if (nmax > 0)
        {
            // gx(irys,0,1) = c00(irys) * gx(irys,0,0)
            // gx(irys,0,n+1) = c00(irys)*gx(irys,0,n)
            // + n*b10(irys)*gx(irys,0,n-1)
            s0x = gx[i];
            s0y = gy[i];
            s0z = gz[i];
            s1x = c00x * s0x;
            s1y = c00y * s0y;
            s1z = c00z * s0z;
            gx[i + dn] = s1x;
            gy[i + dn] = s1y;
            gz[i + dn] = s1z;
            for (n = 1; n < nmax; ++n)
            {
                s2x = c00x * s1x + n * b10 * s0x;
                s2y = c00y * s1y + n * b10 * s0y;
                s2z = c00z * s1z + n * b10 * s0z;
                gx[i + (n + 1) * dn] = s2x;
                gy[i + (n + 1) * dn] = s2y;
                gz[i + (n + 1) * dn] = s2z;
                s0x = s1x;
                s0y = s1y;
                s0z = s1z;
                s1x = s2x;
                s1y = s2y;
                s1z = s2z;
            }
        }

        if (mmax > 0)
        {
            // gx(irys,1,0) = c0p(irys) * gx(irys,0,0)
            // gx(irys,m+1,0) = c0p(irys)*gx(irys,m,0)
            // + m*b01(irys)*gx(irys,m-1,0)
            s0x = gx[i];
            s0y = gy[i];
            s0z = gz[i];
            s1x = c0px * s0x;
            s1y = c0py * s0y;
            s1z = c0pz * s0z;
            gx[i + dm] = s1x;
            gy[i + dm] = s1y;
            gz[i + dm] = s1z;
            for (m = 1; m < mmax; ++m)
            {
                s2x = c0px * s1x + m * b01 * s0x;
                s2y = c0py * s1y + m * b01 * s0y;
                s2z = c0pz * s1z + m * b01 * s0z;
                gx[i + (m + 1) * dm] = s2x;
                gy[i + (m + 1) * dm] = s2y;
                gz[i + (m + 1) * dm] = s2z;
                s0x = s1x;
                s0y = s1y;
                s0z = s1z;
                s1x = s2x;
                s1y = s2y;
                s1z = s2z;
            }

            if (nmax > 0)
            {
                // gx(irys,1,1) = c0p(irys)*gx(irys,0,1) + b00(irys)*gx(irys,0,0)
                // gx(irys,m+1,1) = c0p(irys)*gx(irys,m,1)
                // + m*b01(irys)*gx(irys,m-1,1)
                // + b00(irys)*gx(irys,m,0)
                s0x = gx[i + dn];
                s0y = gy[i + dn];
                s0z = gz[i + dn];
                s1x = c0px * s0x + b00 * gx[i];
                s1y = c0py * s0y + b00 * gy[i];
                s1z = c0pz * s0z + b00 * gz[i];
                gx[i + dn + dm] = s1x;
                gy[i + dn + dm] = s1y;
                gz[i + dn + dm] = s1z;
                for (m = 1; m < mmax; ++m)
                {
                    s2x = c0px * s1x + m * b01 * s0x + b00 * gx[i + m * dm];
                    s2y = c0py * s1y + m * b01 * s0y + b00 * gy[i + m * dm];
                    s2z = c0pz * s1z + m * b01 * s0z + b00 * gz[i + m * dm];
                    gx[i + dn + (m + 1) * dm] = s2x;
                    gy[i + dn + (m + 1) * dm] = s2y;
                    gz[i + dn + (m + 1) * dm] = s2z;
                    s0x = s1x;
                    s0y = s1y;
                    s0z = s1z;
                    s1x = s2x;
                    s1y = s2y;
                    s1z = s2z;
                }
            }
        }

        // gx(irys,m,n+1) = c00(irys)*gx(irys,m,n)
        // + n*b10(irys)*gx(irys,m,n-1)
        // + m*b00(irys)*gx(irys,m-1,n)
        for (m = 1; m <= mmax; ++m)
        {
            off = m * dm;
            j = off + i;
            s0x = gx[j];
            s0y = gy[j];
            s0z = gz[j];
            s1x = gx[j + dn];
            s1y = gy[j + dn];
            s1z = gz[j + dn];
            for (n = 1; n < nmax; ++n)
            {
                s2x = c00x * s1x + n * b10 * s0x + m * b00 * gx[j + n * dn - dm];
                s2y = c00y * s1y + n * b10 * s0y + m * b00 * gy[j + n * dn - dm];
                s2z = c00z * s1z + n * b10 * s0z + m * b00 * gz[j + n * dn - dm];
                gx[j + (n + 1) * dn] = s2x;
                gy[j + (n + 1) * dn] = s2y;
                gz[j + (n + 1) * dn] = s2z;
                s0x = s1x;
                s0y = s1y;
                s0z = s1z;
                s1x = s2x;
                s1y = s2y;
                s1z = s2z;
            }
        }
    }
}

int rys_root1(double X, double *roots, double *weights)
{
    double Y, F1;

    if (X > 33.)
    {
        weights[0] = sqrt(constants::PI_4 / X);
        roots[0] = 0.5E+00 / (X - 0.5E+00);
        return 0;
    }
    else if (X < 3.e-7)
    {
        weights[0] = 1.0E+00 - X / 3.0E+00;
        roots[0] = 0.5E+00 - X / 5.0E+00;
        return 0;
    }

    double E = exp(-X);
    if (X > 15.)
    {
        Y = 1. / X;
        F1 = (((1.9623264149430E-01 * Y - 4.9695241464490E-01) * Y -
               6.0156581186481E-05) *
                  E +
              sqrt(constants::PI_4 / X) - E) *
             Y;
        F1 *= .5;
    }
    else if (X > 10.)
    {
        Y = 1. / X;
        F1 = ((((-1.8784686463512E-01 * Y + 2.2991849164985E-01) * Y -
                4.9893752514047E-01) *
                   Y -
               2.1916512131607E-05) *
                  E +
              sqrt(constants::PI_4 / X) - E) *
             Y;
        F1 *= .5;
    }
    else if (X > 5.)
    {
        Y = 1. / X;
        F1 = (((((((4.6897511375022E-01 * Y - 6.9955602298985E-01) * Y +
                   5.3689283271887E-01) *
                      Y -
                  3.2883030418398E-01) *
                     Y +
                 2.4645596956002E-01) *
                    Y -
                4.9984072848436E-01) *
                   Y -
               3.1501078774085E-06) *
                  E +
              sqrt(constants::PI_4 / X) - E) *
             Y;
        F1 *= .5;
    }
    else if (X > 3.)
    {
        Y = X - 4.0E+00;
        F1 = ((((((((((-2.62453564772299E-11 * Y + 3.24031041623823E-10) * Y -
                      3.614965656163E-09) *
                         Y +
                     3.760256799971E-08) *
                        Y -
                    3.553558319675E-07) *
                       Y +
                   3.022556449731E-06) *
                      Y -
                  2.290098979647E-05) *
                     Y +
                 1.526537461148E-04) *
                    Y -
                8.81947375894379E-04) *
                   Y +
               4.33207949514611E-03) *
                  Y -
              1.75257821619926E-02) *
                 Y +
             5.28406320615584E-02;
    }
    else if (X > 1.)
    {
        Y = X - 2.0E+00;
        F1 = ((((((((((-1.61702782425558E-10 * Y + 1.96215250865776E-09) * Y -
                      2.14234468198419E-08) *
                         Y +
                     2.17216556336318E-07) *
                        Y -
                    1.98850171329371E-06) *
                       Y +
                   1.62429321438911E-05) *
                      Y -
                  1.16740298039895E-04) *
                     Y +
                 7.24888732052332E-04) *
                    Y -
                3.79490003707156E-03) *
                   Y +
               1.61723488664661E-02) *
                  Y -
              5.29428148329736E-02) *
                 Y +
             1.15702180856167E-01;
    }
    else
    {
        F1 = ((((((((-8.36313918003957E-08 * X + 1.21222603512827E-06) * X -
                    1.15662609053481E-05) *
                       X +
                   9.25197374512647E-05) *
                      X -
                  6.40994113129432E-04) *
                     X +
                 3.78787044215009E-03) *
                    X -
                1.85185172458485E-02) *
                   X +
               7.14285713298222E-02) *
                  X -
              1.99999999997023E-01) *
                 X +
             3.33333333333318E-01;
    }

    double WW1 = 2. * X * F1 + E;
    weights[0] = WW1;
    roots[0] = F1 / (WW1 - F1);
    return 0;
}

int rys_root2(double X, double *roots, double *weights)
{

    double R12, R22, W22;
    double RT1, RT2, WW1, WW2;
    double F1, E, Y;

    R12 = 2.75255128608411E-01;
    R22 = 2.72474487139158E+00;
    W22 = 9.17517095361369E-02;

    if (X < 3.e-7)
    {
        RT1 = 1.30693606237085E-01 - 2.90430236082028E-02 * X;
        RT2 = 2.86930639376291E+00 - 6.37623643058102E-01 * X;
        WW1 = 6.52145154862545E-01 - 1.22713621927067E-01 * X;
        WW2 = 3.47854845137453E-01 - 2.10619711404725E-01 * X;
    }
    else if (X < 1.)
    {
        F1 = ((((((((-8.36313918003957E-08 * X + 1.21222603512827E-06) * X -
                    1.15662609053481E-05) *
                       X +
                   9.25197374512647E-05) *
                      X -
                  6.40994113129432E-04) *
                     X +
                 3.78787044215009E-03) *
                    X -
                1.85185172458485E-02) *
                   X +
               7.14285713298222E-02) *
                  X -
              1.99999999997023E-01) *
                 X +
             3.33333333333318E-01;
        WW1 = (X + X) * F1 + exp(-X);
        RT1 = (((((((-2.35234358048491E-09 * X + 2.49173650389842E-08) * X -
                    4.558315364581E-08) *
                       X -
                   2.447252174587E-06) *
                      X +
                  4.743292959463E-05) *
                     X -
                 5.33184749432408E-04) *
                    X +
                4.44654947116579E-03) *
                   X -
               2.90430236084697E-02) *
                  X +
              1.30693606237085E-01;
        RT2 = (((((((-2.47404902329170E-08 * X + 2.36809910635906E-07) * X +
                    1.835367736310E-06) *
                       X -
                   2.066168802076E-05) *
                      X -
                  1.345693393936E-04) *
                     X -
                 5.88154362858038E-05) *
                    X +
                5.32735082098139E-02) *
                   X -
               6.37623643056745E-01) *
                  X +
              2.86930639376289E+00;
        WW2 = ((F1 - WW1) * RT1 + F1) * (1.0E+00 + RT2) / (RT2 - RT1);
        WW1 = WW1 - WW2;
    }
    else if (X < 3.)
    {
        Y = X - 2.0E+00;
        F1 = ((((((((((-1.61702782425558E-10 * Y + 1.96215250865776E-09) * Y -
                      2.14234468198419E-08) *
                         Y +
                     2.17216556336318E-07) *
                        Y -
                    1.98850171329371E-06) *
                       Y +
                   1.62429321438911E-05) *
                      Y -
                  1.16740298039895E-04) *
                     Y +
                 7.24888732052332E-04) *
                    Y -
                3.79490003707156E-03) *
                   Y +
               1.61723488664661E-02) *
                  Y -
              5.29428148329736E-02) *
                 Y +
             1.15702180856167E-01;
        WW1 = (X + X) * F1 + exp(-X);
        RT1 = (((((((((-6.36859636616415E-12 * Y + 8.47417064776270E-11) * Y -
                      5.152207846962E-10) *
                         Y -
                     3.846389873308E-10) *
                        Y +
                    8.472253388380E-08) *
                       Y -
                   1.85306035634293E-06) *
                      Y +
                  2.47191693238413E-05) *
                     Y -
                 2.49018321709815E-04) *
                    Y +
                2.19173220020161E-03) *
                   Y -
               1.63329339286794E-02) *
                  Y +
              8.68085688285261E-02;
        RT2 = (((((((((1.45331350488343E-10 * Y + 2.07111465297976E-09) * Y -
                      1.878920917404E-08) *
                         Y -
                     1.725838516261E-07) *
                        Y +
                    2.247389642339E-06) *
                       Y +
                   9.76783813082564E-06) *
                      Y -
                  1.93160765581969E-04) *
                     Y -
                 1.58064140671893E-03) *
                    Y +
                4.85928174507904E-02) *
                   Y -
               4.30761584997596E-01) *
                  Y +
              1.80400974537950E+00;
        WW2 = ((F1 - WW1) * RT1 + F1) * (1.0E+00 + RT2) / (RT2 - RT1);
        WW1 = WW1 - WW2;
    }
    else if (X < 5.)
    {
        Y = X - 4.0E+00;
        F1 = ((((((((((-2.62453564772299E-11 * Y + 3.24031041623823E-10) * Y -
                      3.614965656163E-09) *
                         Y +
                     3.760256799971E-08) *
                        Y -
                    3.553558319675E-07) *
                       Y +
                   3.022556449731E-06) *
                      Y -
                  2.290098979647E-05) *
                     Y +
                 1.526537461148E-04) *
                    Y -
                8.81947375894379E-04) *
                   Y +
               4.33207949514611E-03) *
                  Y -
              1.75257821619926E-02) *
                 Y +
             5.28406320615584E-02;
        WW1 = (X + X) * F1 + exp(-X);
        RT1 = ((((((((-4.11560117487296E-12 * Y + 7.10910223886747E-11) * Y -
                     1.73508862390291E-09) *
                        Y +
                    5.93066856324744E-08) *
                       Y -
                   9.76085576741771E-07) *
                      Y +
                  1.08484384385679E-05) *
                     Y -
                 1.12608004981982E-04) *
                    Y +
                1.16210907653515E-03) *
                   Y -
               9.89572595720351E-03) *
                  Y +
              6.12589701086408E-02;
        RT2 = (((((((((-1.80555625241001E-10 * Y + 5.44072475994123E-10) * Y +
                      1.603498045240E-08) *
                         Y -
                     1.497986283037E-07) *
                        Y -
                    7.017002532106E-07) *
                       Y +
                   1.85882653064034E-05) *
                      Y -
                  2.04685420150802E-05) *
                     Y -
                 2.49327728643089E-03) *
                    Y +
                3.56550690684281E-02) *
                   Y -
               2.60417417692375E-01) *
                  Y +
              1.12155283108289E+00;
        WW2 = ((F1 - WW1) * RT1 + F1) * (1.0E+00 + RT2) / (RT2 - RT1);
        WW1 = WW1 - WW2;
    }
    else if (X < 10)
    {
        E = exp(-X);
        WW1 = ((((((4.6897511375022E-01 / X - 6.9955602298985E-01) / X +
                   5.3689283271887E-01) /
                      X -
                  3.2883030418398E-01) /
                     X +
                 2.4645596956002E-01) /
                    X -
                4.9984072848436E-01) /
                   X -
               3.1501078774085E-06) *
                  E +
              sqrt(constants::PI_4 / X);
        F1 = (WW1 - E) / (X + X);
        Y = X - 7.5E+00;
        RT1 = (((((((((((((-1.43632730148572E-16 * Y + 2.38198922570405E-16) *
                              Y +
                          1.358319618800E-14) *
                             Y -
                         7.064522786879E-14) *
                            Y -
                        7.719300212748E-13) *
                           Y +
                       7.802544789997E-12) *
                          Y +
                      6.628721099436E-11) *
                         Y -
                     1.775564159743E-09) *
                        Y +
                    1.713828823990E-08) *
                       Y -
                   1.497500187053E-07) *
                      Y +
                  2.283485114279E-06) *
                     Y -
                 3.76953869614706E-05) *
                    Y +
                4.74791204651451E-04) *
                   Y -
               4.60448960876139E-03) *
                  Y +
              3.72458587837249E-02;
        RT2 = ((((((((((((2.48791622798900E-14 * Y - 1.36113510175724E-13) * Y -
                         2.224334349799E-12) *
                            Y +
                        4.190559455515E-11) *
                           Y -
                       2.222722579924E-10) *
                          Y -
                      2.624183464275E-09) *
                         Y +
                     6.128153450169E-08) *
                        Y -
                    4.383376014528E-07) *
                       Y -
                   2.49952200232910E-06) *
                      Y +
                  1.03236647888320E-04) *
                     Y -
                 1.44614664924989E-03) *
                    Y +
                1.35094294917224E-02) *
                   Y -
               9.53478510453887E-02) *
                  Y +
              5.44765245686790E-01;
        WW2 = ((F1 - WW1) * RT1 + F1) * (1.0E+00 + RT2) / (RT2 - RT1);
        WW1 = WW1 - WW2;
    }
    else if (X < 15)
    {
        E = exp(-X);
        WW1 = (((-1.8784686463512E-01 / X + 2.2991849164985E-01) / X -
                4.9893752514047E-01) /
                   X -
               2.1916512131607E-05) *
                  E +
              sqrt(constants::PI_4 / X);
        F1 = (WW1 - E) / (X + X);
        RT1 = ((((-1.01041157064226E-05 * X + 1.19483054115173E-03) * X -
                 6.73760231824074E-02) *
                    X +
                1.25705571069895E+00) *
                   X +
               (((-8.57609422987199E+03 / X + 5.91005939591842E+03) / X -
                 1.70807677109425E+03) /
                    X +
                2.64536689959503E+02) /
                   X -
               2.38570496490846E+01) *
                  E +
              R12 / (X - R12);
        RT2 = (((3.39024225137123E-04 * X - 9.34976436343509E-02) * X -
                4.22216483306320E+00) *
                   X +
               (((-2.08457050986847E+03 / X -
                  1.04999071905664E+03) /
                     X +
                 3.39891508992661E+02) /
                    X -
                1.56184800325063E+02) /
                   X +
               8.00839033297501E+00) *
                  E +
              R22 / (X - R22);
        WW2 = ((F1 - WW1) * RT1 + F1) * (1.0E+00 + RT2) / (RT2 - RT1);
        WW1 = WW1 - WW2;
    }
    else if (X < 33)
    {
        E = exp(-X);
        WW1 = ((1.9623264149430E-01 / X - 4.9695241464490E-01) / X -
               6.0156581186481E-05) *
                  E +
              sqrt(constants::PI_4 / X);
        F1 = (WW1 - E) / (X + X);
        RT1 = ((((-1.14906395546354E-06 * X + 1.76003409708332E-04) * X -
                 1.71984023644904E-02) *
                    X -
                1.37292644149838E-01) *
                   X +
               (-4.75742064274859E+01 / X + 9.21005186542857E+00) / X -
               2.31080873898939E-02) *
                  E +
              R12 / (X - R12);
        RT2 = (((3.64921633404158E-04 * X - 9.71850973831558E-02) * X -
                4.02886174850252E+00) *
                   X +
               (-1.35831002139173E+02 / X -
                8.66891724287962E+01) /
                   X +
               2.98011277766958E+00) *
                  E +
              R22 / (X - R22);
        WW2 = ((F1 - WW1) * RT1 + F1) * (1.0E+00 + RT2) / (RT2 - RT1);
        WW1 = WW1 - WW2;
    }
    else if (X < 40)
    {
        WW1 = sqrt(constants::PI_4 / X);
        E = exp(-X);
        RT1 = (-8.78947307498880E-01 * X + 1.09243702330261E+01) * E + R12 / (X - R12);
        RT2 = (-9.28903924275977E+00 * X + 8.10642367843811E+01) * E + R22 / (X - R22);
        WW2 = (4.46857389308400E+00 * X - 7.79250653461045E+01) * E + W22 * WW1;
        WW1 = WW1 - WW2;
    }
    else
    {
        WW1 = sqrt(constants::PI_4 / X);
        RT1 = R12 / (X - R12);
        RT2 = R22 / (X - R22);
        WW2 = W22 * WW1;
        WW1 = WW1 - WW2;
    }
    roots[0] = RT1;
    roots[1] = RT2;
    weights[0] = WW1;
    weights[1] = WW2;
    return 0;
}

int rys_root3(double X, double *roots, double *weights)
{

    double R13, R23, W23, R33, W33;
    double RT1, RT2, RT3, WW1, WW2, WW3;
    double F1, F2, E, T1, T2, T3, A1, A2, Y;

    R13 = 1.90163509193487E-01;
    R23 = 1.78449274854325E+00;
    W23 = 1.77231492083829E-01;
    R33 = 5.52534374226326E+00;
    W33 = 5.11156880411248E-03;

    if (X < 3.e-7)
    {
        RT1 = 6.03769246832797E-02 - 9.28875764357368E-03 * X;
        RT2 = 7.76823355931043E-01 - 1.19511285527878E-01 * X;
        RT3 = 6.66279971938567E+00 - 1.02504611068957E+00 * X;
        WW1 = 4.67913934572691E-01 - 5.64876917232519E-02 * X;
        WW2 = 3.60761573048137E-01 - 1.49077186455208E-01 * X;
        WW3 = 1.71324492379169E-01 - 1.27768455150979E-01 * X;
    }
    else if (X < 1.)
    {
        RT1 = ((((((-5.10186691538870E-10 * X + 2.40134415703450E-08) * X -
                   5.01081057744427E-07) *
                      X +
                  7.58291285499256E-06) *
                     X -
                 9.55085533670919E-05) *
                    X +
                1.02893039315878E-03) *
                   X -
               9.28875764374337E-03) *
                  X +
              6.03769246832810E-02;
        RT2 = ((((((-1.29646524960555E-08 * X + 7.74602292865683E-08) * X +
                   1.56022811158727E-06) *
                      X -
                  1.58051990661661E-05) *
                     X -
                 3.30447806384059E-04) *
                    X +
                9.74266885190267E-03) *
                   X -
               1.19511285526388E-01) *
                  X +
              7.76823355931033E-01;
        RT3 = ((((((-9.28536484109606E-09 * X - 3.02786290067014E-07) * X -
                   2.50734477064200E-06) *
                      X -
                  7.32728109752881E-06) *
                     X +
                 2.44217481700129E-04) *
                    X +
                4.94758452357327E-02) *
                   X -
               1.02504611065774E+00) *
                  X +
              6.66279971938553E+00;
        F2 = ((((((((-7.60911486098850E-08 * X + 1.09552870123182E-06) * X -
                    1.03463270693454E-05) *
                       X +
                   8.16324851790106E-05) *
                      X -
                  5.55526624875562E-04) *
                     X +
                 3.20512054753924E-03) *
                    X -
                1.51515139838540E-02) *
                   X +
               5.55555554649585E-02) *
                  X -
              1.42857142854412E-01) *
                 X +
             1.99999999999986E-01;
        E = exp(-X);
        F1 = ((X + X) * F2 + E) / 3.0E+00;
        WW1 = (X + X) * F1 + E;
        T1 = RT1 / (RT1 + 1.0E+00);
        T2 = RT2 / (RT2 + 1.0E+00);
        T3 = RT3 / (RT3 + 1.0E+00);
        A2 = F2 - T1 * F1;
        A1 = F1 - T1 * WW1;
        WW3 = (A2 - T2 * A1) / ((T3 - T2) * (T3 - T1));
        WW2 = (T3 * A1 - A2) / ((T3 - T2) * (T2 - T1));
        WW1 = WW1 - WW2 - WW3;
    }
    else if (X < 3.)
    {
        Y = X - 2.0E+00;
        RT1 = ((((((((1.44687969563318E-12 * Y + 4.85300143926755E-12) * Y -
                     6.55098264095516E-10) *
                        Y +
                    1.56592951656828E-08) *
                       Y -
                   2.60122498274734E-07) *
                      Y +
                  3.86118485517386E-06) *
                     Y -
                 5.13430986707889E-05) *
                    Y +
                6.03194524398109E-04) *
                   Y -
               6.11219349825090E-03) *
                  Y +
              4.52578254679079E-02;
        RT2 = (((((((6.95964248788138E-10 * Y - 5.35281831445517E-09) * Y -
                    6.745205954533E-08) *
                       Y +
                   1.502366784525E-06) *
                      Y +
                  9.923326947376E-07) *
                     Y -
                 3.89147469249594E-04) *
                    Y +
                7.51549330892401E-03) *
                   Y -
               8.48778120363400E-02) *
                  Y +
              5.73928229597613E-01;
        RT3 = ((((((((-2.81496588401439E-10 * Y + 3.61058041895031E-09) * Y +
                     4.53631789436255E-08) *
                        Y -
                    1.40971837780847E-07) *
                       Y -
                   6.05865557561067E-06) *
                      Y -
                  5.15964042227127E-05) *
                     Y +
                 3.34761560498171E-05) *
                    Y +
                5.04871005319119E-02) *
                   Y -
               8.24708946991557E-01) *
                  Y +
              4.81234667357205E+00;
        F2 = ((((((((((-1.48044231072140E-10 * Y + 1.78157031325097E-09) * Y -
                      1.92514145088973E-08) *
                         Y +
                     1.92804632038796E-07) *
                        Y -
                    1.73806555021045E-06) *
                       Y +
                   1.39195169625425E-05) *
                      Y -
                  9.74574633246452E-05) *
                     Y +
                 5.83701488646511E-04) *
                    Y -
                2.89955494844975E-03) *
                   Y +
               1.13847001113810E-02) *
                  Y -
              3.23446977320647E-02) *
                 Y +
             5.29428148329709E-02;
        E = exp(-X);
        F1 = ((X + X) * F2 + E) / 3.0E+00;
        WW1 = (X + X) * F1 + E;
        T1 = RT1 / (RT1 + 1.0E+00);
        T2 = RT2 / (RT2 + 1.0E+00);
        T3 = RT3 / (RT3 + 1.0E+00);
        A2 = F2 - T1 * F1;
        A1 = F1 - T1 * WW1;
        WW3 = (A2 - T2 * A1) / ((T3 - T2) * (T3 - T1));
        WW2 = (T3 * A1 - A2) / ((T3 - T2) * (T2 - T1));
        WW1 = WW1 - WW2 - WW3;
    }
    else if (X < 5.)
    {
        Y = X - 4.0E+00;
        RT1 = (((((((1.44265709189601E-11 * Y - 4.66622033006074E-10) * Y +
                    7.649155832025E-09) *
                       Y -
                   1.229940017368E-07) *
                      Y +
                  2.026002142457E-06) *
                     Y -
                 2.87048671521677E-05) *
                    Y +
                3.70326938096287E-04) *
                   Y -
               4.21006346373634E-03) *
                  Y +
              3.50898470729044E-02;
        RT2 = ((((((((-2.65526039155651E-11 * Y + 1.97549041402552E-10) * Y +
                     2.15971131403034E-09) *
                        Y -
                    7.95045680685193E-08) *
                       Y +
                   5.15021914287057E-07) *
                      Y +
                  1.11788717230514E-05) *
                     Y -
                 3.33739312603632E-04) *
                    Y +
                5.30601428208358E-03) *
                   Y -
               5.93483267268959E-02) *
                  Y +
              4.31180523260239E-01;
        RT3 = ((((((((-3.92833750584041E-10 * Y - 4.16423229782280E-09) * Y +
                     4.42413039572867E-08) *
                        Y +
                    6.40574545989551E-07) *
                       Y -
                   3.05512456576552E-06) *
                      Y -
                  1.05296443527943E-04) *
                     Y -
                 6.14120969315617E-04) *
                    Y +
                4.89665802767005E-02) *
                   Y -
               6.24498381002855E-01) *
                  Y +
              3.36412312243724E+00;
        F2 = ((((((((((-2.36788772599074E-11 * Y + 2.89147476459092E-10) * Y -
                      3.18111322308846E-09) *
                         Y +
                     3.25336816562485E-08) *
                        Y -
                    3.00873821471489E-07) *
                       Y +
                   2.48749160874431E-06) *
                      Y -
                  1.81353179793672E-05) *
                     Y +
                 1.14504948737066E-04) *
                    Y -
                6.10614987696677E-04) *
                   Y +
               2.64584212770942E-03) *
                  Y -
              8.66415899015349E-03) *
                 Y +
             1.75257821619922E-02;
        E = exp(-X);
        F1 = ((X + X) * F2 + E) / 3.0E+00;
        WW1 = (X + X) * F1 + E;
        T1 = RT1 / (RT1 + 1.0E+00);
        T2 = RT2 / (RT2 + 1.0E+00);
        T3 = RT3 / (RT3 + 1.0E+00);
        A2 = F2 - T1 * F1;
        A1 = F1 - T1 * WW1;
        WW3 = (A2 - T2 * A1) / ((T3 - T2) * (T3 - T1));
        WW2 = (T3 * A1 - A2) / ((T3 - T2) * (T2 - T1));
        WW1 = WW1 - WW2 - WW3;
    }
    else if (X < 10)
    {
        E = exp(-X);
        WW1 = ((((((4.6897511375022E-01 / X - 6.9955602298985E-01) / X +
                   5.3689283271887E-01) /
                      X -
                  3.2883030418398E-01) /
                     X +
                 2.4645596956002E-01) /
                    X -
                4.9984072848436E-01) /
                   X -
               3.1501078774085E-06) *
                  E +
              sqrt(constants::PI_4 / X);
        F1 = (WW1 - E) / (X + X);
        F2 = (F1 + F1 + F1 - E) / (X + X);
        Y = X - 7.5E+00;
        RT1 = (((((((((((5.74429401360115E-16 * Y + 7.11884203790984E-16) * Y -
                        6.736701449826E-14) *
                           Y -
                       6.264613873998E-13) *
                          Y +
                      1.315418927040E-11) *
                         Y -
                     4.23879635610964E-11) *
                        Y +
                    1.39032379769474E-09) *
                       Y -
                   4.65449552856856E-08) *
                      Y +
                  7.34609900170759E-07) *
                     Y -
                 1.08656008854077E-05) *
                    Y +
                1.77930381549953E-04) *
                   Y -
               2.39864911618015E-03) *
                  Y +
              2.39112249488821E-02;
        RT2 = (((((((((((1.13464096209120E-14 * Y + 6.99375313934242E-15) * Y -
                        8.595618132088E-13) *
                           Y -
                       5.293620408757E-12) *
                          Y -
                      2.492175211635E-11) *
                         Y +
                     2.73681574882729E-09) *
                        Y -
                    1.06656985608482E-08) *
                       Y -
                   4.40252529648056E-07) *
                      Y +
                  9.68100917793911E-06) *
                     Y -
                 1.68211091755327E-04) *
                    Y +
                2.69443611274173E-03) *
                   Y -
               3.23845035189063E-02) *
                  Y +
              2.75969447451882E-01;
        RT3 = ((((((((((((6.66339416996191E-15 * Y + 1.84955640200794E-13) * Y -
                         1.985141104444E-12) *
                            Y -
                        2.309293727603E-11) *
                           Y +
                       3.917984522103E-10) *
                          Y +
                      1.663165279876E-09) *
                         Y -
                     6.205591993923E-08) *
                        Y +
                    8.769581622041E-09) *
                       Y +
                   8.97224398620038E-06) *
                      Y -
                  3.14232666170796E-05) *
                     Y -
                 1.83917335649633E-03) *
                    Y +
                3.51246831672571E-02) *
                   Y -
               3.22335051270860E-01) *
                  Y +
              1.73582831755430E+00;
        T1 = RT1 / (RT1 + 1.0E+00);
        T2 = RT2 / (RT2 + 1.0E+00);
        T3 = RT3 / (RT3 + 1.0E+00);
        A2 = F2 - T1 * F1;
        A1 = F1 - T1 * WW1;
        WW3 = (A2 - T2 * A1) / ((T3 - T2) * (T3 - T1));
        WW2 = (T3 * A1 - A2) / ((T3 - T2) * (T2 - T1));
        WW1 = WW1 - WW2 - WW3;
    }
    else if (X < 15)
    {
        E = exp(-X);
        WW1 = (((-1.8784686463512E-01 / X + 2.2991849164985E-01) / X -
                4.9893752514047E-01) /
                   X -
               2.1916512131607E-05) *
                  E +
              sqrt(constants::PI_4 / X);
        F1 = (WW1 - E) / (X + X);
        F2 = (F1 + F1 + F1 - E) / (X + X);
        Y = X - 12.5E+00;
        RT1 = (((((((((((4.42133001283090E-16 * Y - 2.77189767070441E-15) * Y -
                        4.084026087887E-14) *
                           Y +
                       5.379885121517E-13) *
                          Y +
                      1.882093066702E-12) *
                         Y -
                     8.67286219861085E-11) *
                        Y +
                    7.11372337079797E-10) *
                       Y -
                   3.55578027040563E-09) *
                      Y +
                  1.29454702851936E-07) *
                     Y -
                 4.14222202791434E-06) *
                    Y +
                8.04427643593792E-05) *
                   Y -
               1.18587782909876E-03) *
                  Y +
              1.53435577063174E-02;
        RT2 = (((((((((((6.85146742119357E-15 * Y - 1.08257654410279E-14) * Y -
                        8.579165965128E-13) *
                           Y +
                       6.642452485783E-12) *
                          Y +
                      4.798806828724E-11) *
                         Y -
                     1.13413908163831E-09) *
                        Y +
                    7.08558457182751E-09) *
                       Y -
                   5.59678576054633E-08) *
                      Y +
                  2.51020389884249E-06) *
                     Y -
                 6.63678914608681E-05) *
                    Y +
                1.11888323089714E-03) *
                   Y -
               1.45361636398178E-02) *
                  Y +
              1.65077877454402E-01;
        RT3 = ((((((((((((3.20622388697743E-15 * Y - 2.73458804864628E-14) * Y -
                         3.157134329361E-13) *
                            Y +
                        8.654129268056E-12) *
                           Y -
                       5.625235879301E-11) *
                          Y -
                      7.718080513708E-10) *
                         Y +
                     2.064664199164E-08) *
                        Y -
                    1.567725007761E-07) *
                       Y -
                   1.57938204115055E-06) *
                      Y +
                  6.27436306915967E-05) *
                     Y -
                 1.01308723606946E-03) *
                    Y +
                1.13901881430697E-02) *
                   Y -
               1.01449652899450E-01) *
                  Y +
              7.77203937334739E-01;
        T1 = RT1 / (RT1 + 1.0E+00);
        T2 = RT2 / (RT2 + 1.0E+00);
        T3 = RT3 / (RT3 + 1.0E+00);
        A2 = F2 - T1 * F1;
        A1 = F1 - T1 * WW1;
        WW3 = (A2 - T2 * A1) / ((T3 - T2) * (T3 - T1));
        WW2 = (T3 * A1 - A2) / ((T3 - T2) * (T2 - T1));
        WW1 = WW1 - WW2 - WW3;
    }
    else if (X < 33)
    {
        E = exp(-X);
        WW1 = ((1.9623264149430E-01 / X - 4.9695241464490E-01) / X -
               6.0156581186481E-05) *
                  E +
              sqrt(constants::PI_4 / X);
        F1 = (WW1 - E) / (X + X);
        F2 = (F1 + F1 + F1 - E) / (X + X);
        if (X < 20)
        {
            RT1 = ((((((-2.43270989903742E-06 * X + 3.57901398988359E-04) * X -
                       2.34112415981143E-02) *
                          X +
                      7.81425144913975E-01) *
                         X -
                     1.73209218219175E+01) *
                        X +
                    2.43517435690398E+02) *
                       X +
                   (-1.97611541576986E+04 / X + 9.82441363463929E+03) / X -
                   2.07970687843258E+03) *
                      E +
                  R13 / (X - R13);
            RT2 = (((((-2.62627010965435E-04 * X + 3.49187925428138E-02) * X -
                      3.09337618731880E+00) *
                         X +
                     1.07037141010778E+02) *
                        X -
                    2.36659637247087E+03) *
                       X +
                   ((-2.91669113681020E+06 / X +
                     1.41129505262758E+06) /
                        X -
                    2.91532335433779E+05) /
                       X +
                   3.35202872835409E+04) *
                      E +
                  R23 / (X - R23);
            RT3 = (((((9.31856404738601E-05 * X - 2.87029400759565E-02) * X -
                      7.83503697918455E-01) *
                         X -
                     1.84338896480695E+01) *
                        X +
                    4.04996712650414E+02) *
                       X +
                   (-1.89829509315154E+05 / X +
                    5.11498390849158E+04) /
                       X -
                   6.88145821789955E+03) *
                      E +
                  R33 / (X - R33);
        }
        else
        {
            RT1 = ((((-4.97561537069643E-04 * X - 5.00929599665316E-02) * X +
                     1.31099142238996E+00) *
                        X -
                    1.88336409225481E+01) *
                       X -
                   6.60344754467191E+02 / X + 1.64931462413877E+02) *
                      E +
                  R13 / (X - R13);
            RT2 = ((((-4.48218898474906E-03 * X - 5.17373211334924E-01) * X +
                     1.13691058739678E+01) *
                        X -
                    1.65426392885291E+02) *
                       X -
                   6.30909125686731E+03 / X + 1.52231757709236E+03) *
                      E +
                  R23 / (X - R23);
            RT3 = ((((-1.38368602394293E-02 * X - 1.77293428863008E+00) * X +
                     1.73639054044562E+01) *
                        X -
                    3.57615122086961E+02) *
                       X -
                   1.45734701095912E+04 / X + 2.69831813951849E+03) *
                      E +
                  R33 / (X - R33);
        }
        T1 = RT1 / (RT1 + 1.0E+00);
        T2 = RT2 / (RT2 + 1.0E+00);
        T3 = RT3 / (RT3 + 1.0E+00);
        A2 = F2 - T1 * F1;
        A1 = F1 - T1 * WW1;
        WW3 = (A2 - T2 * A1) / ((T3 - T2) * (T3 - T1));
        WW2 = (T3 * A1 - A2) / ((T3 - T2) * (T2 - T1));
        WW1 = WW1 - WW2 - WW3;
    }
    else if (X < 47)
    {
        WW1 = sqrt(constants::PI_4 / X);
        E = exp(-X);
        RT1 = ((-7.39058467995275E+00 * X + 3.21318352526305E+02) * X -
               3.99433696473658E+03) *
                  E +
              R13 / (X - R13);
        RT2 = ((-7.38726243906513E+01 * X + 3.13569966333873E+03) * X -
               3.86862867311321E+04) *
                  E +
              R23 / (X - R23);
        RT3 = ((-2.63750565461336E+02 * X + 1.04412168692352E+04) * X -
               1.28094577915394E+05) *
                  E +
              R33 / (X - R33);
        WW3 = (((1.52258947224714E-01 * X - 8.30661900042651E+00) * X +
                1.92977367967984E+02) *
                   X -
               1.67787926005344E+03) *
                  E +
              W33 * WW1;
        WW2 = ((6.15072615497811E+01 * X - 2.91980647450269E+03) * X +
               3.80794303087338E+04) *
                  E +
              W23 * WW1;
        WW1 = WW1 - WW2 - WW3;
    }
    else
    {
        WW1 = sqrt(constants::PI_4 / X);
        RT1 = R13 / (X - R13);
        RT2 = R23 / (X - R23);
        RT3 = R33 / (X - R33);
        WW2 = W23 * WW1;
        WW3 = W33 * WW1;
        WW1 = WW1 - WW2 - WW3;
    }
    roots[0] = RT1;
    roots[1] = RT2;
    roots[2] = RT3;
    weights[0] = WW1;
    weights[1] = WW2;
    weights[2] = WW3;
    return 0;
}

int rys_root4(double X, double *roots, double *weights)
{
    double R14, R24, W24, R34, W34, R44, W44;
    double RT1, RT2, RT3, RT4, WW1, WW2, WW3, WW4;
    double Y, E;

    R14 = 1.45303521503316E-01;
    R24 = 1.33909728812636E+00;
    W24 = 2.34479815323517E-01;
    R34 = 3.92696350135829E+00;
    W34 = 1.92704402415764E-02;
    R44 = 8.58863568901199E+00;
    W44 = 2.25229076750736E-04;

    if (X <= 3.0E-7)
    {
        RT1 = 3.48198973061471E-02 - 4.09645850660395E-03 * X;
        RT2 = 3.81567185080042E-01 - 4.48902570656719E-02 * X;
        RT3 = 1.73730726945891E+00 - 2.04389090547327E-01 * X;
        RT4 = 1.18463056481549E+01 - 1.39368301742312E+00 * X;
        WW1 = 3.62683783378362E-01 - 3.13844305713928E-02 * X;
        WW2 = 3.13706645877886E-01 - 8.98046242557724E-02 * X;
        WW3 = 2.22381034453372E-01 - 1.29314370958973E-01 * X;
        WW4 = 1.01228536290376E-01 - 8.28299075414321E-02 * X;
    }
    else if (X <= 1.0)
    {
        RT1 = ((((((-1.95309614628539E-10 * X + 5.19765728707592E-09) * X -
                   1.01756452250573E-07) *
                      X +
                  1.72365935872131E-06) *
                     X -
                 2.61203523522184E-05) *
                    X +
                3.52921308769880E-04) *
                   X -
               4.09645850658433E-03) *
                  X +
              3.48198973061469E-02;
        RT2 = (((((-1.89554881382342E-08 * X + 3.07583114342365E-07) * X +
                  1.270981734393E-06) *
                     X -
                 1.417298563884E-04) *
                    X +
                3.226979163176E-03) *
                   X -
               4.48902570678178E-02) *
                  X +
              3.81567185080039E-01;
        RT3 = ((((((1.77280535300416E-09 * X + 3.36524958870615E-08) * X -
                   2.58341529013893E-07) *
                      X -
                  1.13644895662320E-05) *
                     X -
                 7.91549618884063E-05) *
                    X +
                1.03825827346828E-02) *
                   X -
               2.04389090525137E-01) *
                  X +
              1.73730726945889E+00;
        RT4 = (((((-5.61188882415248E-08 * X - 2.49480733072460E-07) * X +
                  3.428685057114E-06) *
                     X +
                 1.679007454539E-04) *
                    X +
                4.722855585715E-02) *
                   X -
               1.39368301737828E+00) *
                  X +
              1.18463056481543E+01;
        WW1 = ((((((-1.14649303201279E-08 * X + 1.88015570196787E-07) * X -
                   2.33305875372323E-06) *
                      X +
                  2.68880044371597E-05) *
                     X -
                 2.94268428977387E-04) *
                    X +
                3.06548909776613E-03) *
                   X -
               3.13844305680096E-02) *
                  X +
              3.62683783378335E-01;
        WW2 = ((((((((-4.11720483772634E-09 * X + 6.54963481852134E-08) * X -
                     7.20045285129626E-07) *
                        X +
                    6.93779646721723E-06) *
                       X -
                   6.05367572016373E-05) *
                      X +
                  4.74241566251899E-04) *
                     X -
                 3.26956188125316E-03) *
                    X +
                1.91883866626681E-02) *
                   X -
               8.98046242565811E-02) *
                  X +
              3.13706645877886E-01;
        WW3 = ((((((((-3.41688436990215E-08 * X + 5.07238960340773E-07) * X -
                     5.01675628408220E-06) *
                        X +
                    4.20363420922845E-05) *
                       X -
                   3.08040221166823E-04) *
                      X +
                  1.94431864731239E-03) *
                     X -
                 1.02477820460278E-02) *
                    X +
                4.28670143840073E-02) *
                   X -
               1.29314370962569E-01) *
                  X +
              2.22381034453369E-01;
        WW4 = (((((((((4.99660550769508E-09 * X - 7.94585963310120E-08) * X +
                      8.359072409485E-07) *
                         X -
                     7.422369210610E-06) *
                        X +
                    5.763374308160E-05) *
                       X -
                   3.86645606718233E-04) *
                      X +
                  2.18417516259781E-03) *
                     X -
                 9.99791027771119E-03) *
                    X +
                3.48791097377370E-02) *
                   X -
               8.28299075413889E-02) *
                  X +
              1.01228536290376E-01;
    }
    else if (X <= 5)
    {
        Y = X - 3.0E+00;
        RT1 = (((((((((-1.48570633747284E-15 * Y - 1.33273068108777E-13) * Y +
                      4.068543696670E-12) *
                         Y -
                     9.163164161821E-11) *
                        Y +
                    2.046819017845E-09) *
                       Y -
                   4.03076426299031E-08) *
                      Y +
                  7.29407420660149E-07) *
                     Y -
                 1.23118059980833E-05) *
                    Y +
                1.88796581246938E-04) *
                   Y -
               2.53262912046853E-03) *
                  Y +
              2.51198234505021E-02;
        RT2 = (((((((((1.35830583483312E-13 * Y - 2.29772605964836E-12) * Y -
                      3.821500128045E-12) *
                         Y +
                     6.844424214735E-10) *
                        Y -
                    1.048063352259E-08) *
                       Y +
                   1.50083186233363E-08) *
                      Y +
                  3.48848942324454E-06) *
                     Y -
                 1.08694174399193E-04) *
                    Y +
                2.08048885251999E-03) *
                   Y -
               2.91205805373793E-02) *
                  Y +
              2.72276489515713E-01;
        RT3 = (((((((((5.02799392850289E-13 * Y + 1.07461812944084E-11) * Y -
                      1.482277886411E-10) *
                         Y -
                     2.153585661215E-09) *
                        Y +
                    3.654087802817E-08) *
                       Y +
                   5.15929575830120E-07) *
                      Y -
                  9.52388379435709E-06) *
                     Y -
                 2.16552440036426E-04) *
                    Y +
                9.03551469568320E-03) *
                   Y -
               1.45505469175613E-01) *
                  Y +
              1.21449092319186E+00;
        RT4 = (((((((((-1.08510370291979E-12 * Y + 6.41492397277798E-11) * Y +
                      7.542387436125E-10) *
                         Y -
                     2.213111836647E-09) *
                        Y -
                    1.448228963549E-07) *
                       Y -
                   1.95670833237101E-06) *
                      Y -
                  1.07481314670844E-05) *
                     Y +
                 1.49335941252765E-04) *
                    Y +
                4.87791531990593E-02) *
                   Y -
               1.10559909038653E+00) *
                  Y +
              8.09502028611780E+00;
        WW1 = ((((((((((-4.65801912689961E-14 * Y + 7.58669507106800E-13) * Y -
                       1.186387548048E-11) *
                          Y +
                      1.862334710665E-10) *
                         Y -
                     2.799399389539E-09) *
                        Y +
                    4.148972684255E-08) *
                       Y -
                   5.933568079600E-07) *
                      Y +
                  8.168349266115E-06) *
                     Y -
                 1.08989176177409E-04) *
                    Y +
                1.41357961729531E-03) *
                   Y -
               1.87588361833659E-02) *
                  Y +
              2.89898651436026E-01;
        WW2 = ((((((((((((-1.46345073267549E-14 * Y + 2.25644205432182E-13) * Y -
                         3.116258693847E-12) *
                            Y +
                        4.321908756610E-11) *
                           Y -
                       5.673270062669E-10) *
                          Y +
                      7.006295962960E-09) *
                         Y -
                     8.120186517000E-08) *
                        Y +
                    8.775294645770E-07) *
                       Y -
                   8.77829235749024E-06) *
                      Y +
                  8.04372147732379E-05) *
                     Y -
                 6.64149238804153E-04) *
                    Y +
                4.81181506827225E-03) *
                   Y -
               2.88982669486183E-02) *
                  Y +
              1.56247249979288E-01;
        WW3 = (((((((((((((9.06812118895365E-15 * Y - 1.40541322766087E-13) *
                              Y +
                          1.919270015269E-12) *
                             Y -
                         2.605135739010E-11) *
                            Y +
                        3.299685839012E-10) *
                           Y -
                       3.86354139348735E-09) *
                          Y +
                      4.16265847927498E-08) *
                         Y -
                     4.09462835471470E-07) *
                        Y +
                    3.64018881086111E-06) *
                       Y -
                   2.88665153269386E-05) *
                      Y +
                  2.00515819789028E-04) *
                     Y -
                 1.18791896897934E-03) *
                    Y +
                5.75223633388589E-03) *
                   Y -
               2.09400418772687E-02) *
                  Y +
              4.85368861938873E-02;
        WW4 = ((((((((((((((-9.74835552342257E-16 * Y + 1.57857099317175E-14) *
                               Y -
                           2.249993780112E-13) *
                              Y +
                          3.173422008953E-12) *
                             Y -
                         4.161159459680E-11) *
                            Y +
                        5.021343560166E-10) *
                           Y -
                       5.545047534808E-09) *
                          Y +
                      5.554146993491E-08) *
                         Y -
                     4.99048696190133E-07) *
                        Y +
                    3.96650392371311E-06) *
                       Y -
                   2.73816413291214E-05) *
                      Y +
                  1.60106988333186E-04) *
                     Y -
                 7.64560567879592E-04) *
                    Y +
                2.81330044426892E-03) *
                   Y -
               7.16227030134947E-03) *
                  Y +
              9.66077262223353E-03;
    }
    else if (X <= 10.0)
    {
        Y = X - 7.5E+00;
        RT1 = (((((((((4.64217329776215E-15 * Y - 6.27892383644164E-15) * Y +
                      3.462236347446E-13) *
                         Y -
                     2.927229355350E-11) *
                        Y +
                    5.090355371676E-10) *
                       Y -
                   9.97272656345253E-09) *
                      Y +
                  2.37835295639281E-07) *
                     Y -
                 4.60301761310921E-06) *
                    Y +
                8.42824204233222E-05) *
                   Y -
               1.37983082233081E-03) *
                  Y +
              1.66630865869375E-02;
        RT2 = (((((((((2.93981127919047E-14 * Y + 8.47635639065744E-13) * Y -
                      1.446314544774E-11) *
                         Y -
                     6.149155555753E-12) *
                        Y +
                    8.484275604612E-10) *
                       Y -
                   6.10898827887652E-08) *
                      Y +
                  2.39156093611106E-06) *
                     Y -
                 5.35837089462592E-05) *
                    Y +
                1.00967602595557E-03) *
                   Y -
               1.57769317127372E-02) *
                  Y +
              1.74853819464285E-01;
        RT3 = ((((((((((2.93523563363000E-14 * Y - 6.40041776667020E-14) * Y -
                       2.695740446312E-12) *
                          Y +
                      1.027082960169E-10) *
                         Y -
                     5.822038656780E-10) *
                        Y -
                    3.159991002539E-08) *
                       Y +
                   4.327249251331E-07) *
                      Y +
                  4.856768455119E-06) *
                     Y -
                 2.54617989427762E-04) *
                    Y +
                5.54843378106589E-03) *
                   Y -
               7.95013029486684E-02) *
                  Y +
              7.20206142703162E-01;
        RT4 = (((((((((((-1.62212382394553E-14 * Y + 7.68943641360593E-13) * Y +
                        5.764015756615E-12) *
                           Y -
                       1.380635298784E-10) *
                          Y -
                      1.476849808675E-09) *
                         Y +
                     1.84347052385605E-08) *
                        Y +
                    3.34382940759405E-07) *
                       Y -
                   1.39428366421645E-06) *
                      Y -
                  7.50249313713996E-05) *
                     Y -
                 6.26495899187507E-04) *
                    Y +
                4.69716410901162E-02) *
                   Y -
               6.66871297428209E-01) *
                  Y +
              4.11207530217806E+00;
        WW1 = ((((((((((-1.65995045235997E-15 * Y + 6.91838935879598E-14) * Y -
                       9.131223418888E-13) *
                          Y +
                      1.403341829454E-11) *
                         Y -
                     3.672235069444E-10) *
                        Y +
                    6.366962546990E-09) *
                       Y -
                   1.039220021671E-07) *
                      Y +
                  1.959098751715E-06) *
                     Y -
                 3.33474893152939E-05) *
                    Y +
                5.72164211151013E-04) *
                   Y -
               1.05583210553392E-02) *
                  Y +
              2.26696066029591E-01;
        WW2 = ((((((((((((-3.57248951192047E-16 * Y + 6.25708409149331E-15) * Y -
                         9.657033089714E-14) *
                            Y +
                        1.507864898748E-12) *
                           Y -
                       2.332522256110E-11) *
                          Y +
                      3.428545616603E-10) *
                         Y -
                     4.698730937661E-09) *
                        Y +
                    6.219977635130E-08) *
                       Y -
                   7.83008889613661E-07) *
                      Y +
                  9.08621687041567E-06) *
                     Y -
                 9.86368311253873E-05) *
                    Y +
                9.69632496710088E-04) *
                   Y -
               8.14594214284187E-03) *
                  Y +
              8.50218447733457E-02;
        WW3 = (((((((((((((1.64742458534277E-16 * Y - 2.68512265928410E-15) *
                              Y +
                          3.788890667676E-14) *
                             Y -
                         5.508918529823E-13) *
                            Y +
                        7.555896810069E-12) *
                           Y -
                       9.69039768312637E-11) *
                          Y +
                      1.16034263529672E-09) *
                         Y -
                     1.28771698573873E-08) *
                        Y +
                    1.31949431805798E-07) *
                       Y -
                   1.23673915616005E-06) *
                      Y +
                  1.04189803544936E-05) *
                     Y -
                 7.79566003744742E-05) *
                    Y +
                5.03162624754434E-04) *
                   Y -
               2.55138844587555E-03) *
                  Y +
              1.13250730954014E-02;
        WW4 = ((((((((((((((-1.55714130075679E-17 * Y + 2.57193722698891E-16) *
                               Y -
                           3.626606654097E-15) *
                              Y +
                          5.234734676175E-14) *
                             Y -
                         7.067105402134E-13) *
                            Y +
                        8.793512664890E-12) *
                           Y -
                       1.006088923498E-10) *
                          Y +
                      1.050565098393E-09) *
                         Y -
                     9.91517881772662E-09) *
                        Y +
                    8.35835975882941E-08) *
                       Y -
                   6.19785782240693E-07) *
                      Y +
                  3.95841149373135E-06) *
                     Y -
                 2.11366761402403E-05) *
                    Y +
                9.00474771229507E-05) *
                   Y -
               2.78777909813289E-04) *
                  Y +
              5.26543779837487E-04;
    }
    else if (X <= 15)
    {
        Y = X - 12.5E+00;
        RT1 = (((((((((((4.94869622744119E-17 * Y + 8.03568805739160E-16) * Y -
                        5.599125915431E-15) *
                           Y -
                       1.378685560217E-13) *
                          Y +
                      7.006511663249E-13) *
                         Y +
                     1.30391406991118E-11) *
                        Y +
                    8.06987313467541E-11) *
                       Y -
                   5.20644072732933E-09) *
                      Y +
                  7.72794187755457E-08) *
                     Y -
                 1.61512612564194E-06) *
                    Y +
                4.15083811185831E-05) *
                   Y -
               7.87855975560199E-04) *
                  Y +
              1.14189319050009E-02;
        RT2 = (((((((((((4.89224285522336E-16 * Y + 1.06390248099712E-14) * Y -
                        5.446260182933E-14) *
                           Y -
                       1.613630106295E-12) *
                          Y +
                      3.910179118937E-12) *
                         Y +
                     1.90712434258806E-10) *
                        Y +
                    8.78470199094761E-10) *
                       Y -
                   5.97332993206797E-08) *
                      Y +
                  9.25750831481589E-07) *
                     Y -
                 2.02362185197088E-05) *
                    Y +
                4.92341968336776E-04) *
                   Y -
               8.68438439874703E-03) *
                  Y +
              1.15825965127958E-01;
        RT3 = ((((((((((6.12419396208408E-14 * Y + 1.12328861406073E-13) * Y -
                       9.051094103059E-12) *
                          Y -
                      4.781797525341E-11) *
                         Y +
                     1.660828868694E-09) *
                        Y +
                    4.499058798868E-10) *
                       Y -
                   2.519549641933E-07) *
                      Y +
                  4.977444040180E-06) *
                     Y -
                 1.25858350034589E-04) *
                    Y +
                2.70279176970044E-03) *
                   Y -
               3.99327850801083E-02) *
                  Y +
              4.33467200855434E-01;
        RT4 = (((((((((((4.63414725924048E-14 * Y - 4.72757262693062E-14) * Y -
                        1.001926833832E-11) *
                           Y +
                       6.074107718414E-11) *
                          Y +
                      1.576976911942E-09) *
                         Y -
                     2.01186401974027E-08) *
                        Y -
                    1.84530195217118E-07) *
                       Y +
                   5.02333087806827E-06) *
                      Y +
                  9.66961790843006E-06) *
                     Y -
                 1.58522208889528E-03) *
                    Y +
                2.80539673938339E-02) *
                   Y -
               2.78953904330072E-01) *
                  Y +
              1.82835655238235E+00;
        WW4 = (((((((((((((2.90401781000996E-18 * Y - 4.63389683098251E-17) *
                              Y +
                          6.274018198326E-16) *
                             Y -
                         8.936002188168E-15) *
                            Y +
                        1.194719074934E-13) *
                           Y -
                       1.45501321259466E-12) *
                          Y +
                      1.64090830181013E-11) *
                         Y -
                     1.71987745310181E-10) *
                        Y +
                    1.63738403295718E-09) *
                       Y -
                   1.39237504892842E-08) *
                      Y +
                  1.06527318142151E-07) *
                     Y -
                 7.27634957230524E-07) *
                    Y +
                4.12159381310339E-06) *
                   Y -
               1.74648169719173E-05) *
                  Y +
              8.50290130067818E-05;
        WW3 = ((((((((((((-4.19569145459480E-17 * Y + 5.94344180261644E-16) * Y -
                         1.148797566469E-14) *
                            Y +
                        1.881303962576E-13) *
                           Y -
                       2.413554618391E-12) *
                          Y +
                      3.372127423047E-11) *
                         Y -
                     4.933988617784E-10) *
                        Y +
                    6.116545396281E-09) *
                       Y -
                   6.69965691739299E-08) *
                      Y +
                  7.52380085447161E-07) *
                     Y -
                 8.08708393262321E-06) *
                    Y +
                6.88603417296672E-05) *
                   Y -
               4.67067112993427E-04) *
                  Y +
              5.42313365864597E-03;
        WW2 = ((((((((((-6.22272689880615E-15 * Y + 1.04126809657554E-13) * Y -
                       6.842418230913E-13) *
                          Y +
                      1.576841731919E-11) *
                         Y -
                     4.203948834175E-10) *
                        Y +
                    6.287255934781E-09) *
                       Y -
                   8.307159819228E-08) *
                      Y +
                  1.356478091922E-06) *
                     Y -
                 2.08065576105639E-05) *
                    Y +
                2.52396730332340E-04) *
                   Y -
               2.94484050194539E-03) *
                  Y +
              6.01396183129168E-02;
        WW1 = (((-1.8784686463512E-01 / X + 2.2991849164985E-01) / X -
                4.9893752514047E-01) /
                   X -
               2.1916512131607E-05) *
                  exp(-X) +
              sqrt(constants::PI_4 / X) - WW4 - WW3 - WW2;
    }
    else if (X <= 20)
    {
        WW1 = sqrt(constants::PI_4 / X);
        Y = X - 17.5E+00;
        RT1 = (((((((((((4.36701759531398E-17 * Y - 1.12860600219889E-16) * Y -
                        6.149849164164E-15) *
                           Y +
                       5.820231579541E-14) *
                          Y +
                      4.396602872143E-13) *
                         Y -
                     1.24330365320172E-11) *
                        Y +
                    6.71083474044549E-11) *
                       Y +
                   2.43865205376067E-10) *
                      Y +
                  1.67559587099969E-08) *
                     Y -
                 9.32738632357572E-07) *
                    Y +
                2.39030487004977E-05) *
                   Y -
               4.68648206591515E-04) *
                  Y +
              8.34977776583956E-03;
        RT2 = (((((((((((4.98913142288158E-16 * Y - 2.60732537093612E-16) * Y -
                        7.775156445127E-14) *
                           Y +
                       5.766105220086E-13) *
                          Y +
                      6.432696729600E-12) *
                         Y -
                     1.39571683725792E-10) *
                        Y +
                    5.95451479522191E-10) *
                       Y +
                   2.42471442836205E-09) *
                      Y +
                  2.47485710143120E-07) *
                     Y -
                 1.14710398652091E-05) *
                    Y +
                2.71252453754519E-04) *
                   Y -
               4.96812745851408E-03) *
                  Y +
              8.26020602026780E-02;
        RT3 = (((((((((((1.91498302509009E-15 * Y + 1.48840394311115E-14) * Y -
                        4.316925145767E-13) *
                           Y +
                       1.186495793471E-12) *
                          Y +
                      4.615806713055E-11) *
                         Y -
                     5.54336148667141E-10) *
                        Y +
                    3.48789978951367E-10) *
                       Y -
                   2.79188977451042E-09) *
                      Y +
                  2.09563208958551E-06) *
                     Y -
                 6.76512715080324E-05) *
                    Y +
                1.32129867629062E-03) *
                   Y -
               2.05062147771513E-02) *
                  Y +
              2.88068671894324E-01;
        RT4 = (((((((((((-5.43697691672942E-15 * Y - 1.12483395714468E-13) * Y +
                        2.826607936174E-12) *
                           Y -
                       1.266734493280E-11) *
                          Y -
                      4.258722866437E-10) *
                         Y +
                     9.45486578503261E-09) *
                        Y -
                    5.86635622821309E-08) *
                       Y -
                   1.28835028104639E-06) *
                      Y +
                  4.41413815691885E-05) *
                     Y -
                 7.61738385590776E-04) *
                    Y +
                9.66090902985550E-03) *
                   Y -
               1.01410568057649E-01) *
                  Y +
              9.54714798156712E-01;
        WW4 = ((((((((((((-7.56882223582704E-19 * Y + 7.53541779268175E-18) * Y -
                         1.157318032236E-16) *
                            Y +
                        2.411195002314E-15) *
                           Y -
                       3.601794386996E-14) *
                          Y +
                      4.082150659615E-13) *
                         Y -
                     4.289542980767E-12) *
                        Y +
                    5.086829642731E-11) *
                       Y -
                   6.35435561050807E-10) *
                      Y +
                  6.82309323251123E-09) *
                     Y -
                 5.63374555753167E-08) *
                    Y +
                3.57005361100431E-07) *
                   Y -
               2.40050045173721E-06) *
                  Y +
              4.94171300536397E-05;
        WW3 = (((((((((((-5.54451040921657E-17 * Y + 2.68748367250999E-16) * Y +
                        1.349020069254E-14) *
                           Y -
                       2.507452792892E-13) *
                          Y +
                      1.944339743818E-12) *
                         Y -
                     1.29816917658823E-11) *
                        Y +
                    3.49977768819641E-10) *
                       Y -
                   8.67270669346398E-09) *
                      Y +
                  1.31381116840118E-07) *
                     Y -
                 1.36790720600822E-06) *
                    Y +
                1.19210697673160E-05) *
                   Y -
               1.42181943986587E-04) *
                  Y +
              4.12615396191829E-03;
        WW2 = (((((((((((-1.86506057729700E-16 * Y + 1.16661114435809E-15) * Y +
                        2.563712856363E-14) *
                           Y -
                       4.498350984631E-13) *
                          Y +
                      1.765194089338E-12) *
                         Y +
                     9.04483676345625E-12) *
                        Y +
                    4.98930345609785E-10) *
                       Y -
                   2.11964170928181E-08) *
                      Y +
                  3.98295476005614E-07) *
                     Y -
                 5.49390160829409E-06) *
                    Y +
                7.74065155353262E-05) *
                   Y -
               1.48201933009105E-03) *
                  Y +
              4.97836392625268E-02;
        WW1 = ((1.9623264149430E-01 / X - 4.9695241464490E-01) / X -
               6.0156581186481E-05) *
                  exp(-X) +
              WW1 - WW2 - WW3 - WW4;
    }
    else if (X <= 35)
    {
        WW1 = sqrt(constants::PI_4 / X);
        E = exp(-X);
        RT1 = ((((((-4.45711399441838E-05 * X + 1.27267770241379E-03) * X -
                   2.36954961381262E-01) *
                      X +
                  1.54330657903756E+01) *
                     X -
                 5.22799159267808E+02) *
                    X +
                1.05951216669313E+04) *
                   X +
               (-2.51177235556236E+06 / X + 8.72975373557709E+05) / X -
               1.29194382386499E+05) *
                  E +
              R14 / (X - R14);
        RT2 = (((((-7.85617372254488E-02 * X + 6.35653573484868E+00) * X -
                  3.38296938763990E+02) *
                     X +
                 1.25120495802096E+04) *
                    X -
                3.16847570511637E+05) *
                   X +
               ((-1.02427466127427E+09 / X +
                 3.70104713293016E+08) /
                    X -
                5.87119005093822E+07) /
                   X +
               5.38614211391604E+06) *
                  E +
              R24 / (X - R24);
        RT3 = (((((-2.37900485051067E-01 * X + 1.84122184400896E+01) * X -
                  1.00200731304146E+03) *
                     X +
                 3.75151841595736E+04) *
                    X -
                9.50626663390130E+05) *
                   X +
               ((-2.88139014651985E+09 / X +
                 1.06625915044526E+09) /
                    X -
                1.72465289687396E+08) /
                   X +
               1.60419390230055E+07) *
                  E +
              R34 / (X - R34);
        RT4 = ((((((-6.00691586407385E-04 * X - 3.64479545338439E-01) * X +
                   1.57496131755179E+01) *
                      X -
                  6.54944248734901E+02) *
                     X +
                 1.70830039597097E+04) *
                    X -
                2.90517939780207E+05) *
                   X +
               (3.49059698304732E+07 / X - 1.64944522586065E+07) / X +
               2.96817940164703E+06) *
                  E +
              R44 / (X - R44);
        if (X <= 25)
            WW4 = (((((((2.33766206773151E-07 * X -
                         3.81542906607063E-05) *
                            X +
                        3.51416601267000E-03) *
                           X -
                       1.66538571864728E-01) *
                          X +
                      4.80006136831847E+00) *
                         X -
                     8.73165934223603E+01) *
                        X +
                    9.77683627474638E+02) *
                       X +
                   1.66000945117640E+04 / X - 6.14479071209961E+03) *
                      E +
                  W44 * WW1;
        else
            WW4 = ((((((5.74245945342286E-06 * X -
                        7.58735928102351E-05) *
                           X +
                       2.35072857922892E-04) *
                          X -
                      3.78812134013125E-03) *
                         X +
                     3.09871652785805E-01) *
                        X -
                    7.11108633061306E+00) *
                       X +
                   5.55297573149528E+01) *
                      E +
                  W44 * WW1;
        WW3 = ((((((2.36392855180768E-04 * X - 9.16785337967013E-03) * X +
                   4.62186525041313E-01) *
                      X -
                  1.96943786006540E+01) *
                     X +
                 4.99169195295559E+02) *
                    X -
                6.21419845845090E+03) *
                   X +
               ((+5.21445053212414E+07 / X - 1.34113464389309E+07) / X +
                1.13673298305631E+06) /
                   X -
               2.81501182042707E+03) *
                  E +
              W34 * WW1;
        WW2 = ((((((7.29841848989391E-04 * X - 3.53899555749875E-02) * X +
                   2.07797425718513E+00) *
                      X -
                  1.00464709786287E+02) *
                     X +
                 3.15206108877819E+03) *
                    X -
                6.27054715090012E+04) *
                   X +
               (+1.54721246264919E+07 / X - 5.26074391316381E+06) / X +
               7.67135400969617E+05) *
                  E +
              W24 * WW1;
        WW1 = ((1.9623264149430E-01 / X - 4.9695241464490E-01) / X -
               6.0156581186481E-05) *
                  E +
              WW1 - WW2 - WW3 - WW4;
    }
    else if (X <= 53)
    {
        WW1 = sqrt(constants::PI_4 / X);
        E = exp(-X) * pow(X, 4);
        RT4 = ((-2.19135070169653E-03 * X - 1.19108256987623E-01) * X -
               7.50238795695573E-01) *
                  E +
              R44 / (X - R44);
        RT3 = ((-9.65842534508637E-04 * X - 4.49822013469279E-02) * X +
               6.08784033347757E-01) *
                  E +
              R34 / (X - R34);
        RT2 = ((-3.62569791162153E-04 * X - 9.09231717268466E-03) * X +
               1.84336760556262E-01) *
                  E +
              R24 / (X - R24);
        RT1 = ((-4.07557525914600E-05 * X - 6.88846864931685E-04) * X +
               1.74725309199384E-02) *
                  E +
              R14 / (X - R14);
        WW4 = ((5.76631982000990E-06 * X - 7.89187283804890E-05) * X +
               3.28297971853126E-04) *
                  E +
              W44 * WW1;
        WW3 = ((2.08294969857230E-04 * X - 3.77489954837361E-03) * X +
               2.09857151617436E-02) *
                  E +
              W34 * WW1;
        WW2 = ((6.16374517326469E-04 * X - 1.26711744680092E-02) * X +
               8.14504890732155E-02) *
                  E +
              W24 * WW1;
        WW1 = WW1 - WW2 - WW3 - WW4;
    }
    else
    {
        WW1 = sqrt(constants::PI_4 / X);
        RT1 = R14 / (X - R14);
        RT2 = R24 / (X - R24);
        RT3 = R34 / (X - R34);
        RT4 = R44 / (X - R44);
        WW4 = W44 * WW1;
        WW3 = W34 * WW1;
        WW2 = W24 * WW1;
        WW1 = WW1 - WW2 - WW3 - WW4;
    }
    roots[0] = RT1;
    roots[1] = RT2;
    roots[2] = RT3;
    roots[3] = RT4;
    weights[0] = WW1;
    weights[1] = WW2;
    weights[2] = WW3;
    weights[3] = WW4;
    return 0;
}

int rys_root5(double X, double *roots, double *weights)
{
    double R15, R25, W25, R35, W35, R45, W45, R55, W55;
    double RT1, RT2, RT3, RT4, RT5, WW1, WW2, WW3, WW4, WW5;
    double Y, E, XXX;

    R15 = 1.17581320211778E-01;
    R25 = 1.07456201243690E+00;
    W25 = 2.70967405960535E-01;
    R35 = 3.08593744371754E+00;
    W35 = 3.82231610015404E-02;
    R45 = 6.41472973366203E+00;
    W45 = 1.51614186862443E-03;
    R55 = 1.18071894899717E+01;
    W55 = 8.62130526143657E-06;

    if (X < 3.e-7)
    {
        RT1 = 2.26659266316985E-02 - 2.15865967920897E-03 * X;
        RT2 = 2.31271692140903E-01 - 2.20258754389745E-02 * X;
        RT3 = 8.57346024118836E-01 - 8.16520023025515E-02 * X;
        RT4 = 2.97353038120346E+00 - 2.83193369647137E-01 * X;
        RT5 = 1.84151859759051E+01 - 1.75382723579439E+00 * X;
        WW1 = 2.95524224714752E-01 - 1.96867576909777E-02 * X;
        WW2 = 2.69266719309995E-01 - 5.61737590184721E-02 * X;
        WW3 = 2.19086362515981E-01 - 9.71152726793658E-02 * X;
        WW4 = 1.49451349150580E-01 - 1.02979262193565E-01 * X;
        WW5 = 6.66713443086877E-02 - 5.73782817488315E-02 * X;
    }
    else if (X < 1.0)
    {
        RT1 = ((((((-4.46679165328413E-11 * X + 1.21879111988031E-09) * X -
                   2.62975022612104E-08) *
                      X +
                  5.15106194905897E-07) *
                     X -
                 9.27933625824749E-06) *
                    X +
                1.51794097682482E-04) *
                   X -
               2.15865967920301E-03) *
                  X +
              2.26659266316985E-02;
        RT2 = ((((((1.93117331714174E-10 * X - 4.57267589660699E-09) * X +
                   2.48339908218932E-08) *
                      X +
                  1.50716729438474E-06) *
                     X -
                 6.07268757707381E-05) *
                    X +
                1.37506939145643E-03) *
                   X -
               2.20258754419939E-02) *
                  X +
              2.31271692140905E-01;
        RT3 = (((((4.84989776180094E-09 * X + 1.31538893944284E-07) * X -
                  2.766753852879E-06) *
                     X -
                 7.651163510626E-05) *
                    X +
                4.033058545972E-03) *
                   X -
               8.16520022916145E-02) *
                  X +
              8.57346024118779E-01;
        RT4 = ((((-2.48581772214623E-07 * X - 4.34482635782585E-06) * X -
                 7.46018257987630E-07) *
                    X +
                1.01210776517279E-02) *
                   X -
               2.83193369640005E-01) *
                  X +
              2.97353038120345E+00;
        RT5 = (((((-8.92432153868554E-09 * X + 1.77288899268988E-08) * X +
                  3.040754680666E-06) *
                     X +
                 1.058229325071E-04) *
                    X +
                4.596379534985E-02) *
                   X -
               1.75382723579114E+00) *
                  X +
              1.84151859759049E+01;
        WW1 = ((((((-2.03822632771791E-09 * X + 3.89110229133810E-08) * X -
                   5.84914787904823E-07) *
                      X +
                  8.30316168666696E-06) *
                     X -
                 1.13218402310546E-04) *
                    X +
                1.49128888586790E-03) *
                   X -
               1.96867576904816E-02) *
                  X +
              2.95524224714749E-01;
        WW2 = (((((((8.62848118397570E-09 * X - 1.38975551148989E-07) * X +
                    1.602894068228E-06) *
                       X -
                   1.646364300836E-05) *
                      X +
                  1.538445806778E-04) *
                     X -
                 1.28848868034502E-03) *
                    X +
                9.38866933338584E-03) *
                   X -
               5.61737590178812E-02) *
                  X +
              2.69266719309991E-01;
        WW3 = ((((((((-9.41953204205665E-09 * X + 1.47452251067755E-07) * X -
                     1.57456991199322E-06) *
                        X +
                    1.45098401798393E-05) *
                       X -
                   1.18858834181513E-04) *
                      X +
                  8.53697675984210E-04) *
                     X -
                 5.22877807397165E-03) *
                    X +
                2.60854524809786E-02) *
                   X -
               9.71152726809059E-02) *
                  X +
              2.19086362515979E-01;
        WW4 = ((((((((-3.84961617022042E-08 * X + 5.66595396544470E-07) * X -
                     5.52351805403748E-06) *
                        X +
                    4.53160377546073E-05) *
                       X -
                   3.22542784865557E-04) *
                      X +
                  1.95682017370967E-03) *
                     X -
                 9.77232537679229E-03) *
                    X +
                3.79455945268632E-02) *
                   X -
               1.02979262192227E-01) *
                  X +
              1.49451349150573E-01;
        WW5 = (((((((((4.09594812521430E-09 * X - 6.47097874264417E-08) * X +
                      6.743541482689E-07) *
                         X -
                     5.917993920224E-06) *
                        X +
                    4.531969237381E-05) *
                       X -
                   2.99102856679638E-04) *
                      X +
                  1.65695765202643E-03) *
                     X -
                 7.40671222520653E-03) *
                    X +
                2.50889946832192E-02) *
                   X -
               5.73782817487958E-02) *
                  X +
              6.66713443086877E-02;
    }
    else if (X < 5.0)
    {
        Y = X - 3.0E+00;
        RT1 = ((((((((-2.58163897135138E-14 * Y + 8.14127461488273E-13) * Y -
                     2.11414838976129E-11) *
                        Y +
                    5.09822003260014E-10) *
                       Y -
                   1.16002134438663E-08) *
                      Y +
                  2.46810694414540E-07) *
                     Y -
                 4.92556826124502E-06) *
                    Y +
                9.02580687971053E-05) *
                   Y -
               1.45190025120726E-03) *
                  Y +
              1.73416786387475E-02;
        RT2 = (((((((((1.04525287289788E-14 * Y + 5.44611782010773E-14) * Y -
                      4.831059411392E-12) *
                         Y +
                     1.136643908832E-10) *
                        Y -
                    1.104373076913E-09) *
                       Y -
                   2.35346740649916E-08) *
                      Y +
                  1.43772622028764E-06) *
                     Y -
                 4.23405023015273E-05) *
                    Y +
                9.12034574793379E-04) *
                   Y -
               1.52479441718739E-02) *
                  Y +
              1.76055265928744E-01;
        RT3 = (((((((((-6.89693150857911E-14 * Y + 5.92064260918861E-13) * Y +
                      1.847170956043E-11) *
                         Y -
                     3.390752744265E-10) *
                        Y -
                    2.995532064116E-09) *
                       Y +
                   1.57456141058535E-07) *
                      Y -
                  3.95859409711346E-07) *
                     Y -
                 9.58924580919747E-05) *
                    Y +
                3.23551502557785E-03) *
                   Y -
               5.97587007636479E-02) *
                  Y +
              6.46432853383057E-01;
        RT4 = ((((((((-3.61293809667763E-12 * Y - 2.70803518291085E-11) * Y +
                     8.83758848468769E-10) *
                        Y +
                    1.59166632851267E-08) *
                       Y -
                   1.32581997983422E-07) *
                      Y -
                  7.60223407443995E-06) *
                     Y -
                 7.41019244900952E-05) *
                    Y +
                9.81432631743423E-03) *
                   Y -
               2.23055570487771E-01) *
                  Y +
              2.21460798080643E+00;
        RT5 = (((((((((7.12332088345321E-13 * Y + 3.16578501501894E-12) * Y -
                      8.776668218053E-11) *
                         Y -
                     2.342817613343E-09) *
                        Y -
                    3.496962018025E-08) *
                       Y -
                   3.03172870136802E-07) *
                      Y +
                  1.50511293969805E-06) *
                     Y +
                 1.37704919387696E-04) *
                    Y +
                4.70723869619745E-02) *
                   Y -
               1.47486623003693E+00) *
                  Y +
              1.35704792175847E+01;
        WW1 = (((((((((1.04348658616398E-13 * Y - 1.94147461891055E-12) * Y +
                      3.485512360993E-11) *
                         Y -
                     6.277497362235E-10) *
                        Y +
                    1.100758247388E-08) *
                       Y -
                   1.88329804969573E-07) *
                      Y +
                  3.12338120839468E-06) *
                     Y -
                 5.04404167403568E-05) *
                    Y +
                8.00338056610995E-04) *
                   Y -
               1.30892406559521E-02) *
                  Y +
              2.47383140241103E-01;
        WW2 = (((((((((((3.23496149760478E-14 * Y - 5.24314473469311E-13) * Y +
                        7.743219385056E-12) *
                           Y -
                       1.146022750992E-10) *
                          Y +
                      1.615238462197E-09) *
                         Y -
                     2.15479017572233E-08) *
                        Y +
                    2.70933462557631E-07) *
                       Y -
                   3.18750295288531E-06) *
                      Y +
                  3.47425221210099E-05) *
                     Y -
                 3.45558237388223E-04) *
                    Y +
                3.05779768191621E-03) *
                   Y -
               2.29118251223003E-02) *
                  Y +
              1.59834227924213E-01;
        WW3 = ((((((((((((-3.42790561802876E-14 * Y + 5.26475736681542E-13) * Y -
                         7.184330797139E-12) *
                            Y +
                        9.763932908544E-11) *
                           Y -
                       1.244014559219E-09) *
                          Y +
                      1.472744068942E-08) *
                         Y -
                     1.611749975234E-07) *
                        Y +
                    1.616487851917E-06) *
                       Y -
                   1.46852359124154E-05) *
                      Y +
                  1.18900349101069E-04) *
                     Y -
                 8.37562373221756E-04) *
                    Y +
                4.93752683045845E-03) *
                   Y -
               2.25514728915673E-02) *
                  Y +
              6.95211812453929E-02;
        WW4 = (((((((((((((1.04072340345039E-14 * Y - 1.60808044529211E-13) *
                              Y +
                          2.183534866798E-12) *
                             Y -
                         2.939403008391E-11) *
                            Y +
                        3.679254029085E-10) *
                           Y -
                       4.23775673047899E-09) *
                          Y +
                      4.46559231067006E-08) *
                         Y -
                     4.26488836563267E-07) *
                        Y +
                    3.64721335274973E-06) *
                       Y -
                   2.74868382777722E-05) *
                      Y +
                  1.78586118867488E-04) *
                     Y -
                 9.68428981886534E-04) *
                    Y +
                4.16002324339929E-03) *
                   Y -
               1.28290192663141E-02) *
                  Y +
              2.22353727685016E-02;
        WW5 = ((((((((((((((-8.16770412525963E-16 * Y + 1.31376515047977E-14) *
                               Y -
                           1.856950818865E-13) *
                              Y +
                          2.596836515749E-12) *
                             Y -
                         3.372639523006E-11) *
                            Y +
                        4.025371849467E-10) *
                           Y -
                       4.389453269417E-09) *
                          Y +
                      4.332753856271E-08) *
                         Y -
                     3.82673275931962E-07) *
                        Y +
                    2.98006900751543E-06) *
                       Y -
                   2.00718990300052E-05) *
                      Y +
                  1.13876001386361E-04) *
                     Y -
                 5.23627942443563E-04) *
                    Y +
                1.83524565118203E-03) *
                   Y -
               4.37785737450783E-03) *
                  Y +
              5.36963805223095E-03;
    }
    else if (X < 10.0)
    {
        Y = X - 7.5E+00;
        RT1 = ((((((((-1.13825201010775E-14 * Y + 1.89737681670375E-13) * Y -
                     4.81561201185876E-12) *
                        Y +
                    1.56666512163407E-10) *
                       Y -
                   3.73782213255083E-09) *
                      Y +
                  9.15858355075147E-08) *
                     Y -
                 2.13775073585629E-06) *
                    Y +
                4.56547356365536E-05) *
                   Y -
               8.68003909323740E-04) *
                  Y +
              1.22703754069176E-02;
        RT2 = (((((((((-3.67160504428358E-15 * Y + 1.27876280158297E-14) * Y -
                      1.296476623788E-12) *
                         Y +
                     1.477175434354E-11) *
                        Y +
                    5.464102147892E-10) *
                       Y -
                   2.42538340602723E-08) *
                      Y +
                  8.20460740637617E-07) *
                     Y -
                 2.20379304598661E-05) *
                    Y +
                4.90295372978785E-04) *
                   Y -
               9.14294111576119E-03) *
                  Y +
              1.22590403403690E-01;
        RT3 = (((((((((1.39017367502123E-14 * Y - 6.96391385426890E-13) * Y +
                      1.176946020731E-12) *
                         Y +
                     1.725627235645E-10) *
                        Y -
                    3.686383856300E-09) *
                       Y +
                   2.87495324207095E-08) *
                      Y +
                  1.71307311000282E-06) *
                     Y -
                 7.94273603184629E-05) *
                    Y +
                2.00938064965897E-03) *
                   Y -
               3.63329491677178E-02) *
                  Y +
              4.34393683888443E-01;
        RT4 = ((((((((((-1.27815158195209E-14 * Y + 1.99910415869821E-14) * Y +
                       3.753542914426E-12) *
                          Y -
                      2.708018219579E-11) *
                         Y -
                     1.190574776587E-09) *
                        Y +
                    1.106696436509E-08) *
                       Y +
                   3.954955671326E-07) *
                      Y -
                  4.398596059588E-06) *
                     Y -
                 2.01087998907735E-04) *
                    Y +
                7.89092425542937E-03) *
                   Y -
               1.42056749162695E-01) *
                  Y +
              1.39964149420683E+00;
        RT5 = ((((((((((-1.19442341030461E-13 * Y - 2.34074833275956E-12) * Y +
                       6.861649627426E-12) *
                          Y +
                      6.082671496226E-10) *
                         Y +
                     5.381160105420E-09) *
                        Y -
                    6.253297138700E-08) *
                       Y -
                   2.135966835050E-06) *
                      Y -
                  2.373394341886E-05) *
                     Y +
                 2.88711171412814E-06) *
                    Y +
                4.85221195290753E-02) *
                   Y -
               1.04346091985269E+00) *
                  Y +
              7.89901551676692E+00;
        WW1 = (((((((((7.95526040108997E-15 * Y - 2.48593096128045E-13) * Y +
                      4.761246208720E-12) *
                         Y -
                     9.535763686605E-11) *
                        Y +
                    2.225273630974E-09) *
                       Y -
                   4.49796778054865E-08) *
                      Y +
                  9.17812870287386E-07) *
                     Y -
                 1.86764236490502E-05) *
                    Y +
                3.76807779068053E-04) *
                   Y -
               8.10456360143408E-03) *
                  Y +
              2.01097936411496E-01;
        WW2 = (((((((((((1.25678686624734E-15 * Y - 2.34266248891173E-14) * Y +
                        3.973252415832E-13) *
                           Y -
                       6.830539401049E-12) *
                          Y +
                      1.140771033372E-10) *
                         Y -
                     1.82546185762009E-09) *
                        Y +
                    2.77209637550134E-08) *
                       Y -
                   4.01726946190383E-07) *
                      Y +
                  5.48227244014763E-06) *
                     Y -
                 6.95676245982121E-05) *
                    Y +
                8.05193921815776E-04) *
                   Y -
               8.15528438784469E-03) *
                  Y +
              9.71769901268114E-02;
        WW3 = ((((((((((((-8.20929494859896E-16 * Y + 1.37356038393016E-14) * Y -
                         2.022863065220E-13) *
                            Y +
                        3.058055403795E-12) *
                           Y -
                       4.387890955243E-11) *
                          Y +
                      5.923946274445E-10) *
                         Y -
                     7.503659964159E-09) *
                        Y +
                    8.851599803902E-08) *
                       Y -
                   9.65561998415038E-07) *
                      Y +
                  9.60884622778092E-06) *
                     Y -
                 8.56551787594404E-05) *
                    Y +
                6.66057194311179E-04) *
                   Y -
               4.17753183902198E-03) *
                  Y +
              2.25443826852447E-02;
        WW4 = ((((((((((((((-1.08764612488790E-17 * Y + 1.85299909689937E-16) *
                               Y -
                           2.730195628655E-15) *
                              Y +
                          4.127368817265E-14) *
                             Y -
                         5.881379088074E-13) *
                            Y +
                        7.805245193391E-12) *
                           Y -
                       9.632707991704E-11) *
                          Y +
                      1.099047050624E-09) *
                         Y -
                     1.15042731790748E-08) *
                        Y +
                    1.09415155268932E-07) *
                       Y -
                   9.33687124875935E-07) *
                      Y +
                  7.02338477986218E-06) *
                     Y -
                 4.53759748787756E-05) *
                    Y +
                2.41722511389146E-04) *
                   Y -
               9.75935943447037E-04) *
                  Y +
              2.57520532789644E-03;
        WW5 = (((((((((((((((7.28996979748849E-19 * Y - 1.26518146195173E-17) * Y + 1.886145834486E-16) * Y - 2.876728287383E-15) * Y +
                          4.114588668138E-14) *
                             Y -
                         5.44436631413933E-13) *
                            Y +
                        6.64976446790959E-12) *
                           Y -
                       7.44560069974940E-11) *
                          Y +
                      7.57553198166848E-10) *
                         Y -
                     6.92956101109829E-09) *
                        Y +
                    5.62222859033624E-08) *
                       Y -
                   3.97500114084351E-07) *
                      Y +
                  2.39039126138140E-06) *
                     Y -
                 1.18023950002105E-05) *
                    Y +
                4.52254031046244E-05) *
                   Y -
               1.21113782150370E-04) *
                  Y +
              1.75013126731224E-04;
    }
    else if (X < 15.0)
    {
        Y = X - 12.5E+00;
        RT1 = ((((((((((-4.16387977337393E-17 * Y + 7.20872997373860E-16) * Y +
                       1.395993802064E-14) *
                          Y +
                      3.660484641252E-14) *
                         Y -
                     4.154857548139E-12) *
                        Y +
                    2.301379846544E-11) *
                       Y -
                   1.033307012866E-09) *
                      Y +
                  3.997777641049E-08) *
                     Y -
                 9.35118186333939E-07) *
                    Y +
                2.38589932752937E-05) *
                   Y -
               5.35185183652937E-04) *
                  Y +
              8.85218988709735E-03;
        RT2 = ((((((((((-4.56279214732217E-16 * Y + 6.24941647247927E-15) * Y +
                       1.737896339191E-13) *
                          Y +
                      8.964205979517E-14) *
                         Y -
                     3.538906780633E-11) *
                        Y +
                    9.561341254948E-11) *
                       Y -
                   9.772831891310E-09) *
                      Y +
                  4.240340194620E-07) *
                     Y -
                 1.02384302866534E-05) *
                    Y +
                2.57987709704822E-04) *
                   Y -
               5.54735977651677E-03) *
                  Y +
              8.68245143991948E-02;
        RT3 = ((((((((((-2.52879337929239E-15 * Y + 2.13925810087833E-14) * Y +
                       7.884307667104E-13) *
                          Y -
                      9.023398159510E-13) *
                         Y -
                     5.814101544957E-11) *
                        Y -
                    1.333480437968E-09) *
                       Y -
                   2.217064940373E-08) *
                      Y +
                  1.643290788086E-06) *
                     Y -
                 4.39602147345028E-05) *
                    Y +
                1.08648982748911E-03) *
                   Y -
               2.13014521653498E-02) *
                  Y +
              2.94150684465425E-01;
        RT4 = ((((((((((-6.42391438038888E-15 * Y + 5.37848223438815E-15) * Y +
                       8.960828117859E-13) *
                          Y +
                      5.214153461337E-11) *
                         Y -
                     1.106601744067E-10) *
                        Y -
                    2.007890743962E-08) *
                       Y +
                   1.543764346501E-07) *
                      Y +
                  4.520749076914E-06) *
                     Y -
                 1.88893338587047E-04) *
                    Y +
                4.73264487389288E-03) *
                   Y -
               7.91197893350253E-02) *
                  Y +
              8.60057928514554E-01;
        RT5 = (((((((((((-2.24366166957225E-14 * Y + 4.87224967526081E-14) * Y +
                        5.587369053655E-12) *
                           Y -
                       3.045253104617E-12) *
                          Y -
                      1.223983883080E-09) *
                         Y -
                     2.05603889396319E-09) *
                        Y +
                    2.58604071603561E-07) *
                       Y +
                   1.34240904266268E-06) *
                      Y -
                  5.72877569731162E-05) *
                     Y -
                 9.56275105032191E-04) *
                    Y +
                4.23367010370921E-02) *
                   Y -
               5.76800927133412E-01) *
                  Y +
              3.87328263873381E+00;
        WW1 = (((((((((8.98007931950169E-15 * Y + 7.25673623859497E-14) * Y +
                      5.851494250405E-14) *
                         Y -
                     4.234204823846E-11) *
                        Y +
                    3.911507312679E-10) *
                       Y -
                   9.65094802088511E-09) *
                      Y +
                  3.42197444235714E-07) *
                     Y -
                 7.51821178144509E-06) *
                    Y +
                1.94218051498662E-04) *
                   Y -
               5.38533819142287E-03) *
                  Y +
              1.68122596736809E-01;
        WW2 = ((((((((((-1.05490525395105E-15 * Y + 1.96855386549388E-14) * Y -
                       5.500330153548E-13) *
                          Y +
                      1.003849567976E-11) *
                         Y -
                     1.720997242621E-10) *
                        Y +
                    3.533277061402E-09) *
                       Y -
                   6.389171736029E-08) *
                      Y +
                  1.046236652393E-06) *
                     Y -
                 1.73148206795827E-05) *
                    Y +
                2.57820531617185E-04) *
                   Y -
               3.46188265338350E-03) *
                  Y +
              7.03302497508176E-02;
        WW3 = (((((((((((3.60020423754545E-16 * Y - 6.24245825017148E-15) * Y +
                        9.945311467434E-14) *
                           Y -
                       1.749051512721E-12) *
                          Y +
                      2.768503957853E-11) *
                         Y -
                     4.08688551136506E-10) *
                        Y +
                    6.04189063303610E-09) *
                       Y -
                   8.23540111024147E-08) *
                      Y +
                  1.01503783870262E-06) *
                     Y -
                 1.20490761741576E-05) *
                    Y +
                1.26928442448148E-04) *
                   Y -
               1.05539461930597E-03) *
                  Y +
              1.15543698537013E-02;
        WW4 = (((((((((((((2.51163533058925E-18 * Y - 4.31723745510697E-17) *
                              Y +
                          6.557620865832E-16) *
                             Y -
                         1.016528519495E-14) *
                            Y +
                        1.491302084832E-13) *
                           Y -
                       2.06638666222265E-12) *
                          Y +
                      2.67958697789258E-11) *
                         Y -
                     3.23322654638336E-10) *
                        Y +
                    3.63722952167779E-09) *
                       Y -
                   3.75484943783021E-08) *
                      Y +
                  3.49164261987184E-07) *
                     Y -
                 2.92658670674908E-06) *
                    Y +
                2.12937256719543E-05) *
                   Y -
               1.19434130620929E-04) *
                  Y +
              6.45524336158384E-04;
        WW5 = ((((((((((((((-1.29043630202811E-19 * Y + 2.16234952241296E-18) *
                               Y -
                           3.107631557965E-17) *
                              Y +
                          4.570804313173E-16) *
                             Y -
                         6.301348858104E-15) *
                            Y +
                        8.031304476153E-14) *
                           Y -
                       9.446196472547E-13) *
                          Y +
                      1.018245804339E-11) *
                         Y -
                     9.96995451348129E-11) *
                        Y +
                    8.77489010276305E-10) *
                       Y -
                   6.84655877575364E-09) *
                      Y +
                  4.64460857084983E-08) *
                     Y -
                 2.66924538268397E-07) *
                    Y +
                1.24621276265907E-06) *
                   Y -
               4.30868944351523E-06) *
                  Y +
              9.94307982432868E-06;
    }
    else if (X < 20.0)
    {
        Y = X - 17.5E+00;
        RT1 = ((((((((((1.91875764545740E-16 * Y + 7.8357401095707E-16) * Y -
                       3.260875931644E-14) *
                          Y -
                      1.186752035569E-13) *
                         Y +
                     4.275180095653E-12) *
                        Y +
                    3.357056136731E-11) *
                       Y -
                   1.123776903884E-09) *
                      Y +
                  1.231203269887E-08) *
                     Y -
                 3.99851421361031E-07) *
                    Y +
                1.45418822817771E-05) *
                   Y -
               3.49912254976317E-04) *
                  Y +
              6.67768703938812E-03;
        RT2 = ((((((((((2.02778478673555E-15 * Y + 1.01640716785099E-14) * Y -
                       3.385363492036E-13) *
                          Y -
                      1.615655871159E-12) *
                         Y +
                     4.527419140333E-11) *
                        Y +
                    3.853670706486E-10) *
                       Y -
                   1.184607130107E-08) *
                      Y +
                  1.347873288827E-07) *
                     Y -
                 4.47788241748377E-06) *
                    Y +
                1.54942754358273E-04) *
                   Y -
               3.55524254280266E-03) *
                  Y +
              6.44912219301603E-02;
        RT3 = ((((((((((7.79850771456444E-15 * Y + 6.00464406395001E-14) * Y -
                       1.249779730869E-12) *
                          Y -
                      1.020720636353E-11) *
                         Y +
                     1.814709816693E-10) *
                        Y +
                    1.766397336977E-09) *
                       Y -
                   4.603559449010E-08) *
                      Y +
                  5.863956443581E-07) *
                     Y -
                 2.03797212506691E-05) *
                    Y +
                6.31405161185185E-04) *
                   Y -
               1.30102750145071E-02) *
                  Y +
              2.10244289044705E-01;
        RT4 = (((((((((((-2.92397030777912E-15 * Y + 1.94152129078465E-14) * Y +
                        4.859447665850E-13) *
                           Y -
                       3.217227223463E-12) *
                          Y -
                      7.484522135512E-11) *
                         Y +
                     7.19101516047753E-10) *
                        Y +
                    6.88409355245582E-09) *
                       Y -
                   1.44374545515769E-07) *
                      Y +
                  2.74941013315834E-06) *
                     Y -
                 1.02790452049013E-04) *
                    Y +
                2.59924221372643E-03) *
                   Y -
               4.35712368303551E-02) *
                  Y +
              5.62170709585029E-01;
        RT5 = (((((((((((1.17976126840060E-14 * Y + 1.24156229350669E-13) * Y -
                        3.892741622280E-12) *
                           Y -
                       7.755793199043E-12) *
                          Y +
                      9.492190032313E-10) *
                         Y -
                     4.98680128123353E-09) *
                        Y -
                    1.81502268782664E-07) *
                       Y +
                   2.69463269394888E-06) *
                      Y +
                  2.50032154421640E-05) *
                     Y -
                 1.33684303917681E-03) *
                    Y +
                2.29121951862538E-02) *
                   Y -
               2.45653725061323E-01) *
                  Y +
              1.89999883453047E+00;
        WW1 = ((((((((((1.74841995087592E-15 * Y - 6.95671892641256E-16) * Y -
                       3.000659497257E-13) *
                          Y +
                      2.021279817961E-13) *
                         Y +
                     3.853596935400E-11) *
                        Y +
                    1.461418533652E-10) *
                       Y -
                   1.014517563435E-08) *
                      Y +
                  1.132736008979E-07) *
                     Y -
                 2.86605475073259E-06) *
                    Y +
                1.21958354908768E-04) *
                   Y -
               3.86293751153466E-03) *
                  Y +
              1.45298342081522E-01;
        WW2 = ((((((((((-1.11199320525573E-15 * Y + 1.85007587796671E-15) * Y +
                       1.220613939709E-13) *
                          Y +
                      1.275068098526E-12) *
                         Y -
                     5.341838883262E-11) *
                        Y +
                    6.161037256669E-10) *
                       Y -
                   1.009147879750E-08) *
                      Y +
                  2.907862965346E-07) *
                     Y -
                 6.12300038720919E-06) *
                    Y +
                1.00104454489518E-04) *
                   Y -
               1.80677298502757E-03) *
                  Y +
              5.78009914536630E-02;
        WW3 = ((((((((((-9.49816486853687E-16 * Y + 6.67922080354234E-15) * Y +
                       2.606163540537E-15) *
                          Y +
                      1.983799950150E-12) *
                         Y -
                     5.400548574357E-11) *
                        Y +
                    6.638043374114E-10) *
                       Y -
                   8.799518866802E-09) *
                      Y +
                  1.791418482685E-07) *
                     Y -
                 2.96075397351101E-06) *
                    Y +
                3.38028206156144E-05) *
                   Y -
               3.58426847857878E-04) *
                  Y +
              8.39213709428516E-03;
        WW4 = (((((((((((1.33829971060180E-17 * Y - 3.44841877844140E-16) * Y +
                        4.745009557656E-15) *
                           Y -
                       6.033814209875E-14) *
                          Y +
                      1.049256040808E-12) *
                         Y -
                     1.70859789556117E-11) *
                        Y +
                    2.15219425727959E-10) *
                       Y -
                   2.52746574206884E-09) *
                      Y +
                  3.27761714422960E-08) *
                     Y -
                 3.90387662925193E-07) *
                    Y +
                3.46340204593870E-06) *
                   Y -
               2.43236345136782E-05) *
                  Y +
              3.54846978585226E-04;
        WW5 = (((((((((((((2.69412277020887E-20 * Y - 4.24837886165685E-19) *
                              Y +
                          6.030500065438E-18) *
                             Y -
                         9.069722758289E-17) *
                            Y +
                        1.246599177672E-15) *
                           Y -
                       1.56872999797549E-14) *
                          Y +
                      1.87305099552692E-13) *
                         Y -
                     2.09498886675861E-12) *
                        Y +
                    2.11630022068394E-11) *
                       Y -
                   1.92566242323525E-10) *
                      Y +
                  1.62012436344069E-09) *
                     Y -
                 1.23621614171556E-08) *
                    Y +
                7.72165684563049E-08) *
                   Y -
               3.59858901591047E-07) *
                  Y +
              2.43682618601000E-06;
    }
    else if (X < 25.0)
    {
        Y = X - 22.5E+00;
        RT1 = (((((((((-1.13927848238726E-15 * Y + 7.39404133595713E-15) * Y +
                      1.445982921243E-13) *
                         Y -
                     2.676703245252E-12) *
                        Y +
                    5.823521627177E-12) *
                       Y +
                   2.17264723874381E-10) *
                      Y +
                  3.56242145897468E-09) *
                     Y -
                 3.03763737404491E-07) *
                    Y +
                9.46859114120901E-06) *
                   Y -
               2.30896753853196E-04) *
                  Y +
              5.24663913001114E-03;
        RT2 = ((((((((((2.89872355524581E-16 * Y - 1.22296292045864E-14) * Y +
                       6.184065097200E-14) *
                          Y +
                      1.649846591230E-12) *
                         Y -
                     2.729713905266E-11) *
                        Y +
                    3.709913790650E-11) *
                       Y +
                   2.216486288382E-09) *
                      Y +
                  4.616160236414E-08) *
                     Y -
                 3.32380270861364E-06) *
                    Y +
                9.84635072633776E-05) *
                   Y -
               2.30092118015697E-03) *
                  Y +
              5.00845183695073E-02;
        RT3 = ((((((((((1.97068646590923E-15 * Y - 4.89419270626800E-14) * Y +
                       1.136466605916E-13) *
                          Y +
                      7.546203883874E-12) *
                         Y -
                     9.635646767455E-11) *
                        Y -
                    8.295965491209E-11) *
                       Y +
                   7.534109114453E-09) *
                      Y +
                  2.699970652707E-07) *
                     Y -
                 1.42982334217081E-05) *
                    Y +
                3.78290946669264E-04) *
                   Y -
               8.03133015084373E-03) *
                  Y +
              1.58689469640791E-01;
        RT4 = ((((((((((1.33642069941389E-14 * Y - 1.55850612605745E-13) * Y -
                       7.522712577474E-13) *
                          Y +
                      3.209520801187E-11) *
                         Y -
                     2.075594313618E-10) *
                        Y -
                    2.070575894402E-09) *
                       Y +
                   7.323046997451E-09) *
                      Y +
                  1.851491550417E-06) *
                     Y -
                 6.37524802411383E-05) *
                    Y +
                1.36795464918785E-03) *
                   Y -
               2.42051126993146E-02) *
                  Y +
              3.97847167557815E-01;
        RT5 = ((((((((((-6.07053986130526E-14 * Y + 1.04447493138843E-12) * Y -
                       4.286617818951E-13) *
                          Y -
                      2.632066100073E-10) *
                         Y +
                     4.804518986559E-09) *
                        Y -
                    1.835675889421E-08) *
                       Y -
                   1.068175391334E-06) *
                      Y +
                  3.292234974141E-05) *
                     Y -
                 5.94805357558251E-04) *
                    Y +
                8.29382168612791E-03) *
                   Y -
               9.93122509049447E-02) *
                  Y +
              1.09857804755042E+00;
        WW1 = (((((((((-9.10338640266542E-15 * Y + 1.00438927627833E-13) * Y +
                      7.817349237071E-13) *
                         Y -
                     2.547619474232E-11) *
                        Y +
                    1.479321506529E-10) *
                       Y +
                   1.52314028857627E-09) *
                      Y +
                  9.20072040917242E-09) *
                     Y -
                 2.19427111221848E-06) *
                    Y +
                8.65797782880311E-05) *
                   Y -
               2.82718629312875E-03) *
                  Y +
              1.28718310443295E-01;
        WW2 = (((((((((5.52380927618760E-15 * Y - 6.43424400204124E-14) * Y -
                      2.358734508092E-13) *
                         Y +
                     8.261326648131E-12) *
                        Y +
                    9.229645304956E-11) *
                       Y -
                   5.68108973828949E-09) *
                      Y +
                  1.22477891136278E-07) *
                     Y -
                 2.11919643127927E-06) *
                    Y +
                4.23605032368922E-05) *
                   Y -
               1.14423444576221E-03) *
                  Y +
              5.06607252890186E-02;
        WW3 = (((((((((3.99457454087556E-15 * Y - 5.11826702824182E-14) * Y -
                      4.157593182747E-14) *
                         Y +
                     4.214670817758E-12) *
                        Y +
                    6.705582751532E-11) *
                       Y -
                   3.36086411698418E-09) *
                      Y +
                  6.07453633298986E-08) *
                     Y -
                 7.40736211041247E-07) *
                    Y +
                8.84176371665149E-06) *
                   Y -
               1.72559275066834E-04) *
                  Y +
              7.16639814253567E-03;
        WW4 = (((((((((((-2.14649508112234E-18 * Y - 2.45525846412281E-18) * Y +
                        6.126212599772E-16) *
                           Y -
                       8.526651626939E-15) *
                          Y +
                      4.826636065733E-14) *
                         Y -
                     3.39554163649740E-13) *
                        Y +
                    1.67070784862985E-11) *
                       Y -
                   4.42671979311163E-10) *
                      Y +
                  6.77368055908400E-09) *
                     Y -
                 7.03520999708859E-08) *
                    Y +
                6.04993294708874E-07) *
                   Y -
               7.80555094280483E-06) *
                  Y +
              2.85954806605017E-04;
        WW5 = ((((((((((((-5.63938733073804E-21 * Y + 6.92182516324628E-20) * Y -
                         1.586937691507E-18) *
                            Y +
                        3.357639744582E-17) *
                           Y -
                       4.810285046442E-16) *
                          Y +
                      5.386312669975E-15) *
                         Y -
                     6.117895297439E-14) *
                        Y +
                    8.441808227634E-13) *
                       Y -
                   1.18527596836592E-11) *
                      Y +
                  1.36296870441445E-10) *
                     Y -
                 1.17842611094141E-09) *
                    Y +
                7.80430641995926E-09) *
                   Y -
               5.97767417400540E-08) *
                  Y +
              1.65186146094969E-06;
    }
    else if (X < 40)
    {
        WW1 = sqrt(constants::PI_4 / X);
        E = exp(-X);
        RT1 = ((((((((-1.73363958895356E-06 * X + 1.19921331441483E-04) * X -
                     1.59437614121125E-02) *
                        X +
                    1.13467897349442E+00) *
                       X -
                   4.47216460864586E+01) *
                      X +
                  1.06251216612604E+03) *
                     X -
                 1.52073917378512E+04) *
                    X +
                1.20662887111273E+05) *
                   X -
               4.07186366852475E+05) *
                  E +
              R15 / (X - R15);
        RT2 = ((((((((-1.60102542621710E-05 * X + 1.10331262112395E-03) * X -
                     1.50043662589017E-01) *
                        X +
                    1.05563640866077E+01) *
                       X -
                   4.10468817024806E+02) *
                      X +
                  9.62604416506819E+03) *
                     X -
                 1.35888069838270E+05) *
                    X +
                1.06107577038340E+06) *
                   X -
               3.51190792816119E+06) *
                  E +
              R25 / (X - R25);
        RT3 = ((((((((-4.48880032128422E-05 * X + 2.69025112122177E-03) * X -
                     4.01048115525954E-01) *
                        X +
                    2.78360021977405E+01) *
                       X -
                   1.04891729356965E+03) *
                      X +
                  2.36985942687423E+04) *
                     X -
                 3.19504627257548E+05) *
                    X +
                2.34879693563358E+06) *
                   X -
               7.16341568174085E+06) *
                  E +
              R35 / (X - R35);
        RT4 = ((((((((-6.38526371092582E-05 * X - 2.29263585792626E-03) * X -
                     7.65735935499627E-02) *
                        X +
                    9.12692349152792E+00) *
                       X -
                   2.32077034386717E+02) *
                      X +
                  2.81839578728845E+02) *
                     X +
                 9.59529683876419E+04) *
                    X -
                1.77638956809518E+06) *
                   X +
               1.02489759645410E+07) *
                  E +
              R45 / (X - R45);
        RT5 = ((((((((-3.59049364231569E-05 * X - 2.25963977930044E-02) * X +
                     1.12594870794668E+00) *
                        X -
                    4.56752462103909E+01) *
                       X +
                   1.05804526830637E+03) *
                      X -
                  1.16003199605875E+04) *
                     X -
                 4.07297627297272E+04) *
                    X +
                2.22215528319857E+06) *
                   X -
               1.61196455032613E+07) *
                  E +
              R55 / (X - R55);
        WW5 = (((((((((-4.61100906133970E-10 * X + 1.43069932644286E-07) * X -
                      1.63960915431080E-05) *
                         X +
                     1.15791154612838E-03) *
                        X -
                    5.30573476742071E-02) *
                       X +
                   1.61156533367153E+00) *
                      X -
                  3.23248143316007E+01) *
                     X +
                 4.12007318109157E+02) *
                    X -
                3.02260070158372E+03) *
                   X +
               9.71575094154768E+03) *
                  E +
              W55 * WW1;
        WW4 = (((((((((-2.40799435809950E-08 * X + 8.12621667601546E-06) * X -
                      9.04491430884113E-04) *
                         X +
                     6.37686375770059E-02) *
                        X -
                    2.96135703135647E+00) *
                       X +
                   9.15142356996330E+01) *
                      X -
                  1.86971865249111E+03) *
                     X +
                 2.42945528916947E+04) *
                    X -
                1.81852473229081E+05) *
                   X +
               5.96854758661427E+05) *
                  E +
              W45 * WW1;
        WW3 = ((((((((1.83574464457207E-05 * X - 1.54837969489927E-03) * X +
                     1.18520453711586E-01) *
                        X -
                    6.69649981309161E+00) *
                       X +
                   2.44789386487321E+02) *
                      X -
                  5.68832664556359E+03) *
                     X +
                 8.14507604229357E+04) *
                    X -
                6.55181056671474E+05) *
                   X +
               2.26410896607237E+06) *
                  E +
              W35 * WW1;
        WW2 = ((((((((2.77778345870650E-05 * X - 2.22835017655890E-03) * X +
                     1.61077633475573E-01) *
                        X -
                    8.96743743396132E+00) *
                       X +
                   3.28062687293374E+02) *
                      X -
                  7.65722701219557E+03) *
                     X +
                 1.10255055017664E+05) *
                    X -
                8.92528122219324E+05) *
                   X +
               3.10638627744347E+06) *
                  E +
              W25 * WW1;
        WW1 = WW1 - 0.01962E+00 * E - WW2 - WW3 - WW4 - WW5;
    }
    else if (X < 59.0)
    {
        WW1 = sqrt(constants::PI_4 / X);
        XXX = X * X * X;
        E = XXX * exp(-X);
        RT1 = (((-2.43758528330205E-02 * X + 2.07301567989771E+00) * X -
                6.45964225381113E+01) *
                   X +
               7.14160088655470E+02) *
                  E +
              R15 / (X - R15);
        RT2 = (((-2.28861955413636E-01 * X + 1.93190784733691E+01) * X -
                5.99774730340912E+02) *
                   X +
               6.61844165304871E+03) *
                  E +
              R25 / (X - R25);
        RT3 = (((-6.95053039285586E-01 * X + 5.76874090316016E+01) * X -
                1.77704143225520E+03) *
                   X +
               1.95366082947811E+04) *
                  E +
              R35 / (X - R35);
        RT4 = (((-1.58072809087018E+00 * X + 1.27050801091948E+02) * X -
                3.86687350914280E+03) *
                   X +
               4.23024828121420E+04) *
                  E +
              R45 / (X - R45);
        RT5 = (((-3.33963830405396E+00 * X + 2.51830424600204E+02) * X -
                7.57728527654961E+03) *
                   X +
               8.21966816595690E+04) *
                  E +
              R55 / (X - R55);
        E = XXX * E;
        WW5 = ((1.35482430510942E-08 * X - 3.27722199212781E-07) * X +
               2.41522703684296E-06) *
                  E +
              W55 * WW1;
        WW4 = ((1.23464092261605E-06 * X - 3.55224564275590E-05) * X +
               3.03274662192286E-04) *
                  E +
              W45 * WW1;
        WW3 = ((1.34547929260279E-05 * X - 4.19389884772726E-04) * X +
               3.87706687610809E-03) *
                  E +
              W35 * WW1;
        WW2 = ((2.09539509123135E-05 * X - 6.87646614786982E-04) * X +
               6.68743788585688E-03) *
                  E +
              W25 * WW1;
        WW1 = WW1 - WW2 - WW3 - WW4 - WW5;
    }
    else
    {
        WW1 = sqrt(constants::PI_4 / X);
        RT1 = R15 / (X - R15);
        RT2 = R25 / (X - R25);
        RT3 = R35 / (X - R35);
        RT4 = R45 / (X - R45);
        RT5 = R55 / (X - R55);
        WW2 = W25 * WW1;
        WW3 = W35 * WW1;
        WW4 = W45 * WW1;
        WW5 = W55 * WW1;
        WW1 = WW1 - WW2 - WW3 - WW4 - WW5;
    }
    roots[0] = RT1;
    roots[1] = RT2;
    roots[2] = RT3;
    roots[3] = RT4;
    roots[4] = RT5;
    weights[0] = WW1;
    weights[1] = WW2;
    weights[2] = WW3;
    weights[3] = WW4;
    weights[4] = WW5;
    return 0;
}

typedef int QuadratureFunction(int n, double x, double lower, double *roots, double *weights);
int segment_solve(int n, double x, double lower, double *u, double *w,
                  double breakpoint, QuadratureFunction fn1, QuadratureFunction fn2)
{
    int error;
    if (x <= breakpoint)
    {
        error = fn1(n, x, lower, u, w);
    }
    else
    {
        error = fn2(n, x, lower, u, w);
    }
    if (error)
    {
        error = lrys_schmidt(n, x, lower, u, w);
    }
    return error;
}

#define POLYNOMIAL_VALUE1(p, a, order, x) \
    p = a[order];                         \
    for (i = 1; i <= order; i++)          \
    {                                     \
        p = p * x + a[order - i];         \
    }

#define SET_ZERO(a, n, start)   \
    for (k = start; k < n; ++k) \
    {                           \
        for (i = 0; i < n; ++i) \
        {                       \
            a[i + k * n] = 0;   \
        }                       \
    }

int R_dnode(double *a, double *roots, int order)
{
    const double accrt = 1e-15;
    double x0, x1, xi, x1init, p0, p1, pi, p1init;
    int i, m, n;

    x1init = 0;
    p1init = a[0];
    for (m = 0; m < order; ++m)
    {
        x0 = x1init;
        p0 = p1init;
        x1init = roots[m];
        POLYNOMIAL_VALUE1(p1init, a, order, x1init);

        // When all coefficients a are 0, short-circuit the rest code to
        // ensure the roots from the lower order polynomials are preserved
        if (p1init == 0)
        {
            // roots[m] = x1init;
            continue;
        }
        if (p0 * p1init > 0)
        {
            fprintf(stderr, "ROOT NUMBER %d WAS NOT FOUND FOR POLYNOMIAL OF ORDER %d\n",
                    m, order);
            return 1;
        }
        if (x0 <= x1init)
        {
            x1 = x1init;
            p1 = p1init;
        }
        else
        {
            x1 = x0;
            p1 = p0;
            x0 = x1init;
            p0 = p1init;
        }
        // interpolate/extrapolate between [x0,x1]
        if (p1 == 0)
        {
            roots[m] = x1;
            continue;
        }
        else if (p0 == 0)
        {
            roots[m] = x0;
            continue;
        }
        else
        {
            xi = x0 + (x0 - x1) / (p1 - p0) * p0;
        }
        n = 0;
        while (fabs(x1 - x0) > x1 * accrt)
        {
            n++;
            if (n > 200)
            {
                fprintf(stderr, "libcint::rys_roots NO CONV. IN R_dnode\n");
                return 1;
            }
            POLYNOMIAL_VALUE1(pi, a, order, xi);
            if (pi == 0)
            {
                break;
            }
            else if (p0 * pi <= 0)
            {
                x1 = xi;
                p1 = pi;
                xi = x0 * .25 + xi * .75;
            }
            else
            {
                x0 = xi;
                p0 = pi;
                xi = xi * .75 + x1 * .25;
            }
            POLYNOMIAL_VALUE1(pi, a, order, xi);
            if (pi == 0)
            {
                break;
            }
            else if (p0 * pi <= 0)
            {
                x1 = xi;
                p1 = pi;
            }
            else
            {
                x0 = xi;
                p0 = pi;
            }

            xi = x0 + (x0 - x1) / (p1 - p0) * p0;
        }
        roots[m] = xi;
    }
    return 0;
}

void _qr_step(double *A, int nroots, int n0, int n1, double shift)
{
    int m1 = n0 + 1;
    int j, k, m3, j1, j2;
    double c = A[n0 * nroots + n0] - shift;
    double s = A[m1 * nroots + n0];
    double v = sqrt(c * c + s * s);
    double x, y;

    if (v == 0)
    {
        v = 1;
        c = 1;
        s = 0;
    }
    v = 1. / v;
    c *= v;
    s *= v;

    for (k = n0; k < nroots; k++)
    {
        // apply givens rotation from the left
        x = A[n0 * nroots + k];
        y = A[m1 * nroots + k];
        A[n0 * nroots + k] = c * x + s * y;
        A[m1 * nroots + k] = c * y - s * x;
    }

    m3 = std::min(n1, n0 + 3);
    for (k = 0; k < m3; k++)
    {
        // apply givens rotation from the right
        x = A[k * nroots + n0];
        y = A[k * nroots + m1];
        A[k * nroots + n0] = c * x + s * y;
        A[k * nroots + m1] = c * y - s * x;
    }

    for (j = n0; j < n1 - 2; j++)
    {
        j1 = j + 1;
        j2 = j + 2;
        // calculate givens rotation
        c = A[j1 * nroots + j];
        s = A[j2 * nroots + j];
        v = sqrt(c * c + s * s);
        A[j1 * nroots + j] = v;
        A[j2 * nroots + j] = 0;

        if (v == 0)
        {
            v = 1;
            c = 1;
            s = 0;
        }
        v = 1. / v;
        c *= v;
        s *= v;

        for (k = j1; k < nroots; k++)
        {
            // apply givens rotation from the left
            x = A[j1 * nroots + k];
            y = A[j2 * nroots + k];
            A[j1 * nroots + k] = c * x + s * y;
            A[j2 * nroots + k] = c * y - s * x;
        }
        m3 = std::min(n1, j + 4);
        for (k = 0; k < m3; k++)
        {
            // apply givens rotation from the right
            x = A[k * nroots + j1];
            y = A[k * nroots + j2];
            A[k * nroots + j1] = c * x + s * y;
            A[k * nroots + j2] = c * y - s * x;
        }
    }
}

int _hessenberg_qr(double *A, int nroots)
{
    double eps = 1e-15;
    int maxits = 30;
    int n0 = 0;
    int n1 = nroots;
    int its = 0;
    int k, ic, k1;
    for (ic = 0; ic < nroots * maxits; ic++)
    {
        k = n0;
        while (k + 1 < n1)
        {
            double s = fabs(A[k * nroots + k]) + fabs(A[(k + 1) * nroots + k + 1]);
            if (fabs(A[(k + 1) * nroots + k]) < eps * s)
            {
                break;
            }
            k += 1;
        }

        k1 = k + 1;
        if (k1 < n1)
        {
            // deflation found at position (k+1, k)
            A[k1 * nroots + k] = 0;
            n0 = k1;
            its = 0;

            if (n0 + 1 >= n1)
            {
                // block of size at most two has converged
                n0 = 0;
                n1 = k1;
                if (n1 < 2)
                {
                    // QR algorithm has converged
                    return 0;
                }
            }
        }
        else
        {
            int m1 = n1 - 1;
            int m2 = n1 - 2;
            double a11 = A[m1 * nroots + m1];
            double a22 = A[m2 * nroots + m2];
            double shift;
            double t = a11 + a22;
            double s = pow(a11 - a22, 2);
            s += 4 * A[m1 * nroots + m2] * A[m2 * nroots + m1];
            if (s > 0)
            {
                s = sqrt(s);
                double a = (t + s) * .5;
                double b = (t - s) * .5;
                if (fabs(a11 - a) > fabs(a11 - b))
                {
                    shift = b;
                }
                else
                {
                    shift = a;
                }
            }
            else
            {
                if (n1 == 2)
                {
                    fprintf(stderr, "hessenberg_qr: failed to find real roots\n");
                    return 1;
                }
                shift = t * .5;
            }
            its += 1;
            _qr_step(A, nroots, n0, n1, shift);
            if (its > maxits)
            {
                fprintf(stderr, "hessenberg_qr: failed to converge after %d steps\n", its);
                return 1;
            }
        }
    }
    fprintf(stderr, "hessenberg_qr failed\n");
    return 1;
}

int polynomial_roots(double *roots, double *cs, int nroots)
{
    if (nroots == 1)
    {
        roots[0] = -cs[2] / cs[3];
        return 0;
    }
    else if (nroots == 2)
    {
        double dum = sqrt(pow(cs[2 * 3 + 1], 2) - 4 * cs[2 * 3 + 0] * cs[2 * 3 + 2]);
        roots[0] = (-cs[2 * 3 + 1] - dum) / cs[2 * 3 + 2] / 2;
        roots[1] = (-cs[2 * 3 + 1] + dum) / cs[2 * 3 + 2] / 2;
        return 0;
    }

    double A[32 * 32];
    int nroots1 = nroots + 1;
    // reuse the buffer in coefficients
    int i;
    double fac = -1. / cs[nroots * nroots1 + nroots];
    for (i = 0; i < nroots; i++)
    {
        A[nroots - 1 - i] = cs[nroots * nroots1 + i] * fac;
    }
    for (i = nroots; i < nroots * nroots; i++)
    {
        A[i] = 0;
    }
    for (i = 0; i < nroots - 1; i++)
    {
        A[(i + 1) * nroots + i] = 1.;
    }
    int err = _hessenberg_qr(A, nroots);
    if (err == 0)
    {
        for (i = 0; i < nroots; i++)
        {
            roots[nroots - 1 - i] = A[i * nroots + i];
        }
    }
    else
    {
        int k, order;
        double *a;
        double dum = sqrt(cs[2 * nroots1 + 1] * cs[2 * nroots1 + 1] - 4 * cs[2 * nroots1 + 0] * cs[2 * nroots1 + 2]);
        roots[0] = .5 * (-cs[2 * nroots1 + 1] - dum) / cs[2 * nroots1 + 2];
        roots[1] = .5 * (-cs[2 * nroots1 + 1] + dum) / cs[2 * nroots1 + 2];
        for (i = 2; i < nroots; i++)
        {
            roots[i] = 1;
        }
        for (k = 2; k < nroots; ++k)
        {
            order = k + 1;
            a = cs + order * nroots1;
            err = R_dnode(a, roots, order);
            if (err)
            {
                break;
            }
        }
    }
    return err;
}

void fmt1_lgamma_inc_like(long double *f, long double t, int m)
{
    long double b = m + 0.5l;
    long double bi;
    long double e = .5l * expl(-t);
    long double x = e;
    long double s = e;
    long double tol = 2e-20 * e;
    int i;
    for (bi = b + 1.; x > tol; bi += 1.)
    {
        x *= t / bi;
        s += x;
    }
    f[m] = s / b;
    for (i = m; i > 0; i--)
    {
        b -= 1;
        f[i - 1] = (e + t * f[i]) / b;
    }
}

void lgamma_inc_like(long double *f, long double t, int m)
{
    if (t < constants::TURNOVER_POINT[m])
    {
        fmt1_lgamma_inc_like(f, t, m);
    }
    else
    {
        int i;
        long double tt = sqrtl(t);
        f[0] = constants::SQRTPI_4L / tt * erfl(tt);
        if (m > 0)
        {
            long double e = expl(-t);
            long double b = .5l / t;
            for (i = 1; i <= m; i++)
                f[i] = b * ((2 * i - 1) * f[i - 1] - e);
        }
    }
}

void fmt1_lerfc_like(long double *f, long double t, long double lower, int m)
{
    int i;
    long double lower2 = lower * lower;
    long double b = m + 0.5l;
    long double bi;
    long double e = .5l * expl(-t);
    long double e1 = .5l * expl(-t * lower2) * lower;
    e1 *= powl(lower2, m);
    long double x = e;
    long double x1 = e1;
    long double s = e - e1;
    long double div = 1.l;
    long double delta = s;
    long double tol = 2e-20 * fabsl(delta);
    for (bi = b + 1.l; fabsl(delta) > tol; bi += 1.l)
    {
        div *= t / bi;
        x1 *= lower2;
        delta = (x - x1) * div;
        s += delta;
    }
    long double val = s / b;
    f[m] = val;
    for (i = m; i > 0; i--)
    {
        b -= 1.l;
        e1 /= lower2;
        val = (e - e1 + t * val) / b;
        f[i - 1] = val;
    }
}

void fmt_lerfc_like(long double *f, long double t, long double lower, int m)
{
    if (lower == 0)
    {
        return lgamma_inc_like(f, t, m);
    }

    int i;
    long double lower2 = lower * lower;
    // F[m] < .5*sqrt(pi/t) * erfc(low*tt)
    if (t * lower2 > 200.l)
    {
        for (i = 0; i <= m; i++)
        {
            f[i] = 0;
        }
        return;
    }

    if (t < constants::TURNOVER_POINT[m])
    {
        fmt1_lerfc_like(f, t, lower, m);
    }
    else
    {
        long double tt = sqrtl(t);
        // erfc(a) - erfc(b) is more accurate than erf(b) - erf(a)
        long double val = constants::SQRTPI_4L / tt * (erfcl(lower * tt) - erfcl(tt));
        f[0] = val;
        if (m > 0)
        {
            long double e = expl(-t);
            long double e1 = expl(-t * lower2) * lower;
            long double b = .5l / t;
            for (i = 0; i < m; i++)
            {
                val = b * ((2 * i + 1) * val - e + e1);
                e1 *= lower2;
                f[i + 1] = val;
            }
        }
    }
}

int lrys_schmidt(int nroots, double x, double lower, double *roots, double *weights)
{
    int i, k, j, order, error;
    int nroots1 = nroots + 1;
    long double fmt_ints[32 * 2 + 32 * 32];
    long double *qcs = fmt_ints + nroots1 * 2;
    double rt[32 + 32 * 32];
    double *cs = rt + nroots;
    double *a;
    double root, poly, dum, dum0;

    if (lower == 0)
    {
        lgamma_inc_like(fmt_ints, x, nroots * 2);
    }
    else
    {
        fmt_lerfc_like(fmt_ints, x, lower, nroots * 2);
    }

    if (fmt_ints[0] == 0)
    {
        for (k = 0; k < nroots; ++k)
        {
            roots[k] = 0;
            weights[k] = 0;
        }
        return 0;
    }

    if (nroots == 1)
    {
        rt[0] = fmt_ints[1] / fmt_ints[0];
    }
    else
    {
        error = R_lsmit(qcs, fmt_ints, nroots1);
        if (error)
        {
            return error;
        }
        for (k = 1; k < nroots1; k++)
        {
            for (i = 0; i <= k; i++)
            {
                cs[k * nroots1 + i] = qcs[k * nroots1 + i];
            }
        }

        error = polynomial_roots(rt, cs, nroots);
        if (error)
        {
            return error;
        }
    }

    dum0 = 1 / fmt_ints[0];
    for (k = 0; k < nroots; ++k)
    {
        root = rt[k];
        if (root == 1)
        {
            roots[k] = 0;
            weights[k] = 0;
            continue;
        }

        dum = dum0;
        for (j = 1; j < nroots; ++j)
        {
            order = j;
            a = cs + j * nroots1;
            // poly = poly_value1(cs[:order+1,j], order, root);
            POLYNOMIAL_VALUE1(poly, a, order, root);
            dum += poly * poly;
        }
        roots[k] = root / (1 - root);
        weights[k] = 1 / dum;
    }
    return 0;
}

int R_dsmit(double *cs, double *fmt_ints, int n)
{
    int i, j, k;
    double fac, dot, tmp;
    double v[32];

    fac = -fmt_ints[1] / fmt_ints[0];
    tmp = fmt_ints[2] + fac * fmt_ints[1];
    if (tmp <= 0)
    {
        fprintf(stderr, "libcint::rys_roots negative value in sqrt for roots %d (j=1)\n", n - 1);
        SET_ZERO(cs, n, 1);
        return 1;
    }
    tmp = 1 / sqrt(tmp);
    cs[0 + 0 * n] = 1 / sqrt(fmt_ints[0]);
    cs[0 + 1 * n] = fac * tmp;
    cs[1 + 1 * n] = tmp;

    for (j = 2; j < n; ++j)
    {
        for (k = 0; k < j; ++k)
        {
            v[k] = 0;
        }
        fac = fmt_ints[j + j];
        for (k = 0; k < j; ++k)
        {
            dot = 0;
            for (i = 0; i <= k; ++i)
            {
                dot += cs[i + k * n] * fmt_ints[i + j];
            }
            for (i = 0; i <= k; ++i)
            {
                v[i] -= dot * cs[i + k * n];
            }
            fac -= dot * dot;
        }

        if (fac <= 0)
        {
            // set rest coefficients to 0
            SET_ZERO(cs, n, j);
            if (fac == 0)
            {
                return 0;
            }
            fprintf(stderr, "libcint::rys_roots negative value in sqrt for roots %d (j=%d)\n", n - 1, j);
            return j;
        }
        fac = 1 / sqrt(fac);
        cs[j + j * n] = fac;
        for (k = 0; k < j; ++k)
        {
            cs[k + j * n] = fac * v[k];
        }
    }
    return 0;
}

int _rdk_rys_roots(int nroots, double *fmt_ints,
                   double *roots, double *weights)
{
    int i, k, j, order;
    int nroots1 = nroots + 1;
    double rt[32 + 32 * 32];
    double *cs = rt + nroots1;
    double *a;
    double root, poly, dum;

    // to avoid numerical instability for very small fmt integrals
    if (fmt_ints[0] == 0)
    {
        for (k = 0; k < nroots; ++k)
        {
            roots[k] = 0;
            weights[k] = 0;
        }
        return 0;
    }
    if (nroots == 1)
    {
        roots[0] = fmt_ints[1] / (fmt_ints[0] - fmt_ints[1]);
        weights[0] = fmt_ints[0];
        return 0;
    }

    int error = R_dsmit(cs, fmt_ints, nroots1);
    if (error)
    {
        return 1;
    }
    error = polynomial_roots(rt, cs, nroots);
    if (error)
    {
        return error;
    }

    for (k = 0; k < nroots; ++k)
    {
        root = rt[k];
        // When singularity was caught in R_dsmit, they are typically
        // caused by high order Rys polynomials. We assume the contributions
        // from these high order Rys polynomials are negligible. Only the
        // roots obtained from lower polynomials are used.
        if (root == 1)
        {
            roots[k] = 0;
            weights[k] = 0;
            continue;
        }

        dum = 1 / fmt_ints[0];
        for (j = 1; j < nroots; ++j)
        {
            order = j;
            a = cs + j * nroots1;
            // poly = poly_value1(cs[:order+1,j], order, root);
            POLYNOMIAL_VALUE1(poly, a, order, root);
            dum += poly * poly;
        }
        roots[k] = root / (1 - root);
        weights[k] = 1 / dum;
    }
    return 0;
}

int R_lsmit(long double *cs, long double *fmt_ints, int n)
{
    int i, j, k;
    long double fac, dot, tmp;
    long double v[32];

    fac = -fmt_ints[1] / fmt_ints[0];
    tmp = fmt_ints[2] + fac * fmt_ints[1];
    if (tmp <= 0)
    {
        fprintf(stderr, "libcint::rys_roots negative value in sqrtl for roots %d (j=1)\n", n - 1);
        SET_ZERO(cs, n, 1);
        return 1;
    }
    tmp = 1 / sqrtl(tmp);
    cs[0 + 0 * n] = 1 / sqrtl(fmt_ints[0]);
    cs[0 + 1 * n] = fac * tmp;
    cs[1 + 1 * n] = tmp;

    for (j = 2; j < n; ++j)
    {
        for (k = 0; k < j; ++k)
        {
            v[k] = 0;
        }
        fac = fmt_ints[j + j];
        for (k = 0; k < j; ++k)
        {
            dot = 0;
            for (i = 0; i <= k; ++i)
            {
                dot += cs[i + k * n] * fmt_ints[i + j];
            }
            for (i = 0; i <= k; ++i)
            {
                v[i] -= dot * cs[i + k * n];
            }
            fac -= dot * dot;
        }

        if (fac <= 0)
        {
            // set rest coefficients to 0
            SET_ZERO(cs, n, j);
            if (fac == 0)
            {
                return 0;
            }
            fprintf(stderr, "libcint::rys_roots negative value in sqrtl for roots %d (j=%d)\n", n - 1, j);
            return j;
        }
        fac = 1 / sqrtl(fac);
        cs[j + j * n] = fac;
        for (k = 0; k < j; ++k)
        {
            cs[k + j * n] = fac * v[k];
        }
    }
    return 0;
}

void fmt1_gamma_inc_like(double *f, double t, int m)
{
    int i;
    double b = m + 0.5;
    double bi;
    double e = .5 * exp(-t);
    double x = e;
    double s = e;
    double tol = .5 * 2.2204460492503131e-016 * e;
    for (bi = b + 1.; x > tol; bi += 1.)
    {
        x *= t / bi;
        s += x;
    }
    f[m] = s / b;
    for (i = m; i > 0; i--)
    {
        b -= 1.;
        f[i - 1] = (e + t * f[i]) / b;
    }
}

void gamma_inc_like(double *f, double t, int m)
{
    if (t < constants::TURNOVER_POINT[m])
    {
        fmt1_gamma_inc_like(f, t, m);
    }
    else
    {
        int i;
        double tt = sqrt(t);
        f[0] = constants::SQRTPI_4L / tt * erf(tt);
        if (m > 0)
        {
            double e = exp(-t);
            double b = .5 / t;
            for (i = 1; i <= m; i++)
                f[i] = b * ((2 * i - 1) * f[i - 1] - e);
        }
    }
}

void fmt1_erfc_like(double *f, double t, double lower, int m)
{
    int i;
    double lower2 = lower * lower;
    double b = m + 0.5;
    double bi;
    double e = .5 * exp(-t);
    double e1 = .5 * exp(-t * lower2) * lower;
    e1 *= pow(lower2, m);
    double x = e;
    double x1 = e1;
    double s = e - e1;
    double div = 1.;
    double delta = s;
    double tol = .5 * 2.2204460492503131e-016 * fabs(delta);
    for (bi = b + 1.; fabs(delta) > tol; bi += 1.)
    {
        div *= t / bi;
        x1 *= lower2;
        delta = (x - x1) * div;
        s += delta;
    }
    double val = s / b;
    f[m] = val;
    for (i = m; i > 0; i--)
    {
        b -= 1.;
        e1 /= lower2;
        val = (e - e1 + t * val) / b;
        f[i - 1] = val;
    }
}

void fmt_erfc_like(double *f, double t, double lower, int m)
{
    if (lower == 0)
    {
        return gamma_inc_like(f, t, m);
    }

    int i;
    double lower2 = lower * lower;
    // F[m] < .5*sqrt(pi/t) * erfc(low*tt)
    if (t * lower2 > 200)
    {
        for (i = 0; i <= m; i++)
        {
            f[i] = 0;
        }
        return;
    }

    if (t < constants::TURNOVER_POINT[m])
    {
        fmt1_erfc_like(f, t, lower, m);
    }
    else
    {
        double tt = sqrt(t);
        // erfc(a) - erfc(b) is more accurate than erf(b) - erf(a)
        double val = constants::SQRTPI_4L / tt * (erfc(lower * tt) - erfc(tt));
        f[0] = val;
        if (m > 0)
        {
            double e = exp(-t);
            double e1 = exp(-t * lower2) * lower;
            double b = .5 / t;
            for (i = 0; i < m; i++)
            {
                val = b * ((2 * i + 1) * val - e + e1);
                e1 *= lower2;
                f[i + 1] = val;
            }
        }
    }
}

int rys_schmidt(int nroots, double x, double lower, double *roots, double *weights)
{
    double fmt_ints[32 * 2];
    if (lower == 0)
    {
        gamma_inc_like(fmt_ints, x, nroots * 2);
    }
    else
    {
        fmt_erfc_like(fmt_ints, x, lower, nroots * 2);
    }
    return _rdk_rys_roots(nroots, fmt_ints, roots, weights);
}

void wheeler_recursion(int n, const double *alpha, const double *beta, double *moments,
                       double *a, double *b)
{
    int i, j, nc;
    double a0 = alpha[0] + moments[1] / moments[0];
    double b0 = 0;
    double a1, b1;
    a[0] = a0;
    b[0] = b0;
    double buf[32 * 4];
    double *s0 = moments;
    double *sm = buf;
    double *sk = buf + n * 2;
    double *swap;
    for (i = 2; i < n * 2; i++)
    {
        sm[i] = 0.;
    }

    for (i = 1; i < n; i++)
    {
        nc = 2 * (n - i);
        for (j = 0; j < nc; j++)
        {
            sk[j] = (s0[2 + j] - (a0 - alpha[i + j]) * s0[1 + j] - b0 * sm[2 + j] + beta[i + j] * s0[j]);
        }
        a1 = alpha[i] - s0[1] / s0[0] + sk[1] / sk[0];
        b1 = sk[0] / s0[0];
        a[i] = a1;
        b[i] = b1;
        a0 = a1;
        b0 = b1;
        swap = sm;
        sm = s0;
        s0 = sk;
        sk = swap;
    }
}

int _dlaev2(double *eig, double *vec, double *diag, double *diag_off1)
{
    double a = diag[0];
    double b = diag_off1[0];
    double c = diag[1];
    double df, cs, ct, tb, sm, tn, rt, tmp;
    double rt1, rt2, cs1, sn1;
    int sgn1, sgn2;

    sm = a + c;
    df = a - c;
    tb = b + b;

    rt = sqrt(tb * tb + df * df);

    if (sm > 0.)
    {
        rt1 = (sm + rt) * .5;
        sgn1 = 1;
        rt2 = (a * c - b * b) / rt1;
    }
    else if (sm < 0.)
    {
        rt1 = (sm - rt) * .5;
        sgn1 = -1;
        rt2 = (a * c - b * b) / rt1;
    }
    else
    {
        rt1 = rt * .5;
        rt2 = rt * -.5;
        sgn1 = 1;
    }

    /*     Compute the eigenvector */

    if (df >= 0.)
    {
        cs = df + rt;
        sgn2 = 1;
    }
    else
    {
        cs = df - rt;
        sgn2 = -1;
    }

    if (fabs(cs) > fabs(tb))
    {
        ct = -tb / cs;
        sn1 = 1. / sqrt(ct * ct + 1.);
        cs1 = ct * sn1;
    }
    else
    {
        if (b == 0.)
        {
            cs1 = 1.;
            sn1 = 0.;
        }
        else
        {
            tn = -cs / tb;
            cs1 = 1. / sqrt(tn * tn + 1.);
            sn1 = tn * cs1;
        }
    }

    if (sgn1 == sgn2)
    {
        tmp = cs1;
        cs1 = -sn1;
        sn1 = tmp;
    }

    eig[0] = rt2;
    eig[1] = rt1;
    vec[0] = -sn1;
    vec[1] = cs1;
    vec[2] = cs1;
    vec[3] = sn1;
    return 0;
}

int diagonalize(int n, double *diag, double *diag_off1, double *eig, double *vec)
{
    int M = 0;
    int ISUPPZ[32 * 2];
    int TRYRAC = 1;
    /*
        int matrix_layout, char jobz, char range,
        lapack_int n, double* d, double* e, double vl,
        double vu, lapack_int il, lapack_int iu,
        lapack_int* m, double* w, double* z, lapack_int ldz,
        lapack_int nzc, lapack_int* isuppz,
        lapack_logical* tryrac
        */
#if has_RAS == 1
    int INFO = LAPACKE_dstemr(LAPACK_ROW_MAJOR, 'V', 'A', n, diag, diag_off1, 0.0, 0.0, 0, 0, &M, // WE HAVE TO CHECK IF IT IS ACTUALLY ROW MAJOR!!!!
                              eig, vec, n, n, ISUPPZ, &TRYRAC);
    return INFO;
#endif
    return 1;
}

int rys_wheeler_partial(int n, const double *alpha, const double *beta, double *moments,
                        double *roots, double *weights)
{
    double a[32 + 32 + 32 * 32]; // scartch array for LAPACKE
    double *b = a + n;           // Offset by maximum number of roots in scratch
    double *c0 = b + n;          // offset by another maximum number of roots in scratch with size for square matrix
    double mu0 = moments[0];
    int first_seen = 1;
    int i;

    wheeler_recursion(n, alpha, beta, moments, a, b);

    for (i = 1; i < n; i++)
    {
        if (b[i] < 1e-14)
        {
            // very likely we will get numerical issues
            if (!first_seen || b[i] < 0.)
            {
                fprintf(stderr, "libcint rys_wheeler singular value n=%d i=%d b=%g\n",
                        n, i, b[i]);
                for (i = 0; i < n; i++)
                {
                    roots[i] = 0;
                    weights[i] = 0;
                }
                return i;
            }
            first_seen = 0;
        }
        b[i] = sqrt(b[i]);
    }

    int error = diagonalize(n, a, b + 1, roots, c0);

    for (i = 0; i < n; i++)
    {
        roots[i] = roots[i] / (1 - roots[i]);
        weights[i] = c0[i * n] * c0[i * n] * mu0;
    }
    return error;
}

void naive_jacobi_moments(int n, double t, double lower, double *mus)
{
    int i, j, k;
    double s;
    double fmt[32 * 2];
    const double *coef;
    const int *order;

    fmt_erfc_like(fmt, t, lower, n - 1);

    for (i = 0; i < n; i++)
    {
        coef = constants::JACOBI_COEF + i * (i + 1) / 2;
        order = constants::JACOBI_COEF_ORDER + i * (i + 1) / 2;
        s = 0;
        for (j = 0; j <= i; j++)
        {
            k = order[j];
            s += coef[k] * fmt[k];
        }
        mus[i] = s;
    }
}

// Flocke's recipe JCP, 131, 064107
void flocke_jacobi_moments(int n, double t, double *mus)
{
    if (t < 3e-7)
    {
        return naive_jacobi_moments(n, t, 0., mus);
    }

    double t_inv = .5 / t;
    double mu1 = 1.; // DBL_EPSILON; // can be arbitrary number != 0
    double mu2 = 0.;
    double mu0 = 0, rn = 0;
    int i = 0;

    // Miller algorithm
    for (i = n - 1 + 20; i >= n; i--)
    {
        rn = (2 * i + 3) * t_inv + constants::JACOBI_RN_PART2[i];
        mu0 = (mu2 - rn * mu1) / constants::JACOBI_SN[i];
        mu2 = mu1;
        mu1 = mu0;
    }
    for (; i >= 0; i--)
    {
        rn = (2 * i + 3) * t_inv + constants::JACOBI_RN_PART2[i];
        mu0 = (mu2 - rn * mu1) / constants::JACOBI_SN[i];
        mus[i] = mu0;
        mu2 = mu1;
        mu1 = mu0;
    }

    double tt = sqrt(t);
    double norm = constants::SQRTPI_4L * erf(tt) / tt / mu0; // fmt[0]/mu0
    for (i = 0; i < n; i++)
    {
        mus[i] *= norm;
    }
}

// Flocke's recipe JCP, 131, 064107
void lflocke_jacobi_moments(int n, double t, long double *mus)
{
    if (t < 3e-7)
    {
        lnaive_jacobi_moments(n, t, 0., mus);
    }

    long double t_inv = .5l / t;
    long double mu1 = 1.l; // DBL_EPSILON; // can be arbitrary number != 0
    long double mu2 = 0.l;
    long double mu0, rn;
    int i;

    // Miller algorithm
    for (i = n - 1 + 24; i >= n; i--)
    {
        rn = (2 * i + 3) * t_inv + constants::lJACOBI_RN_PART2[i];
        mu0 = (mu2 - rn * mu1) / constants::lJACOBI_SN[i];
        mu2 = mu1;
        mu1 = mu0;
    }
    for (; i >= 0; i--)
    {
        rn = (2 * i + 3) * t_inv + constants::lJACOBI_RN_PART2[i];
        mu0 = (mu2 - rn * mu1) / constants::lJACOBI_SN[i];
        mus[i] = mu0;
        mu2 = mu1;
        mu1 = mu0;
    }

    long double tt = sqrtl(t);
    long double norm = constants::SQRTPI_4L * erfl(tt) / tt / mu0; // fmt[0]/mu0
    for (i = 0; i < n; i++)
    {
        mus[i] *= norm;
    }
}

int rys_jacobi(int n, double x, double lower, double *roots, double *weights)
{
    double moments[32 * 2];
    const double *alpha = constants::JACOBI_ALPHA;
    const double *beta = constants::JACOBI_BETA;

    if (lower == 0)
    {
        flocke_jacobi_moments(n * 2, x, moments);
    }
    else
    {
        naive_jacobi_moments(n * 2, x, lower, moments);
    }
    return rys_wheeler_partial(n, alpha, beta, moments, roots, weights);
}

void rys_roots(int nroots, double x, double *u, double *w)
{
    if (x <= 3e-7)
    {
        int off = nroots * (nroots - 1) / 2;
        int i;
        for (i = 0; i < nroots; i++)
        {
            u[i] = constants::POLY_SMALLX_R0[off + i] + constants::POLY_SMALLX_R1[off + i] * x;
            w[i] = constants::POLY_SMALLX_W0[off + i] + constants::POLY_SMALLX_W1[off + i] * x;
        }
        return;
    }
    else if (x >= 35 + nroots * 5)
    {
        int off = nroots * (nroots - 1) / 2;
        int i;
        double rt;
        double t = sqrt(constants::PI_4 / x);
        for (i = 0; i < nroots; i++)
        {
            rt = constants::POLY_LARGEX_RT[off + i];
            u[i] = rt / (x - rt);
            w[i] = constants::POLY_LARGEX_WW[off + i] * t;
        }
        return;
    }

    int err;
    switch (nroots)
    {
    case 1:
        err = rys_root1(x, u, w);
        break;
    case 2:
        err = rys_root2(x, u, w);
        break;
    case 3:
        err = rys_root3(x, u, w);
        break;
    case 4:
        err = rys_root4(x, u, w);
        break;
    case 5:
        err = rys_root5(x, u, w);
        break;
    case 6:
    case 7:
        err = segment_solve(nroots, x, 0., u, w, 11, rys_jacobi, rys_schmidt);
        break;
    case 8:
        err = segment_solve(nroots, x, 0., u, w, 11, rys_jacobi, lrys_schmidt);
        break;
    case 9:
        err = segment_solve(nroots, x, 0., u, w, 10, lrys_jacobi, lrys_laguerre);
        break;
    case 10:
    case 11:
        err = segment_solve(nroots, x, 0., u, w, 18, lrys_jacobi, lrys_laguerre);
        break;
    case 12:
        err = segment_solve(nroots, x, 0., u, w, 22, lrys_jacobi, lrys_laguerre);
        break;
    default:
        err_not_impl_f("rys_roots: nroots > 8", std::cout);
        err = 0;
        break;
    }
    err_checkf(err == 0, "rys_roots: nroots problem!", std::cout);
}

int segment_solve1(int n, double x, double lower, double *u, double *w,
                   double lower_bp1, double lower_bp2, double breakpoint,
                   QuadratureFunction fn1, QuadratureFunction fn2, QuadratureFunction fn3)
{
    int error;
    if (lower < lower_bp1)
    {
        if (x <= breakpoint)
        {
            error = fn1(n, x, lower, u, w);
        }
        else
        {
            error = fn2(n, x, lower, u, w);
        }
    }
    else if (lower < lower_bp2)
    {
        error = fn3(n, x, lower, u, w);
    }
    else
    {
        return 1;
    }
    if (error)
    {
        error = lrys_schmidt(n, x, lower, u, w);
    }
    return error;
}

void lwheeler_recursion(int n, const long double *alpha, const long double *beta, long double *moments,
                        long double *a, long double *b)
{
    int i, j, nc;
    long double a0 = alpha[0] + moments[1] / moments[0];
    long double b0 = 0;
    long double a1, b1;
    a[0] = a0;
    b[0] = b0;
    long double buf[32 * 4];
    long double *s0 = moments;
    long double *sm = buf;
    long double *sk = buf + n * 2;
    long double *swap;
    for (i = 2; i < n * 2; i++)
    {
        sm[i] = 0.;
    }

    for (i = 1; i < n; i++)
    {
        nc = 2 * (n - i);
        for (j = 0; j < nc; j++)
        {
            sk[j] = (s0[2 + j] - (a0 - alpha[i + j]) * s0[1 + j] - b0 * sm[2 + j] + beta[i + j] * s0[j]);
        }
        a1 = alpha[i] - s0[1] / s0[0] + sk[1] / sk[0];
        b1 = sk[0] / s0[0];
        a[i] = a1;
        b[i] = b1;
        a0 = a1;
        b0 = b1;
        swap = sm;
        sm = s0;
        s0 = sk;
        sk = swap;
    }
}

int lrys_wheeler_partial(int n, const long double *alpha, const long double *beta, long double *moments,
                         double *roots, double *weights)
{
    long double a[32 + 32];
    long double *b = a + n;
    double da[32 + 32 + 32 * 32];
    double *db = da + n;
    double *c0 = db + n;
    double mu0 = moments[0];
    int first_seen = 1;
    int i;

    lwheeler_recursion(n, alpha, beta, moments, a, b);

    da[0] = a[0];
    for (i = 1; i < n; i++)
    {
        if (b[i] < 1e-19)
        {
            // very likely we will get numerical issues
            if (!first_seen || b[i] < 0.)
            {
                fprintf(stderr, "libcint lrys_wheeler singular value n=%d i=%d b=%g\n",
                        n, i, (double)b[i]);
                for (i = 0; i < n; i++)
                {
                    roots[i] = 0;
                    weights[i] = 0;
                }
                return i;
            }
            first_seen = 0;
        }
        da[i] = a[i];
        db[i] = sqrtl(b[i]);
    }

    int error = diagonalize(n, da, db + 1, roots, c0);

    for (i = 0; i < n; i++)
    {
        roots[i] = roots[i] / (1 - roots[i]);
        weights[i] = c0[i * n] * c0[i * n] * mu0;
    }
    return error;
}

void lnaive_jacobi_moments(int n, double t, double lower, long double *mus)
{
    int i, j, k;
    long double s;
    long double fmt[32 * 2];
    const long double *coef;
    const int *order;

    fmt_lerfc_like(fmt, t, lower, n - 1);

    for (i = 0; i < n; i++)
    {
        coef = constants::lJACOBI_COEF + i * (i + 1) / 2;
        order = constants::JACOBI_COEF_ORDER + i * (i + 1) / 2;
        s = 0;
        for (j = 0; j <= i; j++)
        {
            k = order[j];
            s += coef[k] * fmt[k];
        }
        mus[i] = s;
    }
}

int lrys_jacobi(int n, double x, double lower, double *roots, double *weights)
{
    long double moments[32 * 2];
    const long double *alpha = constants::lJACOBI_ALPHA;
    const long double *beta = constants::lJACOBI_BETA;

    if (lower == 0)
    {
        lflocke_jacobi_moments(n * 2, x, moments);
    }
    else
    {
        lnaive_jacobi_moments(n * 2, x, lower, moments);
    }
    return lrys_wheeler_partial(n, alpha, beta, moments, roots, weights);
}

void llaguerre_moments(int n, double t, double lower,
                       long double *alpha, long double *beta, long double *moments)
{
    int i;
    long double tt = sqrtl(t);
    long double t_inv = .5l / t;
    long double t2_inv = .5l / (t * t);
    long double e0 = expl(-t) * t_inv;
    long double l00 = 0.l;
    long double l01 = 1.l;
    long double fac0, fac1, l02;

    alpha[0] = t_inv;
    beta[0] = 0;
    if (lower == 0)
    {
        moments[0] = constants::SQRTPI_4L / tt * erfl(tt);
        moments[1] = -l01 * e0;
        for (i = 1; i < n - 1; i++)
        {
            alpha[i] = (i * 4 + 1) * t_inv;
            beta[i] = i * (i * 2 - 1) * t2_inv;
            fac0 = (i * 4 - 1) * t_inv;
            fac1 = (i - 1) * (i * 2 - 1) * t2_inv;
            l02 = (1. - fac0) * l01 - fac1 * l00;
            l00 = l01;
            l01 = l02;
            moments[i + 1] = -l01 * e0;
        }
    }
    else
    {
        long double lower2 = lower * lower;
        long double l10 = 0.l;
        long double l11 = 1.l;
        long double l12;
        long double et = expl(-t * lower2) * lower * t_inv;
        moments[0] = constants::SQRTPI_4L / tt * (erfcl(lower * tt) - erfcl(tt));
        moments[1] = l11 * et - l01 * e0;
        for (i = 1; i < n - 1; i++)
        {
            alpha[i] = (i * 4 + 1) * t_inv;
            beta[i] = i * (i * 2 - 1) * t2_inv;
            fac0 = (i * 4 - 1) * t_inv;
            fac1 = (i - 1) * (i * 2 - 1) * t2_inv;
            l12 = (lower2 - fac0) * l11 - fac1 * l10;
            l10 = l11;
            l11 = l12;
            l02 = (1. - fac0) * l01 - fac1 * l00;
            l00 = l01;
            l01 = l02;
            moments[i + 1] = l11 * et - l01 * e0;
        }
    }
}

int lrys_laguerre(int n, double x, double lower, double *roots, double *weights)
{
    long double moments[32 * 6];
    long double *alpha = moments + n * 2;
    long double *beta = alpha + n * 2;

    llaguerre_moments(n * 2, x, lower, alpha, beta, moments);

    return lrys_wheeler_partial(n, alpha, beta, moments, roots, weights);
}

void sr_rys_roots(int nroots, double x, double lower, double *u, double *w)
{
    int err = 1;
    switch (nroots)
    {
    case 1:
        err = rys_schmidt(nroots, x, lower, u, w);
        break;
    case 2:
        if (lower < 0.99)
        {
            err = rys_schmidt(nroots, x, lower, u, w);
        }
        else
        {
            err = lrys_jacobi(nroots, x, lower, u, w);
        }
        break;
    case 3:
        if (lower < 0.93)
        {
            err = rys_schmidt(nroots, x, lower, u, w);
        }
        else if (lower < 0.97)
        {
            err = segment_solve(nroots, x, lower, u, w, 10, lrys_jacobi, lrys_laguerre);
        }
        else
        {
            err = lrys_jacobi(nroots, x, lower, u, w);
        }
        break;
    case 4:
        if (lower < 0.8)
        {
            err = rys_schmidt(nroots, x, lower, u, w);
        }
        else if (lower < 0.9)
        {
            err = segment_solve(nroots, x, lower, u, w, 10, lrys_jacobi, lrys_laguerre);
        }
        else
        {
            err = lrys_jacobi(nroots, x, lower, u, w);
        }
        break;
    case 5:
        if (lower < 0.4)
        {
            err = segment_solve(nroots, x, lower, u, w, 50, rys_schmidt, lrys_laguerre);
        }
        else if (lower < 0.8)
        {
            err = segment_solve(nroots, x, lower, u, w, 10, lrys_jacobi, lrys_laguerre);
        }
        else
        {
            err = lrys_jacobi(nroots, x, lower, u, w);
        }
        break;
    case 6:
        if (lower < 0.25)
        {
            err = segment_solve(nroots, x, lower, u, w, 60, rys_schmidt, lrys_laguerre);
        }
        else if (lower < 0.8)
        {
            err = segment_solve(nroots, x, lower, u, w, 10, lrys_jacobi, lrys_laguerre);
        }
        else
        {
            err = lrys_jacobi(nroots, x, lower, u, w);
        }
        break;
    case 7:
        err = segment_solve1(nroots, x, lower, u, w, 0.5, 1., 60, lrys_jacobi, lrys_laguerre, lrys_jacobi);
        break;
    case 8:
    case 9:
    case 10:
        // CINTqrys_jacobi(nroots, x, lower, u, w);
        err = segment_solve1(nroots, x, lower, u, w, 0.15, 1., 60, lrys_jacobi, lrys_laguerre, lrys_jacobi);
        break;
    case 11:
    case 12:
        err = segment_solve1(nroots, x, lower, u, w, 0.15, 1., 60, lrys_jacobi, lrys_laguerre, lrys_jacobi);
        break;
    case 13:
    case 14:
        err = segment_solve1(nroots, x, lower, u, w, 0.25, 1., 60, lrys_jacobi, lrys_laguerre, lrys_jacobi);
        break;
    case 15:
    case 16:
        err = segment_solve1(nroots, x, lower, u, w, 0.25, 0.75, 60, lrys_jacobi, lrys_laguerre, lrys_jacobi);
        break;
    case 17:
        err = segment_solve1(nroots, x, lower, u, w, 0.25, 0.65, 60, lrys_jacobi, lrys_laguerre, lrys_jacobi);
        break;
    case 18:
        err = segment_solve1(nroots, x, lower, u, w, 0.15, 0.65, 60, lrys_jacobi, lrys_laguerre, lrys_jacobi);
        break;
    case 19:
        err = segment_solve1(nroots, x, lower, u, w, 0.15, 0.55, 60, lrys_jacobi, lrys_laguerre, lrys_jacobi);
        break;
    case 20:
    case 21:
        err = segment_solve1(nroots, x, lower, u, w, 0.25, 0.45, 60, lrys_jacobi, lrys_laguerre, lrys_jacobi);
        break;
    case 22:
    case 23:
    case 24:
        err = segment_solve1(nroots, x, lower, u, w, 0.25, 0.35, 60, lrys_jacobi, lrys_laguerre, lrys_jacobi);
        break;
    default:
        fprintf(stderr, "libcint SR-rys_roots does not support nroots=%d\n", nroots);
    }
    if (err)
    {
        fprintf(stderr, "sr_rys_roots fails: nroots=%d x=%.15g lower=%.15g\n",
                nroots, x, lower);
    }
}

int g0_2e(double *g, double *rij, double *rkl, double cutoff, Env *envs)
{
    int irys;
    int nroots = envs->nrys_roots;
    double aij = envs->ai[0] + envs->aj[0];
    double akl = envs->ak[0] + envs->al[0];
    double a0, a1, fac1, x;
    double u[32];
    double *w = g + envs->g_size * 2; // ~ gz
    double xij_kl = rij[0] - rkl[0];
    double yij_kl = rij[1] - rkl[1];
    double zij_kl = rij[2] - rkl[2];
    double rr = xij_kl * xij_kl + yij_kl * yij_kl + zij_kl * zij_kl;

    a1 = aij * akl;
    a0 = a1 / (aij + akl);
    fac1 = sqrt(a0 / (a1 * a1 * a1)) * envs->fac[0];
    x = a0 * rr;
    const double omega = envs->env[8];
    double theta = 0;
    if (omega == 0.)
    {
        rys_roots(nroots, x, u, w);
    }
    else if (omega < 0.)
    {
        // short-range part of range-separated Coulomb
        theta = omega * omega / (omega * omega + a0);
        // very small erfc() leads to ~0 weights. They can cause
        // numerical issue in sr_rys_roots.
        if (theta * x > cutoff || theta * x > 40)
        {
            return 0;
        }
        int rorder = envs->rys_order;
        if (rorder == nroots)
        {
            sr_rys_roots(nroots, x, sqrt(theta), u, w);
        }
        else
        {
            double sqrt_theta = -sqrt(theta);
            rys_roots(rorder, x, u, w);
            rys_roots(rorder, theta * x, u + rorder, w + rorder);
            if (envs->g_size == 2)
            {
                g[0] = 1;
                g[1] = 1;
                g[2] = 1;
                g[3] = 1;
                g[4] *= fac1;
                g[5] *= fac1 * sqrt_theta;
                return 1;
            }
            for (irys = rorder; irys < nroots; irys++)
            {
                double ut = u[irys] * theta;
                u[irys] = ut / (u[irys] + 1. - ut);
                w[irys] *= sqrt_theta;
            }
        }
    }
    else
    { // omega > 0.
        // long-range part of range-separated Coulomb
        theta = omega * omega / (omega * omega + a0);
        x *= theta;
        fac1 *= sqrt(theta);
        rys_roots(nroots, x, u, w);
        /* u[:] = tau^2 / (1 - tau^2)
         * omega^2u^2 = a0 * tau^2 / (theta^-1 - tau^2)
         * transform u[:] to theta^-1 tau^2 / (theta^-1 - tau^2)
         * so the rest code can be reused.
         */
        for (irys = 0; irys < nroots; irys++)
        {
            double ut = u[irys] * theta;
            u[irys] = ut / (u[irys] + 1. - ut);
        }
    }
    if (envs->g_size == 1)
    {
        g[0] = 1;
        g[1] = 1;
        g[2] *= fac1;
        return 1;
    }

    double u2, tmp1, tmp2, tmp3, tmp4, tmp5;
    double rijrx = rij[0] - envs->rx_in_rijrx[0];
    double rijry = rij[1] - envs->rx_in_rijrx[1];
    double rijrz = rij[2] - envs->rx_in_rijrx[2];
    double rklrx = rkl[0] - envs->rx_in_rklrx[0];
    double rklry = rkl[1] - envs->rx_in_rklrx[1];
    double rklrz = rkl[2] - envs->rx_in_rklrx[2];
    Rys2eT bc;
    double *b00 = bc.b00;
    double *b10 = bc.b10;
    double *b01 = bc.b01;
    double *c00x = bc.c00x;
    double *c00y = bc.c00y;
    double *c00z = bc.c00z;
    double *c0px = bc.c0px;
    double *c0py = bc.c0py;
    double *c0pz = bc.c0pz;

    for (irys = 0; irys < nroots; irys++)
    {
        /*
         *u(irys) = t2/(1-t2)
         *t2 = u(irys)/(1+u(irys))
         *u2 = aij*akl/(aij+akl)*t2/(1-t2)
         */
        u2 = a0 * u[irys];
        tmp4 = .5 / (u2 * (aij + akl) + a1);
        tmp5 = u2 * tmp4;
        tmp1 = 2. * tmp5;
        tmp2 = tmp1 * akl;
        tmp3 = tmp1 * aij;
        b00[irys] = tmp5;
        b10[irys] = tmp5 + tmp4 * akl;
        b01[irys] = tmp5 + tmp4 * aij;
        c00x[irys] = rijrx - tmp2 * xij_kl;
        c00y[irys] = rijry - tmp2 * yij_kl;
        c00z[irys] = rijrz - tmp2 * zij_kl;
        c0px[irys] = rklrx + tmp3 * xij_kl;
        c0py[irys] = rklry + tmp3 * yij_kl;
        c0pz[irys] = rklrz + tmp3 * zij_kl;
        w[irys] *= fac1;
    }

    (*envs->f_g0_2d4d)(g, &bc, envs);

    return 1;
}

void init_int2c2e_Env(Env *envs, int *ng, int *shls, int *atm, int natm, int *bas, int nbas, double *env)
{
    envs->natm = natm;
    envs->nbas = nbas;
    envs->atm = atm;
    envs->bas = bas;
    envs->env = env;
    envs->shls = shls;

    const int i_sh = shls[0];
    const int k_sh = shls[1];
    envs->i_l = bas(1, i_sh);
    envs->j_l = 0;
    envs->k_l = bas(1, k_sh);
    envs->l_l = 0;
    envs->x_ctr[0] = bas(3, i_sh);
    envs->x_ctr[1] = bas(3, k_sh);
    envs->x_ctr[2] = 1;
    envs->x_ctr[3] = 1;
    envs->nfi = (envs->i_l + 1) * (envs->i_l + 2) / 2;
    envs->nfj = 1;
    envs->nfk = (envs->k_l + 1) * (envs->k_l + 2) / 2;
    envs->nfl = 1;
    envs->nf = envs->nfi * envs->nfk;

    envs->ri = env + atm(1, bas(0, i_sh));
    envs->rk = env + atm(1, bas(0, k_sh));

    envs->common_factor = (constants::PI3) * 2 / constants::sqr_pi * common_fac_sp(envs->i_l) * common_fac_sp(envs->k_l);
    if (env[0] == 0)
    {
        envs->expcutoff = 60.;
    }
    else
    {
        envs->expcutoff = std::max(40., env[0]);
    }

    envs->gbits = ng[4];
    envs->ncomp_e1 = ng[5];
    envs->ncomp_e2 = ng[6];
    envs->ncomp_tensor = ng[7];

    envs->li_ceil = envs->i_l + ng[0];
    envs->lj_ceil = 0;
    envs->lk_ceil = envs->k_l + ng[2];
    envs->ll_ceil = 0;
    int rys_order = (envs->li_ceil + envs->lk_ceil) / 2 + 1;
    int nrys_roots = rys_order;
    double omega = env[8];
    if (omega < 0 && rys_order <= 3)
    {
        nrys_roots *= 2;
    }
    envs->rys_order = rys_order;
    envs->nrys_roots = nrys_roots;

    int dli = envs->li_ceil + 1;
    int dlk = envs->lk_ceil + 1;
    envs->g_stride_i = nrys_roots;
    envs->g_stride_k = nrys_roots * dli;
    envs->g_stride_l = envs->g_stride_k;
    envs->g_size = nrys_roots * dli * dlk;

    envs->aj[0] = 0;
    envs->al[0] = 0;
    envs->rij[0] = envs->ri[0];
    envs->rij[1] = envs->ri[1];
    envs->rij[2] = envs->ri[2];
    envs->rkl[0] = envs->rk[0];
    envs->rkl[1] = envs->rk[1];
    envs->rkl[2] = envs->rk[2];
    envs->g2d_ijmax = envs->g_stride_i;
    envs->g2d_klmax = envs->g_stride_k;
    envs->rkrl[0] = envs->rk[0];
    envs->rkrl[1] = envs->rk[1];
    envs->rkrl[2] = envs->rk[2];
    envs->rirj[0] = envs->ri[0];
    envs->rirj[1] = envs->ri[1];
    envs->rirj[2] = envs->ri[2];
    envs->rx_in_rklrx = envs->rk;
    envs->rx_in_rijrx = envs->ri;

    if (rys_order <= 2)
    {
        envs->f_g0_2d4d = &g0_2e_2d4d_unrolled;
        if (rys_order != nrys_roots)
        {
            envs->f_g0_2d4d = &srg0_2e_2d4d_unrolled;
        }
    }
    else
    {
        envs->f_g0_2d4d = &g0_2e_2d;
    }
    envs->f_g0_2e = &g0_2e;

    // initialize j_l, j_ctr, nfj because they are used in c2s_sph_1e and
    // CINTg1e_index_xyz
    envs->j_l = envs->k_l;
    envs->nfj = envs->nfk;
    envs->g_stride_j = envs->g_stride_k;
}

void gout2e(double *gout, double *g, int *idx,
            Env *envs, int gout_empty)
{
    int nf = envs->nf;
    int i, ix, iy, iz, n;
    double s;

    if (gout_empty)
    {
        switch (envs->nrys_roots)
        {
        case 1:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] = g[ix] * g[iy] * g[iz];
            }
            break;
        case 2:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] = g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1];
            }
            break;
        case 3:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] = g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1] + g[ix + 2] * g[iy + 2] * g[iz + 2];
            }
            break;
        case 4:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] = g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1] + g[ix + 2] * g[iy + 2] * g[iz + 2] + g[ix + 3] * g[iy + 3] * g[iz + 3];
            }
            break;
        case 5:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] = g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1] + g[ix + 2] * g[iy + 2] * g[iz + 2] + g[ix + 3] * g[iy + 3] * g[iz + 3] + g[ix + 4] * g[iy + 4] * g[iz + 4];
            }
            break;
        case 6:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] = g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1] + g[ix + 2] * g[iy + 2] * g[iz + 2] + g[ix + 3] * g[iy + 3] * g[iz + 3] + g[ix + 4] * g[iy + 4] * g[iz + 4] + g[ix + 5] * g[iy + 5] * g[iz + 5];
            }
            break;
        case 7:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] = g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1] + g[ix + 2] * g[iy + 2] * g[iz + 2] + g[ix + 3] * g[iy + 3] * g[iz + 3] + g[ix + 4] * g[iy + 4] * g[iz + 4] + g[ix + 5] * g[iy + 5] * g[iz + 5] + g[ix + 6] * g[iy + 6] * g[iz + 6];
            }
            break;
        case 8:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] = g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1] + g[ix + 2] * g[iy + 2] * g[iz + 2] + g[ix + 3] * g[iy + 3] * g[iz + 3] + g[ix + 4] * g[iy + 4] * g[iz + 4] + g[ix + 5] * g[iy + 5] * g[iz + 5] + g[ix + 6] * g[iy + 6] * g[iz + 6] + g[ix + 7] * g[iy + 7] * g[iz + 7];
            }
            break;
        default:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                s = 0;
                for (i = 0; i < envs->nrys_roots; i++)
                {
                    s += g[ix + i] * g[iy + i] * g[iz + i];
                }
                gout[n] = s;
            }
            break;
        } // end switch nroots
    }
    else
    { // not flag_acc
        switch (envs->nrys_roots)
        {
        case 1:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] += g[ix] * g[iy] * g[iz];
            }
            break;
        case 2:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] += g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1];
            }
            break;
        case 3:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] += g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1] + g[ix + 2] * g[iy + 2] * g[iz + 2];
            }
            break;
        case 4:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] += g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1] + g[ix + 2] * g[iy + 2] * g[iz + 2] + g[ix + 3] * g[iy + 3] * g[iz + 3];
            }
            break;
        case 5:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] += g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1] + g[ix + 2] * g[iy + 2] * g[iz + 2] + g[ix + 3] * g[iy + 3] * g[iz + 3] + g[ix + 4] * g[iy + 4] * g[iz + 4];
            }
            break;
        case 6:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] += g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1] + g[ix + 2] * g[iy + 2] * g[iz + 2] + g[ix + 3] * g[iy + 3] * g[iz + 3] + g[ix + 4] * g[iy + 4] * g[iz + 4] + g[ix + 5] * g[iy + 5] * g[iz + 5];
            }
            break;
        case 7:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] += g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1] + g[ix + 2] * g[iy + 2] * g[iz + 2] + g[ix + 3] * g[iy + 3] * g[iz + 3] + g[ix + 4] * g[iy + 4] * g[iz + 4] + g[ix + 5] * g[iy + 5] * g[iz + 5] + g[ix + 6] * g[iy + 6] * g[iz + 6];
            }
            break;
        case 8:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] += g[ix] * g[iy] * g[iz] + g[ix + 1] * g[iy + 1] * g[iz + 1] + g[ix + 2] * g[iy + 2] * g[iz + 2] + g[ix + 3] * g[iy + 3] * g[iz + 3] + g[ix + 4] * g[iy + 4] * g[iz + 4] + g[ix + 5] * g[iy + 5] * g[iz + 5] + g[ix + 6] * g[iy + 6] * g[iz + 6] + g[ix + 7] * g[iy + 7] * g[iz + 7];
            }
            break;
        default:
            for (n = 0; n < nf; n++, idx += 3)
            {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                s = 0;
                for (i = 0; i < envs->nrys_roots; i++)
                {
                    s += g[ix + i] * g[iy + i] * g[iz + i];
                }
                gout[n] += s;
            }
            break;
        } // end switch nroots
    }
}

template <typename T>
inline void MALLOC_INSTACK(T*& var, const int& n, double*& cache)
{
    var = (T*)(((intptr_t)cache + 7) & (-(intptr_t)8));
    cache = (double*)(var + (n));
}

int int1e_cache_size(Env *envs)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    int i_prim = bas(2, shls[0]);
    int j_prim = bas(2, shls[1]);
    int *x_ctr = envs->x_ctr;
    int nc = envs->nf * x_ctr[0] * x_ctr[1];
    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
    int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
    int lenj = envs->nf * nc * n_comp;
    int leni = envs->nf * x_ctr[0] * n_comp;
    int len0 = envs->nf * n_comp;
    int pdata_size = (i_prim * j_prim * 5 + i_prim * x_ctr[0] + j_prim * x_ctr[1] + (i_prim + j_prim) * 2 + envs->nf * 3);
    int cache_size = std::max(nc * n_comp + leng + lenj + leni + len0 + pdata_size,
                              nc * n_comp + envs->nf * 8 * 2);
    return cache_size;
}

#define gctrg gout
#define gctrm gctr
#define mempty empty
#define m_ctr n_comp
#define ALIAS_ADDR_IF_EQUAL(x, y) \
    if (y##_ctr == 1)             \
    {                             \
        gctr##x = gctr##y;        \
        x##empty = y##empty;      \
    }                             \
    else                          \
    {                             \
        gctr##x = g1;             \
        g1 += len##x;             \
    }

void prim_to_ctr_0(double *gc, double *gp, double *coeff, size_t nf,
                   int nprim, int nctr, int non0ctr, int *sortedidx)
{
    int i;
    size_t n;
    double c0;

    for (i = 0; i < nctr; i++)
    {
        c0 = coeff[nprim * i];
        for (n = 0; n < nf; n++)
        {
            gc[nf * i + n] = c0 * gp[n];
        }
    }
}

void prim_to_ctr_1(double *gc, double *gp, double *coeff, size_t nf,
                   int nprim, int nctr, int non0ctr, int *sortedidx)
{
    int i, j;
    size_t n;
    double c0;

    for (i = 0; i < non0ctr; i++)
    {
        c0 = coeff[nprim * sortedidx[i]];
        j = sortedidx[i];
        for (n = 0; n < nf; n++)
        {
            gc[nf * j + n] += c0 * gp[n];
        }
    }
}

#define PRIM2CTR(ctrsymb, gp, ngp)                                        \
    if (ctrsymb##_ctr > 1)                                                \
    {                                                                     \
        if (*ctrsymb##empty)                                              \
        {                                                                 \
            prim_to_ctr_0(gctr##ctrsymb, gp, c##ctrsymb + ctrsymb##p,     \
                          ngp, ctrsymb##_prim, ctrsymb##_ctr,             \
                          non0ctr##ctrsymb[ctrsymb##p],                   \
                          non0idx##ctrsymb + ctrsymb##p * ctrsymb##_ctr); \
        }                                                                 \
        else                                                              \
        {                                                                 \
            prim_to_ctr_1(gctr##ctrsymb, gp, c##ctrsymb + ctrsymb##p,     \
                          ngp, ctrsymb##_prim, ctrsymb##_ctr,             \
                          non0ctr##ctrsymb[ctrsymb##p],                   \
                          non0idx##ctrsymb + ctrsymb##p * ctrsymb##_ctr); \
        }                                                                 \
    }                                                                     \
    *ctrsymb##empty = 0

void dmat_transpose(double *a_t, double *a, int m, int n)
{
    int i, j;

    for (j = 0; j < n - 3; j += 4)
    {
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[(j + 0) * m + i] = a[i * n + j + 0];
            a_t[(j + 1) * m + i] = a[i * n + j + 1];
            a_t[(j + 2) * m + i] = a[i * n + j + 2];
            a_t[(j + 3) * m + i] = a[i * n + j + 3];
        }
    }

    switch (n - j)
    {
    case 1:
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[j * m + i] = a[i * n + j];
        }
        break;
    case 2:
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[(j + 0) * m + i] = a[i * n + j + 0];
            a_t[(j + 1) * m + i] = a[i * n + j + 1];
        }
        break;
    case 3:
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[(j + 0) * m + i] = a[i * n + j + 0];
            a_t[(j + 1) * m + i] = a[i * n + j + 1];
            a_t[(j + 2) * m + i] = a[i * n + j + 2];
        }
        break;
    }
}

/*
 * a_t[n,m] += a[m,n]
 */
void dplus_transpose(double *a_t, double *a, int m, int n)
{
    int i, j;

    for (j = 0; j < n - 3; j += 4)
    {
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[(j + 0) * m + i] += a[i * n + j + 0];
            a_t[(j + 1) * m + i] += a[i * n + j + 1];
            a_t[(j + 2) * m + i] += a[i * n + j + 2];
            a_t[(j + 3) * m + i] += a[i * n + j + 3];
        }
    }

    switch (n - j)
    {
    case 1:
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[j * m + i] += a[i * n + j];
        }
        break;
    case 2:
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[(j + 0) * m + i] += a[i * n + j + 0];
            a_t[(j + 1) * m + i] += a[i * n + j + 1];
        }
        break;
    case 3:
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[(j + 0) * m + i] += a[i * n + j + 0];
            a_t[(j + 1) * m + i] += a[i * n + j + 1];
            a_t[(j + 2) * m + i] += a[i * n + j + 2];
        }
        break;
    }
}

#define TRANSPOSE(a)                              \
    if (*empty)                                   \
    {                                             \
        dmat_transpose(gctr, a, nf *nc, n_comp);  \
    }                                             \
    else                                          \
    {                                             \
        dplus_transpose(gctr, a, nf *nc, n_comp); \
    }                                             \
    *empty = 0;

/*
 * GTO = x^{nx}y^{ny}z^{nz}e^{-ar^2}
 */
void cart_comp(int *nx, int *ny, int *nz, const int lmax)
{
    int inc = 0;
    int lx, ly, lz;

    for (lx = lmax; lx >= 0; lx--)
    {
        for (ly = lmax - lx; ly >= 0; ly--)
        {
            lz = lmax - lx - ly;
            nx[inc] = lx;
            ny[inc] = ly;
            nz[inc] = lz;
            inc++;
        }
    }
}

void g1e_index_xyz(int *idx, const Env *envs)
{
    const int i_l = envs->i_l;
    const int j_l = envs->j_l;
    const int nfi = envs->nfi;
    const int nfj = envs->nfj;
    const int di = envs->g_stride_i;
    const int dj = envs->g_stride_j;
    int i, j, n;
    int ofx, ofjx;
    int ofy, ofjy;
    int ofz, ofjz;
    int i_nx[136], i_ny[136], i_nz[136];
    int j_nx[136], j_ny[136], j_nz[136];

    cart_comp(i_nx, i_ny, i_nz, i_l);
    cart_comp(j_nx, j_ny, j_nz, j_l);

    ofx = 0;
    ofy = envs->g_size;
    ofz = envs->g_size * 2;
    n = 0;
    for (j = 0; j < nfj; j++)
    {
        ofjx = ofx + dj * j_nx[j];
        ofjy = ofy + dj * j_ny[j];
        ofjz = ofz + dj * j_nz[j];
        for (i = 0; i < nfi; i++)
        {
            idx[n + 0] = ofjx + di * i_nx[i];
            idx[n + 1] = ofjy + di * i_ny[i];
            idx[n + 2] = ofjz + di * i_nz[i];
            n += 3;
        }
    }
}

void Opt_non0coeff_byshell(int *sortedidx, int *non0ctr, double *ci,
                           int iprim, int ictr)
{
    int ip, j, k, kp;
    int *zeroidx = new int[ictr];
    for (ip = 0; ip < iprim; ip++)
    {
        for (j = 0, k = 0, kp = 0; j < ictr; j++)
        {
            if (ci[iprim * j + ip] != 0)
            {
                sortedidx[k] = j;
                k++;
            }
            else
            {
                zeroidx[kp] = j;
                kp++;
            }
        }
        // Append the index of zero-coeff to sortedidx for function prim_to_ctr_0
        for (j = 0; j < kp; j++)
        {
            sortedidx[k + j] = zeroidx[j];
        }
        non0ctr[ip] = k;
        sortedidx += ictr;
    }
    free(zeroidx);
}

int _2c2e_loop_nopt(double *gctr, Env *envs, double *cache, int *empty)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    double *env = envs->env;
    int i_sh = shls[0];
    int k_sh = shls[1];
    int i_ctr = envs->x_ctr[0];
    int k_ctr = envs->x_ctr[1];
    int i_prim = bas(2, i_sh);
    int k_prim = bas(2, k_sh);
    double *ai = env + bas(5, i_sh);
    double *ak = env + bas(5, k_sh);
    double *ci = env + bas(6, i_sh);
    double *ck = env + bas(6, k_sh);
    double expcutoff = envs->expcutoff;
    double *ri = envs->ri;
    double *rk = envs->rk;
    int n_comp = envs->ncomp_tensor;
    double fac1i, fac1k;
    int ip, kp;
    int _empty[3] = {1, 1, 1};
    int *iempty = _empty + 0;
    int *kempty = _empty + 1;
    int *gempty = _empty + 2;
    int nf = envs->nf;
    const int nc = i_ctr * k_ctr;
    const int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
    const int lenk = nf * nc * n_comp;    // gctrk
    const int leni = nf * i_ctr * n_comp; // gctri
    const int len0 = nf * n_comp;         // gout
    const int len = leng + lenk + leni + len0;
    double *g = NULL;
    MALLOC_INSTACK(g, len, cache);
    double *g1 = g + leng;
    double *gout, *gctri, *gctrk;

    ALIAS_ADDR_IF_EQUAL(k, m);
    ALIAS_ADDR_IF_EQUAL(i, k);
    ALIAS_ADDR_IF_EQUAL(g, i);

    int *idx = NULL;
    MALLOC_INSTACK(idx, envs->nf * 3, cache);
    g1e_index_xyz(idx, envs);

    int *non0ctri = NULL, *non0ctrk = NULL;
    int *non0idxi = NULL, *non0idxk = NULL;
    MALLOC_INSTACK(non0ctri, i_prim + k_prim + i_prim * i_ctr + k_prim * k_ctr, cache);
    non0ctrk = non0ctri + i_prim;
    non0idxi = non0ctrk + k_prim;
    non0idxk = non0idxi + i_prim * i_ctr;
    if (i_ctr > 1)
    {
        Opt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
    }
    if (k_ctr > 1)
    {
        Opt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
    }

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        envs->al[0] = 0; // to use CINTg0_2e
        if (k_ctr == 1)
        {
            fac1k = envs->common_factor * ck[kp];
        }
        else
        {
            fac1k = envs->common_factor;
            *iempty = 1;
        }
        for (ip = 0; ip < i_prim; ip++)
        {
            envs->ai[0] = ai[ip];
            envs->aj[0] = 0;
            if (i_ctr == 1)
            {
                fac1i = fac1k * ci[ip];
            }
            else
            {
                fac1i = fac1k;
            }
            envs->fac[0] = fac1i;
            if ((*envs->f_g0_2e)(g, ri, rk, expcutoff, envs))
            {
                (*envs->f_gout)(gout, g, idx, envs, *gempty);
                PRIM2CTR(i, gout, len0);
            }
        } // end loop i_prim
        if (!*iempty)
        {
            PRIM2CTR(k, gctri, leni);
        }
    } // end loop k_prim

    if (n_comp > 1 && !*kempty)
    {
        TRANSPOSE(gctrk);
    }
    return !*empty;
}

int _2c2e_loop(double *gctr, Env *envs, double *cache, int *empty)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    double *env = envs->env;
    int i_sh = shls[0];
    int k_sh = shls[1];
    int i_ctr = envs->x_ctr[0];
    int k_ctr = envs->x_ctr[1];
    int i_prim = bas(2, i_sh);
    int k_prim = bas(2, k_sh);
    double *ai = env + bas(5, i_sh);
    double *ak = env + bas(5, k_sh);
    double *ci = env + bas(6, i_sh);
    double *ck = env + bas(6, k_sh);
    double expcutoff = envs->expcutoff;
    double *ri = envs->ri;
    double *rk = envs->rk;
    int n_comp = envs->ncomp_tensor;
    double fac1i, fac1k;
    int ip, kp;
    int _empty[3] = {1, 1, 1};
    int *iempty = _empty + 0;
    int *kempty = _empty + 1;
    int *gempty = _empty + 2;
    int *non0ctri = NULL, *non0ctrk = NULL;
    int *non0idxi = NULL, *non0idxk = NULL;
    MALLOC_INSTACK(non0ctri, i_prim + k_prim + i_prim * i_ctr + k_prim * k_ctr, cache);
    non0ctrk = non0ctri + i_prim;
    non0idxi = non0ctrk + k_prim;
    non0idxk = non0idxi + i_prim * i_ctr;
    if (i_ctr > 1)
    {
        Opt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
    }
    if (k_ctr > 1)
    {
        Opt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
    }

    int *idx = envs->opt->index_xyz_array[envs->i_l * 16 + envs->k_l];

    int nf = envs->nf;
    const int nc = i_ctr * k_ctr;
    const int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
    const int lenk = nf * nc * n_comp;    // gctrk
    const int leni = nf * i_ctr * n_comp; // gctri
    const int len0 = nf * n_comp;         // gout
    const int len = leng + lenk + leni + len0;
    double *g = NULL;
    MALLOC_INSTACK(g, len, cache);
    double *g1 = g + leng;
    double *gout, *gctri, *gctrk;

    ALIAS_ADDR_IF_EQUAL(k, m);
    ALIAS_ADDR_IF_EQUAL(i, k);
    ALIAS_ADDR_IF_EQUAL(g, i);

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        if (k_ctr == 1)
        {
            fac1k = envs->common_factor * ck[kp];
        }
        else
        {
            fac1k = envs->common_factor;
            *iempty = 1;
        }
        for (ip = 0; ip < i_prim; ip++)
        {
            envs->ai[0] = ai[ip];
            if (i_ctr == 1)
            {
                fac1i = fac1k * ci[ip];
            }
            else
            {
                fac1i = fac1k;
            }
            envs->fac[0] = fac1i;
            if ((*envs->f_g0_2e)(g, ri, rk, expcutoff, envs))
            {
                (*envs->f_gout)(gout, g, idx, envs, *gempty);
                PRIM2CTR(i, gout, len0);
            }
        } // end loop i_prim
        if (!*iempty)
        {
            PRIM2CTR(k, gctri, leni);
        }
    } // end loop k_prim

    if (n_comp > 1 && !*kempty)
    {
        TRANSPOSE(gctrk);
    }
    return !*empty;
}

static int _len_cart[] = {
    1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120, 136};

struct cart2sp_t
{
    double *cart2sph;
    double *cart2j_lt_lR; // j < kappa, l > 0
    double *cart2j_lt_lI; // j < kappa, l > 0
    double *cart2j_gt_lR; // j > kappa, l < 0
    double *cart2j_gt_lI; // j > kappa, l < 0
};

void g2e_index_xyz(int *idx, const Env *envs)
{
    const int i_l = envs->i_l;
    const int j_l = envs->j_l;
    const int k_l = envs->k_l;
    const int l_l = envs->l_l;
    const int nfi = envs->nfi;
    const int nfj = envs->nfj;
    const int nfk = envs->nfk;
    const int nfl = envs->nfl;
    const int di = envs->g_stride_i;
    const int dk = envs->g_stride_k;
    const int dl = envs->g_stride_l;
    const int dj = envs->g_stride_j;
    int i, j, k, l, n;
    int ofx, ofkx, oflx;
    int ofy, ofky, ofly;
    int ofz, ofkz, oflz;
    int i_nx[136], i_ny[136], i_nz[136];
    int j_nx[136], j_ny[136], j_nz[136];
    int k_nx[136], k_ny[136], k_nz[136];
    int l_nx[136], l_ny[136], l_nz[136];

    cart_comp(i_nx, i_ny, i_nz, i_l);
    cart_comp(j_nx, j_ny, j_nz, j_l);
    cart_comp(k_nx, k_ny, k_nz, k_l);
    cart_comp(l_nx, l_ny, l_nz, l_l);

    ofx = 0;
    ofy = envs->g_size;
    ofz = envs->g_size * 2;
    n = 0;
    for (j = 0; j < nfj; j++)
    {
        for (l = 0; l < nfl; l++)
        {
            oflx = ofx + dj * j_nx[j] + dl * l_nx[l];
            ofly = ofy + dj * j_ny[j] + dl * l_ny[l];
            oflz = ofz + dj * j_nz[j] + dl * l_nz[l];
            for (k = 0; k < nfk; k++)
            {
                ofkx = oflx + dk * k_nx[k];
                ofky = ofly + dk * k_ny[k];
                ofkz = oflz + dk * k_nz[k];
                switch (i_l)
                {
                case 0:
                    idx[n + 0] = ofkx;
                    idx[n + 1] = ofky;
                    idx[n + 2] = ofkz;
                    n += 3;
                    break;
                case 1:
                    idx[n + 0] = ofkx + di;
                    idx[n + 1] = ofky;
                    idx[n + 2] = ofkz;
                    idx[n + 3] = ofkx;
                    idx[n + 4] = ofky + di;
                    idx[n + 5] = ofkz;
                    idx[n + 6] = ofkx;
                    idx[n + 7] = ofky;
                    idx[n + 8] = ofkz + di;
                    n += 9;
                    break;
                case 2:
                    idx[n + 0] = ofkx + di * 2;
                    idx[n + 1] = ofky;
                    idx[n + 2] = ofkz;
                    idx[n + 3] = ofkx + di;
                    idx[n + 4] = ofky + di;
                    idx[n + 5] = ofkz;
                    idx[n + 6] = ofkx + di;
                    idx[n + 7] = ofky;
                    idx[n + 8] = ofkz + di;
                    idx[n + 9] = ofkx;
                    idx[n + 10] = ofky + di * 2;
                    idx[n + 11] = ofkz;
                    idx[n + 12] = ofkx;
                    idx[n + 13] = ofky + di;
                    idx[n + 14] = ofkz + di;
                    idx[n + 15] = ofkx;
                    idx[n + 16] = ofky;
                    idx[n + 17] = ofkz + di * 2;
                    n += 18;
                    break;
                default:
                    for (i = 0; i < nfi; i++)
                    {
                        idx[n + 0] = ofkx + di * i_nx[i]; //(:,ix,kx,lx,jx,1)
                        idx[n + 1] = ofky + di * i_ny[i]; //(:,iy,ky,ly,jy,2)
                        idx[n + 2] = ofkz + di * i_nz[i]; //(:,iz,kz,lz,jz,3)
                        n += 3;
                    } // i
                }
            } // k
        } // l
    } // j
}

double g_trans_cart2sph[] = {
    1, /* factors of s and p are moved to CINTcommon_fac_sp */
       // px
#ifdef PYPZPX
    // py
    0,
    1,
    0,
    // pz
    0,
    0,
    1,
    // px
    1,
    0,
    0,
#else
    // by default, p orbitals are ordered px, py, pz
    // px
    1,
    0,
    0,
    // py
    0,
    1,
    0,
    // pz
    0,
    0,
    1,
#endif
    // dxy
    0,
    1.092548430592079070,
    0,
    0,
    0,
    0,
    // dyz
    0,
    0,
    0,
    0,
    1.092548430592079070,
    0,
    // dz2
    -0.315391565252520002,
    0,
    0,
    -0.315391565252520002,
    0,
    0.630783130505040012,
    // dxz
    0,
    0,
    1.092548430592079070,
    0,
    0,
    0,
    // dy2
    0.546274215296039535,
    0,
    0,
    -0.546274215296039535,
    0,
    0,
    // f-3 ~ fyx2
    0,
    1.770130769779930531,
    0,
    0,
    0,
    0,
    -0.590043589926643510,
    0,
    0,
    0,
    // f-2 ~ fxyz
    0,
    0,
    0,
    0,
    2.890611442640554055,
    0,
    0,
    0,
    0,
    0,
    // f-1 ~ fyz2
    0,
    -0.457045799464465739,
    0,
    0,
    0,
    0,
    -0.457045799464465739,
    0,
    1.828183197857862944,
    0,
    // f0 ~ fz3
    0,
    0,
    -1.119528997770346170,
    0,
    0,
    0,
    0,
    -1.119528997770346170,
    0,
    0.746352665180230782,
    // f1 ~ fxz2
    -0.457045799464465739,
    0,
    0,
    -0.457045799464465739,
    0,
    1.828183197857862944,
    0,
    0,
    0,
    0,
    // f2 ~ fzx2
    0,
    0,
    1.445305721320277020,
    0,
    0,
    0,
    0,
    -1.445305721320277020,
    0,
    0,
    // f3 ~ fx3
    0.590043589926643510,
    0,
    0,
    -1.770130769779930530,
    0,
    0,
    0,
    0,
    0,
    0,
    // g-4 ~ gyx3
    0,
    2.503342941796704538,
    0,
    0,
    0,
    0,
    -2.503342941796704530,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // g-3 ~ gx2yz
    0,
    0,
    0,
    0,
    5.310392309339791593,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.770130769779930530,
    0,
    0,
    0,
    // g-2 ~ gxyz2
    0,
    -0.946174695757560014,
    0,
    0,
    0,
    0,
    -0.946174695757560014,
    0,
    5.677048174545360108,
    0,
    0,
    0,
    0,
    0,
    0,
    // g-1 ~ gyz3
    0,
    0,
    0,
    0,
    -2.007139630671867500,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.007139630671867500,
    0,
    2.676186174229156671,
    0,
    // g0 ~ gz4
    0.317356640745612911,
    0,
    0,
    0.634713281491225822,
    0,
    -2.538853125964903290,
    0,
    0,
    0,
    0,
    0.317356640745612911,
    0,
    -2.538853125964903290,
    0,
    0.846284375321634430,
    // g1 ~ gxz3
    0,
    0,
    -2.007139630671867500,
    0,
    0,
    0,
    0,
    -2.007139630671867500,
    0,
    2.676186174229156671,
    0,
    0,
    0,
    0,
    0,
    // g2 ~ gx2z2
    -0.473087347878780002,
    0,
    0,
    0,
    0,
    2.838524087272680054,
    0,
    0,
    0,
    0,
    0.473087347878780009,
    0,
    -2.838524087272680050,
    0,
    0,
    // g3 ~ gzx3
    0,
    0,
    1.770130769779930531,
    0,
    0,
    0,
    0,
    -5.310392309339791590,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // g4 ~ gy4
    0.625835735449176134,
    0,
    0,
    -3.755014412695056800,
    0,
    0,
    0,
    0,
    0,
    0,
    0.625835735449176134,
    0,
    0,
    0,
    0,
    // h-5 ~ hyx4
    0,
    3.281910284200850514,
    0,
    0,
    0,
    0,
    -6.563820568401701020,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.656382056840170102,
    0,
    0,
    0,
    0,
    0,
    // h-4 ~ hx3yz
    0,
    0,
    0,
    0,
    8.302649259524165115,
    0,
    0,
    0,
    0,
    0,
    0,
    -8.302649259524165110,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // h-3 ~ hyx2z2
    0,
    -1.467714898305751160,
    0,
    0,
    0,
    0,
    -0.978476598870500779,
    0,
    11.741719186446009300,
    0,
    0,
    0,
    0,
    0,
    0,
    0.489238299435250387,
    0,
    -3.913906395482003100,
    0,
    0,
    0,
    // h-2 ~ hxyz3
    0,
    0,
    0,
    0,
    -4.793536784973323750,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.793536784973323750,
    0,
    9.587073569946647510,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // h-1 ~ hyz4
    0,
    0.452946651195696921,
    0,
    0,
    0,
    0,
    0.905893302391393842,
    0,
    -5.435359814348363050,
    0,
    0,
    0,
    0,
    0,
    0,
    0.452946651195696921,
    0,
    -5.435359814348363050,
    0,
    3.623573209565575370,
    0,
    // h0 ~ hx2y2z
    0,
    0,
    1.754254836801353946,
    0,
    0,
    0,
    0,
    3.508509673602707893,
    0,
    -4.678012898136943850,
    0,
    0,
    0,
    0,
    0,
    0,
    1.754254836801353946,
    0,
    -4.678012898136943850,
    0,
    0.935602579627388771,
    // h1 ~ xz4
    0.452946651195696921,
    0,
    0,
    0.905893302391393842,
    0,
    -5.435359814348363050,
    0,
    0,
    0,
    0,
    0.452946651195696921,
    0,
    -5.435359814348363050,
    0,
    3.623573209565575370,
    0,
    0,
    0,
    0,
    0,
    0,
    // h2 ~ hx2z3
    0,
    0,
    -2.396768392486661870,
    0,
    0,
    0,
    0,
    0,
    0,
    4.793536784973323755,
    0,
    0,
    0,
    0,
    0,
    0,
    2.396768392486661877,
    0,
    -4.793536784973323750,
    0,
    0,
    // h3 ~ hx3z2
    -0.489238299435250389,
    0,
    0,
    0.978476598870500775,
    0,
    3.913906395482003101,
    0,
    0,
    0,
    0,
    1.467714898305751163,
    0,
    -11.741719186446009300,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // h4 ~ hzy4
    0,
    0,
    2.075662314881041278,
    0,
    0,
    0,
    0,
    -12.453973889286247600,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.075662314881041278,
    0,
    0,
    0,
    0,
    // h5 ~ hxy4
    0.656382056840170102,
    0,
    0,
    -6.563820568401701020,
    0,
    0,
    0,
    0,
    0,
    0,
    3.281910284200850514,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // i-6
    0,
    4.0991046311514863,
    0,
    0,
    0,
    0,
    -13.6636821038382887,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.0991046311514863,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // i-5
    0,
    0,
    0,
    0,
    11.8330958111587634,
    0,
    0,
    0,
    0,
    0,
    0,
    -23.6661916223175268,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.3666191622317525,
    0,
    0,
    0,
    0,
    0,
    // i-4
    0,
    -2.0182596029148963,
    0,
    0,
    0,
    0,
    0,
    0,
    20.1825960291489679,
    0,
    0,
    0,
    0,
    0,
    0,
    2.0182596029148963,
    0,
    -20.1825960291489679,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // i-3
    0,
    0,
    0,
    0,
    -8.2908473356343109,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.5272315570895412,
    0,
    22.1089262283581647,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.7636157785447706,
    0,
    -7.3696420761193888,
    0,
    0,
    0,
    // i-2
    0,
    0.9212052595149236,
    0,
    0,
    0,
    0,
    1.8424105190298472,
    0,
    -14.7392841522387776,
    0,
    0,
    0,
    0,
    0,
    0,
    0.9212052595149236,
    0,
    -14.7392841522387776,
    0,
    14.7392841522387776,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // i-1
    0,
    0,
    0,
    0,
    2.9131068125936568,
    0,
    0,
    0,
    0,
    0,
    0,
    5.8262136251873136,
    0,
    -11.6524272503746271,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.9131068125936568,
    0,
    -11.6524272503746271,
    0,
    4.6609709001498505,
    0,
    // i0
    -0.3178460113381421,
    0,
    0,
    -0.9535380340144264,
    0,
    5.7212282040865583,
    0,
    0,
    0,
    0,
    -0.9535380340144264,
    0,
    11.4424564081731166,
    0,
    -7.6283042721154111,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.3178460113381421,
    0,
    5.7212282040865583,
    0,
    -7.6283042721154111,
    0,
    1.0171072362820548,
    // i1
    0,
    0,
    2.9131068125936568,
    0,
    0,
    0,
    0,
    5.8262136251873136,
    0,
    -11.6524272503746271,
    0,
    0,
    0,
    0,
    0,
    0,
    2.9131068125936568,
    0,
    -11.6524272503746271,
    0,
    4.6609709001498505,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // i2
    0.4606026297574618,
    0,
    0,
    0.4606026297574618,
    0,
    -7.3696420761193888,
    0,
    0,
    0,
    0,
    -0.4606026297574618,
    0,
    0,
    0,
    7.3696420761193888,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4606026297574618,
    0,
    7.3696420761193888,
    0,
    -7.3696420761193888,
    0,
    0,
    // i3
    0,
    0,
    -2.7636157785447706,
    0,
    0,
    0,
    0,
    5.5272315570895412,
    0,
    7.3696420761193888,
    0,
    0,
    0,
    0,
    0,
    0,
    8.2908473356343109,
    0,
    -22.1089262283581647,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // i4
    -0.5045649007287241,
    0,
    0,
    2.5228245036436201,
    0,
    5.0456490072872420,
    0,
    0,
    0,
    0,
    2.5228245036436201,
    0,
    -30.2738940437234518,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5045649007287241,
    0,
    5.0456490072872420,
    0,
    0,
    0,
    0,
    // i5
    0,
    0,
    2.3666191622317525,
    0,
    0,
    0,
    0,
    -23.6661916223175268,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    11.8330958111587634,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // i6
    0.6831841051919144,
    0,
    0,
    -10.2477615778787161,
    0,
    0,
    0,
    0,
    0,
    0,
    10.2477615778787161,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.6831841051919144,
    0,
    0,
    0,
    0,
    0,
    0,
    // j-7
    0,
    4.9501391276721742,
    0,
    0,
    0,
    0,
    -24.7506956383608703,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    14.8504173830165218,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.7071627325245963,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // j-6
    0,
    0,
    0,
    0,
    15.8757639708114002,
    0,
    0,
    0,
    0,
    0,
    0,
    -52.9192132360380043,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    15.8757639708114002,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // j-5
    0,
    -2.5945778936013020,
    0,
    0,
    0,
    0,
    2.5945778936013020,
    0,
    31.1349347232156219,
    0,
    0,
    0,
    0,
    0,
    0,
    4.6702402084823440,
    0,
    -62.2698694464312439,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5189155787202604,
    0,
    6.2269869446431247,
    0,
    0,
    0,
    0,
    0,
    // j-4
    0,
    0,
    0,
    0,
    -12.4539738892862495,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    41.5132462976208316,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    12.4539738892862495,
    0,
    -41.5132462976208316,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // j-3
    0,
    1.4081304047606462,
    0,
    0,
    0,
    0,
    2.3468840079344107,
    0,
    -28.1626080952129243,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4693768015868821,
    0,
    -18.7750720634752817,
    0,
    37.5501441269505705,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4693768015868821,
    0,
    9.3875360317376408,
    0,
    -12.5167147089835229,
    0,
    0,
    0,
    // j-2
    0,
    0,
    0,
    0,
    6.6379903866747414,
    0,
    0,
    0,
    0,
    0,
    0,
    13.2759807733494828,
    0,
    -35.4026153955986160,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.6379903866747414,
    0,
    -35.4026153955986160,
    0,
    21.2415692373591725,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // j-1
    0,
    -0.4516580379125866,
    0,
    0,
    0,
    0,
    -1.3549741137377600,
    0,
    10.8397929099020782,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.3549741137377600,
    0,
    21.6795858198041564,
    0,
    -21.6795858198041564,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4516580379125866,
    0,
    10.8397929099020782,
    0,
    -21.6795858198041564,
    0,
    5.7812228852811094,
    0,
    // j0
    0,
    0,
    -2.3899496919201728,
    0,
    0,
    0,
    0,
    -7.1698490757605189,
    0,
    14.3396981515210360,
    0,
    0,
    0,
    0,
    0,
    0,
    -7.1698490757605189,
    0,
    28.6793963030420720,
    0,
    -11.4717585212168292,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.3899496919201728,
    0,
    14.3396981515210360,
    0,
    -11.4717585212168292,
    0,
    1.0925484305920790,
    // j1
    -0.4516580379125866,
    0,
    0,
    -1.3549741137377600,
    0,
    10.8397929099020782,
    0,
    0,
    0,
    0,
    -1.3549741137377600,
    0,
    21.6795858198041564,
    0,
    -21.6795858198041564,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4516580379125866,
    0,
    10.8397929099020782,
    0,
    -21.6795858198041564,
    0,
    5.7812228852811094,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // j2
    0,
    0,
    3.3189951933373707,
    0,
    0,
    0,
    0,
    3.3189951933373707,
    0,
    -17.7013076977993080,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.3189951933373707,
    0,
    0,
    0,
    10.6207846186795862,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.3189951933373707,
    0,
    17.7013076977993080,
    0,
    -10.6207846186795862,
    0,
    0,
    // j3
    0.4693768015868821,
    0,
    0,
    -0.4693768015868821,
    0,
    -9.3875360317376408,
    0,
    0,
    0,
    0,
    -2.3468840079344107,
    0,
    18.7750720634752817,
    0,
    12.5167147089835229,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.4081304047606462,
    0,
    28.1626080952129243,
    0,
    -37.5501441269505705,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // j4
    0,
    0,
    -3.1134934723215624,
    0,
    0,
    0,
    0,
    15.5674673616078110,
    0,
    10.3783115744052079,
    0,
    0,
    0,
    0,
    0,
    0,
    15.5674673616078110,
    0,
    -62.2698694464312439,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.1134934723215624,
    0,
    10.3783115744052079,
    0,
    0,
    0,
    0,
    // j5
    -0.5189155787202604,
    0,
    0,
    4.6702402084823440,
    0,
    6.2269869446431247,
    0,
    0,
    0,
    0,
    2.5945778936013020,
    0,
    -62.2698694464312439,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.5945778936013020,
    0,
    31.1349347232156219,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // j6
    0,
    0,
    2.6459606618019000,
    0,
    0,
    0,
    0,
    -39.6894099270284997,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    39.6894099270284997,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.6459606618019000,
    0,
    0,
    0,
    0,
    0,
    0,
    // j7
    0.7071627325245963,
    0,
    0,
    -14.8504173830165218,
    0,
    0,
    0,
    0,
    0,
    0,
    24.7506956383608703,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.9501391276721742,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // l = 8
    0,
    5.83141328139864,
    0,
    0,
    0,
    0,
    -40.81989296979048,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    40.81989296979048,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.83141328139864,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    20.40994648489524,
    0,
    0,
    0,
    0,
    0,
    0,
    -102.0497324244762,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    61.22983945468572,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.91570664069932,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.193996596357255,
    0,
    0,
    0,
    0,
    7.452658724833595,
    0,
    44.71595234900157,
    0,
    0,
    0,
    0,
    0,
    0,
    7.452658724833595,
    0,
    -149.0531744966719,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.193996596357255,
    0,
    44.71595234900157,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -17.24955311049054,
    0,
    0,
    0,
    0,
    0,
    0,
    17.24955311049054,
    0,
    68.99821244196217,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    31.04919559888297,
    0,
    -137.9964248839243,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.449910622098108,
    0,
    13.79964248839243,
    0,
    0,
    0,
    0,
    0,
    0,
    1.913666099037323,
    0,
    0,
    0,
    0,
    1.913666099037323,
    0,
    -45.92798637689575,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.913666099037323,
    0,
    0,
    0,
    76.54664396149292,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.913666099037323,
    0,
    45.92798637689575,
    0,
    -76.54664396149292,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    11.1173953976599,
    0,
    0,
    0,
    0,
    0,
    0,
    18.52899232943316,
    0,
    -74.11596931773265,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.705798465886632,
    0,
    -49.41064621182176,
    0,
    59.29277545418611,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.705798465886632,
    0,
    24.70532310591088,
    0,
    -19.7642584847287,
    0,
    0,
    0,
    0,
    -0.9123045168698189,
    0,
    0,
    0,
    0,
    -2.736913550609457,
    0,
    27.36913550609457,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.736913550609457,
    0,
    54.73827101218914,
    0,
    -72.98436134958553,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.9123045168698189,
    0,
    27.36913550609457,
    0,
    -72.98436134958553,
    0,
    29.19374453983421,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.8164436064573,
    0,
    0,
    0,
    0,
    0,
    0,
    -11.4493308193719,
    0,
    30.5315488516584,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -11.4493308193719,
    0,
    61.06309770331679,
    0,
    -36.63785862199007,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.8164436064573,
    0,
    30.5315488516584,
    0,
    -36.63785862199007,
    0,
    6.978639737521918,
    0,
    0.3180369672047749,
    0,
    0,
    1.272147868819099,
    0,
    -10.1771829505528,
    0,
    0,
    0,
    0,
    1.908221803228649,
    0,
    -30.53154885165839,
    0,
    30.53154885165839,
    0,
    0,
    0,
    0,
    0,
    0,
    1.272147868819099,
    0,
    -30.53154885165839,
    0,
    61.06309770331677,
    0,
    -16.28349272088447,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.3180369672047749,
    0,
    -10.1771829505528,
    0,
    30.53154885165839,
    0,
    -16.28349272088447,
    0,
    1.16310662292032,
    0,
    0,
    -3.8164436064573,
    0,
    0,
    0,
    0,
    -11.4493308193719,
    0,
    30.5315488516584,
    0,
    0,
    0,
    0,
    0,
    0,
    -11.4493308193719,
    0,
    61.06309770331679,
    0,
    -36.63785862199007,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.8164436064573,
    0,
    30.5315488516584,
    0,
    -36.63785862199007,
    0,
    6.978639737521918,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4561522584349095,
    0,
    0,
    -0.9123045168698189,
    0,
    13.68456775304729,
    0,
    0,
    0,
    0,
    0,
    0,
    13.68456775304729,
    0,
    -36.49218067479276,
    0,
    0,
    0,
    0,
    0,
    0,
    0.9123045168698189,
    0,
    -13.68456775304729,
    0,
    0,
    0,
    14.5968722699171,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4561522584349095,
    0,
    -13.68456775304729,
    0,
    36.49218067479276,
    0,
    -14.5968722699171,
    0,
    0,
    0,
    0,
    3.705798465886632,
    0,
    0,
    0,
    0,
    -3.705798465886632,
    0,
    -24.70532310591088,
    0,
    0,
    0,
    0,
    0,
    0,
    -18.52899232943316,
    0,
    49.41064621182176,
    0,
    19.7642584847287,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -11.1173953976599,
    0,
    74.11596931773265,
    0,
    -59.29277545418611,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4784165247593308,
    0,
    0,
    -1.913666099037323,
    0,
    -11.48199659422394,
    0,
    0,
    0,
    0,
    -4.784165247593307,
    0,
    57.40998297111968,
    0,
    19.13666099037323,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.913666099037323,
    0,
    57.40998297111968,
    0,
    -114.8199659422394,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4784165247593308,
    0,
    -11.48199659422394,
    0,
    19.13666099037323,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.449910622098108,
    0,
    0,
    0,
    0,
    31.04919559888297,
    0,
    13.79964248839243,
    0,
    0,
    0,
    0,
    0,
    0,
    17.24955311049054,
    0,
    -137.9964248839243,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -17.24955311049054,
    0,
    68.99821244196217,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5323327660595425,
    0,
    0,
    7.452658724833595,
    0,
    7.452658724833595,
    0,
    0,
    0,
    0,
    0,
    0,
    -111.7898808725039,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -7.452658724833595,
    0,
    111.7898808725039,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.5323327660595425,
    0,
    -7.452658724833595,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.91570664069932,
    0,
    0,
    0,
    0,
    -61.22983945468572,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    102.0497324244762,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -20.40994648489524,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.72892666017483,
    0,
    0,
    -20.40994648489524,
    0,
    0,
    0,
    0,
    0,
    0,
    51.0248662122381,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -20.40994648489524,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.72892666017483,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // l = 9
    0,
    6.740108566678694,
    0,
    0,
    0,
    0,
    -62.9076799556678,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    94.36151993350171,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -26.96043426671477,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.7489009518531882,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    25.41854119163758,
    0,
    0,
    0,
    0,
    0,
    0,
    -177.9297883414631,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    177.9297883414631,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -25.41854119163758,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.814338369408373,
    0,
    0,
    0,
    0,
    15.25735347763349,
    0,
    61.02941391053396,
    0,
    0,
    0,
    0,
    0,
    0,
    7.628676738816745,
    0,
    -305.1470695526698,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -10.89810962688107,
    0,
    183.0882417316019,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.5449054813440533,
    0,
    -8.718487701504852,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -22.65129549625621,
    0,
    0,
    0,
    0,
    0,
    0,
    52.85302282459782,
    0,
    105.7060456491956,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    52.85302282459782,
    0,
    -352.3534854973187,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -22.65129549625621,
    0,
    105.7060456491956,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.436891395195093,
    0,
    0,
    0,
    0,
    0,
    0,
    -68.23295906546261,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.82329590654626,
    0,
    68.23295906546261,
    0,
    136.4659181309252,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.899026232312149,
    0,
    122.8193263178327,
    0,
    -272.9318362618504,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4873782790390186,
    0,
    -13.64659181309252,
    0,
    27.29318362618504,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    16.31079695491669,
    0,
    0,
    0,
    0,
    0,
    0,
    16.31079695491669,
    0,
    -130.4863756393335,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -16.31079695491669,
    0,
    0,
    0,
    130.4863756393335,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -16.31079695491669,
    0,
    130.4863756393335,
    0,
    -130.4863756393335,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.385125560048583,
    0,
    0,
    0,
    0,
    -3.693668160129556,
    0,
    49.864520161749,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.770251120097167,
    0,
    83.107533602915,
    0,
    -166.21506720583,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    16.621506720583,
    0,
    -110.8100448038867,
    0,
    88.64803584310934,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4617085200161945,
    0,
    -16.621506720583,
    0,
    55.40502240194333,
    0,
    -29.54934528103645,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -8.46325696792098,
    0,
    0,
    0,
    0,
    0,
    0,
    -25.38977090376294,
    0,
    84.63256967920979,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -25.38977090376294,
    0,
    169.2651393584196,
    0,
    -135.4121114867357,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -8.46325696792098,
    0,
    84.63256967920979,
    0,
    -135.4121114867357,
    0,
    38.68917471049591,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.451093112065591,
    0,
    0,
    0,
    0,
    1.804372448262364,
    0,
    -18.04372448262364,
    0,
    0,
    0,
    0,
    0,
    0,
    2.706558672393546,
    0,
    -54.13117344787092,
    0,
    72.17489793049457,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.804372448262364,
    0,
    -54.13117344787092,
    0,
    144.3497958609891,
    0,
    -57.73991834439565,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.451093112065591,
    0,
    -18.04372448262364,
    0,
    72.17489793049457,
    0,
    -57.73991834439565,
    0,
    8.248559763485094,
    0,
    0,
    0,
    3.026024588281776,
    0,
    0,
    0,
    0,
    12.1040983531271,
    0,
    -32.27759560833895,
    0,
    0,
    0,
    0,
    0,
    0,
    18.15614752969066,
    0,
    -96.83278682501685,
    0,
    58.0996720950101,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    12.1040983531271,
    0,
    -96.83278682501685,
    0,
    116.1993441900202,
    0,
    -22.1332084171467,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.026024588281776,
    0,
    -32.27759560833895,
    0,
    58.0996720950101,
    0,
    -22.1332084171467,
    0,
    1.229622689841484,
    0.451093112065591,
    0,
    0,
    1.804372448262364,
    0,
    -18.04372448262364,
    0,
    0,
    0,
    0,
    2.706558672393546,
    0,
    -54.13117344787092,
    0,
    72.17489793049457,
    0,
    0,
    0,
    0,
    0,
    0,
    1.804372448262364,
    0,
    -54.13117344787092,
    0,
    144.3497958609891,
    0,
    -57.73991834439565,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.451093112065591,
    0,
    -18.04372448262364,
    0,
    72.17489793049457,
    0,
    -57.73991834439565,
    0,
    8.248559763485094,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.23162848396049,
    0,
    0,
    0,
    0,
    -8.46325696792098,
    0,
    42.3162848396049,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    42.3162848396049,
    0,
    -67.70605574336784,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    8.46325696792098,
    0,
    -42.3162848396049,
    0,
    0,
    0,
    19.34458735524795,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.23162848396049,
    0,
    -42.3162848396049,
    0,
    67.70605574336784,
    0,
    -19.34458735524795,
    0,
    0,
    -0.4617085200161945,
    0,
    0,
    0,
    0,
    16.621506720583,
    0,
    0,
    0,
    0,
    2.770251120097167,
    0,
    -16.621506720583,
    0,
    -55.40502240194333,
    0,
    0,
    0,
    0,
    0,
    0,
    3.693668160129556,
    0,
    -83.107533602915,
    0,
    110.8100448038867,
    0,
    29.54934528103645,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.385125560048583,
    0,
    -49.864520161749,
    0,
    166.21506720583,
    0,
    -88.64803584310934,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.077699238729173,
    0,
    0,
    0,
    0,
    -16.31079695491669,
    0,
    -32.62159390983339,
    0,
    0,
    0,
    0,
    0,
    0,
    -40.77699238729173,
    0,
    163.1079695491669,
    0,
    32.62159390983339,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -16.31079695491669,
    0,
    163.1079695491669,
    0,
    -195.7295634590003,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.077699238729173,
    0,
    -32.62159390983339,
    0,
    32.62159390983339,
    0,
    0,
    0,
    0,
    0.4873782790390186,
    0,
    0,
    -3.899026232312149,
    0,
    -13.64659181309252,
    0,
    0,
    0,
    0,
    -6.82329590654626,
    0,
    122.8193263178327,
    0,
    27.29318362618504,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    68.23295906546261,
    0,
    -272.9318362618504,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.436891395195093,
    0,
    -68.23295906546261,
    0,
    136.4659181309252,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.775215916042701,
    0,
    0,
    0,
    0,
    52.85302282459782,
    0,
    17.61767427486594,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -264.2651141229891,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -52.85302282459782,
    0,
    264.2651141229891,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.775215916042701,
    0,
    -17.61767427486594,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5449054813440533,
    0,
    0,
    10.89810962688107,
    0,
    8.718487701504852,
    0,
    0,
    0,
    0,
    -7.628676738816745,
    0,
    -183.0882417316019,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -15.25735347763349,
    0,
    305.1470695526698,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.814338369408373,
    0,
    -61.02941391053396,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.177317648954698,
    0,
    0,
    0,
    0,
    -88.96489417073154,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    222.4122354268289,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -88.96489417073154,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.177317648954698,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.7489009518531882,
    0,
    0,
    -26.96043426671477,
    0,
    0,
    0,
    0,
    0,
    0,
    94.36151993350171,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -62.9076799556678,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.740108566678694,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // l = 10
    0,
    7.673951182219901,
    0,
    0,
    0,
    0,
    -92.08741418663881,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    193.3835697919415,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -92.08741418663881,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    7.673951182219901,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    30.88705769902543,
    0,
    0,
    0,
    0,
    0,
    0,
    -288.2792051909041,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    432.4188077863561,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -123.5482307961017,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.431895299891715,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.453815461763347,
    0,
    0,
    0,
    0,
    26.72289277058008,
    0,
    80.16867831174027,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -561.1807481821819,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -26.72289277058008,
    0,
    561.1807481821819,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.453815461763347,
    0,
    -80.16867831174027,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -28.63763513582592,
    0,
    0,
    0,
    0,
    0,
    0,
    114.5505405433037,
    0,
    152.7340540577382,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    57.27527027165184,
    0,
    -763.6702702886912,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -81.82181467378834,
    0,
    458.2021621732147,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.091090733689417,
    0,
    -21.81915057967689,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.976705744527138,
    0,
    0,
    0,
    0,
    -3.968940992702851,
    0,
    -95.25458382486842,
    0,
    0,
    0,
    0,
    0,
    0,
    -13.89129347445998,
    0,
    222.2606955913596,
    0,
    222.2606955913597,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.968940992702851,
    0,
    222.2606955913596,
    0,
    -740.8689853045323,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.976705744527138,
    0,
    -95.25458382486842,
    0,
    222.2606955913597,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    22.18705464592268,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -207.0791766952783,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -62.12375300858349,
    0,
    207.0791766952783,
    0,
    248.495012034334,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -35.49928743347628,
    0,
    372.742518051501,
    0,
    -496.990024068668,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.437410929184535,
    0,
    -41.41583533905566,
    0,
    49.6990024068668,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.870976726712969,
    0,
    0,
    0,
    0,
    -3.741953453425937,
    0,
    78.58102252194469,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    78.58102252194469,
    0,
    -314.3240900877788,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.741953453425937,
    0,
    -78.58102252194469,
    0,
    0,
    0,
    209.5493933918525,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.870976726712969,
    0,
    -78.58102252194469,
    0,
    314.3240900877788,
    0,
    -209.5493933918525,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -13.89129347445998,
    0,
    0,
    0,
    0,
    0,
    0,
    -37.04344926522661,
    0,
    166.6955216935197,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -27.78258694891996,
    0,
    277.8258694891996,
    0,
    -333.3910433870395,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    55.56517389783991,
    0,
    -222.2606955913596,
    0,
    127.0061117664912,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.630431158153326,
    0,
    -55.56517389783991,
    0,
    111.1303477956798,
    0,
    -42.33537058883041,
    0,
    0,
    0,
    0,
    0.9081022627604556,
    0,
    0,
    0,
    0,
    3.632409051041822,
    0,
    -43.58890861250187,
    0,
    0,
    0,
    0,
    0,
    0,
    5.448613576562733,
    0,
    -130.7667258375056,
    0,
    217.9445430625093,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.632409051041822,
    0,
    -130.7667258375056,
    0,
    435.8890861250187,
    0,
    -232.4741792666766,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.9081022627604556,
    0,
    -43.58890861250187,
    0,
    217.9445430625093,
    0,
    -232.4741792666766,
    0,
    49.815895557145,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.718637772708116,
    0,
    0,
    0,
    0,
    0,
    0,
    18.87455109083247,
    0,
    -62.91517030277488,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    28.3118266362487,
    0,
    -188.7455109083247,
    0,
    150.9964087266597,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    18.87455109083247,
    0,
    -188.7455109083247,
    0,
    301.9928174533194,
    0,
    -86.28366212951984,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.718637772708116,
    0,
    -62.91517030277488,
    0,
    150.9964087266597,
    0,
    -86.28366212951984,
    0,
    9.587073569946648,
    0,
    -0.3181304937373671,
    0,
    0,
    -1.590652468686835,
    0,
    15.90652468686835,
    0,
    0,
    0,
    0,
    -3.181304937373671,
    0,
    63.62609874747341,
    0,
    -84.83479832996456,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.181304937373671,
    0,
    95.43914812121012,
    0,
    -254.5043949898937,
    0,
    101.8017579959575,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.590652468686835,
    0,
    63.62609874747341,
    0,
    -254.5043949898937,
    0,
    203.6035159919149,
    0,
    -29.08621657027356,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.3181304937373671,
    0,
    15.90652468686835,
    0,
    -84.83479832996456,
    0,
    101.8017579959575,
    0,
    -29.08621657027356,
    0,
    1.292720736456603,
    0,
    0,
    4.718637772708116,
    0,
    0,
    0,
    0,
    18.87455109083247,
    0,
    -62.91517030277488,
    0,
    0,
    0,
    0,
    0,
    0,
    28.3118266362487,
    0,
    -188.7455109083247,
    0,
    150.9964087266597,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    18.87455109083247,
    0,
    -188.7455109083247,
    0,
    301.9928174533194,
    0,
    -86.28366212951984,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.718637772708116,
    0,
    -62.91517030277488,
    0,
    150.9964087266597,
    0,
    -86.28366212951984,
    0,
    9.587073569946648,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4540511313802278,
    0,
    0,
    1.362153394140683,
    0,
    -21.79445430625093,
    0,
    0,
    0,
    0,
    0.9081022627604556,
    0,
    -43.58890861250187,
    0,
    108.9722715312547,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.9081022627604556,
    0,
    0,
    0,
    108.9722715312547,
    0,
    -116.2370896333383,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.362153394140683,
    0,
    43.58890861250187,
    0,
    -108.9722715312547,
    0,
    0,
    0,
    24.9079477785725,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4540511313802278,
    0,
    21.79445430625093,
    0,
    -108.9722715312547,
    0,
    116.2370896333383,
    0,
    -24.9079477785725,
    0,
    0,
    0,
    0,
    -4.630431158153326,
    0,
    0,
    0,
    0,
    0,
    0,
    55.56517389783991,
    0,
    0,
    0,
    0,
    0,
    0,
    27.78258694891996,
    0,
    -55.56517389783991,
    0,
    -111.1303477956798,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    37.04344926522661,
    0,
    -277.8258694891996,
    0,
    222.2606955913596,
    0,
    42.33537058883041,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    13.89129347445998,
    0,
    -166.6955216935197,
    0,
    333.3910433870395,
    0,
    -127.0061117664912,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4677441816782422,
    0,
    0,
    1.403232545034726,
    0,
    19.64525563048617,
    0,
    0,
    0,
    0,
    6.548418543495391,
    0,
    -78.58102252194469,
    0,
    -78.58102252194469,
    0,
    0,
    0,
    0,
    0,
    0,
    6.548418543495391,
    0,
    -196.4525563048617,
    0,
    392.9051126097235,
    0,
    52.38734834796313,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.403232545034726,
    0,
    -78.58102252194469,
    0,
    392.9051126097235,
    0,
    -314.3240900877788,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4677441816782422,
    0,
    19.64525563048617,
    0,
    -78.58102252194469,
    0,
    52.38734834796313,
    0,
    0,
    0,
    0,
    0,
    0,
    4.437410929184535,
    0,
    0,
    0,
    0,
    -35.49928743347628,
    0,
    -41.41583533905566,
    0,
    0,
    0,
    0,
    0,
    0,
    -62.12375300858349,
    0,
    372.742518051501,
    0,
    49.6990024068668,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    207.0791766952783,
    0,
    -496.990024068668,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    22.18705464592268,
    0,
    -207.0791766952783,
    0,
    248.495012034334,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4961176240878564,
    0,
    0,
    -6.449529113142133,
    0,
    -15.8757639708114,
    0,
    0,
    0,
    0,
    -6.945646737229989,
    0,
    222.2606955913596,
    0,
    37.04344926522661,
    0,
    0,
    0,
    0,
    0,
    0,
    6.945646737229989,
    0,
    0,
    0,
    -555.6517389783992,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.449529113142133,
    0,
    -222.2606955913596,
    0,
    555.6517389783992,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4961176240878564,
    0,
    15.8757639708114,
    0,
    -37.04344926522661,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.091090733689417,
    0,
    0,
    0,
    0,
    81.82181467378834,
    0,
    21.81915057967689,
    0,
    0,
    0,
    0,
    0,
    0,
    -57.27527027165184,
    0,
    -458.2021621732147,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -114.5505405433037,
    0,
    763.6702702886912,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    28.63763513582592,
    0,
    -152.7340540577382,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5567269327204184,
    0,
    0,
    15.0316271834513,
    0,
    10.02108478896753,
    0,
    0,
    0,
    0,
    -23.38253117425757,
    0,
    -280.590374091091,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -23.38253117425757,
    0,
    701.4759352277273,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    15.0316271834513,
    0,
    -280.590374091091,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5567269327204184,
    0,
    10.02108478896753,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.431895299891715,
    0,
    0,
    0,
    0,
    -123.5482307961017,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    432.4188077863561,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -288.2792051909041,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    30.88705769902543,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.7673951182219901,
    0,
    0,
    -34.53278031998956,
    0,
    0,
    0,
    0,
    0,
    0,
    161.1529748266179,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -161.1529748266179,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    34.53278031998956,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.7673951182219901,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // l = 11
    0,
    8.631063163659167,
    0,
    0,
    0,
    0,
    -129.4659474548875,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    362.504652873685,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -258.9318949097751,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    43.15531581829584,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.7846421057871971,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    36.80297698805311,
    0,
    0,
    0,
    0,
    0,
    0,
    -441.6357238566373,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    927.4350200989384,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -441.6357238566373,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    36.80297698805311,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.110940374050938,
    0,
    0,
    0,
    0,
    42.59116978375781,
    0,
    102.2188074810188,
    0,
    0,
    0,
    0,
    0,
    0,
    -23.85105507890438,
    0,
    -954.0422031561751,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -51.10940374050938,
    0,
    1431.063304734263,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    19.87587923242031,
    0,
    -408.875229924075,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5678822637834375,
    0,
    11.35764527566875,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -35.19037680383713,
    0,
    0,
    0,
    0,
    0,
    0,
    211.1422608230228,
    0,
    211.1422608230228,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1477.995825761159,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -211.1422608230228,
    0,
    1477.995825761159,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    35.19037680383713,
    0,
    -211.1422608230228,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.532036427339827,
    0,
    0,
    0,
    0,
    -10.59610928201948,
    0,
    -127.1533113842337,
    0,
    0,
    0,
    0,
    0,
    0,
    -21.19221856403896,
    0,
    508.613245536935,
    0,
    339.0754970246234,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.027459794862709,
    0,
    254.3066227684675,
    0,
    -1695.377485123117,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    9.586956017065244,
    0,
    -363.295175383525,
    0,
    1017.22649107387,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5045766324771181,
    0,
    18.16475876917625,
    0,
    -48.43935671780334,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    28.72100542905686,
    0,
    0,
    0,
    0,
    0,
    0,
    -38.29467390540915,
    0,
    -306.3573912432732,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -134.031358668932,
    0,
    714.8339129009709,
    0,
    428.9003477405824,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -38.29467390540915,
    0,
    714.8339129009709,
    0,
    -1429.667825801941,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    28.72100542905686,
    0,
    -306.3573912432732,
    0,
    428.9003477405824,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.369836079783365,
    0,
    0,
    0,
    0,
    -2.369836079783365,
    0,
    113.7521318296015,
    0,
    0,
    0,
    0,
    0,
    0,
    6.63554102339342,
    0,
    0,
    0,
    -530.8432818714737,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    10.4272787510468,
    0,
    -318.5059691228842,
    0,
    530.8432818714737,
    0,
    424.674625497179,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.31777051169671,
    0,
    -182.0034109273624,
    0,
    955.5179073686526,
    0,
    -849.349250994358,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4739672159566729,
    0,
    22.7504263659203,
    0,
    -106.1686563742947,
    0,
    84.9349250994358,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -20.06399012830402,
    0,
    0,
    0,
    0,
    0,
    0,
    -40.12798025660804,
    0,
    280.8958617962563,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    280.8958617962563,
    0,
    -674.150068311015,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    40.12798025660804,
    0,
    -280.8958617962563,
    0,
    0,
    0,
    321.0238420528643,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    20.06399012830402,
    0,
    -280.8958617962563,
    0,
    674.150068311015,
    0,
    -321.0238420528643,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.373687498354136,
    0,
    0,
    0,
    0,
    5.036854160631831,
    0,
    -76.92649990783158,
    0,
    0,
    0,
    0,
    0,
    0,
    6.410541658985967,
    0,
    -205.1373330875509,
    0,
    461.5589994469895,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.747374996708271,
    0,
    -153.8529998156632,
    0,
    769.2649990783159,
    0,
    -615.4119992626527,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4578958327847119,
    0,
    0,
    0,
    153.8529998156632,
    0,
    -410.2746661751018,
    0,
    175.8319997893294,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4578958327847119,
    0,
    25.64216663594386,
    0,
    -153.8529998156632,
    0,
    205.1373330875509,
    0,
    -58.61066659644312,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    10.27973595067153,
    0,
    0,
    0,
    0,
    0,
    0,
    41.11894380268614,
    0,
    -164.4757752107446,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    61.67841570402921,
    0,
    -493.4273256322336,
    0,
    493.4273256322337,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    41.11894380268614,
    0,
    -493.4273256322336,
    0,
    986.8546512644674,
    0,
    -375.9446290531304,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    10.27973595067153,
    0,
    -164.4757752107446,
    0,
    493.4273256322337,
    0,
    -375.9446290531304,
    0,
    62.65743817552173,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4507962425947618,
    0,
    0,
    0,
    0,
    -2.253981212973809,
    0,
    27.04777455568571,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.507962425947618,
    0,
    108.1910982227429,
    0,
    -180.3184970379047,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.507962425947618,
    0,
    162.2866473341143,
    0,
    -540.9554911137142,
    0,
    288.5095952606476,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.253981212973809,
    0,
    108.1910982227429,
    0,
    -540.9554911137142,
    0,
    577.0191905212952,
    0,
    -123.6469693974204,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4507962425947618,
    0,
    27.04777455568571,
    0,
    -180.3184970379047,
    0,
    288.5095952606476,
    0,
    -123.6469693974204,
    0,
    10.99084172421515,
    0,
    0,
    0,
    -3.662285987505434,
    0,
    0,
    0,
    0,
    -18.31142993752717,
    0,
    61.03809979175723,
    0,
    0,
    0,
    0,
    0,
    0,
    -36.62285987505434,
    0,
    244.1523991670289,
    0,
    -195.3219193336232,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -36.62285987505434,
    0,
    366.2285987505434,
    0,
    -585.9657580008695,
    0,
    167.4187880002484,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -18.31142993752717,
    0,
    244.1523991670289,
    0,
    -585.9657580008695,
    0,
    334.8375760004968,
    0,
    -37.20417511116631,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.662285987505434,
    0,
    61.03809979175723,
    0,
    -195.3219193336232,
    0,
    167.4187880002484,
    0,
    -37.20417511116631,
    0,
    1.352879094951502,
    -0.4507962425947618,
    0,
    0,
    -2.253981212973809,
    0,
    27.04777455568571,
    0,
    0,
    0,
    0,
    -4.507962425947618,
    0,
    108.1910982227429,
    0,
    -180.3184970379047,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.507962425947618,
    0,
    162.2866473341143,
    0,
    -540.9554911137142,
    0,
    288.5095952606476,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.253981212973809,
    0,
    108.1910982227429,
    0,
    -540.9554911137142,
    0,
    577.0191905212952,
    0,
    -123.6469693974204,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4507962425947618,
    0,
    27.04777455568571,
    0,
    -180.3184970379047,
    0,
    288.5095952606476,
    0,
    -123.6469693974204,
    0,
    10.99084172421515,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.139867975335767,
    0,
    0,
    0,
    0,
    15.4196039260073,
    0,
    -82.23788760537228,
    0,
    0,
    0,
    0,
    0,
    0,
    10.27973595067153,
    0,
    -164.4757752107446,
    0,
    246.7136628161169,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -10.27973595067153,
    0,
    0,
    0,
    246.7136628161169,
    0,
    -187.9723145265652,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -15.4196039260073,
    0,
    164.4757752107446,
    0,
    -246.7136628161169,
    0,
    0,
    0,
    31.32871908776087,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.139867975335767,
    0,
    82.23788760537228,
    0,
    -246.7136628161169,
    0,
    187.9723145265652,
    0,
    -31.32871908776087,
    0,
    0,
    0.4578958327847119,
    0,
    0,
    0.4578958327847119,
    0,
    -25.64216663594386,
    0,
    0,
    0,
    0,
    -2.747374996708271,
    0,
    0,
    0,
    153.8529998156632,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.410541658985967,
    0,
    153.8529998156632,
    0,
    -153.8529998156632,
    0,
    -205.1373330875509,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.036854160631831,
    0,
    205.1373330875509,
    0,
    -769.2649990783159,
    0,
    410.2746661751018,
    0,
    58.61066659644312,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.373687498354136,
    0,
    76.92649990783158,
    0,
    -461.5589994469895,
    0,
    615.4119992626527,
    0,
    -175.8319997893294,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.015997532076005,
    0,
    0,
    0,
    0,
    15.04799259622802,
    0,
    70.22396544906408,
    0,
    0,
    0,
    0,
    0,
    0,
    70.22396544906407,
    0,
    -280.8958617962563,
    0,
    -168.5375170777538,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    70.22396544906407,
    0,
    -702.2396544906408,
    0,
    842.6875853887689,
    0,
    80.25596051321608,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    15.04799259622802,
    0,
    -280.8958617962563,
    0,
    842.6875853887689,
    0,
    -481.5357630792965,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.015997532076005,
    0,
    70.22396544906408,
    0,
    -168.5375170777538,
    0,
    80.25596051321608,
    0,
    0,
    0,
    0,
    -0.4739672159566729,
    0,
    0,
    3.31777051169671,
    0,
    22.7504263659203,
    0,
    0,
    0,
    0,
    10.4272787510468,
    0,
    -182.0034109273624,
    0,
    -106.1686563742947,
    0,
    0,
    0,
    0,
    0,
    0,
    6.63554102339342,
    0,
    -318.5059691228842,
    0,
    955.5179073686526,
    0,
    84.9349250994358,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.369836079783365,
    0,
    0,
    0,
    530.8432818714737,
    0,
    -849.349250994358,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.369836079783365,
    0,
    113.7521318296015,
    0,
    -530.8432818714737,
    0,
    424.674625497179,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.786834238176144,
    0,
    0,
    0,
    0,
    -62.22884509628987,
    0,
    -51.0595652072122,
    0,
    0,
    0,
    0,
    0,
    0,
    -67.01567933446601,
    0,
    714.8339129009709,
    0,
    71.48339129009707,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    67.01567933446601,
    0,
    0,
    0,
    -1072.250869351456,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    62.22884509628987,
    0,
    -714.8339129009709,
    0,
    1072.250869351456,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.786834238176144,
    0,
    51.0595652072122,
    0,
    -71.48339129009707,
    0,
    0,
    0,
    0,
    0,
    0,
    0.5045766324771181,
    0,
    0,
    -9.586956017065244,
    0,
    -18.16475876917625,
    0,
    0,
    0,
    0,
    -3.027459794862709,
    0,
    363.295175383525,
    0,
    48.43935671780334,
    0,
    0,
    0,
    0,
    0,
    0,
    21.19221856403896,
    0,
    -254.3066227684675,
    0,
    -1017.22649107387,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    10.59610928201948,
    0,
    -508.613245536935,
    0,
    1695.377485123117,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.532036427339827,
    0,
    127.1533113842337,
    0,
    -339.0754970246234,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.398797100479641,
    0,
    0,
    0,
    0,
    118.7675217129503,
    0,
    26.39278260287784,
    0,
    0,
    0,
    0,
    0,
    0,
    -184.7494782201449,
    0,
    -738.9979128805796,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -184.7494782201449,
    0,
    1847.494782201449,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    118.7675217129503,
    0,
    -738.9979128805796,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.398797100479641,
    0,
    26.39278260287784,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5678822637834375,
    0,
    0,
    19.87587923242031,
    0,
    11.35764527566875,
    0,
    0,
    0,
    0,
    -51.10940374050938,
    0,
    -408.875229924075,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -23.85105507890438,
    0,
    1431.063304734263,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    42.59116978375781,
    0,
    -954.0422031561751,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.110940374050938,
    0,
    102.2188074810188,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.680297698805311,
    0,
    0,
    0,
    0,
    -165.613396446239,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    772.8625167491152,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -772.8625167491152,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    165.613396446239,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.680297698805311,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.7846421057871971,
    0,
    0,
    -43.15531581829584,
    0,
    0,
    0,
    0,
    0,
    0,
    258.9318949097751,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -362.504652873685,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    129.4659474548875,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -8.631063163659167,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // l = 12
    0,
    9.609863949407661,
    0,
    0,
    0,
    0,
    -176.1808390724738,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    634.2510206609056,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -634.2510206609056,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    176.1808390724738,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -9.609863949407661,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    43.15531581829583,
    0,
    0,
    0,
    0,
    0,
    0,
    -647.3297372744373,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1812.523264368425,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1294.659474548875,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    215.7765790914791,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.923210528935984,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.784458347938102,
    0,
    0,
    0,
    0,
    63.62904182731912,
    0,
    127.2580836546383,
    0,
    0,
    0,
    0,
    0,
    0,
    -76.35485019278295,
    0,
    -1527.097003855659,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -76.35485019278295,
    0,
    3206.903708096884,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    63.62904182731912,
    0,
    -1527.097003855659,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.784458347938102,
    0,
    127.2580836546383,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -42.2938455917996,
    0,
    0,
    0,
    0,
    0,
    0,
    352.4487132649967,
    0,
    281.9589706119974,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -197.3712794283981,
    0,
    -2631.617059045309,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -422.938455917996,
    0,
    3947.425588567963,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    164.4760661903318,
    0,
    -1127.835882447989,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.699316176866622,
    0,
    31.32877451244415,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.101899446670816,
    0,
    0,
    0,
    0,
    -20.50949723335408,
    0,
    -164.0759778668327,
    0,
    0,
    0,
    0,
    0,
    0,
    -24.6113966800249,
    0,
    984.455867200996,
    0,
    492.2279336004979,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    24.6113966800249,
    0,
    0,
    0,
    -3445.595535203485,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    20.50949723335408,
    0,
    -984.455867200996,
    0,
    3445.595535203485,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.101899446670816,
    0,
    164.0759778668327,
    0,
    -492.2279336004979,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    35.89162015836965,
    0,
    0,
    0,
    0,
    0,
    0,
    -107.6748604751089,
    0,
    -430.6994419004357,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -215.3497209502179,
    0,
    1722.797767601743,
    0,
    689.1191070406971,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    30.76424585003112,
    0,
    861.3988838008713,
    0,
    -3445.595535203485,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    97.42011185843189,
    0,
    -1230.569834001245,
    0,
    2067.357321122091,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.127374308338521,
    0,
    61.52849170006224,
    0,
    -98.44558672009958,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.881335616715016,
    0,
    0,
    0,
    0,
    0.9604452055716719,
    0,
    155.5921233026108,
    0,
    0,
    0,
    0,
    0,
    0,
    17.28801370029009,
    0,
    -207.4561644034811,
    0,
    -829.8246576139245,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    17.28801370029009,
    0,
    -726.0965754121839,
    0,
    1936.257534432491,
    0,
    774.5030137729962,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.9604452055716719,
    0,
    -207.4561644034811,
    0,
    1936.257534432491,
    0,
    -2581.676712576654,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.881335616715016,
    0,
    155.5921233026108,
    0,
    -829.8246576139245,
    0,
    774.5030137729962,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -26.95242673514147,
    0,
    0,
    0,
    0,
    0,
    0,
    -26.95242673514147,
    0,
    431.2388277622634,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    75.4667948583961,
    0,
    0,
    0,
    -1207.468717734338,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    118.5906776346225,
    0,
    -1207.468717734338,
    0,
    1207.468717734338,
    0,
    689.9821244196215,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    37.73339742919805,
    0,
    -689.9821244196215,
    0,
    2173.443691921808,
    0,
    -1379.964248839243,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.390485347028293,
    0,
    86.24776555245269,
    0,
    -241.4937435468676,
    0,
    137.9964248839243,
    0,
    0,
    0,
    0,
    0,
    0,
    1.848921220493557,
    0,
    0,
    0,
    0,
    5.54676366148067,
    0,
    -118.3309581115876,
    0,
    0,
    0,
    0,
    0,
    0,
    3.697842440987113,
    0,
    -236.6619162231752,
    0,
    828.3167067811135,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.697842440987113,
    0,
    0,
    0,
    828.3167067811135,
    0,
    -1325.306730849781,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.54676366148067,
    0,
    236.6619162231752,
    0,
    -828.3167067811135,
    0,
    0,
    0,
    473.3238324463505,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.848921220493557,
    0,
    118.3309581115876,
    0,
    -828.3167067811135,
    0,
    1325.306730849781,
    0,
    -473.3238324463505,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    16.64029098444201,
    0,
    0,
    0,
    0,
    0,
    0,
    61.01440027628737,
    0,
    -310.6187650429175,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    77.65469126072938,
    0,
    -828.3167067811133,
    0,
    1118.227554154503,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    33.28058196888402,
    0,
    -621.237530085835,
    0,
    1863.712590257505,
    0,
    -1064.978623004289,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.54676366148067,
    0,
    0,
    0,
    372.742518051501,
    0,
    -709.9857486695257,
    0,
    236.6619162231752,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.54676366148067,
    0,
    103.5395883476392,
    0,
    -372.742518051501,
    0,
    354.9928743347629,
    0,
    -78.88730540772508,
    0,
    0,
    0,
    0,
    -0.9057827129626244,
    0,
    0,
    0,
    0,
    -4.528913564813122,
    0,
    63.4047899073837,
    0,
    0,
    0,
    0,
    0,
    0,
    -9.057827129626244,
    0,
    253.6191596295348,
    0,
    -507.2383192590696,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -9.057827129626244,
    0,
    380.4287394443022,
    0,
    -1521.714957777209,
    0,
    1014.476638518139,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.528913564813122,
    0,
    253.6191596295348,
    0,
    -1521.714957777209,
    0,
    2028.953277036278,
    0,
    -579.7009362960796,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.9057827129626244,
    0,
    63.4047899073837,
    0,
    -507.2383192590696,
    0,
    1014.476638518139,
    0,
    -579.7009362960796,
    0,
    77.2934581728106,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.620233931023189,
    0,
    0,
    0,
    0,
    0,
    0,
    -28.10116965511595,
    0,
    112.4046786204638,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -56.20233931023189,
    0,
    449.6187144818551,
    0,
    -449.6187144818551,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -56.20233931023189,
    0,
    674.4280717227828,
    0,
    -1348.856143445566,
    0,
    513.8499594078344,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -28.10116965511595,
    0,
    449.6187144818551,
    0,
    -1348.856143445566,
    0,
    1027.699918815669,
    0,
    -171.2833198026115,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.620233931023189,
    0,
    112.4046786204638,
    0,
    -449.6187144818551,
    0,
    513.8499594078344,
    0,
    -171.2833198026115,
    0,
    12.4569687129172,
    0,
    0.318183090330888,
    0,
    0,
    1.909098541985328,
    0,
    -22.90918250382393,
    0,
    0,
    0,
    0,
    4.77274635496332,
    0,
    -114.5459125191197,
    0,
    190.9098541985328,
    0,
    0,
    0,
    0,
    0,
    0,
    6.36366180661776,
    0,
    -229.0918250382393,
    0,
    763.6394167941311,
    0,
    -407.2743556235366,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.77274635496332,
    0,
    -229.0918250382393,
    0,
    1145.459125191197,
    0,
    -1221.82306687061,
    0,
    261.8192286151307,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.909098541985328,
    0,
    -114.5459125191197,
    0,
    763.6394167941311,
    0,
    -1221.82306687061,
    0,
    523.6384572302613,
    0,
    -46.5456406426899,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.318183090330888,
    0,
    -22.90918250382393,
    0,
    190.9098541985328,
    0,
    -407.2743556235366,
    0,
    261.8192286151307,
    0,
    -46.5456406426899,
    0,
    1.410473958869391,
    0,
    0,
    -5.620233931023189,
    0,
    0,
    0,
    0,
    -28.10116965511595,
    0,
    112.4046786204638,
    0,
    0,
    0,
    0,
    0,
    0,
    -56.20233931023189,
    0,
    449.6187144818551,
    0,
    -449.6187144818551,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -56.20233931023189,
    0,
    674.4280717227828,
    0,
    -1348.856143445566,
    0,
    513.8499594078344,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -28.10116965511595,
    0,
    449.6187144818551,
    0,
    -1348.856143445566,
    0,
    1027.699918815669,
    0,
    -171.2833198026115,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.620233931023189,
    0,
    112.4046786204638,
    0,
    -449.6187144818551,
    0,
    513.8499594078344,
    0,
    -171.2833198026115,
    0,
    12.4569687129172,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4528913564813122,
    0,
    0,
    -1.811565425925249,
    0,
    31.70239495369185,
    0,
    0,
    0,
    0,
    -2.264456782406561,
    0,
    95.10718486107555,
    0,
    -253.6191596295348,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    63.4047899073837,
    0,
    -507.2383192590696,
    0,
    507.2383192590696,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.264456782406561,
    0,
    -63.4047899073837,
    0,
    0,
    0,
    507.2383192590696,
    0,
    -289.8504681480398,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.811565425925249,
    0,
    -95.10718486107555,
    0,
    507.2383192590696,
    0,
    -507.2383192590696,
    0,
    0,
    0,
    38.6467290864053,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4528913564813122,
    0,
    -31.70239495369185,
    0,
    253.6191596295348,
    0,
    -507.2383192590696,
    0,
    289.8504681480398,
    0,
    -38.6467290864053,
    0,
    0,
    0,
    0,
    5.54676366148067,
    0,
    0,
    0,
    0,
    5.54676366148067,
    0,
    -103.5395883476392,
    0,
    0,
    0,
    0,
    0,
    0,
    -33.28058196888402,
    0,
    0,
    0,
    372.742518051501,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -77.65469126072938,
    0,
    621.237530085835,
    0,
    -372.742518051501,
    0,
    -354.9928743347629,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -61.01440027628737,
    0,
    828.3167067811133,
    0,
    -1863.712590257505,
    0,
    709.9857486695257,
    0,
    78.88730540772508,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -16.64029098444201,
    0,
    310.6187650429175,
    0,
    -1118.227554154503,
    0,
    1064.978623004289,
    0,
    -236.6619162231752,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4622303051233891,
    0,
    0,
    -0.9244606102467783,
    0,
    -29.5827395278969,
    0,
    0,
    0,
    0,
    -7.857915187097616,
    0,
    88.74821858369071,
    0,
    207.0791766952784,
    0,
    0,
    0,
    0,
    0,
    0,
    -12.9424485434549,
    0,
    414.1583533905567,
    0,
    -828.3167067811135,
    0,
    -331.3266827124453,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -7.857915187097616,
    0,
    414.1583533905567,
    0,
    -2070.791766952784,
    0,
    1656.633413562227,
    0,
    118.3309581115876,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.9244606102467783,
    0,
    88.74821858369071,
    0,
    -828.3167067811135,
    0,
    1656.633413562227,
    0,
    -709.9857486695257,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4622303051233891,
    0,
    -29.5827395278969,
    0,
    207.0791766952784,
    0,
    -331.3266827124453,
    0,
    118.3309581115876,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.390485347028293,
    0,
    0,
    0,
    0,
    37.73339742919805,
    0,
    86.24776555245269,
    0,
    0,
    0,
    0,
    0,
    0,
    118.5906776346225,
    0,
    -689.9821244196215,
    0,
    -241.4937435468676,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    75.4667948583961,
    0,
    -1207.468717734338,
    0,
    2173.443691921808,
    0,
    137.9964248839243,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -26.95242673514147,
    0,
    0,
    0,
    1207.468717734338,
    0,
    -1379.964248839243,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -26.95242673514147,
    0,
    431.2388277622634,
    0,
    -1207.468717734338,
    0,
    689.9821244196215,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.480222602785836,
    0,
    0,
    5.762671233430032,
    0,
    25.93202055043514,
    0,
    0,
    0,
    0,
    12.96601027521757,
    0,
    -337.1162671556568,
    0,
    -138.3041096023208,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -363.048287706092,
    0,
    1936.257534432491,
    0,
    129.0838356288327,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -12.96601027521757,
    0,
    363.048287706092,
    0,
    0,
    0,
    -1936.257534432491,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.762671233430032,
    0,
    337.1162671556568,
    0,
    -1936.257534432491,
    0,
    1936.257534432491,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.480222602785836,
    0,
    -25.93202055043514,
    0,
    138.3041096023208,
    0,
    -129.0838356288327,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.127374308338521,
    0,
    0,
    0,
    0,
    -97.42011185843189,
    0,
    -61.52849170006224,
    0,
    0,
    0,
    0,
    0,
    0,
    -30.76424585003112,
    0,
    1230.569834001245,
    0,
    98.44558672009958,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    215.3497209502179,
    0,
    -861.3988838008713,
    0,
    -2067.357321122091,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    107.6748604751089,
    0,
    -1722.797767601743,
    0,
    3445.595535203485,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -35.89162015836965,
    0,
    430.6994419004357,
    0,
    -689.1191070406971,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.512737430833852,
    0,
    0,
    -13.33117320168015,
    0,
    -20.50949723335408,
    0,
    0,
    0,
    0,
    7.691061462507781,
    0,
    553.7564253005602,
    0,
    61.52849170006224,
    0,
    0,
    0,
    0,
    0,
    0,
    43.06994419004357,
    0,
    -861.3988838008714,
    0,
    -1722.797767601743,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    7.691061462507781,
    0,
    -861.3988838008714,
    0,
    4306.994419004357,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -13.33117320168015,
    0,
    553.7564253005602,
    0,
    -1722.797767601743,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.512737430833852,
    0,
    -20.50949723335408,
    0,
    61.52849170006224,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.699316176866622,
    0,
    0,
    0,
    0,
    164.4760661903318,
    0,
    31.32877451244415,
    0,
    0,
    0,
    0,
    0,
    0,
    -422.938455917996,
    0,
    -1127.835882447989,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -197.3712794283981,
    0,
    3947.425588567963,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    352.4487132649967,
    0,
    -2631.617059045309,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -42.2938455917996,
    0,
    281.9589706119974,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5784458347938102,
    0,
    0,
    25.45161673092765,
    0,
    12.72580836546383,
    0,
    0,
    0,
    0,
    -95.44356274097868,
    0,
    -572.6613764458722,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2672.419756747403,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    95.44356274097868,
    0,
    -2672.419756747403,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -25.45161673092765,
    0,
    572.6613764458722,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.5784458347938102,
    0,
    -12.72580836546383,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.923210528935984,
    0,
    0,
    0,
    0,
    -215.7765790914791,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1294.659474548875,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1812.523264368425,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    647.3297372744373,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -43.15531581829583,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.8008219957839717,
    0,
    0,
    -52.85425172174213,
    0,
    0,
    0,
    0,
    0,
    0,
    396.406887913066,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -739.9595241043899,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    396.406887913066,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -52.85425172174213,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.8008219957839717,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // l = 13
    0,
    10.60900254488917,
    0,
    0,
    0,
    0,
    -233.3980559875617,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1050.291251944028,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1400.38833592537,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    583.4951399689042,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -63.65401526933501,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.8160771188376283,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    49.93431784259574,
    0,
    0,
    0,
    0,
    0,
    0,
    -915.4624937809218,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3295.664977611319,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3295.664977611319,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    915.4624937809218,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -49.93431784259574,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.473297372744374,
    0,
    0,
    0,
    0,
    90.62616321842124,
    0,
    155.359136945865,
    0,
    0,
    0,
    0,
    0,
    0,
    -174.7790290640981,
    0,
    -2330.387054187975,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -77.67956847293249,
    0,
    6525.083751726329,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    161.8324343186094,
    0,
    -4660.774108375949,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -31.77800528438147,
    0,
    776.7956847293249,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.5884815793403977,
    0,
    -14.12355790416954,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -49.93431784259574,
    0,
    0,
    0,
    0,
    0,
    0,
    549.2774962685531,
    0,
    366.1849975123687,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -659.1329955222637,
    0,
    -4394.219970148425,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -659.1329955222637,
    0,
    9227.861937311692,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    549.2774962685531,
    0,
    -4394.219970148425,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -49.93431784259574,
    0,
    366.1849975123687,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.685411261829863,
    0,
    0,
    0,
    0,
    -34.35968258675233,
    0,
    -206.158095520514,
    0,
    0,
    0,
    0,
    0,
    0,
    -17.17984129337616,
    0,
    1717.984129337616,
    0,
    687.1936517350465,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    68.71936517350465,
    0,
    -962.0711124290651,
    0,
    -6413.807416193768,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    28.63306882229361,
    0,
    -2061.58095520514,
    0,
    9620.711124290651,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -17.70044254469059,
    0,
    801.7259270242209,
    0,
    -2748.774606940186,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.5206012513144292,
    0,
    -22.90645505783488,
    0,
    76.35485019278295,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    43.68089589976209,
    0,
    0,
    0,
    0,
    0,
    0,
    -218.4044794988104,
    0,
    -582.4119453301611,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -262.0853753985725,
    0,
    3494.471671980967,
    0,
    1048.34150159429,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    262.0853753985725,
    0,
    0,
    0,
    -7338.39051116003,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    218.4044794988104,
    0,
    -3494.471671980967,
    0,
    7338.39051116003,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -43.68089589976209,
    0,
    582.4119453301611,
    0,
    -1048.34150159429,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.404978058421841,
    0,
    0,
    0,
    0,
    6.809956116843682,
    0,
    204.2986835053105,
    0,
    0,
    0,
    0,
    0,
    0,
    30.64480252579657,
    0,
    -612.8960505159314,
    0,
    -1225.792101031863,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    17.51131572902661,
    0,
    -1225.792101031863,
    0,
    4903.168404127451,
    0,
    1307.511574433987,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -12.16063592293515,
    0,
    175.1131572902661,
    0,
    2451.584202063726,
    0,
    -6537.557872169935,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -8.755657864513306,
    0,
    554.5249980858427,
    0,
    -3502.263145805322,
    0,
    3922.534723301961,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4864254369174059,
    0,
    -29.18552621504435,
    0,
    175.1131572902661,
    0,
    -186.7873677762839,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -34.53278031998955,
    0,
    0,
    0,
    0,
    0,
    0,
    11.51092677332985,
    0,
    621.5900457598119,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    207.1966819199373,
    0,
    -828.7867276797492,
    0,
    -1989.088146431398,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    207.1966819199373,
    0,
    -2900.753546879122,
    0,
    4641.205675006596,
    0,
    1326.058764287599,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    11.51092677332985,
    0,
    -828.7867276797492,
    0,
    4641.205675006596,
    0,
    -4420.195880958662,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -34.53278031998955,
    0,
    621.5900457598119,
    0,
    -1989.088146431398,
    0,
    1326.058764287599,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.334148624627139,
    0,
    0,
    0,
    0,
    4.668297249254278,
    0,
    -168.058700973154,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.20146752432885,
    0,
    -168.058700973154,
    0,
    1344.469607785232,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -16.8058700973154,
    0,
    470.5643627248312,
    0,
    0,
    0,
    -2509.6766011991,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -13.53806202283741,
    0,
    739.4582842818776,
    0,
    -3764.51490179865,
    0,
    2509.6766011991,
    0,
    1075.575686228186,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.800978349552567,
    0,
    235.2821813624156,
    0,
    -2151.151372456371,
    0,
    4517.41788215838,
    0,
    -2151.151372456371,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4668297249254278,
    0,
    -33.6117401946308,
    0,
    268.8939215570464,
    0,
    -501.93532023982,
    0,
    215.1151372456371,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    23.76708941910389,
    0,
    0,
    0,
    0,
    0,
    0,
    71.30126825731166,
    0,
    -507.0312409408829,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    47.53417883820777,
    0,
    -1014.062481881766,
    0,
    2129.531211951708,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -47.53417883820777,
    0,
    0,
    0,
    2129.531211951708,
    0,
    -2433.749956516238,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -71.30126825731166,
    0,
    1014.062481881766,
    0,
    -2129.531211951708,
    0,
    0,
    0,
    676.0416545878439,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -23.76708941910389,
    0,
    507.0312409408829,
    0,
    -2129.531211951708,
    0,
    2433.749956516238,
    0,
    -676.0416545878439,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.36713941034431,
    0,
    0,
    0,
    0,
    -6.379983914940114,
    0,
    109.3711528275448,
    0,
    0,
    0,
    0,
    0,
    0,
    -11.39282841953592,
    0,
    401.0275603676643,
    0,
    -1020.797426390418,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -9.114262735628734,
    0,
    510.3987131952091,
    0,
    -2722.126470374449,
    0,
    2449.913823337004,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.278565683907184,
    0,
    218.7423056550896,
    0,
    -2041.594852780836,
    0,
    4083.189705561673,
    0,
    -1749.938445240717,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.9114262735628734,
    0,
    -36.45705094251494,
    0,
    0,
    0,
    816.6379411123346,
    0,
    -1166.625630160478,
    0,
    311.1001680427941,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4557131367814367,
    0,
    -36.45705094251494,
    0,
    340.2658087968061,
    0,
    -816.6379411123346,
    0,
    583.312815080239,
    0,
    -103.7000560142647,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -12.09143589391947,
    0,
    0,
    0,
    0,
    0,
    0,
    -60.45717946959737,
    0,
    282.1335041914544,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -120.9143589391947,
    0,
    1128.534016765818,
    0,
    -1354.240820118981,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -120.9143589391947,
    0,
    1692.801025148726,
    0,
    -4062.722460356943,
    0,
    1934.629743027116,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -60.45717946959737,
    0,
    1128.534016765818,
    0,
    -4062.722460356943,
    0,
    3869.259486054231,
    0,
    -859.8354413453848,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -12.09143589391947,
    0,
    282.1335041914544,
    0,
    -1354.240820118981,
    0,
    1934.629743027116,
    0,
    -859.8354413453848,
    0,
    93.80022996495107,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4506212100730813,
    0,
    0,
    0,
    0,
    2.703727260438488,
    0,
    -37.85218164613883,
    0,
    0,
    0,
    0,
    0,
    0,
    6.75931815109622,
    0,
    -189.2609082306941,
    0,
    378.5218164613883,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    9.012424201461626,
    0,
    -378.5218164613883,
    0,
    1514.087265845553,
    0,
    -1009.391510563702,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.75931815109622,
    0,
    -378.5218164613883,
    0,
    2271.13089876833,
    0,
    -3028.174531691106,
    0,
    865.1927233403161,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.703727260438488,
    0,
    -189.2609082306941,
    0,
    1514.087265845553,
    0,
    -3028.174531691106,
    0,
    1730.385446680632,
    0,
    -230.7180595574176,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4506212100730813,
    0,
    -37.85218164613883,
    0,
    378.5218164613883,
    0,
    -1009.391510563702,
    0,
    865.1927233403161,
    0,
    -230.7180595574176,
    0,
    13.98291270044955,
    0,
    0,
    0,
    4.298652372786529,
    0,
    0,
    0,
    0,
    25.79191423671917,
    0,
    -103.1676569468767,
    0,
    0,
    0,
    0,
    0,
    0,
    64.47978559179793,
    0,
    -515.8382847343835,
    0,
    515.8382847343835,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    85.97304745573058,
    0,
    -1031.676569468767,
    0,
    2063.353138937534,
    0,
    -786.0392910238224,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    64.47978559179793,
    0,
    -1031.676569468767,
    0,
    3095.029708406301,
    0,
    -2358.117873071467,
    0,
    393.0196455119112,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    25.79191423671917,
    0,
    -515.8382847343835,
    0,
    2063.353138937534,
    0,
    -2358.117873071467,
    0,
    786.0392910238224,
    0,
    -57.16649389264163,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.298652372786529,
    0,
    -103.1676569468767,
    0,
    515.8382847343835,
    0,
    -786.0392910238224,
    0,
    393.0196455119112,
    0,
    -57.16649389264163,
    0,
    1.46580753570876,
    0.4506212100730813,
    0,
    0,
    2.703727260438488,
    0,
    -37.85218164613883,
    0,
    0,
    0,
    0,
    6.75931815109622,
    0,
    -189.2609082306941,
    0,
    378.5218164613883,
    0,
    0,
    0,
    0,
    0,
    0,
    9.012424201461626,
    0,
    -378.5218164613883,
    0,
    1514.087265845553,
    0,
    -1009.391510563702,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.75931815109622,
    0,
    -378.5218164613883,
    0,
    2271.13089876833,
    0,
    -3028.174531691106,
    0,
    865.1927233403161,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.703727260438488,
    0,
    -189.2609082306941,
    0,
    1514.087265845553,
    0,
    -3028.174531691106,
    0,
    1730.385446680632,
    0,
    -230.7180595574176,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4506212100730813,
    0,
    -37.85218164613883,
    0,
    378.5218164613883,
    0,
    -1009.391510563702,
    0,
    865.1927233403161,
    0,
    -230.7180595574176,
    0,
    13.98291270044955,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.045717946959737,
    0,
    0,
    0,
    0,
    -24.18287178783895,
    0,
    141.0667520957272,
    0,
    0,
    0,
    0,
    0,
    0,
    -30.22858973479868,
    0,
    423.2002562871816,
    0,
    -677.1204100594905,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    282.1335041914544,
    0,
    -1354.240820118981,
    0,
    967.3148715135579,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    30.22858973479868,
    0,
    -282.1335041914544,
    0,
    0,
    0,
    967.3148715135579,
    0,
    -429.9177206726924,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    24.18287178783895,
    0,
    -423.2002562871816,
    0,
    1354.240820118981,
    0,
    -967.3148715135579,
    0,
    0,
    0,
    46.90011498247553,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.045717946959737,
    0,
    -141.0667520957272,
    0,
    677.1204100594905,
    0,
    -967.3148715135579,
    0,
    429.9177206726924,
    0,
    -46.90011498247553,
    0,
    0,
    -0.4557131367814367,
    0,
    0,
    -0.9114262735628734,
    0,
    36.45705094251494,
    0,
    0,
    0,
    0,
    2.278565683907184,
    0,
    36.45705094251494,
    0,
    -340.2658087968061,
    0,
    0,
    0,
    0,
    0,
    0,
    9.114262735628734,
    0,
    -218.7423056550896,
    0,
    0,
    0,
    816.6379411123346,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    11.39282841953592,
    0,
    -510.3987131952091,
    0,
    2041.594852780836,
    0,
    -816.6379411123346,
    0,
    -583.312815080239,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.379983914940114,
    0,
    -401.0275603676643,
    0,
    2722.126470374449,
    0,
    -4083.189705561673,
    0,
    1166.625630160478,
    0,
    103.7000560142647,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.36713941034431,
    0,
    -109.3711528275448,
    0,
    1020.797426390418,
    0,
    -2449.913823337004,
    0,
    1749.938445240717,
    0,
    -311.1001680427941,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.941772354775972,
    0,
    0,
    0,
    0,
    -11.88354470955194,
    0,
    -126.7578102352207,
    0,
    0,
    0,
    0,
    0,
    0,
    -101.0101300311915,
    0,
    380.2734307056622,
    0,
    532.3828029879271,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -166.3696259337272,
    0,
    1774.60934329309,
    0,
    -2129.531211951708,
    0,
    -608.4374891290595,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -101.0101300311915,
    0,
    1774.60934329309,
    0,
    -5323.828029879271,
    0,
    3042.187445645297,
    0,
    169.010413646961,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -11.88354470955194,
    0,
    380.2734307056622,
    0,
    -2129.531211951708,
    0,
    3042.187445645297,
    0,
    -1014.062481881766,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.941772354775972,
    0,
    -126.7578102352207,
    0,
    532.3828029879271,
    0,
    -608.4374891290595,
    0,
    169.010413646961,
    0,
    0,
    0,
    0,
    0.4668297249254278,
    0,
    0,
    -2.800978349552567,
    0,
    -33.6117401946308,
    0,
    0,
    0,
    0,
    -13.53806202283741,
    0,
    235.2821813624156,
    0,
    268.8939215570464,
    0,
    0,
    0,
    0,
    0,
    0,
    -16.8058700973154,
    0,
    739.4582842818776,
    0,
    -2151.151372456371,
    0,
    -501.93532023982,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.20146752432885,
    0,
    470.5643627248312,
    0,
    -3764.51490179865,
    0,
    4517.41788215838,
    0,
    215.1151372456371,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.668297249254278,
    0,
    -168.058700973154,
    0,
    0,
    0,
    2509.6766011991,
    0,
    -2151.151372456371,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.334148624627139,
    0,
    -168.058700973154,
    0,
    1344.469607785232,
    0,
    -2509.6766011991,
    0,
    1075.575686228186,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.755463386664925,
    0,
    0,
    0,
    0,
    69.0655606399791,
    0,
    103.5983409599687,
    0,
    0,
    0,
    0,
    0,
    0,
    155.397511439953,
    0,
    -1346.778432479592,
    0,
    -331.5146910718997,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1450.376773439561,
    0,
    4641.205675006596,
    0,
    221.0097940479331,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -155.397511439953,
    0,
    1450.376773439561,
    0,
    0,
    0,
    -3315.146910718997,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -69.0655606399791,
    0,
    1346.778432479592,
    0,
    -4641.205675006596,
    0,
    3315.146910718997,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.755463386664925,
    0,
    -103.5983409599687,
    0,
    331.5146910718997,
    0,
    -221.0097940479331,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4864254369174059,
    0,
    0,
    8.755657864513306,
    0,
    29.18552621504435,
    0,
    0,
    0,
    0,
    12.16063592293515,
    0,
    -554.5249980858427,
    0,
    -175.1131572902661,
    0,
    0,
    0,
    0,
    0,
    0,
    -17.51131572902661,
    0,
    -175.1131572902661,
    0,
    3502.263145805322,
    0,
    186.7873677762839,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -30.64480252579657,
    0,
    1225.792101031863,
    0,
    -2451.584202063726,
    0,
    -3922.534723301961,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.809956116843682,
    0,
    612.8960505159314,
    0,
    -4903.168404127451,
    0,
    6537.557872169935,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.404978058421841,
    0,
    -204.2986835053105,
    0,
    1225.792101031863,
    0,
    -1307.511574433987,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.460111987470261,
    0,
    0,
    0,
    0,
    -141.9629116742268,
    0,
    -72.80149316627014,
    0,
    0,
    0,
    0,
    0,
    0,
    81.90167981205391,
    0,
    1965.640315489294,
    0,
    131.0426876992863,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    458.6494069475019,
    0,
    -3057.662712983346,
    0,
    -3669.195255580015,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    81.90167981205391,
    0,
    -3057.662712983346,
    0,
    9172.988138950038,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -141.9629116742268,
    0,
    1965.640315489294,
    0,
    -3669.195255580015,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.460111987470261,
    0,
    -72.80149316627014,
    0,
    131.0426876992863,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.5206012513144292,
    0,
    0,
    -17.70044254469059,
    0,
    -22.90645505783488,
    0,
    0,
    0,
    0,
    28.63306882229361,
    0,
    801.7259270242209,
    0,
    76.35485019278295,
    0,
    0,
    0,
    0,
    0,
    0,
    68.71936517350465,
    0,
    -2061.58095520514,
    0,
    -2748.774606940186,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -17.17984129337616,
    0,
    -962.0711124290651,
    0,
    9620.711124290651,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -34.35968258675233,
    0,
    1717.984129337616,
    0,
    -6413.807416193768,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.685411261829863,
    0,
    -206.158095520514,
    0,
    687.1936517350465,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.993431784259574,
    0,
    0,
    0,
    0,
    219.7109985074212,
    0,
    36.61849975123687,
    0,
    0,
    0,
    0,
    0,
    0,
    -823.9162444028297,
    0,
    -1647.832488805659,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    7689.884947759744,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    823.9162444028297,
    0,
    -7689.884947759744,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -219.7109985074212,
    0,
    1647.832488805659,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.993431784259574,
    0,
    -36.61849975123687,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5884815793403977,
    0,
    0,
    31.77800528438147,
    0,
    14.12355790416954,
    0,
    0,
    0,
    0,
    -161.8324343186094,
    0,
    -776.7956847293249,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    77.67956847293249,
    0,
    4660.774108375949,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    174.7790290640981,
    0,
    -6525.083751726329,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -90.62616321842124,
    0,
    2330.387054187975,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.473297372744374,
    0,
    -155.359136945865,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.161193153549645,
    0,
    0,
    0,
    0,
    -274.6387481342766,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2059.790611007074,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3844.942473879872,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2059.790611007074,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -274.6387481342766,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.161193153549645,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.8160771188376283,
    0,
    0,
    -63.65401526933501,
    0,
    0,
    0,
    0,
    0,
    0,
    583.4951399689042,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1400.38833592537,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1050.291251944028,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -233.3980559875617,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    10.60900254488917,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // l = 14
    0,
    11.62730916290334,
    0,
    0,
    0,
    0,
    -302.3100382354867,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1662.705210295177,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2850.351789077446,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1662.705210295177,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -302.3100382354867,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    11.62730916290334,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    57.13122714353754,
    0,
    0,
    0,
    0,
    0,
    0,
    -1256.886997157826,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5655.991487210216,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -7541.321982946955,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3142.217492894565,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -342.7873628612252,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.394709780272118,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -7.176531019523646,
    0,
    0,
    0,
    0,
    124.3932043384099,
    0,
    186.5898065076148,
    0,
    0,
    0,
    0,
    0,
    0,
    -342.0813119306271,
    0,
    -3420.813119306271,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    12314.92722950258,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    342.0813119306271,
    0,
    -12314.92722950258,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -124.3932043384099,
    0,
    3420.813119306271,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    7.176531019523646,
    0,
    -186.5898065076148,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -58.09962199636542,
    0,
    0,
    0,
    0,
    0,
    0,
    813.3947079491158,
    0,
    464.7969759709233,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1568.689793901866,
    0,
    -6971.95463956385,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -697.195463956385,
    0,
    19521.47299077878,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1452.490549909135,
    0,
    -13943.9092791277,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -285.2163261639757,
    0,
    2323.984879854617,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.281783817851402,
    0,
    -42.25427054281121,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.281783817851402,
    0,
    0,
    0,
    0,
    -52.81783817851402,
    0,
    -253.5256232568673,
    0,
    0,
    0,
    0,
    0,
    0,
    11.61992439927308,
    0,
    2788.78185582554,
    0,
    929.5939519418467,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    139.439092791277,
    0,
    -3346.538226990648,
    0,
    -11155.12742330216,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    11.61992439927308,
    0,
    -3346.538226990648,
    0,
    23425.76758893454,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -52.81783817851402,
    0,
    2788.78185582554,
    0,
    -11155.12742330216,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.281783817851402,
    0,
    -253.5256232568673,
    0,
    929.5939519418467,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    52.07313853625346,
    0,
    0,
    0,
    0,
    0,
    0,
    -381.8696825991921,
    0,
    -763.7393651983841,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -190.934841299596,
    0,
    6364.494709986534,
    0,
    1527.478730396768,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    763.7393651983841,
    0,
    -3564.117037592459,
    0,
    -14256.46815036984,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    318.2247354993267,
    0,
    -7637.393651983841,
    0,
    21384.70222555475,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -196.720745581402,
    0,
    2970.097531327049,
    0,
    -6109.914921587073,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.78590428180594,
    0,
    -84.85992946648712,
    0,
    169.7198589329742,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.940231044985851,
    0,
    0,
    0,
    0,
    15.7609241799434,
    0,
    260.0552489690662,
    0,
    0,
    0,
    0,
    0,
    0,
    43.34254149484436,
    0,
    -1300.276244845331,
    0,
    -1733.701659793774,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1560.331493814397,
    0,
    10402.20995876265,
    0,
    2080.441991752529,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -43.34254149484436,
    0,
    1560.331493814397,
    0,
    0,
    0,
    -14563.0939422677,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -15.7609241799434,
    0,
    1300.276244845331,
    0,
    -10402.20995876265,
    0,
    14563.0939422677,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.940231044985851,
    0,
    -260.0552489690662,
    0,
    1733.701659793774,
    0,
    -2080.441991752529,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -42.78485868831644,
    0,
    0,
    0,
    0,
    0,
    0,
    85.56971737663287,
    0,
    855.6971737663287,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    385.0637281948479,
    0,
    -2567.091521298986,
    0,
    -3080.509825558783,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    220.0364161113417,
    0,
    -5134.183042597972,
    0,
    12322.03930223513,
    0,
    2347.055105187644,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -152.8030667439873,
    0,
    733.4547203711389,
    0,
    6161.019651117567,
    0,
    -11735.27552593822,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -110.0182080556708,
    0,
    2322.606614508606,
    0,
    -8801.456644453667,
    0,
    7041.165315562933,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.112122669759491,
    0,
    -122.2424533951898,
    0,
    440.0728322226833,
    0,
    -335.2935864553778,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.829363009969403,
    0,
    0,
    0,
    0,
    1.886242006646268,
    0,
    -226.3490407975522,
    0,
    0,
    0,
    0,
    0,
    0,
    -17.91929906313955,
    0,
    75.44968026585074,
    0,
    2037.14136717797,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -33.95235611963283,
    0,
    1358.094244785313,
    0,
    -2716.188489570627,
    0,
    -4345.901583313003,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -17.91929906313955,
    0,
    1358.094244785313,
    0,
    -9506.659713497193,
    0,
    10140.43702773034,
    0,
    2172.950791656501,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.886242006646268,
    0,
    75.44968026585074,
    0,
    -2716.188489570627,
    0,
    10140.43702773034,
    0,
    -7243.169305521671,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.829363009969403,
    0,
    -226.3490407975522,
    0,
    2037.14136717797,
    0,
    -4345.901583313003,
    0,
    2172.950791656501,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    31.633240116575,
    0,
    0,
    0,
    0,
    0,
    0,
    63.26648023315,
    0,
    -759.1977627977999,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -56.939832209835,
    0,
    -759.1977627977999,
    0,
    3644.14926142944,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -227.75932883934,
    0,
    2125.75373583384,
    0,
    0,
    0,
    -4858.86568190592,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -183.472792676135,
    0,
    3340.47015631032,
    0,
    -10203.61793200243,
    0,
    4858.86568190592,
    0,
    1619.62189396864,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -37.95988813989,
    0,
    1062.87686791692,
    0,
    -5830.638818287104,
    0,
    8745.958227430655,
    0,
    -3239.24378793728,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.326648023315,
    0,
    -151.83955255956,
    0,
    728.8298522858879,
    0,
    -971.7731363811839,
    0,
    323.924378793728,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.835933153488193,
    0,
    0,
    0,
    0,
    -7.343732613952774,
    0,
    165.2339838139374,
    0,
    0,
    0,
    0,
    0,
    0,
    -9.179665767440967,
    0,
    495.7019514418122,
    0,
    -1762.495827348666,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    330.4679676278748,
    0,
    -3524.991654697331,
    0,
    4934.988316576264,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    9.179665767440967,
    0,
    -330.4679676278748,
    0,
    0,
    0,
    4934.988316576264,
    0,
    -4229.989985636798,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    7.343732613952774,
    0,
    -495.7019514418122,
    0,
    3524.991654697331,
    0,
    -4934.988316576264,
    0,
    0,
    0,
    939.997774585955,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.835933153488193,
    0,
    -165.2339838139374,
    0,
    1762.495827348666,
    0,
    -4934.988316576264,
    0,
    4229.989985636798,
    0,
    -939.997774585955,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -19.37540204348254,
    0,
    0,
    0,
    0,
    0,
    0,
    -90.41854286958517,
    0,
    516.677387826201,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -161.4616836956878,
    0,
    1894.483755362737,
    0,
    -2893.393371826726,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -129.1693469565502,
    0,
    2411.161143188938,
    0,
    -7715.715658204601,
    0,
    4960.10292313153,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -32.29233673913756,
    0,
    1033.354775652402,
    0,
    -5786.786743653451,
    0,
    8266.838205219216,
    0,
    -2755.612735073072,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    12.91693469565502,
    0,
    -172.225795942067,
    0,
    0,
    0,
    1653.367641043843,
    0,
    -1837.075156715381,
    0,
    400.8163978288105,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.458467347827512,
    0,
    -172.225795942067,
    0,
    964.4644572755752,
    0,
    -1653.367641043843,
    0,
    918.5375783576907,
    0,
    -133.6054659429368,
    0,
    0,
    0,
    0,
    0.9043663200508067,
    0,
    0,
    0,
    0,
    5.42619792030484,
    0,
    -86.81916672487744,
    0,
    0,
    0,
    0,
    0,
    0,
    13.5654948007621,
    0,
    -434.0958336243872,
    0,
    1012.890278456903,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    18.08732640101613,
    0,
    -868.1916672487744,
    0,
    4051.561113827614,
    0,
    -3241.248891062091,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    13.5654948007621,
    0,
    -868.1916672487744,
    0,
    6077.341670741421,
    0,
    -9723.746673186273,
    0,
    3472.766668995098,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.42619792030484,
    0,
    -434.0958336243872,
    0,
    4051.561113827614,
    0,
    -9723.746673186273,
    0,
    6945.533337990195,
    0,
    -1234.761482309368,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.9043663200508067,
    0,
    -86.81916672487744,
    0,
    1012.890278456903,
    0,
    -3241.248891062091,
    0,
    3472.766668995098,
    0,
    -1234.761482309368,
    0,
    112.2510438463062,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.521478277491721,
    0,
    0,
    0,
    0,
    0,
    0,
    39.12886966495032,
    0,
    -182.6013917697682,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    97.82217416237581,
    0,
    -913.0069588488409,
    0,
    1095.608350618609,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    130.4295655498344,
    0,
    -1826.013917697682,
    0,
    4382.433402474436,
    0,
    -2086.873048797351,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    97.82217416237581,
    0,
    -1826.013917697682,
    0,
    6573.650103711654,
    0,
    -6260.619146392052,
    0,
    1391.248699198234,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    39.12886966495032,
    0,
    -913.0069588488409,
    0,
    4382.433402474436,
    0,
    -6260.619146392052,
    0,
    2782.497398396467,
    0,
    -303.5451707341601,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.521478277491721,
    0,
    -182.6013917697682,
    0,
    1095.608350618609,
    0,
    -2086.873048797351,
    0,
    1391.248699198234,
    0,
    -303.5451707341601,
    0,
    15.56641901200821,
    0,
    -0.3182155563368222,
    0,
    0,
    -2.227508894357756,
    0,
    31.18512452100858,
    0,
    0,
    0,
    0,
    -6.682526683073267,
    0,
    187.1107471260515,
    0,
    -374.2214942521029,
    0,
    0,
    0,
    0,
    0,
    0,
    -11.13754447178878,
    0,
    467.7768678151287,
    0,
    -1871.107471260515,
    0,
    1247.404980840343,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -11.13754447178878,
    0,
    623.7024904201716,
    0,
    -3742.214942521029,
    0,
    4989.619923361373,
    0,
    -1425.605692388964,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.682526683073267,
    0,
    467.7768678151287,
    0,
    -3742.214942521029,
    0,
    7484.429885042059,
    0,
    -4276.817077166891,
    0,
    570.2422769555854,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.227508894357756,
    0,
    187.1107471260515,
    0,
    -1871.107471260515,
    0,
    4989.619923361373,
    0,
    -4276.817077166891,
    0,
    1140.484553911171,
    0,
    -69.12027599461642,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.3182155563368222,
    0,
    31.18512452100858,
    0,
    -374.2214942521029,
    0,
    1247.404980840343,
    0,
    -1425.605692388964,
    0,
    570.2422769555854,
    0,
    -69.12027599461642,
    0,
    1.519126944936625,
    0,
    0,
    6.521478277491721,
    0,
    0,
    0,
    0,
    39.12886966495032,
    0,
    -182.6013917697682,
    0,
    0,
    0,
    0,
    0,
    0,
    97.82217416237581,
    0,
    -913.0069588488409,
    0,
    1095.608350618609,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    130.4295655498344,
    0,
    -1826.013917697682,
    0,
    4382.433402474436,
    0,
    -2086.873048797351,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    97.82217416237581,
    0,
    -1826.013917697682,
    0,
    6573.650103711654,
    0,
    -6260.619146392052,
    0,
    1391.248699198234,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    39.12886966495032,
    0,
    -913.0069588488409,
    0,
    4382.433402474436,
    0,
    -6260.619146392052,
    0,
    2782.497398396467,
    0,
    -303.5451707341601,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.521478277491721,
    0,
    -182.6013917697682,
    0,
    1095.608350618609,
    0,
    -2086.873048797351,
    0,
    1391.248699198234,
    0,
    -303.5451707341601,
    0,
    15.56641901200821,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4521831600254033,
    0,
    0,
    2.260915800127017,
    0,
    -43.40958336243872,
    0,
    0,
    0,
    0,
    4.06964844022863,
    0,
    -173.6383334497549,
    0,
    506.4451392284517,
    0,
    0,
    0,
    0,
    0,
    0,
    2.260915800127017,
    0,
    -217.0479168121936,
    0,
    1519.335417685355,
    0,
    -1620.624445531046,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.260915800127017,
    0,
    0,
    0,
    1012.890278456903,
    0,
    -3241.248891062091,
    0,
    1736.383334497549,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.06964844022863,
    0,
    217.0479168121936,
    0,
    -1012.890278456903,
    0,
    0,
    0,
    1736.383334497549,
    0,
    -617.380741154684,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.260915800127017,
    0,
    173.6383334497549,
    0,
    -1519.335417685355,
    0,
    3241.248891062091,
    0,
    -1736.383334497549,
    0,
    0,
    0,
    56.12552192315309,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4521831600254033,
    0,
    43.40958336243872,
    0,
    -506.4451392284517,
    0,
    1620.624445531046,
    0,
    -1736.383334497549,
    0,
    617.380741154684,
    0,
    -56.12552192315309,
    0,
    0,
    0,
    0,
    -6.458467347827512,
    0,
    0,
    0,
    0,
    -12.91693469565502,
    0,
    172.225795942067,
    0,
    0,
    0,
    0,
    0,
    0,
    32.29233673913756,
    0,
    172.225795942067,
    0,
    -964.4644572755752,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    129.1693469565502,
    0,
    -1033.354775652402,
    0,
    0,
    0,
    1653.367641043843,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    161.4616836956878,
    0,
    -2411.161143188938,
    0,
    5786.786743653451,
    0,
    -1653.367641043843,
    0,
    -918.5375783576907,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    90.41854286958517,
    0,
    -1894.483755362737,
    0,
    7715.715658204601,
    0,
    -8266.838205219216,
    0,
    1837.075156715381,
    0,
    133.6054659429368,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    19.37540204348254,
    0,
    -516.677387826201,
    0,
    2893.393371826726,
    0,
    -4960.10292313153,
    0,
    2755.612735073072,
    0,
    -400.8163978288105,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4589832883720484,
    0,
    0,
    0.4589832883720484,
    0,
    41.30849595348435,
    0,
    0,
    0,
    0,
    8.720682479068919,
    0,
    -82.6169919069687,
    0,
    -440.6239568371664,
    0,
    0,
    0,
    0,
    0,
    0,
    20.65424797674218,
    0,
    -702.244431209234,
    0,
    1321.871870511499,
    0,
    1233.747079144066,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    20.65424797674218,
    0,
    -1156.637886697562,
    0,
    6168.73539572033,
    0,
    -4934.988316576264,
    0,
    -1057.497496409199,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    8.720682479068919,
    0,
    -702.244431209234,
    0,
    6168.73539572033,
    0,
    -12337.47079144066,
    0,
    5287.487482045997,
    0,
    234.9994436464888,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4589832883720484,
    0,
    -82.6169919069687,
    0,
    1321.871870511499,
    0,
    -4934.988316576264,
    0,
    5287.487482045997,
    0,
    -1409.996661878933,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4589832883720484,
    0,
    41.30849595348435,
    0,
    -440.6239568371664,
    0,
    1233.747079144066,
    0,
    -1057.497496409199,
    0,
    234.9994436464888,
    0,
    0,
    0,
    0,
    0,
    0,
    6.326648023315,
    0,
    0,
    0,
    0,
    -37.95988813989,
    0,
    -151.83955255956,
    0,
    0,
    0,
    0,
    0,
    0,
    -183.472792676135,
    0,
    1062.87686791692,
    0,
    728.8298522858879,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -227.75932883934,
    0,
    3340.47015631032,
    0,
    -5830.638818287104,
    0,
    -971.7731363811839,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -56.939832209835,
    0,
    2125.75373583384,
    0,
    -10203.61793200243,
    0,
    8745.958227430655,
    0,
    323.924378793728,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    63.26648023315,
    0,
    -759.1977627977999,
    0,
    0,
    0,
    4858.86568190592,
    0,
    -3239.24378793728,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    31.633240116575,
    0,
    -759.1977627977999,
    0,
    3644.14926142944,
    0,
    -4858.86568190592,
    0,
    1619.62189396864,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4715605016615671,
    0,
    0,
    -5.187165518277238,
    0,
    -37.72484013292537,
    0,
    0,
    0,
    0,
    -18.39085956480112,
    0,
    452.6980815951044,
    0,
    339.5235611963283,
    0,
    0,
    0,
    0,
    0,
    0,
    -12.73213354486231,
    0,
    1018.570683588985,
    0,
    -4413.806295552268,
    0,
    -724.3169305521671,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    12.73213354486231,
    0,
    0,
    0,
    -4753.329856748596,
    0,
    10140.43702773034,
    0,
    362.1584652760835,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    18.39085956480112,
    0,
    -1018.570683588985,
    0,
    4753.329856748596,
    0,
    0,
    0,
    -5432.376979141253,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.187165518277238,
    0,
    -452.6980815951044,
    0,
    4413.806295552268,
    0,
    -10140.43702773034,
    0,
    5432.376979141253,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4715605016615671,
    0,
    37.72484013292537,
    0,
    -339.5235611963283,
    0,
    724.3169305521671,
    0,
    -362.1584652760835,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.112122669759491,
    0,
    0,
    0,
    0,
    110.0182080556708,
    0,
    122.2424533951898,
    0,
    0,
    0,
    0,
    0,
    0,
    152.8030667439873,
    0,
    -2322.606614508606,
    0,
    -440.0728322226833,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -220.0364161113417,
    0,
    -733.4547203711389,
    0,
    8801.456644453667,
    0,
    335.2935864553778,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -385.0637281948479,
    0,
    5134.183042597972,
    0,
    -6161.019651117567,
    0,
    -7041.165315562933,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -85.56971737663287,
    0,
    2567.091521298986,
    0,
    -12322.03930223513,
    0,
    11735.27552593822,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    42.78485868831644,
    0,
    -855.6971737663287,
    0,
    3080.509825558783,
    0,
    -2347.055105187644,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4925288806232314,
    0,
    0,
    12.31322201558078,
    0,
    32.50690612113327,
    0,
    0,
    0,
    0,
    5.417817686855545,
    0,
    -845.179559149465,
    0,
    -216.7127074742218,
    0,
    0,
    0,
    0,
    0,
    0,
    -48.7603591816999,
    0,
    487.603591816999,
    0,
    5851.243101803988,
    0,
    260.0552489690662,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -48.7603591816999,
    0,
    2730.580114175195,
    0,
    -9101.933713917315,
    0,
    -7281.546971133852,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.417817686855545,
    0,
    487.603591816999,
    0,
    -9101.933713917315,
    0,
    18203.86742783463,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    12.31322201558078,
    0,
    -845.179559149465,
    0,
    5851.243101803988,
    0,
    -7281.546971133852,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4925288806232314,
    0,
    32.50690612113327,
    0,
    -216.7127074742218,
    0,
    260.0552489690662,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.78590428180594,
    0,
    0,
    0,
    0,
    -196.720745581402,
    0,
    -84.85992946648712,
    0,
    0,
    0,
    0,
    0,
    0,
    318.2247354993267,
    0,
    2970.097531327049,
    0,
    169.7198589329742,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    763.7393651983841,
    0,
    -7637.393651983841,
    0,
    -6109.914921587073,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -190.934841299596,
    0,
    -3564.117037592459,
    0,
    21384.70222555475,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -381.8696825991921,
    0,
    6364.494709986534,
    0,
    -14256.46815036984,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    52.07313853625346,
    0,
    -763.7393651983841,
    0,
    1527.478730396768,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.5281783817851402,
    0,
    0,
    -22.71167041676103,
    0,
    -25.35256232568673,
    0,
    0,
    0,
    0,
    63.90958419600196,
    0,
    1115.512742330216,
    0,
    92.95939519418467,
    0,
    0,
    0,
    0,
    0,
    0,
    87.14943299454813,
    0,
    -4183.17278373831,
    0,
    -4183.17278373831,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -87.14943299454813,
    0,
    0,
    0,
    19521.47299077878,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -63.90958419600196,
    0,
    4183.17278373831,
    0,
    -19521.47299077878,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    22.71167041676103,
    0,
    -1115.512742330216,
    0,
    4183.17278373831,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5281783817851402,
    0,
    25.35256232568673,
    0,
    -92.95939519418467,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.281783817851402,
    0,
    0,
    0,
    0,
    285.2163261639757,
    0,
    42.25427054281121,
    0,
    0,
    0,
    0,
    0,
    0,
    -1452.490549909135,
    0,
    -2323.984879854617,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    697.195463956385,
    0,
    13943.9092791277,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1568.689793901866,
    0,
    -19521.47299077878,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -813.3947079491158,
    0,
    6971.95463956385,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    58.09962199636542,
    0,
    -464.7969759709233,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5980442516269705,
    0,
    0,
    38.87287635575308,
    0,
    15.54915054230123,
    0,
    0,
    0,
    0,
    -256.5609839479703,
    0,
    -1026.243935791881,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    256.5609839479703,
    0,
    7696.82951843911,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    256.5609839479703,
    0,
    -14367.41510108634,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -256.5609839479703,
    0,
    7696.82951843911,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    38.87287635575308,
    0,
    -1026.243935791881,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5980442516269705,
    0,
    15.54915054230123,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.394709780272118,
    0,
    0,
    0,
    0,
    -342.7873628612252,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3142.217492894565,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -7541.321982946955,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5655.991487210216,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1256.886997157826,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    57.13122714353754,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.830522083064524,
    0,
    0,
    -75.57750955887168,
    0,
    0,
    0,
    0,
    0,
    0,
    831.3526051475885,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2494.057815442766,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2494.057815442766,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -831.3526051475885,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    75.57750955887168,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.830522083064524,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    // l = 15
    0,
    12.66375976286059,
    0,
    0,
    0,
    0,
    -384.1340461401045,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2535.28470452469,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5432.752938267193,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4225.47450754115,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1152.402138420314,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    88.64631834002413,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.8442506508573726,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    64.73811759282017,
    0,
    0,
    0,
    0,
    0,
    0,
    -1683.191057413324,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    9257.550815773284,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -15870.0871127542,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    9257.550815773284,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1683.191057413324,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    64.73811759282017,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -7.893350484654575,
    0,
    0,
    0,
    0,
    165.7603601777461,
    0,
    221.0138135703281,
    0,
    0,
    0,
    0,
    0,
    0,
    -607.7879873184023,
    0,
    -4862.303898547218,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    260.480565993601,
    0,
    21880.36754346248,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    607.7879873184023,
    0,
    -29173.82339128331,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -386.7741737480742,
    0,
    12155.75974636805,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    46.75292210141556,
    0,
    -1326.082881421969,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.6071808065118904,
    0,
    17.00106258233293,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -66.77884815276176,
    0,
    0,
    0,
    0,
    0,
    0,
    1157.500034647871,
    0,
    578.7500173239353,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3183.125095281644,
    0,
    -10610.41698427215,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    38197.50114337973,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3183.125095281644,
    0,
    -38197.50114337973,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1157.500034647871,
    0,
    10610.41698427215,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    66.77884815276176,
    0,
    -578.7500173239353,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5.890314521234328,
    0,
    0,
    0,
    0,
    -76.57408877604626,
    0,
    -306.2963551041851,
    0,
    0,
    0,
    0,
    0,
    0,
    76.57408877604626,
    0,
    4288.148971458591,
    0,
    1225.18542041674,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    229.7222663281388,
    0,
    -8270.001587812996,
    0,
    -18377.7813062511,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -76.57408877604626,
    0,
    -3675.556261250221,
    0,
    51457.78765750309,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -118.3417735629806,
    0,
    7657.408877604626,
    0,
    -36755.56261250221,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    28.38060632958358,
    0,
    -1503.636652329636,
    0,
    6125.927102083701,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.5354831382940298,
    0,
    27.84512319128955,
    0,
    -111.3804927651582,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    61.05447148378159,
    0,
    0,
    0,
    0,
    0,
    0,
    -610.5447148378159,
    0,
    -976.8715437405055,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    134.3198372643195,
    0,
    10745.58698114556,
    0,
    2149.117396229112,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1611.838047171834,
    0,
    -12894.70437737467,
    0,
    -25789.40875474934,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    134.3198372643195,
    0,
    -12894.70437737467,
    0,
    54157.75838497362,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -610.5447148378159,
    0,
    10745.58698114556,
    0,
    -25789.40875474934,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    61.05447148378159,
    0,
    -976.8715437405055,
    0,
    2149.117396229112,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.486569049517132,
    0,
    0,
    0,
    0,
    28.4149373136085,
    0,
    323.0329715652335,
    0,
    0,
    0,
    0,
    0,
    0,
    49.35225954468845,
    0,
    -2368.908458145046,
    0,
    -2368.908458145046,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -49.35225954468845,
    0,
    -1184.454229072523,
    0,
    19740.90381787538,
    0,
    3158.544610860061,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -93.22093469552263,
    0,
    4737.816916290091,
    0,
    -11054.90613801021,
    0,
    -29479.74970136057,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -10.46866111553997,
    0,
    1974.090381787538,
    0,
    -23689.08458145046,
    0,
    44219.62455204085,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    16.45075318156282,
    0,
    -1220.34678146866,
    0,
    9212.421781675177,
    0,
    -12634.17844344024,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4985076721685702,
    0,
    35.89255239613705,
    0,
    -263.2120509050051,
    0,
    350.9494012066734,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -51.69118335186266,
    0,
    0,
    0,
    0,
    0,
    0,
    206.7647334074506,
    0,
    1137.206033740979,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    568.6030168704893,
    0,
    -5686.030168704893,
    0,
    -4548.824134963914,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6823.236202445871,
    0,
    27292.94480978348,
    0,
    3898.992115683355,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -568.6030168704893,
    0,
    6823.236202445871,
    0,
    0,
    0,
    -27292.94480978348,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -206.7647334074506,
    0,
    5686.030168704893,
    0,
    -27292.94480978348,
    0,
    27292.94480978348,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    51.69118335186266,
    0,
    -1137.206033740979,
    0,
    4548.824134963914,
    0,
    -3898.992115683355,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.334384020345036,
    0,
    0,
    0,
    0,
    -3.334384020345036,
    0,
    -293.4257937903632,
    0,
    0,
    0,
    0,
    0,
    0,
    -36.6782242237954,
    0,
    586.8515875807264,
    0,
    2934.257937903632,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -47.15771685916551,
    0,
    2640.832144113269,
    0,
    -8802.773813710896,
    0,
    -7042.219050968717,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.239746317685057,
    0,
    1509.046939493296,
    0,
    -17605.54762742179,
    0,
    28168.87620387487,
    0,
    4024.125171982124,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    20.48264469640522,
    0,
    -1047.949263537011,
    0,
    2515.078232488827,
    0,
    14084.43810193743,
    0,
    -20120.62585991062,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    8.097789763695088,
    0,
    -754.5234697466482,
    0,
    7964.414402881287,
    0,
    -20120.62585991062,
    0,
    12072.37551594637,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4763405743350052,
    0,
    41.91797054148046,
    0,
    -419.1797054148046,
    0,
    1006.031292995531,
    0,
    -574.8750245688748,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    40.21623606427654,
    0,
    0,
    0,
    0,
    0,
    0,
    26.81082404285103,
    0,
    -1072.432961714041,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -254.7028284070847,
    0,
    357.4776539046803,
    0,
    5791.137993255822,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -482.5948327713185,
    0,
    6434.597770284246,
    0,
    -7721.517324341095,
    0,
    -8824.591227818395,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -254.7028284070847,
    0,
    6434.597770284246,
    0,
    -27025.31063519383,
    0,
    20590.71286490959,
    0,
    3431.785477484931,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    26.81082404285103,
    0,
    357.4776539046803,
    0,
    -7721.517324341095,
    0,
    20590.71286490959,
    0,
    -11439.28492494977,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    40.21623606427654,
    0,
    -1072.432961714041,
    0,
    5791.137993255822,
    0,
    -8824.591227818395,
    0,
    3431.785477484931,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.312653286194929,
    0,
    0,
    0,
    0,
    -6.937959858584787,
    0,
    231.2653286194929,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4625306572389858,
    0,
    462.5306572389858,
    0,
    -2775.183943433915,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    20.81387957575436,
    0,
    -416.2775915150872,
    0,
    -2775.183943433915,
    0,
    8880.588618988527,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    30.06449272053408,
    0,
    -1665.110366060349,
    0,
    7770.515041614961,
    0,
    0,
    0,
    -8880.588618988527,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    16.1885730033645,
    0,
    -1341.338905993059,
    0,
    12210.80935110922,
    0,
    -24865.64813316788,
    0,
    8880.588618988527,
    0,
    2368.156965063607,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2.312653286194929,
    0,
    -277.5183943433915,
    0,
    3885.257520807481,
    0,
    -14208.94179038164,
    0,
    15985.05951417935,
    0,
    -4736.313930127214,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4625306572389858,
    0,
    46.25306572389858,
    0,
    -555.0367886867829,
    0,
    1776.117723797705,
    0,
    -1776.117723797705,
    0,
    473.6313930127214,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -27.44175328360903,
    0,
    0,
    0,
    0,
    0,
    0,
    -109.7670131344361,
    0,
    823.252598508271,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -137.2087664180452,
    0,
    2469.757795524813,
    0,
    -5268.816630452934,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1646.505197016542,
    0,
    -10537.63326090587,
    0,
    10537.63326090587,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    137.2087664180452,
    0,
    -1646.505197016542,
    0,
    0,
    0,
    10537.63326090587,
    0,
    -7025.088840603912,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    109.7670131344361,
    0,
    -2469.757795524813,
    0,
    10537.63326090587,
    0,
    -10537.63326090587,
    0,
    0,
    0,
    1277.288880109802,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    27.44175328360903,
    0,
    -823.252598508271,
    0,
    5268.816630452934,
    0,
    -10537.63326090587,
    0,
    7025.088840603912,
    0,
    -1277.288880109802,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.363030880952603,
    0,
    0,
    0,
    0,
    7.723841658731416,
    0,
    -147.2073351428811,
    0,
    0,
    0,
    0,
    0,
    0,
    17.71940145238384,
    0,
    -686.9675640001119,
    0,
    1962.764468571748,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    20.44546321428904,
    0,
    -1226.727792857343,
    0,
    7196.803051429743,
    0,
    -7327.654016001193,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    11.35859067460502,
    0,
    -981.3822342858741,
    0,
    9159.567520001492,
    0,
    -19540.41070933652,
    0,
    9421.269449144391,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.363030880952603,
    0,
    -245.3455585714685,
    0,
    3925.528937143496,
    0,
    -14655.30803200239,
    0,
    15702.11574857399,
    0,
    -4187.230866286396,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.363030880952603,
    0,
    98.13822342858741,
    0,
    -654.2548228572494,
    0,
    0,
    0,
    3140.423149714797,
    0,
    -2791.487244190931,
    0,
    507.543135307442,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.454343626984201,
    0,
    49.0691117142937,
    0,
    -654.2548228572494,
    0,
    2442.551338667064,
    0,
    -3140.423149714797,
    0,
    1395.743622095465,
    0,
    -169.1810451024807,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    13.90024211921377,
    0,
    0,
    0,
    0,
    0,
    0,
    83.40145271528264,
    0,
    -444.8077478148407,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    208.5036317882066,
    0,
    -2224.038739074204,
    0,
    3113.654234703885,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    278.0048423842755,
    0,
    -4448.077478148407,
    0,
    12454.61693881554,
    0,
    -7116.923965037452,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    208.5036317882066,
    0,
    -4448.077478148407,
    0,
    18681.92540822331,
    0,
    -21350.77189511235,
    0,
    5930.769970864543,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    83.40145271528264,
    0,
    -2224.038739074204,
    0,
    12454.61693881554,
    0,
    -21350.77189511235,
    0,
    11861.53994172909,
    0,
    -1725.31490061514,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    13.90024211921377,
    0,
    -444.8077478148407,
    0,
    3113.654234703885,
    0,
    -7116.923965037452,
    0,
    5930.769970864543,
    0,
    -1725.31490061514,
    0,
    132.7165308165492,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4505094349975498,
    0,
    0,
    0,
    0,
    -3.153566044982849,
    0,
    50.45705671972558,
    0,
    0,
    0,
    0,
    0,
    0,
    -9.460698134948546,
    0,
    302.7423403183535,
    0,
    -706.3987940761581,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -15.76783022491424,
    0,
    756.8558507958837,
    0,
    -3531.99397038079,
    0,
    2825.595176304632,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -15.76783022491424,
    0,
    1009.141134394512,
    0,
    -7063.987940761581,
    0,
    11302.38070521853,
    0,
    -4036.564537578046,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -9.460698134948546,
    0,
    756.8558507958837,
    0,
    -7063.987940761581,
    0,
    16953.57105782779,
    0,
    -12109.69361273414,
    0,
    2152.834420041625,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.153566044982849,
    0,
    302.7423403183535,
    0,
    -3531.99397038079,
    0,
    11302.38070521853,
    0,
    -12109.69361273414,
    0,
    4305.668840083249,
    0,
    -391.4244400075681,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4505094349975498,
    0,
    50.45705671972558,
    0,
    -706.3987940761581,
    0,
    2825.595176304632,
    0,
    -4036.564537578046,
    0,
    2152.834420041625,
    0,
    -391.4244400075681,
    0,
    17.20546989044255,
    0,
    0,
    0,
    -4.935083598341307,
    0,
    0,
    0,
    0,
    -34.54558518838915,
    0,
    161.2127308791494,
    0,
    0,
    0,
    0,
    0,
    0,
    -103.6367555651675,
    0,
    967.2763852748962,
    0,
    -1160.731662329875,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -172.7279259419458,
    0,
    2418.190963187241,
    0,
    -5803.658311649377,
    0,
    2763.646815071132,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -172.7279259419458,
    0,
    3224.254617582987,
    0,
    -11607.31662329875,
    0,
    11054.58726028453,
    0,
    -2456.574946729895,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -103.6367555651675,
    0,
    2418.190963187241,
    0,
    -11607.31662329875,
    0,
    16581.88089042679,
    0,
    -7369.724840189685,
    0,
    803.9699825661475,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -34.54558518838915,
    0,
    967.2763852748962,
    0,
    -5803.658311649377,
    0,
    11054.58726028453,
    0,
    -7369.724840189685,
    0,
    1607.939965132295,
    0,
    -82.4584597503741,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.935083598341307,
    0,
    161.2127308791494,
    0,
    -1160.731662329875,
    0,
    2763.646815071132,
    0,
    -2456.574946729895,
    0,
    803.9699825661475,
    0,
    -82.4584597503741,
    0,
    1.570637328578554,
    -0.4505094349975498,
    0,
    0,
    -3.153566044982849,
    0,
    50.45705671972558,
    0,
    0,
    0,
    0,
    -9.460698134948546,
    0,
    302.7423403183535,
    0,
    -706.3987940761581,
    0,
    0,
    0,
    0,
    0,
    0,
    -15.76783022491424,
    0,
    756.8558507958837,
    0,
    -3531.99397038079,
    0,
    2825.595176304632,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -15.76783022491424,
    0,
    1009.141134394512,
    0,
    -7063.987940761581,
    0,
    11302.38070521853,
    0,
    -4036.564537578046,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -9.460698134948546,
    0,
    756.8558507958837,
    0,
    -7063.987940761581,
    0,
    16953.57105782779,
    0,
    -12109.69361273414,
    0,
    2152.834420041625,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.153566044982849,
    0,
    302.7423403183535,
    0,
    -3531.99397038079,
    0,
    11302.38070521853,
    0,
    -12109.69361273414,
    0,
    4305.668840083249,
    0,
    -391.4244400075681,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4505094349975498,
    0,
    50.45705671972558,
    0,
    -706.3987940761581,
    0,
    2825.595176304632,
    0,
    -4036.564537578046,
    0,
    2152.834420041625,
    0,
    -391.4244400075681,
    0,
    17.20546989044255,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.950121059606886,
    0,
    0,
    0,
    0,
    34.75060529803443,
    0,
    -222.4038739074204,
    0,
    0,
    0,
    0,
    0,
    0,
    62.55108953646198,
    0,
    -889.6154956296814,
    0,
    1556.827117351943,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    34.75060529803443,
    0,
    -1112.019369537102,
    0,
    4670.481352055828,
    0,
    -3558.461982518726,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -34.75060529803443,
    0,
    0,
    0,
    3113.654234703885,
    0,
    -7116.923965037452,
    0,
    2965.384985432271,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -62.55108953646198,
    0,
    1112.019369537102,
    0,
    -3113.654234703885,
    0,
    0,
    0,
    2965.384985432271,
    0,
    -862.6574503075699,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -34.75060529803443,
    0,
    889.6154956296814,
    0,
    -4670.481352055828,
    0,
    7116.923965037452,
    0,
    -2965.384985432271,
    0,
    0,
    0,
    66.35826540827461,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.950121059606886,
    0,
    222.4038739074204,
    0,
    -1556.827117351943,
    0,
    3558.461982518726,
    0,
    -2965.384985432271,
    0,
    862.6574503075699,
    0,
    -66.35826540827461,
    0,
    0,
    0.454343626984201,
    0,
    0,
    1.363030880952603,
    0,
    -49.0691117142937,
    0,
    0,
    0,
    0,
    -1.363030880952603,
    0,
    -98.13822342858741,
    0,
    654.2548228572494,
    0,
    0,
    0,
    0,
    0,
    0,
    -11.35859067460502,
    0,
    245.3455585714685,
    0,
    654.2548228572494,
    0,
    -2442.551338667064,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -20.44546321428904,
    0,
    981.3822342858741,
    0,
    -3925.528937143496,
    0,
    0,
    0,
    3140.423149714797,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -17.71940145238384,
    0,
    1226.727792857343,
    0,
    -9159.567520001492,
    0,
    14655.30803200239,
    0,
    -3140.423149714797,
    0,
    -1395.743622095465,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -7.723841658731416,
    0,
    686.9675640001119,
    0,
    -7196.803051429743,
    0,
    19540.41070933652,
    0,
    -15702.11574857399,
    0,
    2791.487244190931,
    0,
    169.1810451024807,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1.363030880952603,
    0,
    147.2073351428811,
    0,
    -1962.764468571748,
    0,
    7327.654016001193,
    0,
    -9421.269449144391,
    0,
    4187.230866286396,
    0,
    -507.543135307442,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.860438320902258,
    0,
    0,
    0,
    0,
    6.860438320902258,
    0,
    205.8131496270677,
    0,
    0,
    0,
    0,
    0,
    0,
    130.3483280971429,
    0,
    -411.6262992541355,
    0,
    -1317.204157613234,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    308.7197244406016,
    0,
    -3498.823543660152,
    0,
    3951.612472839701,
    0,
    2634.408315226467,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    308.7197244406016,
    0,
    -5762.768189557897,
    0,
    18440.85820658527,
    0,
    -10537.63326090587,
    0,
    -1756.272210150978,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    130.3483280971429,
    0,
    -3498.823543660152,
    0,
    18440.85820658527,
    0,
    -26344.08315226467,
    0,
    8781.36105075489,
    0,
    319.3222200274506,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.860438320902258,
    0,
    -411.6262992541355,
    0,
    3951.612472839701,
    0,
    -10537.63326090587,
    0,
    8781.36105075489,
    0,
    -1915.933320164703,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.860438320902258,
    0,
    205.8131496270677,
    0,
    -1317.204157613234,
    0,
    2634.408315226467,
    0,
    -1756.272210150978,
    0,
    319.3222200274506,
    0,
    0,
    0,
    0,
    -0.4625306572389858,
    0,
    0,
    2.312653286194929,
    0,
    46.25306572389858,
    0,
    0,
    0,
    0,
    16.1885730033645,
    0,
    -277.5183943433915,
    0,
    -555.0367886867829,
    0,
    0,
    0,
    0,
    0,
    0,
    30.06449272053408,
    0,
    -1341.338905993059,
    0,
    3885.257520807481,
    0,
    1776.117723797705,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    20.81387957575436,
    0,
    -1665.110366060349,
    0,
    12210.80935110922,
    0,
    -14208.94179038164,
    0,
    -1776.117723797705,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4625306572389858,
    0,
    -416.2775915150872,
    0,
    7770.515041614961,
    0,
    -24865.64813316788,
    0,
    15985.05951417935,
    0,
    473.6313930127214,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.937959858584787,
    0,
    462.5306572389858,
    0,
    -2775.183943433915,
    0,
    0,
    0,
    8880.588618988527,
    0,
    -4736.313930127214,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2.312653286194929,
    0,
    231.2653286194929,
    0,
    -2775.183943433915,
    0,
    8880.588618988527,
    0,
    -8880.588618988527,
    0,
    2368.156965063607,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.702706010712756,
    0,
    0,
    0,
    0,
    -73.72976611784032,
    0,
    -178.7388269523402,
    0,
    0,
    0,
    0,
    0,
    0,
    -261.4055344177975,
    0,
    2144.865923428082,
    0,
    965.1896655426369,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -180.9730622892444,
    0,
    4825.948327713185,
    0,
    -12547.46565205428,
    0,
    -1470.765204636399,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    180.9730622892444,
    0,
    0,
    0,
    -13512.65531759692,
    0,
    20590.71286490959,
    0,
    571.9642462474885,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    261.4055344177975,
    0,
    -4825.948327713185,
    0,
    13512.65531759692,
    0,
    0,
    0,
    -8579.463693712328,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    73.72976611784032,
    0,
    -2144.865923428082,
    0,
    12547.46565205428,
    0,
    -20590.71286490959,
    0,
    8579.463693712328,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.702706010712756,
    0,
    178.7388269523402,
    0,
    -965.1896655426369,
    0,
    1470.765204636399,
    0,
    -571.9642462474885,
    0,
    0,
    0,
    0,
    0,
    0,
    0.4763405743350052,
    0,
    0,
    -8.097789763695088,
    0,
    -41.91797054148046,
    0,
    0,
    0,
    0,
    -20.48264469640522,
    0,
    754.5234697466482,
    0,
    419.1797054148046,
    0,
    0,
    0,
    0,
    0,
    0,
    5.239746317685057,
    0,
    1047.949263537011,
    0,
    -7964.414402881287,
    0,
    -1006.031292995531,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    47.15771685916551,
    0,
    -1509.046939493296,
    0,
    -2515.078232488827,
    0,
    20120.62585991062,
    0,
    574.8750245688748,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    36.6782242237954,
    0,
    -2640.832144113269,
    0,
    17605.54762742179,
    0,
    -14084.43810193743,
    0,
    -12072.37551594637,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.334384020345036,
    0,
    -586.8515875807264,
    0,
    8802.773813710896,
    0,
    -28168.87620387487,
    0,
    20120.62585991062,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -3.334384020345036,
    0,
    293.4257937903632,
    0,
    -2934.257937903632,
    0,
    7042.219050968717,
    0,
    -4024.125171982124,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.461397918982832,
    0,
    0,
    0,
    0,
    161.5349479745708,
    0,
    142.1507542176223,
    0,
    0,
    0,
    0,
    0,
    0,
    71.07537710881116,
    0,
    -3695.91960965818,
    0,
    -568.6030168704893,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -639.6783939793004,
    0,
    2132.261313264335,
    0,
    15352.28145550321,
    0,
    487.3740144604194,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -639.6783939793004,
    0,
    11940.66335428027,
    0,
    -23881.32670856055,
    0,
    -13646.47240489174,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    71.07537710881116,
    0,
    2132.261313264335,
    0,
    -23881.32670856055,
    0,
    34116.18101222936,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    161.5349479745708,
    0,
    -3695.91960965818,
    0,
    15352.28145550321,
    0,
    -13646.47240489174,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.461397918982832,
    0,
    142.1507542176223,
    0,
    -568.6030168704893,
    0,
    487.3740144604194,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.4985076721685702,
    0,
    0,
    16.45075318156282,
    0,
    35.89255239613705,
    0,
    0,
    0,
    0,
    -10.46866111553997,
    0,
    -1220.34678146866,
    0,
    -263.2120509050051,
    0,
    0,
    0,
    0,
    0,
    0,
    -93.22093469552263,
    0,
    1974.090381787538,
    0,
    9212.421781675177,
    0,
    350.9494012066734,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -49.35225954468845,
    0,
    4737.816916290091,
    0,
    -23689.08458145046,
    0,
    -12634.17844344024,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    49.35225954468845,
    0,
    -1184.454229072523,
    0,
    -11054.90613801021,
    0,
    44219.62455204085,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    28.4149373136085,
    0,
    -2368.908458145046,
    0,
    19740.90381787538,
    0,
    -29479.74970136057,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.486569049517132,
    0,
    323.0329715652335,
    0,
    -2368.908458145046,
    0,
    3158.544610860061,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    6.105447148378159,
    0,
    0,
    0,
    0,
    -262.5342273802609,
    0,
    -97.68715437405055,
    0,
    0,
    0,
    0,
    0,
    0,
    738.7591049537573,
    0,
    4298.234792458224,
    0,
    214.9117396229112,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1007.398779482396,
    0,
    -16118.38047171834,
    0,
    -9671.028283031004,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1007.398779482396,
    0,
    0,
    0,
    45131.46532081135,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -738.7591049537573,
    0,
    16118.38047171834,
    0,
    -45131.46532081135,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    262.5342273802609,
    0,
    -4298.234792458224,
    0,
    9671.028283031004,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -6.105447148378159,
    0,
    97.68715437405055,
    0,
    -214.9117396229112,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.5354831382940298,
    0,
    0,
    -28.38060632958358,
    0,
    -27.84512319128955,
    0,
    0,
    0,
    0,
    118.3417735629806,
    0,
    1503.636652329636,
    0,
    111.3804927651582,
    0,
    0,
    0,
    0,
    0,
    0,
    76.57408877604626,
    0,
    -7657.408877604626,
    0,
    -6125.927102083701,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -229.7222663281388,
    0,
    3675.556261250221,
    0,
    36755.56261250221,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -76.57408877604626,
    0,
    8270.001587812996,
    0,
    -51457.78765750309,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    76.57408877604626,
    0,
    -4288.148971458591,
    0,
    18377.7813062511,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.890314521234328,
    0,
    306.2963551041851,
    0,
    -1225.18542041674,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.564904012730147,
    0,
    0,
    0,
    0,
    361.7187608274595,
    0,
    48.22916811032794,
    0,
    0,
    0,
    0,
    0,
    0,
    -2387.343821461233,
    0,
    -3183.125095281644,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2387.343821461233,
    0,
    23873.43821461233,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2387.343821461233,
    0,
    -44563.75133394302,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2387.343821461233,
    0,
    23873.43821461233,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    361.7187608274595,
    0,
    -3183.125095281644,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -5.564904012730147,
    0,
    48.22916811032794,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -0.6071808065118904,
    0,
    0,
    46.75292210141556,
    0,
    17.00106258233293,
    0,
    0,
    0,
    0,
    -386.7741737480742,
    0,
    -1326.082881421969,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    607.7879873184023,
    0,
    12155.75974636805,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    260.480565993601,
    0,
    -29173.82339128331,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -607.7879873184023,
    0,
    21880.36754346248,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    165.7603601777461,
    0,
    -4862.303898547218,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -7.893350484654575,
    0,
    221.0138135703281,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4.624151256630012,
    0,
    0,
    0,
    0,
    -420.7977643533311,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    4628.775407886642,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -13886.32622365993,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    13886.32622365993,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4628.775407886642,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    420.7977643533311,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4.624151256630012,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.8442506508573726,
    0,
    0,
    -88.64631834002413,
    0,
    0,
    0,
    0,
    0,
    0,
    1152.402138420314,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -4225.47450754115,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    5432.752938267193,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -2535.28470452469,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    384.1340461401045,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -12.66375976286059,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
};

// [*] = n(n+1)(n+2)(n+3)/4+(n+1)(n+2)(n+3)/6
static struct cart2sp_t g_c2s[] = {
    {g_trans_cart2sph, NULL, NULL, NULL},
    {g_trans_cart2sph + 1, NULL, NULL, NULL},
    {g_trans_cart2sph + 10, NULL, NULL, NULL},
    {g_trans_cart2sph + 40, NULL, NULL, NULL},
    {g_trans_cart2sph + 110, NULL, NULL, NULL},
    {g_trans_cart2sph + 245, NULL, NULL, NULL},
    {g_trans_cart2sph + 476, NULL, NULL, NULL},
    {g_trans_cart2sph + 840, NULL, NULL, NULL},
    {g_trans_cart2sph + 1380, NULL, NULL, NULL},
    {g_trans_cart2sph + 2145, NULL, NULL, NULL},
    {g_trans_cart2sph + 3190, NULL, NULL, NULL},
    {g_trans_cart2sph + 4576, NULL, NULL, NULL},
    {g_trans_cart2sph + 6370, NULL, NULL, NULL},
    {g_trans_cart2sph + 8645, NULL, NULL, NULL, NULL},
    {g_trans_cart2sph + 11480, NULL, NULL, NULL, NULL},
    {g_trans_cart2sph + 14960, NULL, NULL, NULL, NULL},
};

void dgemm_NN1(int m, int n, int k,
               double *a, double *b, double *c, int ldc)
{
    int i, j, kp;
    double bi;
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < m; i++)
        {
            c[i + ldc * j] = 0;
        }
        for (kp = 0; kp < k; kp++)
        {
            bi = b[kp + k * j];
#pragma ivdep
            for (i = 0; i < m; i++)
            {
                c[i + ldc * j] += a[i + m * kp] * bi;
            }
        }
    }
}

void dgemm_NN(int m, int n, int k,
              double *a, double *b, double *c)
{
    dgemm_NN1(m, n, k, a, b, c, m);
}

void dgemm_TN(int m, int n, int k,
              double *a, double *b, double *c)
{
    int i, j, kp;
    double ci;
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < m; i++)
        {
            ci = 0;
#pragma ivdep
            for (kp = 0; kp < k; kp++)
            {
                ci += a[kp + k * i] * b[kp + k * j];
            }
            c[i + m * j] = ci;
        }
    }
}

void dgemm_NT(int m, int n, int k,
              double *a, double *b, double *c)
{
    int i, j, kp;
    double bi;
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < m; i++)
        {
            c[i + m * j] = 0;
        }
        for (kp = 0; kp < k; kp++)
        {
            bi = b[j + n * kp];
#pragma ivdep
            for (i = 0; i < m; i++)
            {
                c[i + m * j] += a[i + m * kp] * bi;
            }
        }
    }
}

// transform integrals from cartesian to spheric
double *a_bra_cart2spheric(double *gsph, int nket, double *gcart, int l)
{
    int nf = _len_cart[l];
    int nd = l * 2 + 1;
    dgemm_TN(nd, nket, nf, g_c2s[l].cart2sph, gcart, gsph);
    return gsph;
}

double *a_ket_cart2spheric(double *gsph, double *gcart,
                           int lds, int nbra, int l)
{
    int nf = _len_cart[l];
    int nd = l * 2 + 1;
    dgemm_NN1(nbra, nd, nf, gcart, g_c2s[l].cart2sph, gsph, lds);
    return gsph;
}

// transform s function from cartesian to spheric
double *s_bra_cart2spheric(double *gsph, int nket, double *gcart, int l)
{
    return gcart;
}
double *s_ket_cart2spheric(double *gsph, double *gcart, int lds, int nbra, int l)
{
    return gcart;
}
double *s_ket_cart2spheric_copy(double *gsph, double *gcart,
                                int lds, int nbra, int l)
{
    int i;
    for (i = 0; i < nbra; i++)
    {
        gsph[i] = gcart[i];
    }
    return gsph;
}

// transform p function from cartesian to spheric
double *p_bra_cart2spheric(double *gsph, int nket, double *gcart, int l)
{
#ifdef PYPZPX
    int i;
    for (i = 0; i < nket; i++)
    {
        gsph[i * 3 + 0] = gcart[i * 3 + 1]; // py
        gsph[i * 3 + 1] = gcart[i * 3 + 2]; // pz
        gsph[i * 3 + 2] = gcart[i * 3 + 0]; // px
    }
    return gsph;
#else
    return gcart;
#endif
}
double *p_ket_cart2spheric(double *gsph, double *gcart,
                           int lds, int nbra, int l)
{
#ifdef PYPZPX
    int i;
    for (i = 0; i < nbra; i++)
    {
        gsph[0 * nbra + i] = gcart[1 * nbra + i]; // py
        gsph[1 * nbra + i] = gcart[2 * nbra + i]; // pz
        gsph[2 * nbra + i] = gcart[0 * nbra + i]; // px
    }
    return gsph;
#else
    return gcart;
#endif
}
double *p_ket_cart2spheric_copy(double *gsph, double *gcart,
                                int lds, int nbra, int l)
{
    int i;
#ifdef PYPZPX
    for (i = 0; i < nbra; i++)
    {
        gsph[0 * nbra + i] = gcart[1 * nbra + i]; // py
        gsph[1 * nbra + i] = gcart[2 * nbra + i]; // pz
        gsph[2 * nbra + i] = gcart[0 * nbra + i]; // px
    }
#else
    for (i = 0; i < nbra; i++)
    {
        gsph[0 * lds + i] = gcart[0 * nbra + i];
        gsph[1 * lds + i] = gcart[1 * nbra + i];
        gsph[2 * lds + i] = gcart[2 * nbra + i];
    }
#endif
    return gsph;
}

// transform d function from cartesian to spheric
double *d_bra_cart2spheric(double *gsph, int nket, double *gcart, int l)
{
    double *coeff_c2s = g_c2s[2].cart2sph;
    double *pgsph = gsph;
    int i;
    for (i = 0; i < nket; i++)
    {
        gsph[0] = coeff_c2s[1] * gcart[1];
        gsph[1] = coeff_c2s[10] * gcart[4];
        gsph[2] = coeff_c2s[12] * gcart[0] + coeff_c2s[15] * gcart[3] + coeff_c2s[17] * gcart[5];
        gsph[3] = coeff_c2s[20] * gcart[2];
        gsph[4] = coeff_c2s[24] * gcart[0] + coeff_c2s[27] * gcart[3];
        gsph += 5;
        gcart += 6;
    }
    return pgsph;
}
double *d_ket_cart2spheric(double *gsph, double *gcart,
                           int lds, int nbra, int l)
{
    double *coeff_c2s = g_c2s[2].cart2sph;
    double *pgsph = gsph;
    int i;
    for (i = 0; i < nbra; i++)
    {
        gsph[0 * lds + i] = coeff_c2s[1] * gcart[1 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[1 * lds + i] = coeff_c2s[10] * gcart[4 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[2 * lds + i] = coeff_c2s[12] * gcart[0 * nbra + i] + coeff_c2s[15] * gcart[3 * nbra + i] + coeff_c2s[17] * gcart[5 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[3 * lds + i] = coeff_c2s[20] * gcart[2 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[4 * lds + i] = coeff_c2s[24] * gcart[0 * nbra + i] + coeff_c2s[27] * gcart[3 * nbra + i];
    }
    return pgsph;
}

// transform f function from cartesian to spheric
double *f_bra_cart2spheric(double *gsph, int nket, double *gcart, int l)
{
    double *coeff_c2s = g_c2s[3].cart2sph;
    double *pgsph = gsph;
    int i;
    for (i = 0; i < nket; i++)
    {
        gsph[0] = coeff_c2s[1] * gcart[1] + coeff_c2s[6] * gcart[6];
        gsph[1] = coeff_c2s[14] * gcart[4];
        gsph[2] = coeff_c2s[21] * gcart[1] + coeff_c2s[26] * gcart[6] + coeff_c2s[28] * gcart[8];
        gsph[3] = coeff_c2s[32] * gcart[2] + coeff_c2s[37] * gcart[7] + coeff_c2s[39] * gcart[9];
        gsph[4] = coeff_c2s[40] * gcart[0] + coeff_c2s[43] * gcart[3] + coeff_c2s[45] * gcart[5];
        gsph[5] = coeff_c2s[52] * gcart[2] + coeff_c2s[57] * gcart[7];
        gsph[6] = coeff_c2s[60] * gcart[0] + coeff_c2s[63] * gcart[3];
        gsph += 7;
        gcart += 10;
    }
    return pgsph;
}
double *f_ket_cart2spheric(double *gsph, double *gcart,
                           int lds, int nbra, int l)
{
    double *coeff_c2s = g_c2s[3].cart2sph;
    double *pgsph = gsph;
    int i;
    for (i = 0; i < nbra; i++)
    {
        gsph[0 * lds + i] = coeff_c2s[1] * gcart[1 * nbra + i] + coeff_c2s[6] * gcart[6 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[1 * lds + i] = coeff_c2s[14] * gcart[4 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[2 * lds + i] = coeff_c2s[21] * gcart[1 * nbra + i] + coeff_c2s[26] * gcart[6 * nbra + i] + coeff_c2s[28] * gcart[8 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[3 * lds + i] = coeff_c2s[32] * gcart[2 * nbra + i] + coeff_c2s[37] * gcart[7 * nbra + i] + coeff_c2s[39] * gcart[9 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[4 * lds + i] = coeff_c2s[40] * gcart[0 * nbra + i] + coeff_c2s[43] * gcart[3 * nbra + i] + coeff_c2s[45] * gcart[5 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[5 * lds + i] = coeff_c2s[52] * gcart[2 * nbra + i] + coeff_c2s[57] * gcart[7 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[6 * lds + i] = coeff_c2s[60] * gcart[0 * nbra + i] + coeff_c2s[63] * gcart[3 * nbra + i];
    }
    return pgsph;
}

// transform g function from cartesian to spheric
double *g_bra_cart2spheric(double *gsph, int nket, double *gcart, int l)
{
    double *coeff_c2s = g_c2s[4].cart2sph;
    double *pgsph = gsph;
    int i;
    for (i = 0; i < nket; i++)
    {
        gsph[0] = coeff_c2s[1] * gcart[1] + coeff_c2s[6] * gcart[6];
        gsph[1] = coeff_c2s[19] * gcart[4] + coeff_c2s[26] * gcart[11];
        gsph[2] = coeff_c2s[31] * gcart[1] + coeff_c2s[36] * gcart[6] + coeff_c2s[38] * gcart[8];
        gsph[3] = coeff_c2s[49] * gcart[4] + coeff_c2s[56] * gcart[11] + coeff_c2s[58] * gcart[13];
        gsph[4] = coeff_c2s[60] * gcart[0] + coeff_c2s[63] * gcart[3] + coeff_c2s[65] * gcart[5] + coeff_c2s[70] * gcart[10] + coeff_c2s[72] * gcart[12] + coeff_c2s[74] * gcart[14];
        gsph[5] = coeff_c2s[77] * gcart[2] + coeff_c2s[82] * gcart[7] + coeff_c2s[84] * gcart[9];
        gsph[6] = coeff_c2s[90] * gcart[0] + coeff_c2s[95] * gcart[5] + coeff_c2s[100] * gcart[10] + coeff_c2s[102] * gcart[12];
        gsph[7] = coeff_c2s[107] * gcart[2] + coeff_c2s[112] * gcart[7];
        gsph[8] = coeff_c2s[120] * gcart[0] + coeff_c2s[123] * gcart[3] + coeff_c2s[130] * gcart[10];
        gsph += 9;
        gcart += 15;
    }
    return pgsph;
}
double *g_ket_cart2spheric(double *gsph, double *gcart,
                           int lds, int nbra, int l)
{
    double *coeff_c2s = g_c2s[4].cart2sph;
    double *pgsph = gsph;
    int i;
    for (i = 0; i < nbra; i++)
    {
        gsph[0 * lds + i] = coeff_c2s[1] * gcart[1 * nbra + i] + coeff_c2s[6] * gcart[6 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[1 * lds + i] = coeff_c2s[19] * gcart[4 * nbra + i] + coeff_c2s[26] * gcart[11 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[2 * lds + i] = coeff_c2s[31] * gcart[1 * nbra + i] + coeff_c2s[36] * gcart[6 * nbra + i] + coeff_c2s[38] * gcart[8 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[3 * lds + i] = coeff_c2s[49] * gcart[4 * nbra + i] + coeff_c2s[56] * gcart[11 * nbra + i] + coeff_c2s[58] * gcart[13 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[4 * lds + i] = coeff_c2s[60] * gcart[0 * nbra + i] + coeff_c2s[63] * gcart[3 * nbra + i] + coeff_c2s[65] * gcart[5 * nbra + i] + coeff_c2s[70] * gcart[10 * nbra + i] + coeff_c2s[72] * gcart[12 * nbra + i] + coeff_c2s[74] * gcart[14 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[5 * lds + i] = coeff_c2s[77] * gcart[2 * nbra + i] + coeff_c2s[82] * gcart[7 * nbra + i] + coeff_c2s[84] * gcart[9 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[6 * lds + i] = coeff_c2s[90] * gcart[0 * nbra + i] + coeff_c2s[95] * gcart[5 * nbra + i] + coeff_c2s[100] * gcart[10 * nbra + i] + coeff_c2s[102] * gcart[12 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[7 * lds + i] = coeff_c2s[107] * gcart[2 * nbra + i] + coeff_c2s[112] * gcart[7 * nbra + i];
    }
    for (i = 0; i < nbra; i++)
    {
        gsph[8 * lds + i] = coeff_c2s[120] * gcart[0 * nbra + i] + coeff_c2s[123] * gcart[3 * nbra + i] + coeff_c2s[130] * gcart[10 * nbra + i];
    }
    return pgsph;
}

double *(*c2s_bra_sph[])(double *, int, double *, int) = {
    s_bra_cart2spheric,
    p_bra_cart2spheric,
    d_bra_cart2spheric,
    f_bra_cart2spheric,
    g_bra_cart2spheric,
    a_bra_cart2spheric,
    a_bra_cart2spheric,
    a_bra_cart2spheric,
    a_bra_cart2spheric,
    a_bra_cart2spheric,
    a_bra_cart2spheric,
    a_bra_cart2spheric,
    a_bra_cart2spheric,
    a_bra_cart2spheric,
    a_bra_cart2spheric,
    a_bra_cart2spheric,
};

double *(*c2s_ket_sph[])(double *, double *, int, int, int) = {
    s_ket_cart2spheric,
    p_ket_cart2spheric,
    d_ket_cart2spheric,
    f_ket_cart2spheric,
    g_ket_cart2spheric,
    a_ket_cart2spheric,
    a_ket_cart2spheric,
    a_ket_cart2spheric,
    a_ket_cart2spheric,
    a_ket_cart2spheric,
    a_ket_cart2spheric,
    a_ket_cart2spheric,
    a_ket_cart2spheric,
    a_ket_cart2spheric,
    a_ket_cart2spheric,
    a_ket_cart2spheric,
};

static void dcopy_ij(double *out, double *gctr,
                     const int ni, const int nj, const int mi, const int mj)
{
    int i, j;

    for (j = 0; j < mj; j++)
    {
        for (i = 0; i < mi; i++)
        {
            out[j * ni + i] = gctr[j * mi + i];
        }
    }
}

void c2s_sph_1e(double *opij, double *gctr, int *dims,
                Env *envs, double *cache)
{
    int i_l = envs->i_l;
    int j_l = envs->j_l;
    int i_ctr = envs->x_ctr[0];
    int j_ctr = envs->x_ctr[1];
    int di = i_l * 2 + 1;
    int dj = j_l * 2 + 1;
    int ni = dims[0];
    int nj = dims[1];
    int ofj = ni * dj;
    int nfi = envs->nfi;
    int nf = envs->nf;
    int ic, jc;
    int buflen = nfi * dj;
    double *buf1 = NULL, *buf2 = NULL;
    MALLOC_INSTACK(buf1, buflen, cache);
    MALLOC_INSTACK(buf2, buflen, cache);
    double *pij = NULL;
    double *tmp1 = NULL;

    for (jc = 0; jc < j_ctr; jc++)
    {
        for (ic = 0; ic < i_ctr; ic++)
        {
            tmp1 = (c2s_ket_sph[j_l])(buf1, gctr, nfi, nfi, j_l);
            tmp1 = (c2s_bra_sph[i_l])(buf2, dj, tmp1, i_l);
            pij = opij + ofj * jc + di * ic;
            dcopy_ij(pij, tmp1, ni, nj, di, dj);
            gctr += nf;
        }
    }
}

void c2s_dset0(double *out, int *dims, int *counts)
{
    int ni = dims[0];
    int nj = dims[1];
    int nk = dims[2];
    size_t nij = ni * nj;
    size_t nijk = nij * nk;
    int i, j, k, l;
    if (dims == counts)
    {
        for (i = 0; i < nijk * counts[3]; i++)
        {
            out[i] = 0;
        }
        return;
    }
    int di = counts[0];
    int dj = counts[1];
    int dk = counts[2];
    int dl = counts[3];
    double *pout;
    for (l = 0; l < dl; l++)
    {
        for (k = 0; k < dk; k++)
        {
            pout = out + k * nij;
            for (j = 0; j < dj; j++)
            {
                for (i = 0; i < di; i++)
                {
                    pout[j * ni + i] = 0;
                }
            }
        }
        out += nijk;
    }
}

int _2c2e_drv(double *out, int *dims, Env *envs, Opt *opt, double *cache, void (*f_c2s)(double *, double *, int *, Env *, double *))
{
    int *x_ctr = envs->x_ctr;
    int nc = envs->nf * x_ctr[0] * x_ctr[1];
    int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
    if (out == NULL)
    {
        return int1e_cache_size(envs);
    }
    double *stack = NULL;
    if (cache == NULL)
    {
        size_t cache_size = int1e_cache_size(envs);
        stack = (double *)malloc(sizeof(double) * cache_size);
        cache = stack;
    }
    double *gctr = NULL;
    MALLOC_INSTACK(gctr, nc * n_comp, cache);

    int n;
    int empty = 1;
    if (opt != NULL)
    {
        envs->opt = opt;
        _2c2e_loop(gctr, envs, cache, &empty);
    }
    else
    {
        _2c2e_loop_nopt(gctr, envs, cache, &empty);
    }

    int counts[4];
    if (f_c2s == &c2s_sph_1e)
    {
        counts[0] = (envs->i_l * 2 + 1) * x_ctr[0];
        counts[1] = (envs->k_l * 2 + 1) * x_ctr[1];
    }
    else
    {
        counts[0] = envs->nfi * x_ctr[0];
        counts[1] = envs->nfk * x_ctr[1];
    }
    counts[2] = 1;
    counts[3] = 1;
    if (dims == NULL)
    {
        dims = counts;
    }
    int nout = dims[0] * dims[1];
    if (!empty)
    {
        for (n = 0; n < n_comp; n++)
        {
            (*f_c2s)(out + nout * n, gctr + nc * n, dims, envs, cache);
        }
    }
    else
    {
        for (n = 0; n < n_comp; n++)
        {
            c2s_dset0(out + nout * n, dims, counts);
        }
    }
    if (stack != NULL)
    {
        free(stack);
    }
    return !empty;
}

int int2c2e_sph(double *out, int *dims, int *shls, int *atms, int natms, int *bas, int nbas, double *env, Opt *opt, double *cache)
{
    int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
    Env envs;
    init_int2c2e_Env(&envs, ng, shls, atms, natms, bas, nbas, env);
    envs.f_gout = &gout2e;
    return _2c2e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
};

int GTOmax_cache_size(int (*intor)(double *, int *, int *, int *, int, int *, int, double *, Opt *, double *), int *shls_slice, int ncenter,
                      int *atm, int natm, int *bas, int nbas, double *env)
{
    int i;
    int i0 = shls_slice[0];
    int i1 = shls_slice[1];
    for (i = 1; i < ncenter; i++)
    {
        i0 = std::min(i0, shls_slice[i * 2]);
        i1 = std::max(i1, shls_slice[i * 2 + 1]);
    }
    int (*f)(double *, int *, int *, int *, int, int *, int, double *, Opt *, double *) = (int (*)(double *, int *, int *, int *, int, int *, int, double *, Opt *, double *))intor;
    int cache_size = 0;
    int n;
    int shls[4]{i0, i0, i0, i0};
    for (i = i0; i < i1; i++)
    {
        shls[0] = i;
        shls[1] = i;
        shls[2] = i;
        shls[3] = i;
        n = (*f)(NULL, NULL, shls, atm, natm, bas, nbas, env, NULL, NULL);
        cache_size = std::max(cache_size, n);
    }
    return cache_size;
}

void GTOint2c(int (*intor)(double *, int *, int *, int *, int, int *, int, double *, Opt *, double *),
              double *mat, int comp, int hermi,
              int *shls_slice, int *ao_loc, Opt *opt,
              int *atm, int natm, int *bas, int nbas, double *env)
{
    const int ish0 = shls_slice[0];
    const int ish1 = shls_slice[1];
    const int jsh0 = shls_slice[2];
    const int jsh1 = shls_slice[3];
    const int nish = ish1 - ish0;
    const int njsh = jsh1 - jsh0;
    const int naoi = ao_loc[ish1] - ao_loc[ish0];
    const int naoj = ao_loc[jsh1] - ao_loc[jsh0];
    const int cache_size = GTOmax_cache_size(intor, shls_slice, 2,
                                             atm, natm, bas, nbas, env);
//#pragma omp parallel
    {
        int dims[] = {naoi, naoj};
        int ish, jsh, ij, i0, j0;
        int shls[2];
        double *cache = (double *)malloc(sizeof(double) * cache_size);
//#pragma omp for schedule(dynamic, 4)
        for (ij = 0; ij < nish * njsh; ij++)
        {
            ish = ij / njsh;
            jsh = ij % njsh;
            ish += ish0;
            jsh += jsh0;
            shls[0] = ish;
            shls[1] = jsh;
            i0 = ao_loc[ish] - ao_loc[ish0];
            j0 = ao_loc[jsh] - ao_loc[jsh0];
            (*intor)(mat + j0 * naoi + i0, dims, shls,
                     atm, natm, bas, nbas, env, opt, cache);
        }
        free(cache);
    }
}

void g0_il2d_4d(double *g, Env *envs)
{
    int lk = envs->lk_ceil;
    int lj = envs->lj_ceil;
    if (lj == 0 && lk == 0)
    {
        return;
    }
    int nmax = envs->li_ceil + envs->lj_ceil;
    int mmax = envs->lk_ceil + envs->ll_ceil;
    // int li = envs->li_ceil;
    int ll = envs->ll_ceil;
    int nroots = envs->nrys_roots;
    int i, j, k, l, ptr, n;
    int di = envs->g_stride_i;
    int dk = envs->g_stride_k;
    int dl = envs->g_stride_l;
    int dj = envs->g_stride_j;
    double *rirj = envs->rirj;
    double *rkrl = envs->rkrl;
    double *gx = g;
    double *gy = gx + envs->g_size;
    double *gz = gy + 2 * envs->g_size;
    double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;
    double rx, ry, rz;

    // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
    rx = rkrl[0];
    ry = rkrl[1];
    rz = rkrl[2];
    p1x = gx - dk;
    p1y = gy - dk;
    p1z = gz - dk;
    p2x = gx - dk + dl;
    p2y = gy - dk + dl;
    p2z = gz - dk + dl;
    for (k = 1; k <= lk; k++)
    {
        for (l = 0; l <= mmax - k; l++)
        {
            for (i = 0; i <= nmax; i++)
            {
                ptr = l * dl + k * dk + i * di;
                for (n = ptr; n < ptr + nroots; n++)
                {
                    gx[n] = rx * p1x[n] + p2x[n];
                    gy[n] = ry * p1y[n] + p2y[n];
                    gz[n] = rz * p1z[n] + p2z[n];
                }
            }
        }
    }

    // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
    rx = rirj[0];
    ry = rirj[1];
    rz = rirj[2];
    p1x = gx - dj;
    p1y = gy - dj;
    p1z = gz - dj;
    p2x = gx - dj + di;
    p2y = gy - dj + di;
    p2z = gz - dj + di;
    for (j = 1; j <= lj; j++)
    {
        for (l = 0; l <= ll; l++)
        {
            for (k = 0; k <= lk; k++)
            {
                ptr = j * dj + l * dl + k * dk;
                for (n = ptr; n < ptr + dk - di * j; n++)
                {
                    gx[n] = rx * p1x[n] + p2x[n];
                    gy[n] = ry * p1y[n] + p2y[n];
                    gz[n] = rz * p1z[n] + p2z[n];
                }
            }
        }
    }
}

void g0_lj2d_4d(double *g, Env *envs)
{
    int li = envs->li_ceil;
    int lk = envs->lk_ceil;
    if (li == 0 && lk == 0)
    {
        return;
    }
    int nmax = envs->li_ceil + envs->lj_ceil;
    int mmax = envs->lk_ceil + envs->ll_ceil;
    // int ll = envs->ll_ceil;
    int lj = envs->lj_ceil;
    int nroots = envs->nrys_roots;
    int i, j, k, l, ptr, n;
    int di = envs->g_stride_i;
    int dk = envs->g_stride_k;
    int dl = envs->g_stride_l;
    int dj = envs->g_stride_j;
    double *rirj = envs->rirj;
    double *rkrl = envs->rkrl;
    double *gx = g;
    double *gy = gx + envs->g_size;
    double *gz = gy + 2 * envs->g_size;
    double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;
    double rx, ry, rz;

    // g(i,...,j) = rirj * g(i-1,...,j) +  g(i-1,...,j+1)
    rx = rirj[0];
    ry = rirj[1];
    rz = rirj[2];
    p1x = gx - di;
    p1y = gy - di;
    p1z = gz - di;
    p2x = gx - di + dj;
    p2y = gy - di + dj;
    p2z = gz - di + dj;
    for (i = 1; i <= li; i++)
    {
        for (j = 0; j <= nmax - i; j++)
        {
            for (l = 0; l <= mmax; l++)
            {
                ptr = j * dj + l * dl + i * di;
                for (n = ptr; n < ptr + nroots; n++)
                {
                    gx[n] = rx * p1x[n] + p2x[n];
                    gy[n] = ry * p1y[n] + p2y[n];
                    gz[n] = rz * p1z[n] + p2z[n];
                }
            }
        }
    }

    // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
    rx = rkrl[0];
    ry = rkrl[1];
    rz = rkrl[2];
    p1x = gx - dk;
    p1y = gy - dk;
    p1z = gz - dk;
    p2x = gx - dk + dl;
    p2y = gy - dk + dl;
    p2z = gz - dk + dl;
    for (j = 0; j <= lj; j++)
    {
        for (k = 1; k <= lk; k++)
        {
            for (l = 0; l <= mmax - k; l++)
            {
                ptr = j * dj + l * dl + k * dk;
                for (n = ptr; n < ptr + dk; n++)
                {
                    gx[n] = rx * p1x[n] + p2x[n];
                    gy[n] = ry * p1y[n] + p2y[n];
                    gz[n] = rz * p1z[n] + p2z[n];
                }
            }
        }
    }
}

void g0_2e_il2d4d(double *g, Rys2eT *bc, Env *envs)
{
    g0_2e_2d(g, bc, envs);
    g0_il2d_4d(g, envs);
}
void g0_2e_lj2d4d(double *g, Rys2eT *bc, Env *envs)
{
    g0_2e_2d(g, bc, envs);
    g0_lj2d_4d(g, envs);
}

void init_int3c2e_Env(Env *envs, int *ng, int *shls,
                      int *atm, int natm, int *bas, int nbas, double *env)
{
    envs->natm = natm;
    envs->nbas = nbas;
    envs->atm = atm;
    envs->bas = bas;
    envs->env = env;
    envs->shls = shls;

    const int i_sh = shls[0];
    const int j_sh = shls[1];
    const int k_sh = shls[2];
    envs->i_l = bas(1, i_sh);
    envs->j_l = bas(1, j_sh);
    envs->k_l = bas(1, k_sh);
    envs->l_l = 0;
    envs->x_ctr[0] = bas(3, i_sh);
    envs->x_ctr[1] = bas(3, j_sh);
    envs->x_ctr[2] = bas(3, k_sh);
    envs->x_ctr[3] = 1;
    envs->nfi = (envs->i_l + 1) * (envs->i_l + 2) / 2;
    envs->nfj = (envs->j_l + 1) * (envs->j_l + 2) / 2;
    envs->nfk = (envs->k_l + 1) * (envs->k_l + 2) / 2;
    envs->nfl = 1;
    envs->nf = envs->nfi * envs->nfk * envs->nfj;

    envs->ri = env + atm(1, bas(0, i_sh));
    envs->rj = env + atm(1, bas(0, j_sh));
    envs->rk = env + atm(1, bas(0, k_sh));

    envs->common_factor = (constants::PI3) * 2 / constants::sqr_pi * common_fac_sp(envs->i_l) * common_fac_sp(envs->j_l) * common_fac_sp(envs->k_l);
    if (env[0] == 0)
    {
        envs->expcutoff = 60.;
    }
    else
    {
        envs->expcutoff = std::max(40., env[0]);
    }

    envs->gbits = ng[4];
    envs->ncomp_e1 = ng[5];
    envs->ncomp_e2 = ng[6];
    envs->ncomp_tensor = ng[7];

    envs->li_ceil = envs->i_l + ng[0];
    envs->lj_ceil = envs->j_l + ng[1];
    envs->lk_ceil = 0; // to reuse CINTg0_2e_lj2d4d
    envs->ll_ceil = envs->k_l + ng[2];
    int rys_order = (envs->li_ceil + envs->lj_ceil + envs->ll_ceil) / 2 + 1;
    int nrys_roots = rys_order;
    double omega = env[8];
    if (omega < 0 && rys_order <= 3)
    {
        nrys_roots *= 2;
    }
    envs->rys_order = rys_order;
    envs->nrys_roots = nrys_roots;

    int dli, dlj, dlk;
    int ibase = envs->li_ceil > envs->lj_ceil;
    if (ibase)
    {
        dli = envs->li_ceil + envs->lj_ceil + 1;
        dlj = envs->lj_ceil + 1;
    }
    else
    {
        dli = envs->li_ceil + 1;
        dlj = envs->li_ceil + envs->lj_ceil + 1;
    }
    dlk = envs->ll_ceil + 1;

    envs->g_stride_i = nrys_roots;
    envs->g_stride_k = nrys_roots * dli;
    envs->g_stride_l = nrys_roots * dli;
    envs->g_stride_j = nrys_roots * dli * dlk;
    envs->g_size = nrys_roots * dli * dlk * dlj;

    envs->al[0] = 0;
    envs->rkl[0] = envs->rk[0];
    envs->rkl[1] = envs->rk[1];
    envs->rkl[2] = envs->rk[2];
    envs->g2d_klmax = envs->g_stride_k;
    envs->rkrl[0] = envs->rk[0];
    envs->rkrl[1] = envs->rk[1];
    envs->rkrl[2] = envs->rk[2];
    // in g0_2d rklrx = rkl - rx = 0 => rkl = rx
    envs->rx_in_rklrx = envs->rk;

    if (ibase)
    {
        envs->g2d_ijmax = envs->g_stride_i;
        envs->rx_in_rijrx = envs->ri;
        envs->rirj[0] = envs->ri[0] - envs->rj[0];
        envs->rirj[1] = envs->ri[1] - envs->rj[1];
        envs->rirj[2] = envs->ri[2] - envs->rj[2];
    }
    else
    {
        envs->g2d_ijmax = envs->g_stride_j;
        envs->rx_in_rijrx = envs->rj;
        envs->rirj[0] = envs->rj[0] - envs->ri[0];
        envs->rirj[1] = envs->rj[1] - envs->ri[1];
        envs->rirj[2] = envs->rj[2] - envs->ri[2];
    }

    if (rys_order <= 2)
    {
        envs->f_g0_2d4d = &g0_2e_2d4d_unrolled;
        if (rys_order != nrys_roots)
        {
            envs->f_g0_2d4d = &srg0_2e_2d4d_unrolled;
        }
    }
    else if (ibase)
    {
        envs->f_g0_2d4d = &g0_2e_il2d4d;
    }
    else
    {
        envs->f_g0_2d4d = &g0_2e_lj2d4d;
    }
    envs->f_g0_2e = &g0_2e;
}

void dcopy_iklj(double *fijkl, const double *gctr,
                const int ni, const int nj, const int nk, const int nl,
                const int mi, const int mj, const int mk, const int ml)
{
    const size_t nij = ni * nj;
    const size_t nijk = nij * nk;
    const size_t mik = mi * mk;
    const size_t mikl = mik * ml;
    int i, j, k, l;
    double *pijkl;
    const double *pgctr;

    switch (mi)
    {
    case 1:
        for (l = 0; l < ml; l++)
        {
            for (k = 0; k < mk; k++)
            {
                pijkl = fijkl + k * nij;
                pgctr = gctr + k * mi;
#pragma ivdep
                for (j = 0; j < mj; j++)
                {
                    pijkl[ni * j] = pgctr[mikl * j];
                }
            }
            fijkl += nijk;
            gctr += mik;
        }
        break;
    case 3:
        for (l = 0; l < ml; l++)
        {
            for (k = 0; k < mk; k++)
            {
                pijkl = fijkl + k * nij;
                pgctr = gctr + k * mi;
#pragma ivdep
                for (j = 0; j < mj; j++)
                {
                    pijkl[ni * j + 0] = pgctr[mikl * j + 0];
                    pijkl[ni * j + 1] = pgctr[mikl * j + 1];
                    pijkl[ni * j + 2] = pgctr[mikl * j + 2];
                }
            }
            fijkl += nijk;
            gctr += mik;
        }
        break;
    case 5:
        for (l = 0; l < ml; l++)
        {
            for (k = 0; k < mk; k++)
            {
                pijkl = fijkl + k * nij;
                pgctr = gctr + k * mi;
#pragma ivdep
                for (j = 0; j < mj; j++)
                {
                    pijkl[ni * j + 0] = pgctr[mikl * j + 0];
                    pijkl[ni * j + 1] = pgctr[mikl * j + 1];
                    pijkl[ni * j + 2] = pgctr[mikl * j + 2];
                    pijkl[ni * j + 3] = pgctr[mikl * j + 3];
                    pijkl[ni * j + 4] = pgctr[mikl * j + 4];
                }
            }
            fijkl += nijk;
            gctr += mik;
        }
        break;
    case 6:
        for (l = 0; l < ml; l++)
        {
            for (k = 0; k < mk; k++)
            {
                pijkl = fijkl + k * nij;
                pgctr = gctr + k * mi;
#pragma ivdep
                for (j = 0; j < mj; j++)
                {
                    pijkl[ni * j + 0] = pgctr[mikl * j + 0];
                    pijkl[ni * j + 1] = pgctr[mikl * j + 1];
                    pijkl[ni * j + 2] = pgctr[mikl * j + 2];
                    pijkl[ni * j + 3] = pgctr[mikl * j + 3];
                    pijkl[ni * j + 4] = pgctr[mikl * j + 4];
                    pijkl[ni * j + 5] = pgctr[mikl * j + 5];
                }
            }
            fijkl += nijk;
            gctr += mik;
        }
        break;
    case 7:
        for (l = 0; l < ml; l++)
        {
            for (k = 0; k < mk; k++)
            {
                pijkl = fijkl + k * nij;
                pgctr = gctr + k * mi;
#pragma ivdep
                for (j = 0; j < mj; j++)
                {
                    pijkl[ni * j + 0] = pgctr[mikl * j + 0];
                    pijkl[ni * j + 1] = pgctr[mikl * j + 1];
                    pijkl[ni * j + 2] = pgctr[mikl * j + 2];
                    pijkl[ni * j + 3] = pgctr[mikl * j + 3];
                    pijkl[ni * j + 4] = pgctr[mikl * j + 4];
                    pijkl[ni * j + 5] = pgctr[mikl * j + 5];
                    pijkl[ni * j + 6] = pgctr[mikl * j + 6];
                }
            }
            fijkl += nijk;
            gctr += mik;
        }
        break;
    default:
        for (l = 0; l < ml; l++)
        {
            for (k = 0; k < mk; k++)
            {
                pijkl = fijkl + k * nij;
                pgctr = gctr + k * mi;
                for (j = 0; j < mj; j++)
                {
#pragma ivdep
                    for (i = 0; i < mi; i++)
                    {
                        pijkl[ni * j + i] = pgctr[mikl * j + i];
                    }
                }
            }
            fijkl += nijk;
            gctr += mik;
        }
    }
}

double *sph2e_inner(double *gsph, double *gcart,
                    int l, int nbra, int ncall, int sizsph, int sizcart)
{
    int n;
    switch (l)
    {
#ifdef PYPZPX
    case 0:
        return gcart;
    case 1:
        for (n = 0; n < ncall; n++)
        {
            p_ket_cart2spheric(gsph + n * sizsph, gcart + n * sizcart, nbra, nbra, l);
        }
        break;
#else
    case 0:
    case 1:
        return gcart;
#endif
    case 2:
        for (n = 0; n < ncall; n++)
        {
            d_ket_cart2spheric(gsph + n * sizsph, gcart + n * sizcart, nbra, nbra, l);
        }
        break;
    case 3:
        for (n = 0; n < ncall; n++)
        {
            f_ket_cart2spheric(gsph + n * sizsph, gcart + n * sizcart, nbra, nbra, l);
        }
        break;
    case 4:
        for (n = 0; n < ncall; n++)
        {
            g_ket_cart2spheric(gsph + n * sizsph, gcart + n * sizcart, nbra, nbra, l);
        }
        break;
    default:
        for (n = 0; n < ncall; n++)
        {
            a_ket_cart2spheric(gsph + n * sizsph, gcart + n * sizcart, nbra, nbra, l);
        }
    }
    return gsph;
}

void c2s_sph_3c2e1(double *bufijk, double *gctr, int *dims,
                   Env *envs, double *cache)
{
    int i_l = envs->i_l;
    int j_l = envs->j_l;
    int k_l = envs->k_l;
    int i_ctr = envs->x_ctr[0];
    int j_ctr = envs->x_ctr[1];
    int k_ctr = envs->x_ctr[2];
    int di = i_l * 2 + 1;
    int dj = j_l * 2 + 1;
    int dk = k_l * 2 + 1;
    int ni = dims[0];
    int nj = dims[1];
    int nk = dims[2];
    int nfi = envs->nfi;
    int nfk = envs->nfk;
    int nf = envs->nf;
    int nfik = nfi * nfk;
    int ofj = ni * dj;
    int ofk = ni * nj * dk;
    int ic, jc, kc;
    int buflen = nfi * nfk * dj;
    double *buf1 = NULL;
    MALLOC_INSTACK(buf1, buflen * 3, cache);
    double *buf2 = buf1 + buflen;
    double *buf3 = buf2 + buflen;
    double *pijk;
    double *tmp1;

    for (kc = 0; kc < k_ctr; kc++)
    {
        for (jc = 0; jc < j_ctr; jc++)
        {
            for (ic = 0; ic < i_ctr; ic++)
            {
                tmp1 = (c2s_ket_sph[j_l])(buf1, gctr, nfik, nfik, j_l);
                tmp1 = sph2e_inner(buf2, tmp1, k_l, nfi, dj, nfi * dk, nfik);
                tmp1 = (c2s_bra_sph[i_l])(buf3, dk * dj, tmp1, i_l);
                pijk = bufijk + ofk * kc + ofj * jc + di * ic;
                dcopy_iklj(pijk, tmp1, ni, nj, nk, 1, di, dj, dk, 1);
                gctr += nf;
            }
        }
    }
}

int set_pairdata(PairData *pairdata, double *ai, double *aj, double *ri, double *rj,
                 double *log_maxci, double *log_maxcj,
                 int li_ceil, int lj_ceil, int iprim, int jprim,
                 double rr_ij, double expcutoff, double *env)
{
    int ip, jp, n;
    double aij, eij, cceij, wj;
    // Normally
    //    (aj*d/sqrt(aij)+1)^li * (ai*d/sqrt(aij)+1)^lj
    //    * pi^1.5/aij^{(li+lj+3)/2} * exp(-ai*aj/aij*rr_ij)
    // is a good approximation for overlap integrals.
    //    <~ (aj*d/aij+1/sqrt(aij))^li * (ai*d/aij+1/sqrt(aij))^lj * (pi/aij)^1.5
    //    <~ (d+1/sqrt(aij))^(li+lj) * (pi/aij)^1.5
    aij = ai[iprim - 1] + aj[jprim - 1];
    double log_rr_ij = 1.7 - 1.5 * log(aij);
    int lij = li_ceil + lj_ceil;
    if (lij > 0)
    {
        double dist_ij = sqrt(rr_ij);
        double omega = env[8];
        if (omega < 0)
        {
            double r_guess = 8.;
            double omega2 = omega * omega;
            double theta = omega2 / (omega2 + aij);
            log_rr_ij += lij * log(dist_ij + theta * r_guess + 1.);
        }
        else
        {
            log_rr_ij += lij * log(dist_ij + 1.);
        }
    }
    PairData *pdata;

    int empty = 1;
    for (n = 0, jp = 0; jp < jprim; jp++)
    {
        for (ip = 0; ip < iprim; ip++, n++)
        {
            aij = 1 / (ai[ip] + aj[jp]);
            eij = rr_ij * ai[ip] * aj[jp] * aij;
            cceij = eij - log_rr_ij - log_maxci[ip] - log_maxcj[jp];
            pdata = pairdata + n;
            pdata->cceij = cceij;
            if (cceij < expcutoff)
            {
                empty = 0;
                wj = aj[jp] * aij;
                pdata->rij[0] = ri[0] + wj * (rj[0] - ri[0]);
                pdata->rij[1] = ri[1] + wj * (rj[1] - ri[1]);
                pdata->rij[2] = ri[2] + wj * (rj[2] - ri[2]);
                pdata->eij = exp(-eij);
            }
            else
            {
                pdata->rij[0] = 1e18;
                pdata->rij[1] = 1e18;
                pdata->rij[2] = 1e18;
                pdata->eij = 0;
            }
        }
    }
    return empty;
}

#define PAIRDATA_NON0IDX_SIZE(ps) \
    int *bas = envs->bas;         \
    int *shls = envs->shls;       \
    int i_prim = bas(2, shls[0]); \
    int j_prim = bas(2, shls[1]); \
    int k_prim = bas(2, shls[2]); \
    int ps = (i_prim * j_prim * 5 + i_prim * x_ctr[0] + j_prim * x_ctr[1] + k_prim * x_ctr[2] + (i_prim + j_prim) * 2 + k_prim + envs->nf * 3 + 16);

int _3c2e_loop(double *gctr, Env *envs, double *cache, int *empty)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    double *env = envs->env;
    int i_sh = shls[0];
    int j_sh = shls[1];
    Opt *opt = envs->opt;
    if (opt->pairdata != NULL && opt->pairdata[i_sh * opt->nbas + j_sh] == ((void *)0xffffffffffffffffuL))
    {
        return 0;
    }
    int k_sh = shls[2];
    int i_ctr = envs->x_ctr[0];
    int j_ctr = envs->x_ctr[1];
    int k_ctr = envs->x_ctr[2];
    int i_prim = bas(2, i_sh);
    int j_prim = bas(2, j_sh);
    int k_prim = bas(2, k_sh);
    double *ai = env + bas(5, i_sh);
    double *aj = env + bas(5, j_sh);
    double *ak = env + bas(5, k_sh);
    double *ci = env + bas(6, i_sh);
    double *cj = env + bas(6, j_sh);
    double *ck = env + bas(6, k_sh);
    double expcutoff = envs->expcutoff;
    double rr_ij = (envs->rirj)[0] * (envs->rirj)[0] + (envs->rirj)[1] * (envs->rirj)[1] + (envs->rirj)[2] * (envs->rirj)[2];
    PairData *pdata_base = NULL, *pdata_ij = NULL;
    if (opt->pairdata != NULL)
    {
        pdata_base = opt->pairdata[i_sh * opt->nbas + j_sh];
    }
    else
    {
        double *log_maxci = opt->log_max_coeff[i_sh];
        double *log_maxcj = opt->log_max_coeff[j_sh];
        MALLOC_INSTACK(pdata_base, i_prim * j_prim, cache);
        if (set_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                         log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                         i_prim, j_prim, rr_ij, expcutoff, env))
        {
            return 0;
        }
    }
    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
    size_t nf = envs->nf;
    double fac1i, fac1j, fac1k;
    int ip, jp, kp;
    int _empty[4] = {1, 1, 1, 1};
    int *iempty = _empty + 0;
    int *jempty = _empty + 1;
    int *kempty = _empty + 2;
    int *gempty = _empty + 3;
    int *non0ctri = opt->non0ctr[i_sh];
    int *non0ctrj = opt->non0ctr[j_sh];
    int *non0idxi = opt->sortedidx[i_sh];
    int *non0idxj = opt->sortedidx[j_sh];
    int *non0ctrk = NULL, *non0idxk = NULL;
    MALLOC_INSTACK(non0ctrk, k_prim + k_prim * k_ctr, cache);
    non0idxk = non0ctrk + k_prim;
    Opt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
    double expij, cutoff;
    double *rij;
    double *rkl = envs->rkl;
    int *idx = opt->index_xyz_array[envs->i_l * 16 * 16 + envs->j_l * 16 + envs->k_l];
    if (idx == NULL)
    {
        MALLOC_INSTACK(idx, nf * 3, cache);
        g2e_index_xyz(idx, envs);
    }
    double omega = env[8];
    if (omega < 0 && envs->rys_order > 1)
    {
        double r_guess = 8.;
        double omega2 = omega * omega;
        int lij = envs->li_ceil + envs->lj_ceil;
        if (lij > 0)
        {
            double dist_ij = sqrt(rr_ij);
            double aij = ai[i_prim - 1] + aj[j_prim - 1];
            double theta = omega2 / (omega2 + aij);
            expcutoff += lij * log(
                                   (dist_ij + theta * r_guess + 1.) / (dist_ij + 1.));
        }
        if (envs->lk_ceil > 0)
        {
            double theta = omega2 / (omega2 + ak[k_prim - 1]);
            expcutoff += envs->lk_ceil * log(theta * r_guess + 1.);
        }
    }
    int nc = i_ctr * j_ctr * k_ctr;
    // (irys,i,j,k,coord,0:1); +1 for nabla-r12
    int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
    int lenk = nf * nc * n_comp;            // gctrk
    int lenj = nf * i_ctr * j_ctr * n_comp; // gctrj
    int leni = nf * i_ctr * n_comp;         // gctri
    int len0 = nf * n_comp;                 // gout
    int len = leng + lenk + lenj + leni + len0;
    double *g = NULL;
    MALLOC_INSTACK(g, len, cache);
    double *g1 = g + leng;
    double *gout, *gctri, *gctrj, *gctrk;

    ALIAS_ADDR_IF_EQUAL(k, m);
    ALIAS_ADDR_IF_EQUAL(j, k);
    ALIAS_ADDR_IF_EQUAL(i, j);
    ALIAS_ADDR_IF_EQUAL(g, i);

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        if (k_ctr == 1)
        {
            fac1k = envs->common_factor * ck[kp];
        }
        else
        {
            fac1k = envs->common_factor;
            *jempty = 1;
        }
        pdata_ij = pdata_base;
        for (jp = 0; jp < j_prim; jp++)
        {
            envs->aj[0] = aj[jp];
            if (j_ctr == 1)
            {
                fac1j = fac1k * cj[jp];
            }
            else
            {
                fac1j = fac1k;
                *iempty = 1;
            }
            for (ip = 0; ip < i_prim; ip++, pdata_ij++)
            {
                if (pdata_ij->cceij > expcutoff)
                {
                    continue;
                }
                envs->ai[0] = ai[ip];
                expij = pdata_ij->eij;
                rij = pdata_ij->rij;
                cutoff = expcutoff - pdata_ij->cceij;
                if (i_ctr == 1)
                {
                    fac1i = fac1j * ci[ip] * expij;
                }
                else
                {
                    fac1i = fac1j * expij;
                }
                envs->fac[0] = fac1i;
                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs))
                {
                    (*envs->f_gout)(gout, g, idx, envs, *gempty);
                    PRIM2CTR(i, gout, len0);
                }
            } // end loop i_prim
            if (!*iempty)
            {
                PRIM2CTR(j, gctri, leni);
            }
        } // end loop j_prim
        if (!*jempty)
        {
            PRIM2CTR(k, gctrj, lenj);
        }
    } // end loop k_prim

    if (n_comp > 1 && !*kempty)
    {
        TRANSPOSE(gctrk);
    }
    return !*empty;
}

int _3c2e_n11_loop(double *gctr, Env *envs, double *cache, int *empty)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    double *env = envs->env;
    int i_sh = shls[0];
    int j_sh = shls[1];
    Opt *opt = envs->opt;
    if (opt->pairdata != NULL && opt->pairdata[i_sh * opt->nbas + j_sh] == ((void *)0xffffffffffffffffuL))
    {
        return 0;
    }
    int k_sh = shls[2];
    int i_ctr = envs->x_ctr[0];
    int j_ctr = envs->x_ctr[1];
    int k_ctr = envs->x_ctr[2];
    int i_prim = bas(2, i_sh);
    int j_prim = bas(2, j_sh);
    int k_prim = bas(2, k_sh);
    double *ai = env + bas(5, i_sh);
    double *aj = env + bas(5, j_sh);
    double *ak = env + bas(5, k_sh);
    double *ci = env + bas(6, i_sh);
    double *cj = env + bas(6, j_sh);
    double *ck = env + bas(6, k_sh);
    double expcutoff = envs->expcutoff;
    double rr_ij = (envs->rirj)[0] * (envs->rirj)[0] + (envs->rirj)[1] * (envs->rirj)[1] + (envs->rirj)[2] * (envs->rirj)[2];
    PairData *pdata_base = NULL, *pdata_ij = NULL;
    if (opt->pairdata != NULL)
    {
        pdata_base = opt->pairdata[i_sh * opt->nbas + j_sh];
    }
    else
    {
        double *log_maxci = opt->log_max_coeff[i_sh];
        double *log_maxcj = opt->log_max_coeff[j_sh];
        MALLOC_INSTACK(pdata_base, i_prim * j_prim, cache);
        if (set_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                         log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                         i_prim, j_prim, rr_ij, expcutoff, env))
        {
            return 0;
        }
    }
    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
    size_t nf = envs->nf;
    double fac1i, fac1j, fac1k;
    int ip, jp, kp;
    int _empty[4] = {1, 1, 1, 1};
    int *iempty = _empty + 0;
    int *jempty = _empty + 1;
    int *kempty = _empty + 2;
    int *gempty = _empty + 3;
    int *non0ctri = opt->non0ctr[i_sh];
    int *non0idxi = opt->sortedidx[i_sh];
    int *non0ctrk = NULL, *non0idxk = NULL;
    MALLOC_INSTACK(non0ctrk, k_prim + k_prim * k_ctr, cache);
    non0idxk = non0ctrk + k_prim;
    Opt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
    double expij, cutoff;
    double *rij;
    double *rkl = envs->rkl;
    int *idx = opt->index_xyz_array[envs->i_l * 16 * 16 + envs->j_l * 16 + envs->k_l];
    if (idx == NULL)
    {
        MALLOC_INSTACK(idx, nf * 3, cache);
        g2e_index_xyz(idx, envs);
    }
    double omega = env[8];
    if (omega < 0 && envs->rys_order > 1)
    {
        double r_guess = 8.;
        double omega2 = omega * omega;
        int lij = envs->li_ceil + envs->lj_ceil;
        if (lij > 0)
        {
            double dist_ij = sqrt(rr_ij);
            double aij = ai[i_prim - 1] + aj[j_prim - 1];
            double theta = omega2 / (omega2 + aij);
            expcutoff += lij * log(
                                   (dist_ij + theta * r_guess + 1.) / (dist_ij + 1.));
        }
        if (envs->lk_ceil > 0)
        {
            double theta = omega2 / (omega2 + ak[k_prim - 1]);
            expcutoff += envs->lk_ceil * log(theta * r_guess + 1.);
        }
    }
    int nc = i_ctr;
    int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
    int leni = nf * i_ctr * n_comp; // gctri
    int len0 = nf * n_comp;         // gout
    int len = leng + leni + len0;
    double *g = NULL;
    MALLOC_INSTACK(g, len, cache);
    double *g1 = g + leng;
    double *gout, *gctri;
    ALIAS_ADDR_IF_EQUAL(i, m);
    gout = g1;

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        fac1k = envs->common_factor * ck[kp];
        pdata_ij = pdata_base;
        for (jp = 0; jp < j_prim; jp++)
        {
            envs->aj[0] = aj[jp];
            fac1j = fac1k * cj[jp];
            for (ip = 0; ip < i_prim; ip++, pdata_ij++)
            {
                if (pdata_ij->cceij > expcutoff)
                {
                    continue;
                }
                envs->ai[0] = ai[ip];
                expij = pdata_ij->eij;
                rij = pdata_ij->rij;
                cutoff = expcutoff - pdata_ij->cceij;
                fac1i = fac1j * expij;
                envs->fac[0] = fac1i;
                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs))
                {
                    (*envs->f_gout)(gout, g, idx, envs, 1);
                    PRIM2CTR(i, gout, len0);
                }
            } // end loop i_prim
        } // end loop j_prim
    } // end loop k_prim

    if (n_comp > 1 && !*iempty)
    {
        TRANSPOSE(gctri);
    }
    return !*empty;
}

int _3c2e_1n1_loop(double *gctr, Env *envs, double *cache, int *empty)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    double *env = envs->env;
    int i_sh = shls[0];
    int j_sh = shls[1];
    Opt *opt = envs->opt;
    if (opt->pairdata != NULL && opt->pairdata[i_sh * opt->nbas + j_sh] == ((void *)0xffffffffffffffffuL))
    {
        return 0;
    }
    int k_sh = shls[2];
    int i_ctr = envs->x_ctr[0];
    int j_ctr = envs->x_ctr[1];
    int k_ctr = envs->x_ctr[2];
    int i_prim = bas(2, i_sh);
    int j_prim = bas(2, j_sh);
    int k_prim = bas(2, k_sh);
    double *ai = env + bas(5, i_sh);
    double *aj = env + bas(5, j_sh);
    double *ak = env + bas(5, k_sh);
    double *ci = env + bas(6, i_sh);
    double *cj = env + bas(6, j_sh);
    double *ck = env + bas(6, k_sh);
    double expcutoff = envs->expcutoff;
    double rr_ij = (envs->rirj)[0] * (envs->rirj)[0] + (envs->rirj)[1] * (envs->rirj)[1] + (envs->rirj)[2] * (envs->rirj)[2];
    PairData *pdata_base = NULL, *pdata_ij = NULL;
    if (opt->pairdata != NULL)
    {
        pdata_base = opt->pairdata[i_sh * opt->nbas + j_sh];
    }
    else
    {
        double *log_maxci = opt->log_max_coeff[i_sh];
        double *log_maxcj = opt->log_max_coeff[j_sh];
        MALLOC_INSTACK(pdata_base, i_prim * j_prim, cache);
        if (set_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                         log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                         i_prim, j_prim, rr_ij, expcutoff, env))
        {
            return 0;
        }
    }
    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
    size_t nf = envs->nf;
    double fac1i, fac1j, fac1k;
    int ip, jp, kp;
    int _empty[4] = {1, 1, 1, 1};
    int *iempty = _empty + 0;
    int *jempty = _empty + 1;
    int *kempty = _empty + 2;
    int *gempty = _empty + 3;
    int *non0ctrj = opt->non0ctr[j_sh];
    int *non0idxj = opt->sortedidx[j_sh];
    int *non0ctrk = NULL, *non0idxk = NULL;
    MALLOC_INSTACK(non0ctrk, k_prim + k_prim * k_ctr, cache);
    non0idxk = non0ctrk + k_prim;
    Opt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
    double expij, cutoff;
    double *rij;
    double *rkl = envs->rkl;
    int *idx = opt->index_xyz_array[envs->i_l * 16 * 16 + envs->j_l * 16 + envs->k_l];
    if (idx == NULL)
    {
        MALLOC_INSTACK(idx, nf * 3, cache);
        g2e_index_xyz(idx, envs);
    }
    double omega = env[8];
    if (omega < 0 && envs->rys_order > 1)
    {
        double r_guess = 8.;
        double omega2 = omega * omega;
        int lij = envs->li_ceil + envs->lj_ceil;
        if (lij > 0)
        {
            double dist_ij = sqrt(rr_ij);
            double aij = ai[i_prim - 1] + aj[j_prim - 1];
            double theta = omega2 / (omega2 + aij);
            expcutoff += lij * log(
                                   (dist_ij + theta * r_guess + 1.) / (dist_ij + 1.));
        }
        if (envs->lk_ceil > 0)
        {
            double theta = omega2 / (omega2 + ak[k_prim - 1]);
            expcutoff += envs->lk_ceil * log(theta * r_guess + 1.);
        }
    }
    int nc = j_ctr;
    int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
    int lenj = nf * j_ctr * n_comp; // gctrj
    int len0 = nf * n_comp;         // gout
    int len = leng + lenj + len0;
    double *g = NULL;
    MALLOC_INSTACK(g, len, cache);
    double *g1 = g + leng;
    double *gout, *gctrj;
    ALIAS_ADDR_IF_EQUAL(j, m);
    gout = g1;

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        fac1k = envs->common_factor * ck[kp];
        pdata_ij = pdata_base;
        for (jp = 0; jp < j_prim; jp++)
        {
            envs->aj[0] = aj[jp];
            fac1j = fac1k;
            *iempty = 1;
            for (ip = 0; ip < i_prim; ip++, pdata_ij++)
            {
                if (pdata_ij->cceij > expcutoff)
                {
                    continue;
                }
                envs->ai[0] = ai[ip];
                expij = pdata_ij->eij;
                rij = pdata_ij->rij;
                cutoff = expcutoff - pdata_ij->cceij;
                fac1i = fac1j * ci[ip] * expij;
                envs->fac[0] = fac1i;
                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs))
                {
                    (*envs->f_gout)(gout, g, idx, envs, *iempty);
                    *iempty = 0;
                }
            } // end loop i_prim
            if (!*iempty)
            {
                PRIM2CTR(j, gout, len0);
            }
        } // end loop j_prim
    } // end loop k_prim

    if (n_comp > 1 && !*jempty)
    {
        TRANSPOSE(gctrj);
    }
    return !*empty;
}

int _3c2e_111_loop(double *gctr, Env *envs, double *cache, int *empty)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    double *env = envs->env;
    int i_sh = shls[0];
    int j_sh = shls[1];
    Opt *opt = envs->opt;
    if (opt->pairdata != NULL && opt->pairdata[i_sh * opt->nbas + j_sh] == ((void *)0xffffffffffffffffuL))
    {
        return 0;
    }
    int k_sh = shls[2];
    int i_ctr = envs->x_ctr[0];
    int j_ctr = envs->x_ctr[1];
    int k_ctr = envs->x_ctr[2];
    int i_prim = bas(2, i_sh);
    int j_prim = bas(2, j_sh);
    int k_prim = bas(2, k_sh);
    double *ai = env + bas(5, i_sh);
    double *aj = env + bas(5, j_sh);
    double *ak = env + bas(5, k_sh);
    double *ci = env + bas(6, i_sh);
    double *cj = env + bas(6, j_sh);
    double *ck = env + bas(6, k_sh);
    double expcutoff = envs->expcutoff;
    double rr_ij = (envs->rirj)[0] * (envs->rirj)[0] + (envs->rirj)[1] * (envs->rirj)[1] + (envs->rirj)[2] * (envs->rirj)[2];
    PairData *pdata_base = NULL, *pdata_ij = NULL;
    if (opt->pairdata != NULL)
    {
        pdata_base = opt->pairdata[i_sh * opt->nbas + j_sh];
    }
    else
    {
        double *log_maxci = opt->log_max_coeff[i_sh];
        double *log_maxcj = opt->log_max_coeff[j_sh];
        MALLOC_INSTACK(pdata_base, i_prim * j_prim, cache);
        if (set_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                         log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                         i_prim, j_prim, rr_ij, expcutoff, env))
        {
            return 0;
        }
    }
    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
    size_t nf = envs->nf;
    double fac1i, fac1j, fac1k;
    int ip, jp, kp;
    int _empty[4] = {1, 1, 1, 1};
    int *iempty = _empty + 0;
    int *jempty = _empty + 1;
    int *kempty = _empty + 2;
    int *gempty = _empty + 3;
    int *non0ctrk = NULL, *non0idxk = NULL;
    MALLOC_INSTACK(non0ctrk, k_prim + k_prim * k_ctr, cache);
    non0idxk = non0ctrk + k_prim;
    Opt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
    double expij, cutoff;
    double *rij;
    double *rkl = envs->rkl;
    int *idx = opt->index_xyz_array[envs->i_l * 16 * 16 + envs->j_l * 16 + envs->k_l];
    if (idx == NULL)
    {
        MALLOC_INSTACK(idx, nf * 3, cache);
        g2e_index_xyz(idx, envs);
    }
    double omega = env[8];
    if (omega < 0 && envs->rys_order > 1)
    {
        double r_guess = 8.;
        double omega2 = omega * omega;
        int lij = envs->li_ceil + envs->lj_ceil;
        if (lij > 0)
        {
            double dist_ij = sqrt(rr_ij);
            double aij = ai[i_prim - 1] + aj[j_prim - 1];
            double theta = omega2 / (omega2 + aij);
            expcutoff += lij * log(
                                   (dist_ij + theta * r_guess + 1.) / (dist_ij + 1.));
        }
        if (envs->lk_ceil > 0)
        {
            double theta = omega2 / (omega2 + ak[k_prim - 1]);
            expcutoff += envs->lk_ceil * log(theta * r_guess + 1.);
        }
    }
    int nc = 1;
    int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
    int len0 = envs->nf * n_comp;
    int len = leng + len0;
    double *g = NULL;
    MALLOC_INSTACK(g, len, cache);
    double *gout;
    if (n_comp == 1)
    {
        gout = gctr;
        gempty = empty;
    }
    else
    {
        gout = g + leng;
    }

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        fac1k = envs->common_factor * ck[kp];

        pdata_ij = pdata_base;
        for (jp = 0; jp < j_prim; jp++)
        {
            envs->aj[0] = aj[jp];
            fac1j = fac1k * cj[jp];
            for (ip = 0; ip < i_prim; ip++, pdata_ij++)
            {
                if (pdata_ij->cceij > expcutoff)
                {
                    continue;
                }
                envs->ai[0] = ai[ip];
                expij = pdata_ij->eij;
                rij = pdata_ij->rij;
                cutoff = expcutoff - pdata_ij->cceij;
                fac1i = fac1j * ci[ip] * expij;
                envs->fac[0] = fac1i;
                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs))
                {
                    (*envs->f_gout)(gout, g, idx, envs, *gempty);
                    *gempty = 0;
                }
            } // end loop i_prim
        } // end loop j_prim
    } // end loop k_prim

    if (n_comp > 1 && !*gempty)
    {
        TRANSPOSE(gout);
    }
    return !*empty;
}

int (*f_3c2e_loop[8])(double *, Env *, double *, int *) = {
    _3c2e_loop,
    _3c2e_loop,
    _3c2e_loop,
    _3c2e_n11_loop,
    _3c2e_loop,
    _3c2e_1n1_loop,
    _3c2e_loop,
    _3c2e_111_loop,
};

void Opt_log_max_pgto_coeff(double *log_maxc, double *coeff, int nprim, int nctr)
{
    int i, ip;
    double maxc;
    for (ip = 0; ip < nprim; ip++)
    {
        maxc = 0;
        for (i = 0; i < nctr; i++)
        {
            maxc = std::max(maxc, fabs(coeff[i * nprim + ip]));
        }
        log_maxc[ip] = log(maxc);
    }
}

int _3c2e_loop_nopt(double *gctr, Env *envs, double *cache, int *empty)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    double *env = envs->env;
    int i_sh = shls[0];
    int j_sh = shls[1];
    int k_sh = shls[2];
    int i_ctr = envs->x_ctr[0];
    int j_ctr = envs->x_ctr[1];
    int k_ctr = envs->x_ctr[2];
    int i_prim = bas(2, i_sh);
    int j_prim = bas(2, j_sh);
    int k_prim = bas(2, k_sh);
    double *ai = env + bas(5, i_sh);
    double *aj = env + bas(5, j_sh);
    double *ak = env + bas(5, k_sh);
    double *ci = env + bas(6, i_sh);
    double *cj = env + bas(6, j_sh);
    double *ck = env + bas(6, k_sh);

    double expcutoff = envs->expcutoff;
    const double rr_ij = (envs->rirj)[0] * (envs->rirj)[0] + (envs->rirj)[1] * (envs->rirj)[1] + (envs->rirj)[2] * (envs->rirj)[2];
    double *log_maxci = NULL, *log_maxcj;
    PairData *pdata_base = NULL, *pdata_ij;
    MALLOC_INSTACK(log_maxci, i_prim + j_prim, cache);
    MALLOC_INSTACK(pdata_base, i_prim * j_prim, cache);
    log_maxcj = log_maxci + i_prim;
    Opt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
    Opt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
    if (set_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                     log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                     i_prim, j_prim, rr_ij, expcutoff, env))
    {
        return 0;
    }

    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
    size_t nf = envs->nf;
    double fac1i, fac1j, fac1k;
    int ip, jp, kp;
    int _empty[4] = {1, 1, 1, 1};
    int *iempty = _empty + 0;
    int *jempty = _empty + 1;
    int *kempty = _empty + 2;
    int *gempty = _empty + 3;

    double expij, cutoff;
    double *rij;
    double *rkl = envs->rk;
    double omega = env[8];
    if (omega < 0 && envs->rys_order > 1)
    {
        double r_guess = 8.;
        double omega2 = omega * omega;
        int lij = envs->li_ceil + envs->lj_ceil;
        if (lij > 0)
        {
            double dist_ij = sqrt(rr_ij);
            double aij = ai[i_prim - 1] + aj[j_prim - 1];
            double theta = omega2 / (omega2 + aij);
            expcutoff += lij * log(
                                   (dist_ij + theta * r_guess + 1.) / (dist_ij + 1.));
        }
        if (envs->lk_ceil > 0)
        {
            double theta = omega2 / (omega2 + ak[k_prim - 1]);
            expcutoff += envs->lk_ceil * log(theta * r_guess + 1.);
        }
    }

    int *idx = NULL;
    MALLOC_INSTACK(idx, nf * 3, cache);
    g2e_index_xyz(idx, envs);

    int *non0ctri = NULL, *non0ctrj, *non0ctrk;
    int *non0idxi = NULL, *non0idxj, *non0idxk;
    MALLOC_INSTACK(non0ctri, i_prim + j_prim + k_prim + i_prim * i_ctr + j_prim * j_ctr + k_prim * k_ctr, cache);
    non0ctrj = non0ctri + i_prim;
    non0ctrk = non0ctrj + j_prim;
    non0idxi = non0ctrk + k_prim;
    non0idxj = non0idxi + i_prim * i_ctr;
    non0idxk = non0idxj + j_prim * j_ctr;
    Opt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
    Opt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);
    Opt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);

    int nc = i_ctr * j_ctr * k_ctr;
    // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
    int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
    int lenk = nf * nc * n_comp;            // gctrk
    int lenj = nf * i_ctr * j_ctr * n_comp; // gctrj
    int leni = nf * i_ctr * n_comp;         // gctri
    int len0 = nf * n_comp;                 // gout
    int len = leng + lenk + lenj + leni + len0;
    double *g = NULL;
    MALLOC_INSTACK(g, len, cache); // must be allocated last in this function
    double *g1 = g + leng;
    double *gout, *gctri, *gctrj, *gctrk;

    ALIAS_ADDR_IF_EQUAL(k, m);
    ALIAS_ADDR_IF_EQUAL(j, k);
    ALIAS_ADDR_IF_EQUAL(i, j);
    ALIAS_ADDR_IF_EQUAL(g, i);

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        if (k_ctr == 1)
        {
            fac1k = envs->common_factor * ck[kp];
        }
        else
        {
            fac1k = envs->common_factor;
            *jempty = 1;
        }

        pdata_ij = pdata_base;
        for (jp = 0; jp < j_prim; jp++)
        {
            envs->aj[0] = aj[jp];
            if (j_ctr == 1)
            {
                fac1j = fac1k * cj[jp];
            }
            else
            {
                fac1j = fac1k;
                *iempty = 1;
            }
            for (ip = 0; ip < i_prim; ip++, pdata_ij++)
            {
                if (pdata_ij->cceij > expcutoff)
                {
                    goto i_contracted;
                }
                envs->ai[0] = ai[ip];
                expij = pdata_ij->eij;
                rij = pdata_ij->rij;
                cutoff = expcutoff - pdata_ij->cceij;
                if (i_ctr == 1)
                {
                    fac1i = fac1j * ci[ip] * expij;
                }
                else
                {
                    fac1i = fac1j * expij;
                }
                envs->fac[0] = fac1i;
                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs))
                {
                    (*envs->f_gout)(gout, g, idx, envs, *gempty);
                    PRIM2CTR(i, gout, len0);
                }
            i_contracted:;
            } // end loop i_prim
            if (!*iempty)
            {
                PRIM2CTR(j, gctri, leni);
            }
        } // end loop j_prim
        if (!*jempty)
        {
            PRIM2CTR(k, gctrj, lenj);
        }
    } // end loop k_prim

    if (n_comp > 1 && !*kempty)
    {
        TRANSPOSE(gctrk);
    }
    return !*empty;
}

int _3c2e_drv(double *out, int *dims, Env *envs, Opt *opt,
              double *cache, void (*f_e1_c2s)(double *, double *, int *, Env *, double *), int is_ssc)
{
    int *x_ctr = envs->x_ctr;
    int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
    if (out == NULL)
    {
        PAIRDATA_NON0IDX_SIZE(pdata_size);
        int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
        int len0 = envs->nf * n_comp;
        int cache_size = std::max(leng + len0 + nc * n_comp * 3 + pdata_size,
                                  nc * n_comp + envs->nf * 3);
        return cache_size;
    }
    double *stack = NULL;
    if (cache == NULL)
    {
        PAIRDATA_NON0IDX_SIZE(pdata_size);
        int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
        int len0 = envs->nf * n_comp;
        int cache_size = std::max(leng + len0 + nc * n_comp * 3 + pdata_size,
                                  nc * n_comp + envs->nf * 3);
        stack = (double *)malloc(sizeof(double) * cache_size);
        cache = stack;
    }
    double *gctr = NULL;
    MALLOC_INSTACK(gctr, nc * n_comp, cache);

    int n;
    int empty = 1;
    if (opt != NULL)
    {
        envs->opt = opt;
        n = ((envs->x_ctr[0] == 1) << 2) + ((envs->x_ctr[1] == 1) << 1) + (envs->x_ctr[2] == 1);
        f_3c2e_loop[n](gctr, envs, cache, &empty);
    }
    else
    {
        _3c2e_loop_nopt(gctr, envs, cache, &empty);
    }

    int counts[4];
    if (f_e1_c2s == &c2s_sph_3c2e1)
    {
        counts[0] = (envs->i_l * 2 + 1) * x_ctr[0];
        counts[1] = (envs->j_l * 2 + 1) * x_ctr[1];
        if (is_ssc)
        {
            counts[2] = envs->nfk * x_ctr[2];
        }
        else
        {
            counts[2] = (envs->k_l * 2 + 1) * x_ctr[2];
        }
    }
    else
    {
        counts[0] = envs->nfi * x_ctr[0];
        counts[1] = envs->nfj * x_ctr[1];
        counts[2] = envs->nfk * x_ctr[2];
    }
    counts[3] = 1;
    if (dims == NULL)
    {
        dims = counts;
    }
    int nout = dims[0] * dims[1] * dims[2];
    if (!empty)
    {
        for (n = 0; n < n_comp; n++)
        {
            (*f_e1_c2s)(out + nout * n, gctr + nc * n, dims, envs, cache);
        }
    }
    else
    {
        for (n = 0; n < n_comp; n++)
        {
            c2s_dset0(out + nout * n, dims, counts);
        }
    }
    if (stack != NULL)
    {
        free(stack);
    }
    return !empty;
}

int int3c2e_sph(double *out, int *dims, int *shls, int *atm, int natm, int *bas, int nbas, double *env, Opt *opt, double *cache)
{
    int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
    Env envs;
    init_int3c2e_Env(&envs, ng, shls, atm, natm, bas, nbas, env);
    envs.f_gout = &gout2e;
    return _3c2e_drv(out, dims, &envs, opt, cache, &c2s_sph_3c2e1, 0);
};

void GTOnr3c_fill_s1(int (*intor)(double *, int *, int *, int *, int, int *, int, double *, Opt *, double *), double *out, double *buf,
                     int comp, int jobid,
                     int *shls_slice, int *ao_loc, Opt *cintopt,
                     int *atm, int natm, int *bas, int nbas, double *env)
{
    const int ish0 = shls_slice[0];
    const int ish1 = shls_slice[1];
    const int jsh0 = shls_slice[2];
    const int jsh1 = shls_slice[3];
    const int ksh0 = shls_slice[4];
    const int ksh1 = shls_slice[5];
    const int nksh = ksh1 - ksh0;

    const int ksh = jobid % nksh + ksh0;
    const int jstart = jobid / nksh * 8 + jsh0;
    const int jend = std::min(jstart + 8, jsh1);
    if (jstart >= jend)
    {
        return;
    }

    const int naoi = ao_loc[ish1] - ao_loc[ish0];
    const int naoj = ao_loc[jsh1] - ao_loc[jsh0];
    const int naok = ao_loc[ksh1] - ao_loc[ksh0];
    int dims[] = {naoi, naoj, naok};

    const int k0 = ao_loc[ksh] - ao_loc[ksh0];
    out += naoi * naoj * k0;

    int ish, jsh, i0, j0;
    int shls[3] = {0, 0, ksh};

    for (jsh = jstart; jsh < jend; jsh++)
    {
        for (ish = ish0; ish < ish1; ish++)
        {
            shls[0] = ish;
            shls[1] = jsh;
            i0 = ao_loc[ish] - ao_loc[ish0];
            j0 = ao_loc[jsh] - ao_loc[jsh0];
            (*intor)(out + j0 * naoi + i0, dims, shls, atm, natm, bas, nbas, env,
                     cintopt, buf);
        }
    }
}

int GTOmax_shell_dim(const int *ao_loc, const int *shls_slice, int ncenter)
{
    int i;
    int i0 = shls_slice[0];
    int i1 = shls_slice[1];
    int di = 0;
    for (i = 1; i < ncenter; i++)
    {
        i0 = std::min(i0, shls_slice[i * 2]);
        i1 = std::max(i1, shls_slice[i * 2 + 1]);
    }
    for (i = i0; i < i1; i++)
    {
        di = std::max(di, ao_loc[i + 1] - ao_loc[i]);
    }
    return di;
}

void GTOnr3c_drv(int (*intor)(double *, int *, int *, int *, int, int *, int, double *, Opt *, double *),
                 void (*fill)(int (*intor)(double *, int *, int *, int *, int, int *, int, double *, Opt *, double *),
                              double *, double *, int, int, int *, int *, Opt *, int *, int, int *, int, double *),
                 double *eri, int comp,
                 int *shls_slice, int *ao_loc, Opt *cintopt,
                 int *atm, int natm, int *bas, int nbas, double *env)
{
    const int ish0 = shls_slice[0];
    const int ish1 = shls_slice[1];
    const int jsh0 = shls_slice[2];
    const int jsh1 = shls_slice[3];
    const int ksh0 = shls_slice[4];
    const int ksh1 = shls_slice[5];
    const int nish = ish1 - ish0;
    const int njsh = jsh1 - jsh0;
    const int nksh = ksh1 - ksh0;
    const int di = GTOmax_shell_dim(ao_loc, shls_slice, 3);
    const int cache_size = GTOmax_cache_size(intor, shls_slice, 3,
                                             atm, natm, bas, nbas, env);
    const int njobs = (std::max(nish, njsh) / 8 + 1) * nksh;

#pragma omp parallel
    {
        int jobid;
        double *buf = (double *)malloc(sizeof(double) * (di * di * di * comp + cache_size));
#pragma omp for nowait schedule(dynamic)
        for (jobid = 0; jobid < njobs; jobid++)
        {
            (*fill)(intor, eri, buf, comp, jobid, shls_slice, ao_loc,
                    cintopt, atm, natm, bas, nbas, env);
        }
        free(buf);
    }
}

void init_2e_optimizer(Opt &opt, int *atm, int natm, int *bas, int nbas, double *env)
{
    opt.index_xyz_array = NULL;
    opt.non0ctr = NULL;
    opt.sortedidx = NULL;
    opt.nbas = nbas;
    opt.log_max_coeff = NULL;
    opt.pairdata = NULL;
}

void Opt_set_log_maxc(Opt &opt, int *atm, int natm,
                      int *bas, int nbas, double *env)
{
    int i, iprim, ictr;
    double *ci;
    size_t tot_prim = 0;
    for (i = 0; i < nbas; i++)
    {
        tot_prim += bas(2, i);
    }
    if (tot_prim == 0)
    {
        return;
    }

    opt.log_max_coeff = (double **)malloc(sizeof(double *) * std::max(nbas, 1));
    double *plog_maxc = (double *)malloc(sizeof(double) * tot_prim);
    opt.log_max_coeff[0] = plog_maxc;
    for (i = 0; i < nbas; i++)
    {
        iprim = bas(2, i);
        ictr = bas(3, i);
        ci = env + bas(6, i);
        opt.log_max_coeff[i] = plog_maxc;
        Opt_log_max_pgto_coeff(plog_maxc, ci, iprim, ictr);
        plog_maxc += iprim;
    }
}

void Opt_setij(Opt &opt, int *ng,
               int *atm, int natm, int *bas, int nbas, double *env)
{
    int i, j, ip, jp;
    int iprim, jprim, li, lj;
    double *ai, *aj, *ri, *rj;
    double expcutoff;
    if (env[0] == 0)
    {
        expcutoff = 60.;
    }
    else
    {
        expcutoff = std::max(40., env[0]);
    }

    if (opt.log_max_coeff == NULL)
    {
        Opt_set_log_maxc(opt, atm, natm, bas, nbas, env);
    }
    double **log_max_coeff = opt.log_max_coeff;
    double *log_maxci, *log_maxcj;

    size_t tot_prim = 0;
    for (i = 0; i < nbas; i++)
    {
        tot_prim += bas(2, i);
    }
    if (tot_prim == 0 || tot_prim > 2048)
    {
        return;
    }
    opt.pairdata = (PairData **)malloc(sizeof(PairData *) * std::max(nbas * nbas, 1));
    PairData *pdata = (PairData *)malloc(sizeof(PairData) * tot_prim * tot_prim);
    opt.pairdata[0] = pdata;

    int ijkl_inc;
    if ((ng[0] + ng[1]) > (ng[2] + ng[3]))
    {
        ijkl_inc = ng[0] + ng[1];
    }
    else
    {
        ijkl_inc = ng[2] + ng[3];
    }

    int empty;
    double rr;
    PairData *pdata0;
    for (i = 0; i < nbas; i++)
    {
        ri = env + atm(1, bas(0, i));
        ai = env + bas(5, i);
        iprim = bas(2, i);
        li = bas(1, i);
        log_maxci = log_max_coeff[i];

        for (j = 0; j <= i; j++)
        {
            rj = env + atm(1, bas(0, j));
            aj = env + bas(5, j);
            jprim = bas(2, j);
            lj = bas(1, j);
            log_maxcj = log_max_coeff[j];
            rr = (ri[0] - rj[0]) * (ri[0] - rj[0]) + (ri[1] - rj[1]) * (ri[1] - rj[1]) + (ri[2] - rj[2]) * (ri[2] - rj[2]);

            empty = set_pairdata(pdata, ai, aj, ri, rj, log_maxci, log_maxcj,
                                 li + ijkl_inc, lj, iprim, jprim, rr, expcutoff, env);
            if (i == 0 && j == 0)
            {
                opt.pairdata[0] = pdata;
                pdata += iprim * jprim;
            }
            else if (!empty)
            {
                opt.pairdata[i * nbas + j] = pdata;
                pdata += iprim * jprim;
                if (i != j)
                {
                    opt.pairdata[j * nbas + i] = pdata;
                    pdata0 = opt.pairdata[i * nbas + j];
                    // transpose pairdata
                    for (ip = 0; ip < iprim; ip++)
                    {
                        for (jp = 0; jp < jprim; jp++, pdata++)
                        {
                            memcpy(pdata, pdata0 + jp * iprim + ip, sizeof(PairData));
                        }
                    }
                }
            }
            else
            {
                opt.pairdata[i * nbas + j] = ((PairData *)0xffffffffffffffffuL);
                opt.pairdata[j * nbas + i] = ((PairData *)0xffffffffffffffffuL);
            }
        }
    }
}

void Opt_set_non0coeff(Opt &opt, int *atm, int natm,
                       int *bas, int nbas, double *env)
{
    int i, iprim, ictr;
    double *ci;
    size_t tot_prim = 0;
    size_t tot_prim_ctr = 0;
    for (i = 0; i < nbas; i++)
    {
        tot_prim += bas(2, i);
        tot_prim_ctr += bas(2, i) * bas(3, i);
    }
    if (tot_prim == 0)
    {
        return;
    }

    opt.non0ctr = (int **)malloc(sizeof(int *) * std::max(nbas, 1));
    opt.sortedidx = (int **)malloc(sizeof(int *) * std::max(nbas, 1));
    int *pnon0ctr = (int *)malloc(sizeof(int) * tot_prim);
    int *psortedidx = (int *)malloc(sizeof(int) * tot_prim_ctr);
    opt.non0ctr[0] = pnon0ctr;
    opt.sortedidx[0] = psortedidx;
    for (i = 0; i < nbas; i++)
    {
        iprim = bas(2, i);
        ictr = bas(3, i);
        ci = env + bas(6, i);
        opt.non0ctr[i] = pnon0ctr;
        opt.sortedidx[i] = psortedidx;
        Opt_non0coeff_byshell(psortedidx, pnon0ctr, ci, iprim, ictr);
        pnon0ctr += iprim;
        psortedidx += iprim * ictr;
    }
}

int _make_fakebas(int *fakebas, int *bas, int nbas, double *env)
{
    int i;
    int max_l = 0;
    for (i = 0; i < nbas; i++)
    {
        max_l = std::max(max_l, bas(1, i));
    }

    int fakenbas = max_l + 1;
    for (i = 0; i < 8 * fakenbas; i++)
    {
        fakebas[i] = 0;
    }
    // fakebas only initializes ANG_OF, since the others does not
    // affect index_xyz
    for (i = 0; i <= max_l; i++)
    {
        fakebas[8 * i + 1] = i;
    }
    return max_l;
}

int *_allocate_index_xyz(Opt &opt, int max_l, int l_allow, int order)
{
    int i;
    int cumcart = (l_allow + 1) * (l_allow + 2) * (l_allow + 3) / 6;
    size_t ll = max_l + 1;
    size_t cc = cumcart;
    for (i = 1; i < order; i++)
    {
        ll *= 16;
        cc *= cumcart;
    }
    int *buf = (int *)malloc(sizeof(int) * cc * 3);
    int **ppbuf = (int **)malloc(sizeof(int *) * ll);
    ppbuf[0] = buf;
    for (i = 1; i < ll; i++)
    {
        ppbuf[i] = NULL;
    }
    opt.index_xyz_array = ppbuf;
    return buf;
}

void gen_idx(Opt &opt,
             void (*finit)(Env *, int *, int *, int *, int, int *, int, double *),
             void (*findex_xyz)(int *, const Env *),
             int order, int l_allow, int *ng,
             int *atm, int natm, int *bas, int nbas, double *env)
{
    int i, j, k, l, ptr;
    int fakebas[8 * 16];
    int max_l = _make_fakebas(fakebas, bas, nbas, env);
    int fakenbas = max_l + 1;
    // index_xyz bufsize may blow up for large max_l
    l_allow = std::min(max_l, l_allow);
    int *buf = _allocate_index_xyz(opt, max_l, l_allow, order);

    Env envs;
    int shls[4] = {
        0,
    };
    if (order == 2)
    {
        for (i = 0; i <= l_allow; i++)
        {
            for (j = 0; j <= l_allow; j++)
            {
                shls[0] = i;
                shls[1] = j;
                (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                ptr = i * 16 + j;
                opt.index_xyz_array[ptr] = buf;
                (*findex_xyz)(buf, &envs);
                buf += envs.nf * 3;
            }
        }
    }
    else if (order == 3)
    {
        for (i = 0; i <= l_allow; i++)
        {
            for (j = 0; j <= l_allow; j++)
            {
                for (k = 0; k <= l_allow; k++)
                {
                    shls[0] = i;
                    shls[1] = j;
                    shls[2] = k;
                    (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                    ptr = i * 16 * 16 + j * 16 + k;
                    opt.index_xyz_array[ptr] = buf;
                    (*findex_xyz)(buf, &envs);
                    buf += envs.nf * 3;
                }
            }
        }
    }
    else
    {
        for (i = 0; i <= l_allow; i++)
        {
            for (j = 0; j <= l_allow; j++)
            {
                for (k = 0; k <= l_allow; k++)
                {
                    for (l = 0; l <= l_allow; l++)
                    {
                        shls[0] = i;
                        shls[1] = j;
                        shls[2] = k;
                        shls[3] = l;
                        (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                        ptr = i * 16 * 16 * 16 + j * 16 * 16 + k * 16 + l;
                        opt.index_xyz_array[ptr] = buf;
                        (*findex_xyz)(buf, &envs);
                        buf += envs.nf * 3;
                    }
                }
            }
        }
    }
}

Opt int3c2e_optimizer(int *atm, int natm, int *bas, int nbas, double *env)
{
    int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
    Opt o;
    init_2e_optimizer(o, atm, natm, bas, nbas, env);
    Opt_setij(o, ng, atm, natm, bas, nbas, env);
    Opt_set_non0coeff(o, atm, natm, bas, nbas, env);
    gen_idx(o, &init_int3c2e_Env, &g2e_index_xyz,
            3, 12, ng, atm, natm, bas, nbas, env);
    // LEAKY! THIS OPTIMIZER ALLOCATES MEMORY; THIS WILL NEED CLEANUP LATER ON!
    return o;
}

Opt int2c2e_optimizer(int *atm, int natm, int *bas, int nbas, double *env)
{
    int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
    Opt o;
    init_2e_optimizer(o, atm, natm, bas, nbas, env);
    Opt_set_log_maxc(o, atm, natm, bas, nbas, env);
    Opt_set_non0coeff(o, atm, natm, bas, nbas, env);
    gen_idx(o, &init_int2c2e_Env, &g1e_index_xyz,
            2, 15, ng, atm, natm, bas, nbas, env);
    // LEAKY! THIS OPTIMIZER ALLOCATES MEMORY; THIS WILL NEED CLEANUP LATER ON!
    return o;
}