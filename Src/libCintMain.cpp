#include "libCintMain.h"


#include "cint_funcs.h"
#include "constants.h"

extern "C" {
    extern CINTOptimizerFunction int3c2e_optimizer;
    extern CINTIntegralFunction int3c2e_sph;

    extern CINTOptimizerFunction int2c2e_optimizer;
    extern CINTIntegralFunction int2c2e_sph;

    extern CINTOptimizerFunction int1e_ovlp_optimizer;
    extern CINTIntegralFunction int1e_ovlp_sph;

    extern CINTOptimizerFunction int3c1e_optimizer;
    extern CINTIntegralFunction int3c1e_sph;
}

#include "nos_math.h"
//
#if defined(__APPLE__)
// On macOS we are using Accelerate for BLAS/LAPACK
#include <Accelerate/Accelerate.h>
#else
// Linux/Windows with oneMKL
#include <mkl.h>
#endif

#define BLKSIZE 8
/*
 * out[naoi,naoj,naok,comp] in F-order
 */
void GTOnr3c_fill_s1(
    int (*intor)(double*, int*, int*, int*, int, int*, int, double*, CINTOpt*, double*),
    double* out, double* buf, int comp, int jobid, int* shls_slice, int* ao_loc, CINTOpt* cintopt,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    const int ish0 = shls_slice[0];
    const int ish1 = shls_slice[1];
    const int jsh0 = shls_slice[2];
    const int jsh1 = shls_slice[3];
    const int ksh0 = shls_slice[4];
    const int ksh1 = shls_slice[5];
    const int nksh = ksh1 - ksh0;

    const int ksh = jobid % nksh + ksh0;
    const int jstart = jobid / nksh * BLKSIZE + jsh0;
    const int jend = std::min(jstart + BLKSIZE, jsh1);
    if (jstart >= jend)
    {
        return;
    }

    const size_t naoi = ao_loc[ish1] - ao_loc[ish0];
    const size_t naoj = ao_loc[jsh1] - ao_loc[jsh0];
    const size_t naok = ao_loc[ksh1] - ao_loc[ksh0];
    int dims[] = { (int)naoi, (int)naoj, (int)naok };

    const int k0 = ao_loc[ksh] - ao_loc[ksh0];
    out += naoi * naoj * k0;

    int ish, jsh, i0, j0;
    int shls[3] = { 0, 0, ksh };

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

size_t GTOmax_cache_size(
    int (*intor)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *),
    int *shls_slice, int ncenter,
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
    int (*f)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *) = (int (*)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *))intor;
    int cache_size = 0;
    int n;
    int shls[4];
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


void GTOnr3c_drv(
    int (*intor)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *),
    void (*fill)(int (*intor)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *), double *, double *, int, int, int *, int *, CINTOpt *, int *, int, int *, int, double *),
    double *eri, int comp, int *shls_slice, int *ao_loc, CINTOpt *cintopt,
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
    const int njobs = (std::max(nish, njsh) / BLKSIZE + 1) * nksh;

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


/*
 * mat(naoi,naoj,comp) in F-order
 */
void GTOint2c(int (*intor)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *), double *mat, int comp, int hermi,
              int *shls_slice, int *ao_loc, CINTOpt *opt,
              int *atm, int natm, int *bas, int nbas, double *env)
{
    const int ish0 = shls_slice[0];
    const int ish1 = shls_slice[1];
    const int jsh0 = shls_slice[2];
    const int jsh1 = shls_slice[3];
    const int nish = ish1 - ish0;
    const int njsh = jsh1 - jsh0;
    const size_t naoi = ao_loc[ish1] - ao_loc[ish0];
    const size_t naoj = ao_loc[jsh1] - ao_loc[jsh0];
    const int cache_size = GTOmax_cache_size(intor, shls_slice, 2,
                                             atm, natm, bas, nbas, env);
#pragma omp parallel
    {
        int dims[] = {(int)naoi, (int)naoj};
        int ish, jsh, ij, i0, j0;
        int shls[2];
        double *cache = (double *)malloc(sizeof(double) * cache_size);
#pragma omp for schedule(dynamic, 4)
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



// Function to compute three-center two-electron integrals (eri3c)
void computeEri3c(Int_Params &param1,
                  Int_Params &param2,
                  vec &eri3c)
{
    int nQM = param1.get_nbas();
    int nAux = param2.get_nbas();

    Int_Params combined(param1, param2);
    // combined.print_data("combined");

    ivec bas = combined.get_bas();
    ivec atm = combined.get_atm();
    vec env = combined.get_env();

    ivec shl_slice = {
        0,
        nQM,
        0,
        nQM,
        nQM,
        nQM + nAux,
    };

    int nat = combined.get_natoms();
    int nbas = combined.get_nbas();

    assert(shl_slice[1] <= nbas);
    assert(shl_slice[3] <= nbas);
    assert(shl_slice[5] <= nbas);

    ivec aoloc = make_loc(bas, nbas);
    int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
    int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];
    int naok = aoloc[shl_slice[5]] - aoloc[shl_slice[4]];

    CINTOpt *opty = nullptr;
    int3c2e_optimizer(&opty, atm.data(), nat, bas.data(), nbas, env.data());

    // Compute integrals
    vec res((size_t)naoi * (size_t)naoj * (size_t)naok, 0.0);
    eri3c.resize((size_t)naoi * (size_t)naoj * (size_t)naok, 0.0);

    GTOnr3c_drv(int3c2e_sph, GTOnr3c_fill_s1, res.data(), 1, shl_slice.data(), aoloc.data(), opty, atm.data(), nat, bas.data(), nbas, env.data());

    // FOR TESTING PURPOSES!!!!
    // GTOnr3c_drv(int3c2e_sph, res.data(), 1, shl_slice.data(), aoloc.data(), NULL, atm.data(), nat, bas.data(), nbas, env.data());

    // res is in fortran order, write the result in regular ordering
    for (int k = 0; k < naok; k++)
    {
        for (int j = 0; j < naoj; j++)
        {
            for (int i = 0; i < naoi; i++)
            {
                std::size_t idx_F = i + j * (size_t)naoi + k * ((size_t)naoi * (size_t)naoj);
                std::size_t idx_C = i * ((size_t)naoj * (size_t)naok) + j * (size_t)naok + k;
                eri3c[idx_C] = res[idx_F];
            }
        }
    }
}


// --- Traits for kernels ------------------------------------------------------
void Coulomb2C::optimizer(CINTOpt*& opt,
    int* atm, int nat, int* bas, int nbas, double* env) {
    int2c2e_optimizer(&opt, atm, nat, bas, nbas, env);
}
void Coulomb2C::drv(double* out, int comp, int* shl_slice, int* aoloc,
    CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env) {
    GTOint2c(int2c2e_sph, out, comp, 0, shl_slice, aoloc, opt, atm, nat, bas, nbas, env);
}

void Overlap2C::optimizer(CINTOpt*& opt,
    const int*, int, const int*, int, const double*) {
    opt = nullptr; // no optimizer needed for overlap
}
void Overlap2C::drv(double* out, int comp, int* shl_slice, int* aoloc,
    CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env) {
    GTOint2c(int1e_ovlp_sph, out, comp, 0, shl_slice, aoloc, opt, atm, nat, bas, nbas, env);
}

template <typename Kernel>
void compute2C(Int_Params& params, vec& ret) {
    ivec bas = params.get_bas();
    ivec atm = params.get_atm();
    vec env = params.get_env();

    int nbas = params.get_nbas();
    int nat = params.get_natoms();

    ivec shl_slice = { 0, nbas, 0, nbas };
    ivec aoloc = make_loc(bas, nbas);

    int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
    int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];

    CINTOpt* opty = nullptr;
    Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());

    // Compute integrals
    vec res((size_t)naoi * (size_t)naoj, 0.0);
    ret.resize((size_t)naoi * (size_t)naoj, 0.0);
    Kernel::drv(res.data(), 1, shl_slice.data(), aoloc.data(), opty, atm.data(), nat, bas.data(), nbas, env.data());

    // res is in fortran order, write the result in regular ordering
    for (int i = 0; i < naoi; i++)
    {
        for (int j = 0; j < naoj; j++)
        {
            ret[(size_t)j * (size_t)naoi + i] = res[(size_t)i * (size_t)naoj + j];
        }
    }
}
template void compute2C<Coulomb2C>(Int_Params& params, vec& ret);
template void compute2C<Overlap2C>(Int_Params& params, vec& ret);

// Function to calculate the number of 3-center 2-electron integrals to compute at once based on the available memory
// naoi = number of basis functions in the first shell
// naoj = number of basis functions in the second shell
// aoloc = Running total of functions to computed based on the order of the basis functions
// nQM = number of basis functions in the QM basis
// nAux = number of basis functions in the auxiliary basis
// max_RAM = maximum available memory in MB
// Returns the number of functions to compute at once
ivec calc_3c_steps(const unsigned long long int naoi, const unsigned long long int naoj, const ivec aoloc, const int nQM, const int nAux, const double max_mem)
{
    ivec steps = {nQM};

    // First check if the maximum memory is enough to compute all integrals at once
    // Calculate minimum memory needed for all integrals
    unsigned long long int naok_end = (size_t)(aoloc[(size_t)nQM + nAux]) - aoloc[nQM];
    double min_mem = static_cast<double>(sizeof(double) * naoi * naoj * naok_end) * 1e-6 + 200; // Small buffer of 200MB for other memory usage
    if (min_mem < max_mem)
    {
        steps.push_back(nQM + nAux);
        return steps;
    }

    // Calculate maximum number of basis functions for every iteration to stay under the memory limit
    int current_step = 0;
    unsigned long long int naok_max = static_cast<unsigned long long int>(max_mem / ((static_cast<double>(sizeof(double) * naoi * naoj)) * 1e-6));
    for (int bas_i = 1; bas_i <= nAux; bas_i++)
    {
        unsigned long long int naok = aoloc[nQM + bas_i] - aoloc[steps[current_step]];
        if (naok > naok_max)
        {
            steps.push_back(nQM + bas_i - 1);
            current_step += 1;
        }
    }
    steps.push_back(nQM + nAux);
    return steps;
}

// --- Traits for kernels ------------------------------------------------------
void Coulomb3C::optimizer(CINTOpt*& opt,
    int* atm, int nat, int* bas, int nbas, double* env) {
    int3c2e_optimizer(&opt, atm, nat, bas, nbas, env);
}
void Coulomb3C::drv(double* out, int comp, int* shl_slice, int* aoloc,
    CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env) {
    GTOnr3c_drv(int3c2e_sph, GTOnr3c_fill_s1, out, comp, shl_slice, aoloc, opt, atm, nat, bas, nbas, env);
}

void Overlap3C::optimizer(CINTOpt*& opt,
    const int*, int, const int*, int, const double*) {
    opt = nullptr;
}
void Overlap3C::drv(double* out, int comp, int* shl_slice, int* aoloc,
    CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env) {
    (void)opt;
    GTOnr3c_drv(int3c1e_sph, GTOnr3c_fill_s1, out, comp, shl_slice, aoloc, nullptr, atm, nat, bas, nbas, env);
}

inline double* aligned_alloc_d(size_t n)
{
    void* p = nullptr;

#if defined(__APPLE__)
    // POSIX (Linux, macOS)
    if (posix_memalign(&p, 64, n * sizeof(double)) != 0)
        throw std::bad_alloc();
#else
    p = MKL_malloc(n * sizeof(double), 64);
    if (!p) throw std::bad_alloc();
#endif
    return reinterpret_cast<double*>(p);
}

inline void aligned_free_d(double* p)
{
#if defined(__APPLE__)
    free(p);
#else
    MKL_free(p);
#endif
}

//This currently scales as O(N^2), if this becomes a bottleneck we can implement cell-lists or k-d trees
vec2 build_distance_list(const std::vector<atom>& atoms)
{
    std::vector<std::pair<int, int>> pairs;
    vec2 pair_list(atoms.size(), vec(atoms.size()));

    const int natoms = atoms.size();
#pragma omp parallel for collapse(2)
    for (int i = 0; i < natoms; i++) {
        for (int j = i+1 ; j < natoms; j++) {
            const double dist = atoms[i].distance_to(atoms[j]);
            pair_list[i][j] = dist;
            pair_list[j][i] = dist;
        }
    }
    return pair_list;
}

vec build_largest_exp(const std::vector<atom>& atoms)
{
    //fill vector with largest number
    vec res(atoms.size(), 2E100 );
    const int natoms = atoms.size();
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < natoms; i++) {
        for (int s = 0; s < atoms[i].get_basis_set_size(); s++) {
           const double exp = atoms[i].get_basis_set_exponent(s);
           res[i] = std::min(res[i], exp);
       }
    }
    return res;
}

//Ivec contains the ao indices for the given wave object
//The list contains (0, n_ao_atom1, last_item + n_ao_atom2, ...)
//So that ao indices for atom i are in [ao_indices_per_atom[i-1], ao_indices_per_atom[i])
ivec generate_bas_indices_per_atom(const Int_Params& params)
{
    const ivec bas = params.get_bas();
    int nbas = params.get_nbas();
    const int natoms = params.get_natoms();
    ivec bas_indices_per_atom(natoms, 0);
    for (int i = 0; i < nbas; i++)
    {
        bas_indices_per_atom[bas(ATOM_OF, i)]++;
    }

    ivec bas_indices_location(natoms + 1, 0);
    std::partial_sum(bas_indices_per_atom.begin(), bas_indices_per_atom.end(), bas_indices_location.begin() + 1);

    return std::move(bas_indices_location);
}


void get_matrix_slice(const dMatrix2& dm, std::pair<int, int> row_range, std::pair<int, int> col_range, vec& dm_slice, const double weight)
{
    const int nrows = row_range.second - row_range.first;
    const int ncols = col_range.second - col_range.first;
    
    // Direct indexing without mdspan overhead
    int idx = 0;
    for (int i = 0; i < nrows; i++)
    {
        const int row_idx = row_range.first + i;
        for (int j = 0; j < ncols; j++)
        {
            dm_slice[idx++] = dm(row_idx, col_range.first + j) * weight;
        }
    }
}

template <typename Kernel>
void computeRho(
    const Int_Params& normal_basis, 
    const Int_Params& aux_basis,
    const dMatrix2& dm,
    vec& rho)
{
    Int_Params combined(normal_basis,aux_basis);

    ivec bas = combined.get_bas();
    ivec atm = combined.get_atm();
    vec  env = combined.get_env();

    const int natoms = normal_basis.get_natoms();
    const int nQM = normal_basis.get_nbas();
    const int nAux = aux_basis.get_nbas();
    const int nat = combined.get_natoms();
    const int nbas = combined.get_nbas();

    ivec aoloc = make_loc(bas, nbas);

    rho.resize(aoloc[nQM + nAux] - aoloc[nQM], 0.0);

    CINTOpt* opty = nullptr;
    Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());

    ivec bas_orbital_indices = generate_bas_indices_per_atom(normal_basis);
    ivec bas_aux_indices = generate_bas_indices_per_atom(aux_basis);

    double exp_cutoff = 0.5*constants::exp_cutoff;


    vec2 distance_list = build_distance_list(normal_basis.get_atoms());
    vec worst_exponents = build_largest_exp(normal_basis.get_atoms());
    bvec2 screened(natoms, bvec(natoms, false));
    // Pre-compute maximum block size to avoid repeated allocations
    int max_block_ij = 0;
    for (int atom_i = 0; atom_i < natoms; atom_i++) {
        for (int atom_j = atom_i; atom_j < natoms; atom_j++) {
            if (-pow(distance_list[atom_i][atom_j], 2) * (worst_exponents[atom_i] + worst_exponents[atom_j]) < exp_cutoff) {
                //std::cout << "Pre-screening: " << atom_i << " : " << atom_j << " Distance: " << distance_list[atom_i][atom_j] << " Criteria: " << -pow(distance_list[atom_i][atom_j], 2) * (worst_exponents[atom_i] + worst_exponents[atom_j]) << " < " << exp_cutoff << std::endl;
                screened[atom_i][atom_j] = true;
                continue;
            }
            const int naoi = aoloc[bas_orbital_indices[atom_i + 1]] - aoloc[bas_orbital_indices[atom_i]];
            const int naoj = aoloc[bas_orbital_indices[atom_j + 1]] - aoloc[bas_orbital_indices[atom_j]];
            max_block_ij = std::max(max_block_ij, naoi * naoj);
        }
    }
//    int skipped = 0;
//#pragma omp parallel for schedule(dynamic)
    for (int atm_idx = 0; atm_idx < natoms; atm_idx++) {
        double* rho_atom = rho.data() + aoloc[nQM + bas_aux_indices[atm_idx]] - aoloc[nQM];
        int shl_slice[6] = {
            0, 0,  // i shells
            0,0,  // j shells
            nQM + bas_aux_indices[atm_idx], nQM + bas_aux_indices[atm_idx + 1] };

        int naok = aoloc[shl_slice[5]] - aoloc[shl_slice[4]];

        vec res(max_block_ij * naok);
        vec dm_slice(max_block_ij);
        for (int atom_i = 0; atom_i < natoms; atom_i++) {
            for (int atom_j = atom_i; atom_j < natoms; atom_j++) {
                //if (screened[atom_i][atom_j]) skipped++;
                if (screened[atom_i][atom_j]) continue;
                // Hoist weight calculation before kernel call
                const double weight = 2.0 - static_cast<double>(atom_i == atom_j);
                
                shl_slice[0] = bas_orbital_indices[atom_i];
                shl_slice[1] = bas_orbital_indices[atom_i + 1];
                shl_slice[2] = bas_orbital_indices[atom_j];
                shl_slice[3] = bas_orbital_indices[atom_j + 1];

                const int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
                const int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];
                const int block_ij = naoi * naoj;

                Kernel::drv(res.data(),
                    1,
                    shl_slice,
                    aoloc.data(),
                    opty,
                    atm.data(), nat,
                    bas.data(), nbas,
                    env.data());

                // Inline optimized matrix slice extraction
                const int row_start = aoloc[shl_slice[2]];
                const int row_end = aoloc[shl_slice[3]];
                const int col_start = aoloc[shl_slice[0]];
                const int col_end = aoloc[shl_slice[1]];
                
                int idx = 0;
                for (int i = row_start; i < row_end; i++) {
                    for (int j = col_start; j < col_end; j++) {
                        dm_slice[idx++] = dm(i, j) * weight;
                    }
                }

                // accumulate into rho[aux_ao0 .. aux_ao0+naok)
                cblas_dgemv(CblasRowMajor,
                    CblasNoTrans,
                    naok,
                    block_ij,
                    1.0,
                    res.data(),
                    block_ij,
                    dm_slice.data(),
                    1,
                    1.0,
                    rho_atom,
                    1);
            }
        }
    }
    if (opty) {
        delete opty;
        opty = nullptr;
    }
//    std::cout << "Skipped " << skipped << " blocks due to distance cutoff." << std::endl;
}
template void computeRho<Coulomb3C>(
    const Int_Params& normal_basis,
    const Int_Params& aux_basis,
    const dMatrix2& dm,
    vec& rho);
template void computeRho<Overlap3C>(
    const Int_Params& normal_basis,
    const Int_Params& aux_basis,
    const dMatrix2& dm,
    vec& rho);


template <typename Kernel>
void compute3C(Int_Params& param1,
    Int_Params& param2,
    vec& eri3c) {
    int nQM = param1.get_nbas();
    int nAux = param2.get_nbas();
    Int_Params combined(param1, param2);
    ivec bas = combined.get_bas();
    ivec atm = combined.get_atm();
    vec env = combined.get_env();
    int nat = combined.get_natoms();
    int nbas = combined.get_nbas();
    ivec aoloc = make_loc(bas, nbas);
    unsigned long long int naoi = aoloc[nQM] - aoloc[0];
    unsigned long long int naoj = aoloc[nQM] - aoloc[0];
    unsigned long long int naok = aoloc[nQM + nAux] - aoloc[nQM];
    eri3c.resize(naoi * naoj * naok, 0.0);
    CINTOpt* opty = nullptr;
    Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());
    ivec shl_slice = { 0, nQM, 0, nQM, nQM, nQM + nAux };
    Kernel::drv(eri3c.data(), 1, shl_slice.data(), aoloc.data(), opty, atm.data(), nat, bas.data(), nbas, env.data());
}
template void compute3C<Coulomb3C>(Int_Params& param1,
    Int_Params& param2,
    vec& eri3c);
template void compute3C<Overlap3C>(Int_Params& param1,
    Int_Params& param2,
    vec& eri3c);