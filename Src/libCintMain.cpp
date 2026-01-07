#include "constants.h"
#include "libCintMain.h"
#include "libCintKernels.h"
//
#if defined(__APPLE__)
// On macOS we are using Accelerate for BLAS/LAPACK
#include <Accelerate/Accelerate.h>
#else
// Linux/Windows with oneMKL
#include <mkl.h>
#endif

// Function to compute three-center two-electron integrals (eri3c)
template <typename Kernel>
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

    ivec aoloc = Kernel::gen_loc(bas, nbas);
    int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
    int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];
    int naok = aoloc[shl_slice[5]] - aoloc[shl_slice[4]];

    CINTOpt* opty = nullptr;
    Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());

    // Compute integrals
    vec res((size_t)naoi * (size_t)naoj * (size_t)naok, 0.0);
    eri3c.resize((size_t)naoi * (size_t)naoj * (size_t)naok, 0.0);

    Kernel::drv(res.data(),
        1,
        shl_slice.data(),
        aoloc.data(),
        opty,
        atm.data(), nat,
        bas.data(), nbas,
        env.data());

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
template void computeEri3c<Coulomb3C_SPH>(Int_Params &param1,
                                        Int_Params &param2,
                                        vec& eri3c);

template <typename Kernel>
void compute2C(Int_Params& params, vec& ret) {
    ivec bas = params.get_bas();
    ivec atm = params.get_atm();
    vec env = params.get_env();

    int nbas = params.get_nbas();
    int nat = params.get_natoms();

    ivec shl_slice = { 0, nbas, 0, nbas };
    ivec aoloc = Kernel::gen_loc(bas, nbas);

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
template void compute2C<Coulomb2C_SPH>(Int_Params& params, vec& ret);
template void compute2C<Overlap2C_SPH>(Int_Params& params, vec& ret);
template void compute2C<Overlap2C_CRT>(Int_Params& params, vec& ret);


void calc_screend_functions_and_max_ij(
    const std::vector<atom>& atoms,
    const ivec& aoloc,
    const ivec& bas_orbital_indices,
    bvec2& screened,
    int& max_ij
) {
    const int natoms = static_cast<int>(atoms.size());
    if (natoms == 0) {
        screened.clear();
        max_ij = 0.0;
        return;
    }

    // Initialize screened matrix (natoms x natoms) only once
    screened.assign(natoms, bvec(natoms, false));
    max_ij = 0.0;

    // --- 1) Precompute worst exponents per atom (parallel) ---

    std::vector<double> worst_exp(natoms);

#pragma omp parallel for schedule(static)
    for (int i = 0; i < natoms; ++i) {
        const int nbas = atoms[i].get_basis_set_size();

        // Use a local variable; std::numeric_limits<double>::infinity() is more "semantic"
        double w = std::numeric_limits<double>::infinity();
        for (int s = 0; s < nbas; ++s) {
            const double e = atoms[i].get_basis_set_exponent(s);
            if (e < w) w = e;
        }
        worst_exp[i] = w;
    }

    const double exp_cutoff = 0.5 * constants::exp_cutoff;

    // --- 2) Loop over atom pairs and do screening + max-block computation (parallel) ---

    std::string output = "";    //Only used for debug output

    int max_block_ij = 0;
#pragma omp parallel
    {
        int local_max = 0;
        std::string local_output = "";

#pragma omp for schedule(dynamic) nowait
        for (int atom_i = 0; atom_i < natoms; ++atom_i) {
            for (int atom_j = atom_i; atom_j < natoms; ++atom_j) {
                const double dist = atoms[atom_i].distance_to(atoms[atom_j]);
                const double dist2 = dist * dist;

                const double crit = -dist2 * (worst_exp[atom_i] + worst_exp[atom_j]);
                if (crit < exp_cutoff) {
                    //local_output += "Screening atom pair (" + std::to_string(atom_i) + ", " + std::to_string(atom_j) + ") with distance " + std::to_string(dist) + " and criterion " + std::to_string(crit) + " < " + std::to_string(exp_cutoff) + "\n";
                    screened[atom_i][atom_j] = true;
                    continue;
                }

                const int bi = bas_orbital_indices[atom_i];
                const int bip1 = bas_orbital_indices[atom_i + 1];
                const int bj = bas_orbital_indices[atom_j];
                const int bjp1 = bas_orbital_indices[atom_j + 1];

                const int naoi = aoloc[bip1] - aoloc[bi];
                const int naoj = aoloc[bjp1] - aoloc[bj];

                const int block_ij = naoi * naoj;
                if (block_ij > local_max) {
                    local_max = block_ij;
                }
            }
        }
        #pragma omp critical
        {
            if (local_max > max_block_ij) {
                max_block_ij = local_max;
            }
            output += local_output;    //Only used for debug output
        }

    }
    max_ij = max_block_ij;    //Only used for debug output

    std::cout << output << std::flush;
    int skipped = std::accumulate(screened.begin(), screened.end(), 0,
        [](int sum, const bvec& row) {
            return sum + std::count(row.begin(), row.end(), true);
        });
    std::cout << "Screened out " << skipped << " atom pairs due to overlap criteria." << std::endl;
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

    //ivec aoloc = make_loc(bas, nbas);
    ivec aoloc = Kernel::gen_loc(bas, nbas);

    rho.resize(aoloc[nQM + nAux] - aoloc[nQM], 0.0);

    ivec bas_orbital_indices = generate_bas_indices_per_atom(normal_basis);
    ivec bas_aux_indices = generate_bas_indices_per_atom(aux_basis);

    bvec2 screened;
    int max_block_ij = 0;
    calc_screend_functions_and_max_ij(
        normal_basis.get_atoms(),
        aoloc,
        bas_orbital_indices,
        screened,
        max_block_ij
    );

    CINTOpt* opty = nullptr;
    Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());

#pragma omp parallel for schedule(dynamic)
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
            shl_slice[0] = bas_orbital_indices[atom_i];
            shl_slice[1] = bas_orbital_indices[atom_i + 1];
            const int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
            for (int atom_j = atom_i; atom_j < natoms; atom_j++) {
                if (screened[atom_i][atom_j]) continue;
                // Hoist weight calculation before kernel call
                const double weight = 2.0 - static_cast<double>(atom_i == atom_j);
                
                shl_slice[2] = bas_orbital_indices[atom_j];
                shl_slice[3] = bas_orbital_indices[atom_j + 1];

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
}
template void computeRho<Coulomb3C_SPH>(
    const Int_Params& normal_basis,
    const Int_Params& aux_basis,
    const dMatrix2& dm,
    vec& rho);
template void computeRho<Overlap3C_SPH>(
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
    //ivec aoloc = make_loc(bas, nbas);
    ivec aoloc = Kernel::gen_loc(bas, nbas);
    unsigned long long int naoi = aoloc[nQM] - aoloc[0];
    unsigned long long int naoj = aoloc[nQM] - aoloc[0];
    unsigned long long int naok = aoloc[nQM + nAux] - aoloc[nQM];
    eri3c.resize(naoi * naoj * naok, 0.0);
    CINTOpt* opty = nullptr;
    Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());
    ivec shl_slice = { 0, nQM, 0, nQM, nQM, nQM + nAux };
    Kernel::drv(eri3c.data(), 1, shl_slice.data(), aoloc.data(), opty, atm.data(), nat, bas.data(), nbas, env.data());
}
template void compute3C<Coulomb3C_SPH>(Int_Params& param1,
    Int_Params& param2,
    vec& eri3c);
template void compute3C<Overlap3C_SPH>(Int_Params& param1,
    Int_Params& param2,
    vec& eri3c);