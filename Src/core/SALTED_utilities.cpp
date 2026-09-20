#include "pch.h"
#include "SALTED_utilities.h"
#include "integration_params.h"
#include "SALTED_io.h"
#include "constants.h"
#include "atoms.h"
#include "cube.h"
#include "wfn_class.h"
#include "metatensor.hpp"
#include "featomic.hpp"
#include "aux_density.h"
#ifdef NOSPHERA2_USE_GPU
#include "aux_density_gpu.h"
#endif

std::vector<cvec2> SALTED_Utils::complex_to_real_transformation(std::vector<int> sizes)
{
    const double sqrt_2 = sqrt(2.0);
    using namespace std;
    vector<cvec2> matrices{};
    for (int i = 0; i < sizes.size(); i++)
    {
        int lval = (sizes[i] - 1) / 2;
        int st = (lval & 1) ? 1 : -1;

        cvec2 transformed_matrix(sizes[i], cvec(sizes[i], 0.0));
        for (int j = 0; j < lval; j++)
        {
            transformed_matrix[j][j] = complex<double>(0.0, 1.0);
            transformed_matrix[j][sizes[i] - j - 1] = complex<double>(0.0, st);
            transformed_matrix[sizes[i] - j - 1][j] = complex<double>(1.0, 0.0);
            transformed_matrix[sizes[i] - j - 1][sizes[i] - j - 1] = complex<double>(-st, 0.0);
            st = -st;
        }
        transformed_matrix[lval][lval] = sqrt_2;
        // Divide each element by sqrt(2.0)
        for (auto& row : transformed_matrix)
        {
            for (auto& elem : row)
            {
                elem /= sqrt_2;
            }
        }
        matrices.push_back(transformed_matrix);
    }
    return matrices;
}

int SALTED_Utils::get_lmax_max(std::unordered_map<std::string, int>& lmax)
{
    int lmax_max = 0;
    for (auto& [key, value] : lmax)
    {
        if (value > lmax_max)
        {
            lmax_max = value;
        }
    }
    return lmax_max;
}

void SALTED_Utils::set_lmax_nmax(std::unordered_map<std::string, int>& lmax, std::unordered_map<std::string, int>& nmax, const std::array<std::vector<primitive>, 118>& basis_set, std::vector<std::string> species)
{
    // lmax = {"C": 5, "H":2,...} with the numbers beeing the maximum angular momentum (type) for the given atom
    // nmax = {C0: 10, C1: 7, ...} with the numbers beeing the maximum number of primitives for the given atom and type

    for (auto& spe : species)
    {
        int atom_index = constants::get_Z_from_label(spe.c_str());
        // get the last element of the basis set for the given atom
        lmax[spe] = basis_set[atom_index].back().get_type();
        // initialize nmax with symbol + type
        for (int i = 0; i < basis_set[atom_index].back().get_type() + 1; i++)
        {
            nmax[spe + std::to_string(i)] = 0;
        }
        // count the number of primitives for the given atom and type
        for (int i = 0; i < basis_set[atom_index].size(); i++)
        {
            nmax[spe + std::to_string(basis_set[atom_index][i].get_type())] += 1;
        }
    }
}


std::vector<char> SALTED_Utils::filter_input(WFN& wavy, options& opt, const SALTEDConfig& config) {
    // Two kinds of atom cannot be predicted, and both are handed to the spherical
    // Thakkar fill instead of guessed at:
    //
    //  1. a species the model was never trained on;
    //  2. an atom with nothing inside the descriptor cutoff.
    //
    // Removing them disturbs nobody: the test is symmetric, so an atom with
    // nothing within rcut is also nobody's neighbour within rcut, and every other
    // atom's environment is exactly what it was.
    //
    // The test here is purely geometric, while featomic additionally ignores
    // neighbours whose species is outside neighspe. An atom with neighbours but no
    // ALLOWED ones therefore still reaches equicomb - the zero guard there catches
    // it, leaves it spherical rather than NaN, and says so in the log.
    const int ncen_in = wavy.get_ncen();
    std::vector<char> use_thakkar(ncen_in, 0);
    for (int a = 0; a < ncen_in; a++)
        if (std::find(config.species.begin(), config.species.end(), std::string(constants::atnr2letter(wavy.get_atom_charge(a)))) == config.species.end())
            use_thakkar[a] = 1;
    const int n_unknown = static_cast<int>(std::count(use_thakkar.begin(), use_thakkar.end(), (char)1));

    const double rcut = std::min(config.rcut1, config.rcut2);
    // rcut is in Angstrom, the coordinates may not be
    const double rcut_internal = wavy.get_isBohr() ? constants::ang2bohr(rcut) : rcut;
    const double cut_sq = rcut_internal * rcut_internal;
    int n_isolated = 0;
#pragma omp parallel for reduction(+ : n_isolated)
    for (int a = 0; a < ncen_in; a++)
    {
        if (use_thakkar[a]) continue;
        bool lonely = true;
        for (int b = 0; b < ncen_in && lonely; b++)
        {
            if (b == a) continue;
            double d_sq = 0.0;
            for (unsigned int ax = 0; ax < 3; ax++)
            {
                const double dx = wavy.get_atom_coordinate(a, ax) - wavy.get_atom_coordinate(b, ax);
                d_sq += dx * dx;
            }
            if (d_sq < cut_sq) lonely = false;
        }
        if (lonely)
        {
            use_thakkar[a] = 1;   // distinct indices, and char so there is no bitfield to race on
            ++n_isolated;
        }
    }

    if (n_unknown + n_isolated > 0)
    {
        if (n_unknown > 0)
        {
            std::cout << "WARNING: Not all species in the structure are known to the model. The following species are not known: ";
            for (int a = 0; a < ncen_in; a++)
            {
                if (std::find(config.species.begin(), config.species.end(), std::string(constants::atnr2letter(wavy.get_atom_charge(a)))) == config.species.end())
                {
                    std::cout << constants::atnr2letter(wavy.get_atom_charge(a)) << " ";
                }
            }
            std::cout << std::endl;
        }
        if (n_isolated > 0)
        {
            std::cout << "WARNING: " << n_isolated << " atom(s) have no neighbour within the "
                << rcut << " A descriptor cutoff, so there is no environment to predict from."
                << " Isolated solvent is the usual cause." << std::endl;
        }
        std::cout << "I will fill out these atoms using spherical Thakkar densities!\n";
        // make a copy of initial wavefunction, to leave the initial one untouched!
        for (int a = ncen_in - 1; a >= 0; a--)
        {
            if (use_thakkar[a])
            {
                wavy.erase_atom(a);
            }
        }
        // remove all known basis sets, to not get problems with newly loaded ones
        for (int a = 0; a < wavy.get_ncen(); a++)
        {
            wavy.clear_atom_basis_set(a);
        }
        //std::filesystem::path new_fn = wavy.get_path().parent_path() / "SALTED_temp.xyz"; //I think this is not actually neccecary....
        //wavy.write_xyz(new_fn);
        //wavy.set_path(new_fn);
        opt.needs_Thakkar_fill = true;
        return use_thakkar;
    }
    return {};
}

std::string SALTED_Utils::FeatomicHyperParameters::to_json() const
{
    std::ostringstream oss;
    oss << "{\n"
        << "  \"cutoff\": {  \
                    \"radius\": " << this->cutoff_radius << " , \"smoothing\": \
                                  {\"type\": \"" << this->cutoff_function.type << "\", \"width\": " << this->cutoff_function.width << "} }, \n"
        << "  \"density\": { \
                    \"type\": \"Gaussian\", \"width\": " << this->atomic_gaussian_width << ", \"center_atom_weight\": " << this->center_atom_weight << "},\n"
        << "  \"basis\": { \
                    \"type\": \"TensorProduct\", \"max_angular\": " << this->max_angular << ", \"radial\": {\"type\":  \"" << this->radial_basis.type << "\", \"max_radial\": " << this->max_radial << "} , \"spline_accuracy\": " << this->radial_basis.spline_accuracy << "}\n"
        << "}";
    return oss.str();
}


// Used to generate metatensor::TensorMap and save the buffer location into the descriptor_buffer
static metatensor::TensorMap get_feats_projs(featomic::SimpleSystem featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters)
{
    // size_t nspe1 = neighspe.size();
    std::vector<std::array<int32_t, 4>> keys_array;
    keys_array.reserve((parameters.max_angular + 1) * parameters.species.size() * parameters.neighspe.size());

    for (int l = 0; l < parameters.max_angular + 1; ++l)
    {
        for (const std::string& center_spe : parameters.species)
        {
            int32_t center_z = constants::get_Z_from_label(center_spe.c_str()) + 1;
            for (const std::string& neigh_spe : parameters.neighspe)
            {
                int32_t neigh_z = constants::get_Z_from_label(neigh_spe.c_str()) + 1;
                // Directly emplace back initializer_lists into keys_array
                keys_array.push_back({ l, 1, center_z, neigh_z});
            }
        }
    }

    // Assuming metatensor::Labels expects a flat sequence of integers for each label
    std::vector<int32_t> flattened_keys;
    for (const auto& subVector : keys_array)
    {
        flattened_keys.insert(flattened_keys.end(), subVector.begin(), subVector.end());
    }

    // Convert keys_array to rascaline::Labels
    std::vector<std::string> names = { "o3_lambda", "o3_sigma", "center_type", "neighbor_type" };
    metatensor::Labels keys_selection(names, flattened_keys.data(), flattened_keys.size() / names.size());


    //create the calculator with its name and parameters
    //Do not ask me, why Featomic expects the max_radial to be one less than the actual number of radial basis functions, but it does, so here we are
    SALTED_Utils::FeatomicHyperParameters modif_param = parameters;
    modif_param.max_radial -= 1;
    //featomic's rayon pool is built on first use from RAYON_NUM_THREADS and otherwise
    //takes every logical core, ignoring -cpus and OMP_NUM_THREADS (48 threads on an
    //8-thread run, 51 s of futex in the test suite)
#ifdef _OPENMP
    if (std::getenv("RAYON_NUM_THREADS") == nullptr) { // Flawfinder: ignore
        const std::string n = std::to_string(omp_get_max_threads());
#ifdef _WIN32
        _putenv_s("RAYON_NUM_THREADS", n.c_str());
#else
        setenv("RAYON_NUM_THREADS", n.c_str(), 1);
#endif
    }
#endif
    auto calculator = featomic::Calculator("spherical_expansion", modif_param.to_json().c_str());

    featomic::CalculationOptions calc_opts;
    calc_opts.selected_keys = keys_selection;
    calc_opts.use_native_system = true;
    // run the calculation
    metatensor::TensorMap descriptor = calculator.compute(featomic_system, calc_opts);

    // The descriptor is a metatensor `TensorMap`, containing multiple blocks.
    // We can transform it to a single block containing a dense representation,
    // with one sample for each atom-centered environment.
    descriptor = descriptor.keys_to_samples("center_type");
    descriptor = descriptor.keys_to_properties("neighbor_type");
    // descriptor.save("spx_pred.npy");

    return descriptor;
}

// Reads the descriptor buffer and fills the expansion coefficients vector
static SALTEDDescriptors get_expansion_coeffs(std::vector<uint8_t> descriptor_buffer, const featomic::SimpleSystem& featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters)
{
    int n_atoms = (int)featomic_system.size();
    int nspe = (int)parameters.species.size();
    metatensor::TensorMap descriptor = metatensor::TensorMap::load_buffer(descriptor_buffer);
    const int nchannels = nspe * parameters.max_radial;
    SALTEDDescriptors omega(n_atoms, nchannels, parameters.max_angular);
    for (int l = 0; l < parameters.max_angular + 1; ++l)
    {
        cvec2 c2r = SALTED_Utils::complex_to_real_transformation({ (2 * l) + 1 })[0];
        metatensor::TensorBlock descriptor_block = descriptor.block_by_id(l);
        metatensor::NDArray<double> descriptor_values = descriptor_block.values();

        // descriptor_values is contiguous in its property (d) dimension. The
        // packed l slab keeps the corresponding output m values contiguous.
        for (int a = 0; a < n_atoms; ++a)
        {
            for (int r = 0; r < 2 * l + 1; ++r)
            {
                for (int d = 0; d < nchannels; ++d)
                {
                    cdouble* output = omega.block(a, d, l);
                    const double value = descriptor_values(a, r, d);
                    for (int c = 0; c < 2 * l + 1; ++c)
                    {
                        output[c] += conj(c2r[r][c]) * value;
                    }
                }
            }
        }
        c2r.clear();
    }

    return omega;
}


SALTEDDescriptors SALTED_Utils::calculate_SALTED_descriptors(const featomic::SimpleSystem& featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters)
{
    metatensor::TensorMap descriptor = get_feats_projs(featomic_system, parameters);
    std::vector<uint8_t> descriptor_buffer = descriptor.save_buffer();
    return get_expansion_coeffs(descriptor_buffer, featomic_system, parameters);
}


namespace
{
    // Constructing a calculator splines the radial integral for every (n, l)
    // pair to `spline_accuracy`, which depends only on the hyperparameters, so
    // it is computed once and kept.
    //
    // Keyed on the parameter JSON, so a caller that changes any hyperparameter
    // gets a new calculator rather than a silently wrong one -- a descriptor of
    // the right shape computed with the wrong settings.
    //
    // `featomic::Calculator` is move-only, hence the indirection, and it is not
    // safe to use one calculator from several threads at once. The only caller
    // of `calculate_SOAP_Powerspectrum` is the serial descriptor path; SALTED
    // builds its own calculator in `get_feats_projs`.
    featomic::Calculator& cached_calculator(const std::string& name, const std::string& json)
    {
        static std::map<std::string, std::unique_ptr<featomic::Calculator>> cache;
        const std::string key = name + '\n' + json;
        auto found = cache.find(key);
        if (found == cache.end())
        {
            found = cache.emplace(key, std::make_unique<featomic::Calculator>(name, json)).first;
        }
        return *found->second;
    }
}

//FEATOMIC POWER Spectrum
metatensor::TensorMap SALTED_Utils::calculate_SOAP_Powerspectrum(featomic::SimpleSystem featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters) {
    // Phase timings, off unless NOSPHERA2_TIME_SOAP is set. Caching the
    // calculator does not remove the fixed per-call cost; it is in what follows.
    const bool time_phases = std::getenv("NOSPHERA2_TIME_SOAP") != nullptr; // Flawfinder: ignore
    auto mark = std::chrono::steady_clock::now();
    auto lap = [&mark](const char* what, bool on) {
        if (!on) return;
        const auto now = std::chrono::steady_clock::now();
        std::cout << "  SOAP_PHASE " << what << " "
                  << std::chrono::duration<double>(now - mark).count() << std::endl;
        mark = now;
    };

    // Built once per set of hyperparameters and kept: see `cached_calculator`.
    auto& calculator = cached_calculator("soap_power_spectrum", parameters.to_json());
    lap("calculator", time_phases);

    std::vector<std::array<int32_t,3>> keys_array;
    for (const std::string& center_type : parameters.species)
    {
        int32_t z_center = constants::get_Z_from_label(center_type.c_str()) + 1;

        for (size_t i = 0; i < parameters.neighspe.size(); ++i)
        {
            int32_t z1 = constants::get_Z_from_label(parameters.neighspe[i].c_str()) + 1;

            for (size_t j = i; j < parameters.neighspe.size(); ++j)
            {
                int32_t z2 = constants::get_Z_from_label(parameters.neighspe[j].c_str()) + 1;

                keys_array.push_back({ z_center, z1, z2 });
            }
        }
    }

    // Assuming metatensor::Labels expects a flat sequence of integers for each label
    std::vector<int32_t> flattened_keys;
    for (const auto& subVector : keys_array)
    {
        flattened_keys.insert(flattened_keys.end(), subVector.begin(), subVector.end());
    }

    // Convert keys_array to rascaline::Labels
    std::vector<std::string> names = {"center_type", "neighbor_1_type", "neighbor_2_type"};
    metatensor::Labels keys_selection(names, flattened_keys.data(), flattened_keys.size() / names.size());

    featomic::CalculationOptions calc_opts;
    calc_opts.use_native_system = true;
    calc_opts.selected_keys = keys_selection;
    // run the calculation
    // Initialize descriptor directly from computation result
    lap("keys", time_phases);
    metatensor::TensorMap descriptor = calculator.compute(featomic_system, calc_opts);
    lap("compute", time_phases);

    // The descriptor is a metatensor `TensorMap`, containing multiple blocks.
    // We can transform it to a single block containing a dense representation,
    // with one sample for each atom-centered environment.
    descriptor = descriptor.keys_to_samples("center_type");
    lap("keys_to_samples", time_phases);
    descriptor = descriptor.keys_to_properties(svec{ "neighbor_1_type" , "neighbor_2_type" });
    lap("keys_to_properties", time_phases);

    return descriptor;
}

aux_density_table::aux_density_table(const std::vector<atom>& atoms)
{
    n_at = (int)atoms.size();
    cx.resize(n_at), cy.resize(n_at), cz.resize(n_at), r2_max.resize(n_at), Z.resize(n_at), sh_start.resize(n_at + 1);
    for (int a = 0; a < n_at; a++) {
        cx[a] = atoms[a].get_coordinate(0), cy[a] = atoms[a].get_coordinate(1), cz[a] = atoms[a].get_coordinate(2), Z[a] = atoms[a].get_charge() - atoms[a].get_ECP_electrons();
        const std::vector<unsigned int> sc = atoms[a].get_shellcount();
        double alpha_min = DBL_MAX;
        sh_start[a] = n_sh;
        int prim = 0;
        for (int s = 0; s < (int)sc.size(); s++) {
            const int l = atoms[a].get_basis_set_type(prim);
            err_checkf(l <= 8, "Aux basis shells above l = 8 are not supported on the grid", std::cout);
            sh_l.push_back(l), pr_start.push_back(n_pr), coef_off.push_back(n_coef), sh_atom.push_back(a);
            
            for (int m = -l; m <= l; ++m)
            {
                coef_shell.push_back(n_sh);
                coef_m.push_back(m);
            }

            vec exponents(sc[s]);
            vec coefficients(sc[s]);
            for (int p = 0; p < sc[s]; ++p, prim++)
            {
                const primitive& pr = atoms[a].get_basis_set_entry(prim).get_primitive();
                exponents[p] = pr.get_exp();
                coefficients[p] = pr.get_coef();

                pr_exp_l32.push_back(std::pow(exponents[p], l + 1.5));
                alpha_min = std::min(alpha_min, exponents[p]);
                pr_exp.push_back(exponents[p]);
                int slot = 0;
                while (slot < (int)uniq_exp.size() && !(uniq_exp[slot] == exponents[p] && uniq_l[slot] == l)) slot++;
                if (slot == (int)uniq_exp.size())
                    uniq_exp.push_back(exponents[p]), uniq_l.push_back(l), uniq_exp_l32.push_back(pr_exp_l32.back());
                pr_uniq.push_back(slot);
            }
            coefficients = Int_Params::normalize_gto(coefficients, exponents, l);

            pr_norm.insert(
                pr_norm.end(),
                coefficients.begin(),
                coefficients.end());

            n_pr += sc[s], n_coef += 2 * l + 1, n_sh++;
        }
        r2_max[a] = 46.0517 / alpha_min;
    }
    sh_start[n_at] = n_sh, pr_start.push_back(n_pr);

    err_checkf(
        static_cast<int>(sh_l.size()) == n_sh &&
        static_cast<int>(sh_atom.size()) == n_sh &&
        static_cast<int>(coef_off.size()) == n_sh &&
        static_cast<int>(pr_start.size()) == n_sh + 1,
        "Invalid auxiliary shell table",
        std::cout);

    err_checkf(
        static_cast<int>(pr_exp.size()) == n_pr &&
        static_cast<int>(pr_norm.size()) == n_pr &&
        static_cast<int>(pr_exp_l32.size()) == n_pr,
        "Invalid auxiliary primitive table",
        std::cout);

    err_checkf(
        static_cast<int>(coef_shell.size()) == n_coef &&
        static_cast<int>(coef_m.size()) == n_coef,
        "Invalid auxiliary coefficient table",
        std::cout);
}

double aux_density_table::operator()(const double x, const double y, const double z, const double* coefs) const
{
    return aux_density::at(x, y, z, n_at, cx.data(), cy.data(), cz.data(), r2_max.data(), sh_start.data(), sh_l.data(), pr_start.data(), coef_off.data(), pr_exp.data(), pr_norm.data(), coefs);
}

double aux_density_table::operator()(const double x, const double y, const double z, const double* coefs, double& gx, double& gy, double& gz) const
{
    return aux_density::at_grad(x, y, z, n_at, cx.data(), cy.data(), cz.data(), r2_max.data(), sh_start.data(), sh_l.data(), pr_start.data(), coef_off.data(), pr_exp.data(), pr_norm.data(), coefs, gx, gy, gz);
}

double aux_density_table::operator()(const double x, const double y, const double z, const double* coefs, double& gx, double& gy, double& gz, double& lap) const
{
    return aux_density::at_lap(x, y, z, n_at, cx.data(), cy.data(), cz.data(), r2_max.data(), sh_start.data(), sh_l.data(), pr_start.data(), coef_off.data(), pr_exp.data(), pr_norm.data(), coefs, gx, gy, gz, lap);
}

double aux_density_table::operator()(const double x, const double y, const double z, const double* coefs, double& gx, double& gy, double& gz, double* H) const
{
    return aux_density::at_hess(x, y, z, n_at, cx.data(), cy.data(), cz.data(), r2_max.data(), sh_start.data(), sh_l.data(), pr_start.data(), coef_off.data(), pr_exp.data(), pr_norm.data(), coefs, gx, gy, gz, H);
}

static inline cdouble apply_i_to_l(
    const int l,
    const double value)
{
    switch (l & 3)
    {
    case 0:
        return cdouble(value, 0.0);

    case 1:
        return cdouble(0.0, value);

    case 2:
        return cdouble(-value, 0.0);

    case 3:
        return cdouble(0.0, -value);
    }

    return constants::cnull;
}

cdouble aux_density_table::fourier_atom(
    const double kx,
    const double ky,
    const double kz,
    const double* coefs,
    const int atom_idx) const
{
    err_checkf(
        atom_idx >= 0 && atom_idx < n_at,
        "Invalid atom index in aux_density_table::fourier_atom",
        std::cout
    );

    const double H2 =
        kx * kx +
        ky * ky +
        kz * kz;

    const double H = std::sqrt(H2);

    double k[4];

    if (H > 0.0) [[likely]]
    {
        k[0] = kx / H;
        k[1] = ky / H;
        k[2] = kz / H;
    }
    else
    {
        k[0] = 0.0;
        k[1] = 0.0;
        k[2] = 1.0;
    }

    k[3] = H;

    cdouble sf = constants::cnull;

    for (int s = sh_start[atom_idx]; s < sh_start[atom_idx + 1]; ++s)
    {
        const int l = sh_l[s];

        if (H == 0.0 && l > 0)
            continue;

        double Hl_over_2l = 1.0;

        for (int i = 0; i < l; ++i)
            Hl_over_2l *= 0.5 * H;

        double radial = 0.0;

        for (int p = pr_start[s]; p < pr_start[s + 1]; ++p)
        {
            radial += pr_norm[p] * Hl_over_2l * std::exp(-H2 / (4.0 * pr_exp[p])) / pr_exp_l32[p];
        }

        const double angular = constants::spherical_harmonic(l, k[0], k[1], k[2], coefs + coef_off[s]);

        sf += apply_i_to_l(l, constants::PI3_2 * radial * angular);
    }

    return sf;
}

double aux_density_table::lap(const double x, const double y, const double z, const double* coefs) const
{
    double gx, gy, gz, lap;
    (*this)(x, y, z, coefs, gx, gy, gz, lap);
    return lap;
}
double aux_density_table::eli(const double x, const double y, const double z, const double* coefs) const
{
    double gx, gy, gz, lap;
    const double rho = (*this)(x, y, z, coefs, gx, gy, gz, lap);
    return aux_density::eli_from_density(rho, gx * gx + gy * gy + gz * gz, lap);
}
double aux_density_table::esp(const double x, const double y, const double z, const double* coefs) const
{
    return aux_density::esp_at(x, y, z, n_at, cx.data(), cy.data(), cz.data(), Z.data(), sh_start.data(), sh_l.data(), pr_start.data(), coef_off.data(), pr_exp.data(), pr_norm.data(), coefs);
}

void calc_density_ML(const aux_density_table& t, const vec& coefficients, const int np, const double* x, const double* y, const double* z, double* rho, double* gx, double* gy, double* gz, double* lap)
{
    err_checkf((int)coefficients.size() == t.n_coef, "Coefficient count does not match the auxiliary basis", std::cout);
#ifdef NOSPHERA2_USE_GPU
    if (aux_density_gpu_enabled() && aux_density_gpu_eval(t.n_at, t.cx.data(), t.cy.data(), t.cz.data(), t.r2_max.data(), t.n_sh, t.sh_start.data(), t.sh_l.data(), t.pr_start.data(), t.coef_off.data(), t.n_pr, t.pr_exp.data(), t.pr_norm.data(), t.n_coef, coefficients.data(), np, x, y, z, rho, gx, gy, gz, lap)) {
        static std::atomic<bool> announced{ false };
        if (!announced.exchange(true) && !constants::hide_gpu_notes)
            std::cout << "GPU in use: fitted density on the grid" << std::endl;
        return;
    }
#endif
    if (gx == nullptr) {
#pragma omp parallel for
        for (int p = 0; p < np; p++) rho[p] = t(x[p], y[p], z[p], coefficients.data());
        return;
    }
    if (lap == nullptr) {
#pragma omp parallel for
        for (int p = 0; p < np; p++) rho[p] = t(x[p], y[p], z[p], coefficients.data(), gx[p], gy[p], gz[p]);
        return;
    }
#pragma omp parallel for
    for (int p = 0; p < np; p++) rho[p] = t(x[p], y[p], z[p], coefficients.data(), gx[p], gy[p], gz[p], lap[p]);
}


double aux_density_table::shell_radial_moment(
    const int shell) const
{
    const int l = sh_l[shell];

    const double prefactor =
        0.5 * std::tgamma(l + 1.5);

    double integral = 0.0;

    for (int p = pr_start[shell];
        p < pr_start[shell + 1];
        ++p)
    {
        // pr_exp_l32[p] = alpha^(l + 3/2)
        integral +=
            pr_norm[p]
            * prefactor
            / pr_exp_l32[p];
    }

    return integral;
}

double aux_density_table::shell_population_integral(
    const int shell) const
{
    err_checkf(
        sh_l[shell] == 0,
        "Population integral requested for non-s auxiliary shell",
        std::cout
    );

    double integral = 0.0;

    for (int p = pr_start[shell];
        p < pr_start[shell + 1];
        ++p)
    {
        // For l=0:
        //
        // ∫ exp(-alpha r²) Y00 d³r
        //
        // = pi / (2 alpha^(3/2))
        integral +=
            pr_norm[p]
            * constants::PI
            / (2.0 * pr_exp_l32[p]);
    }

    return integral;
}

/**
 * Calculates the atomic density for a given list of atoms and coefficients.
 *
 * @param atoms The list of atoms.
 * @param coefs The coefficients used in the calculation.
 * @return The atomic density for each atom.
 */
vec calc_atomic_density(
    const std::vector<atom>& atoms,
    const vec& coefs)
{
    vec atom_elecs(atoms.size(), 0.0);
    int coef_counter = 0;
    for (int a = 0; a < atoms.size(); ++a)
    {
        int prim = 0;
        for (unsigned int shell = 0; shell < atoms[a].get_shellcount().size(); ++shell)
        {
            const int type = atoms[a].get_basis_set_entry(prim).get_type();
            const unsigned int nprim = atoms[a].get_shellcount()[shell];

            // Only s-functions have a non-zero integral
            // over all space.
            if (type != 0)
            {
                coef_counter += 2 * type + 1; prim += nprim;
                continue;
            }

            vec shell_coefs(nprim);
            vec shell_exps(nprim);

            // Collect RAW contraction coefficients and exponents.
            for (unsigned int e = 0; e < nprim; ++e, ++prim)
            {
                const basis_set_entry& bf = atoms[a].get_basis_set_entry(prim);
                shell_coefs[e] = bf.get_coefficient();
                shell_exps[e] = bf.get_exponent();
            }

            //   primitive normalization, contraction normalization
            const vec normalized_coefs = Int_Params::normalize_gto(shell_coefs, shell_exps, 0);

            double radial_integral = 0.0;
            for (unsigned int e = 0; e < nprim; ++e)
            {
                // Integral:
                //
                // ∫ exp(-alpha*r²) Y_00 d³r
                //
                // with Y_00 = 1 / sqrt(4*pi)
                //
                // = pi / (2 * alpha^(3/2))
                radial_integral +=
                    normalized_coefs[e] / (2.0 * std::pow(shell_exps[e], 1.5));
            }

            atom_elecs[a] += radial_integral * coefs[coef_counter] * constants::PI;

            ++coef_counter;
        }

        atom_elecs[a] += atoms[a].get_ECP_electrons();
    }

    return atom_elecs;
}

double apply_charge_constraint(const std::vector<atom>& atoms, vec& coefs,
                               int net_charge, bool spherical_fill_used,
                               int n_filled, double filled_eeq_charge,
                               double applied_fill_charge, std::ostream& file)
{
    // The auxiliary fit does not conserve the integral: measured over this
    // training set the REFERENCE coefficients are already -0.208 % short and
    // the ML prediction -0.235 %, so most of the deficit is inherited from the
    // density fitting rather than learned. The exact target is fixed by the
    // composition, so it can be imposed here.
    //
    // With the spherical Thakkar fill active, `atoms` is already only the
    // ML-predicted subset (the predictor erases the filled atoms from its
    // copy of the wavefunction), and Thakkar densities are exactly neutral,
    // so the sum of nuclear charges over THESE atoms is the right target.
    //
    // Only l=0 functions have a non-zero integral, so only they are scaled.
    // The scaling is GLOBAL on purpose: the atom-centred auxiliary basis is not
    // an atomic-charge partition (measured per-element ratios span 0.899 for Al
    // to 1.150 for B, with ~10 % scatter, which is bonding rather than error),
    // so a per-element correction would distort the chemistry to fix a total.
    vec atom_elecs = calc_atomic_density(atoms, coefs);

    double predicted = 0.0, target = 0.0, ecp = 0.0;
    for (int a = 0; a < static_cast<int>(atoms.size()); a++)
    {
        predicted += atom_elecs[a];
        target += static_cast<double>(atoms[a].get_charge());
        ecp += static_cast<double>(atoms[a].get_ECP_electrons());
    }
    // calc_atomic_density folds the ECP electrons in, but those do not come
    // from the coefficients and must not be rescaled.
    const double from_coefs = predicted - ecp;
    double target_from_coefs = target - ecp;

    // A charged system does not integrate to the sum of the nuclear charges.
    // Ignoring this would silently pull an anion back to neutral: -1 on a
    // 300-electron system is 0.33 %, well inside the sanity guard below, so it
    // would be applied without complaint and be wrong.
    // The spherically filled atoms are no longer assumed neutral: each carries
    // Z - q. Those q electrons have to come from somewhere, and the only place
    // they can come from is the predicted region - so the target moves by
    // exactly the charge the fill took on. Without this the two halves disagree
    // and the SYSTEM total stops being exact, which is the one guarantee this
    // whole mechanism exists to provide.
    if (std::abs(applied_fill_charge) > 1e-12)
    {
        target_from_coefs += applied_fill_charge;
        file << "Charge constraint: " << std::showpos << std::fixed << std::setprecision(3)
             << applied_fill_charge << std::noshowpos
             << " e moved to the spherically filled atom(s), so the predicted"
             << " region is targeted accordingly and the system total stays exact."
             << std::endl;
    }

    if (net_charge != 0)
    {
        if (spherical_fill_used)
        {
            // Some atoms are Thakkar-filled, so `atoms` is only the ML subset.
            // Thakkar densities are strictly neutral, but how a NET charge
            // divides between the ML region and the filled region is not
            // defined here. Refusing beats guessing.
            file << "Charge constraint SKIPPED: net charge " << net_charge
                 << " on a structure that also uses the spherical fill. The split"
                 << " of that charge between the predicted and filled regions is"
                 << " undefined, so the density is left alone." << std::endl;
            return 1.0;
        }
        target_from_coefs -= static_cast<double>(net_charge);
        file << "Charge constraint: net charge " << net_charge
             << " taken into account." << std::endl;
    }

    if (from_coefs <= 0.0 || target_from_coefs <= 0.0)
    {
        file << "Charge constraint skipped: non-positive electron count." << std::endl;
        return 1.0;
    }
    // With the spherical fill in play the target assumes every filled atom is
    // NEUTRAL. That is exactly right for a Thakkar density of a neutral atom
    // and wrong whenever the filled atom is a formal ion. Seen in practice on a
    // calcium salt: Ca filled as neutral (20 e) puts the ML target at 116, the
    // model predicted 116.48 because the organic part is really a dianion, and
    // the constraint then pulled it the WRONG way. The net charge is zero, so
    // no charge check can catch this - say so and let the user judge.
    if (spherical_fill_used)
    {
        // The TOTAL stays exact either way - the filled atoms carry a fixed
        // neutral count and the rest is imposed here - so F(000) is safe. What
        // the neutral assumption gets wrong is WHERE the electrons sit, and an
        // exact total hides that. Quantify it instead.
        if (std::isnan(filled_eeq_charge))
        {
            file << "NOTE: " << n_filled << " atom(s) are spherically filled and their"
                 << " charge could not be estimated, so they are left neutral." << std::endl;
        }
        else
        {
            const double missed = filled_eeq_charge - applied_fill_charge;
            file << "NOTE: " << n_filled << " atom(s) are spherically filled. EEQ puts "
                 << std::showpos << std::fixed << std::setprecision(2) << filled_eeq_charge
                 << std::noshowpos << " e on them, of which " << std::showpos
                 << applied_fill_charge << std::noshowpos << " e is applied." << std::endl;
            if (std::abs(missed) > 0.05)
                file << "      " << std::showpos << std::setprecision(2) << missed
                     << std::noshowpos << " e could not be applied because no ion is"
                     << " tabulated for that element; those atoms stay neutral and that"
                     << " charge remains on the predicted region." << std::endl;
        }
    }

    const double factor = target_from_coefs / from_coefs;
    if (std::abs(factor - 1.0) > 0.05)
    {
        file << "Charge constraint SKIPPED: factor " << factor
             << " is further than 5 % from unity, which means something else is"
             << " wrong - refusing to paper over it." << std::endl;
        return 1.0;
    }

    int coef_counter = 0;
    for (int a = 0; a < static_cast<int>(atoms.size()); a++)
    {
        int prim = 0;
        for (unsigned int shell = 0; shell < atoms[a].get_shellcount().size(); shell++)
        {
            const int type = atoms[a].get_basis_set_entry(prim).get_type();
            if (type != 0)
            {
                coef_counter += (2 * type + 1);
                prim += atoms[a].get_shellcount()[shell];
                continue;
            }
            prim += atoms[a].get_shellcount()[shell];
            coefs[coef_counter] *= factor;
            coef_counter++;
        }
    }

    // On the training chemistry this factor sits at 1.00235 +/- 0.0008. A much
    // larger one means the density shape is being extrapolated, and a fixed
    // integral would otherwise hide that.
    if (std::abs(factor - 1.0) > 0.005)
    {
        file << "NOTE: correction of " << std::fixed << std::setprecision(3)
             << 100.0 * (factor - 1.0) << " % is well outside the ~0.24 % seen on"
             << " training-like systems - treat this prediction with caution."
             << std::endl;
    }

    file << "Charge constraint applied: " << std::fixed << std::setprecision(4)
         << from_coefs << " -> " << target_from_coefs << " electrons over "
         << atoms.size() << (spherical_fill_used ? " ML-predicted" : "")
         << " atoms (factor " << std::setprecision(8) << factor << ")" << std::endl;
    return factor;
}

void calc_cube_ML(const vec& data, WFN& dummy, cube& cube_data, const int& atom_nr)
{
    _time_point start = get_time();

    const int s1 = cube_data.get_size(0), s2 = cube_data.get_size(1), s3 = cube_data.get_size(2), total_size = s1 * s2 * s3;
    std::cout << "Lets go into the loop! There is " << total_size << " points" << std::endl;

    ProgressBar* progress = new ProgressBar(total_size, 60, "=", " ", "Calculating Values");
    vec v1{
        cube_data.get_vector(0, 0),
        cube_data.get_vector(1, 0),
        cube_data.get_vector(2, 0) },
        v2{
            cube_data.get_vector(0, 1),
            cube_data.get_vector(1, 1),
            cube_data.get_vector(2, 1) },
            v3{
                cube_data.get_vector(0, 2),
                cube_data.get_vector(1, 2),
                cube_data.get_vector(2, 2) },
                orig{
                    cube_data.get_origin(0),
                    cube_data.get_origin(1),
                    cube_data.get_origin(2) };

    if (atom_nr != -1)
        std::cout << "Calculation for atom " << atom_nr << std::endl;

    const std::vector<atom> atoms = dummy.get_atoms();
    //atom_nr selects one atom: its own table and the slice of the coefficients that belongs to it
    const aux_density_table full(atoms);
    const aux_density_table t(atom_nr == -1 ? atoms : std::vector<atom>{ atoms[atom_nr] });
    vec coefs = data;
    if (atom_nr != -1) {
        const int off = full.coef_off[full.sh_start[atom_nr]];
        coefs.assign(data.begin() + off, data.begin() + off + t.n_coef);
    }
#pragma omp parallel for schedule(dynamic)
    for (int index = 0; index < total_size; index++)
    {
        int i = index / (s2 * s3);
        int j = (index / s3) % s2;
        int k = index % s3;

        cube_data.set_value(i, j, k, t(
            i * v1[0] + j * v2[0] + k * v3[0] + orig[0],
            i * v1[1] + j * v2[1] + k * v3[1] + orig[1],
            i * v1[2] + j * v2[2] + k * v3[2] + orig[2], coefs.data()));
        progress->update();
    }
    delete (progress);

    using namespace std;
    _time_point end = get_time();
    if (get_sec(start, end) < 60)
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
    cube_data.calc_dv();
    std::cout << "Number of electrons: " << std::fixed << std::setprecision(4) << cube_data.sum() << std::endl;
};

cube calc_cube_ML(const vec& data, WFN& dummy, const int& atom_nr)
{
    properties_options opts;
    readxyzMinMax_fromWFN(dummy, opts);
    cube CubeRho(opts.NbSteps, dummy.get_ncen(), true);
    CubeRho.give_parent_wfn(dummy);

    for (int i = 0; i < 3; i++)
    {
        CubeRho.set_origin(i, opts.MinMax[i]);
        CubeRho.set_vector(i, i, (opts.MinMax[i + 3] - opts.MinMax[i]) / opts.NbSteps[i]);
    }
    CubeRho.set_comment1("Calculated density using NoSpherA2 from ML Data");
    CubeRho.set_comment2("from " + dummy.get_path().string());
    CubeRho.set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_RI_rho.cube");

    calc_cube_ML(data, dummy, CubeRho, atom_nr);

    return CubeRho;
};

#include "integrator.h"
#include "libCintMain.h"
#include "nos_math.h"
#include "npy.h"
void create_SALTED_training_data(const WFN& orbital, const WFN& aux, const options& opts) {
    std::cout << "Calculating density fitting coefficients..." << std::endl;
    DensityFitting::CONFIG config = DensityFitting::config_from_options(opts);
    config.analyze_quality = true;
    //config.restrain_type = DensityFitting::RESTRAINT_TYPE::SIMPLE_AND_TIK;
    //config.charge_scheme = DensityFitting::CHARGE_SCHEME::HIRSHFELD;
    //if (wavy->get_origin() == e_origin::ptb)
    //    config.restraint_strength = 1.0e-4;

    vec coefs = DensityFitting::density_fit(orbital, aux, config);

    vec overlap;
    Int_Params aux_basis(aux);
    compute2C<Overlap2C_SPH>(aux_basis, overlap);
    const int nao_max = aux_basis.get_nao();

    dMatrix1 coefs_vec(coefs.size());
    coefs_vec.container() = coefs;
    dMatrix2 overlap_mat(nao_max, nao_max);
    overlap_mat.container() = overlap;

    npy::write_npy("coefficients.npy", 
        npy::npy_data<double>{
            coefs, 
            { static_cast<unsigned long>(coefs.size()) },
            false}
    );

    npy::write_npy("overlap.npy",
        npy::npy_data<double>{
        overlap,
        { static_cast<unsigned long>(nao_max), static_cast<unsigned long>(nao_max) },
            false}
    );

}
