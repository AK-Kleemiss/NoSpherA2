#include "pch.h"
#include "SALTED_predictor.h"
#include "SALTED_utilities.h"
#include "SALTED_equicomb.h"
#include "nos_math.h"
#include "constants.h"
#include "wfn_class.h"
#include "basis_set.h"
#include <filesystem>


SALTEDPredictor::SALTEDPredictor(const WFN &wavy_in, options &opt_in)
{
    std::filesystem::path _path = opt_in.salted_model_dir;
    SALTED_DIR = opt_in.salted_model_dir;
    debug = opt_in.debug;

    config.salted_filename = find_first_salted_file(opt_in.salted_model_dir);

    if (config.salted_filename == "") {
        if (opt_in.coef_file != "") {
            std::cout << "Using density coefficients found in: " << opt_in.coef_file << std::endl;
            config.dfbasis = "cc-pvqz-jkfit";
            config.salted_filename = "coefficient file";
            coef_file = opt_in.coef_file;
            wavy = wavy_in;
            return;
        }
        std::cout << "No SALTED binary file found in directory: " << opt_in.salted_model_dir << std::endl;
        exit(1);
    }

    if (opt_in.debug) std::cout << "Using SALTED Binary file: " << config.salted_filename << std::endl;
    _path = _path / config.salted_filename;
    SALTED_BINARY_FILE file = SALTED_BINARY_FILE(_path);
    file.populate_config(config);

    // Lambda blocks held at once; 0 would hold every block
    lam_group_limit = 1;

    // Two kinds of atom go to the spherical Thakkar fill instead: a species the model
    // was never trained on, and an atom with nothing inside the descriptor cutoff.
    // The latter has a descriptor whose equivariant lam >= 1 components are all exactly
    // zero (only l = 0 survives), so the 1/sqrt(inner) normalisation is +inf and the
    // prediction is NaN. Dropping them changes no other atom's environment, the
    // neighbour test being symmetric.
    // The test here is purely geometric while featomic also ignores neighbours outside
    // neighspe, so an atom with only disallowed neighbours still reaches equicomb,
    // where the zero guard leaves it spherical rather than NaN.
    const int ncen_in = wavy_in.get_ncen();
    std::vector<char> use_thakkar(ncen_in, 0);
    for (int a = 0; a < ncen_in; a++)
        if (find(config.species.begin(), config.species.end(), std::string(constants::atnr2letter(wavy_in.get_atom_charge(a)))) == config.species.end())
            use_thakkar[a] = 1;
    const int n_unknown = static_cast<int>(std::count(use_thakkar.begin(), use_thakkar.end(), (char)1));

    const double rcut = std::min(config.rcut1, config.rcut2);
    // rcut is in Angstrom, the coordinates may not be
    const double rcut_internal = wavy_in.get_isBohr() ? constants::ang2bohr(rcut) : rcut;
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
                const double dx = wavy_in.get_atom_coordinate(a, ax) - wavy_in.get_atom_coordinate(b, ax);
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
                if (find(config.species.begin(), config.species.end(), std::string(constants::atnr2letter(wavy_in.get_atom_charge(a)))) == config.species.end())
                {
                    std::cout << constants::atnr2letter(wavy_in.get_atom_charge(a)) << " ";
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
        wavy = wavy_in; // make a copy of initial wavefunction, to leave the initial one untouched!
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
        std::filesystem::path new_fn = wavy.get_path().parent_path() / "SALTED_temp.xyz";
        wavy.write_xyz(new_fn);
        wavy.set_path(new_fn);
        opt_in.needs_Thakkar_fill = true;
    }
    else
    {
        wavy = wavy_in;
    }
    wavy.write_xyz("temp_rascaline.xyz");
    natoms = wavy.get_ncen();
    config.predict_filename = "temp_rascaline.xyz";
    if (wavy.get_nmo() != 0)
        wavy.clear_MOs(); // Delete unneccesarry MOs, since we are predicting anyway.

    wavy.delete_basis_set();
    if (file.basis_set_defined()) {
        load_basis_into_WFN(wavy, file.read_basis_set());
        bbasis_set_loaded = true;
    }
}

const std::string SALTEDPredictor::get_dfbasis_name() const
{
    return config.dfbasis;
}

void calculateConjugate(SALTEDDescriptors& v2)
{
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(v2.values().size()); ++i)
    {
        v2.values()[i] = std::conj(v2.values()[i]);
    }
}

void SALTEDPredictor::setup_atomic_environment()
{
    const std::shared_ptr<std::array<std::vector<primitive>, 118>> bs = wavy.get_basis_set_ptr();
    SALTED_Utils::set_lmax_nmax(lmax, nmax, *bs, config.species);

    atomic_symbols.reserve(wavy.get_ncen());
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        std::string label = wavy.get_atom_label(i);
        // Deuterium is hydrogen for the electron density; without this a joint
        // X-ray/neutron structure is refused with "Excluded species: D"
        if (label == "D" || label == "d")
        {
            label = "H";
        }
        atomic_symbols.emplace_back(label);
    }
    // # Define system excluding atoms that belong to species not listed in SALTED input
    atomic_symbols = SALTED_Utils::filter_species(atomic_symbols, config.species);

    // Print all Atomic symbols
    if (debug)
    {
        std::cout << "Atomic symbols: ";
        for (const auto &symbol : atomic_symbols)
        {
            std::cout << symbol << " ";
        }
        std::cout << std::endl;
    }

    natoms = static_cast<int>(atomic_symbols.size());
    for (int i = 0; i < atomic_symbols.size(); i++)
    {
        atom_idx[atomic_symbols[i]].push_back(i);
        natom_dict[atomic_symbols[i]] += 1;
    }


    SALTED_Utils::FeatomicHyperParameters hp{
        .cutoff_radius = config.rcut1,
        .max_radial = config.nrad1,
        .max_angular = config.nang1,
        .atomic_gaussian_width = config.sig1,
        .center_atom_weight = 1.0,
        .species = config.species,
        .neighspe = config.neighspe1,
        .radial_basis = {
            .type = "Gto",
            .spline_accuracy = 1e-6
        },
        .cutoff_function = {
            .type = "ShiftedCosine",
            .width = 0.1
        }
    };

    featomic::SimpleSystem featomic_system = SALTED_Utils::gen_featomic_system(config.predict_filename);
    // RASCALINE (Generate descriptors)
    const auto _t_desc = std::chrono::steady_clock::now();
    v1 = SALTED_Utils::calculate_SALTED_descriptors(featomic_system, hp);

    if ((config.nrad2 != config.nrad1) || (config.nang2 != config.nang1) || (config.sig2 != config.sig1) || (config.rcut2 != config.rcut1) || (config.neighspe2 != config.neighspe1))
    {
        hp.max_radial = config.nrad2;
        hp.max_angular = config.nang2;
        hp.atomic_gaussian_width = config.sig2;
        hp.cutoff_radius = config.rcut2;
        hp.neighspe = config.neighspe2;
        v2 = SALTED_Utils::calculate_SALTED_descriptors(featomic_system, hp);
    }
    else
    {
        // Same hyperparameters: v2 would duplicate v1 exactly, so record that instead.
        // equicomb then reads conj(v1); conjugation only flips the sign of the imaginary
        // part, exact in IEEE, so the result is bit-identical
        v2_is_conj_of_v1 = true;
        v2.clear();
    }

    // Conjugate v2 once here rather than per use in equicomb; skipped when v2 is conj(v1)
    if (ProgressBar::report_counts)
        std::cout << "[stages] featomic descriptors "
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - _t_desc).count()
                  << " s" << std::endl;
    if (!v2_is_conj_of_v1)
        calculateConjugate(v2);
    std::cout << "Descriptor sets " << (v2_is_conj_of_v1 ? "identical: sharing one copy"
                                                        : "differ: two copies held") << std::endl;
    // END RASCALINE
}


void SALTEDPredictor::read_model_data() {
    const auto _t_model = std::chrono::steady_clock::now();
    const std::filesystem::path _SALTEDpath = SALTED_DIR / config.salted_filename;
    // Kept open for the whole prediction; the matrices are fetched lambda by lambda
    model_file = std::make_unique<SALTED_BINARY_FILE>(_SALTEDpath);
    SALTED_BINARY_FILE &file = *model_file;
    if (config.field) {
        err_not_impl_f("Calculations using 'Field = True' are not yet supported", std::cout);
    }
    weights = file.read_weights();
    wigner3j = file.read_wigners();

    if (config.average) av_coefs = file.read_averages();
    if (config.sparsify) vfps = file.read_fps();


    // Only present species need their model data; absent ones are needed for their
    // shape alone, because `weights` is one flat vector laid out over every species
    // the model knows and the absent widths in front shift a present species' offset
    std::unordered_set<std::string> present;
    for (const std::string &spe : config.species)
        if (atom_idx.find(spe) != atom_idx.end()) present.insert(spe);

    // Indexing only, no payload: offset and shape per (species, lambda). The matrices
    // are read in load_model_lambda() and dropped again, each block used once per run
    feat_index = file.index_lambda_based_data("FEATS");
    proj_index = file.index_lambda_based_data("PROJ");
    model_species = present;

    // From the shapes alone: Mspe, the number of sparse environments of a present
    // species, and the projector width of every species the model knows
    for (const auto &[k, ref] : proj_index)
        proj_dims[k] = { ref.rows, ref.cols };
    for (const std::string &spe : present)
    {
        const auto it = feat_index.find(spe + "0");
        if (it != feat_index.end()) Mspe[spe] = static_cast<int>(it->second.rows);
    }

    if (ProgressBar::report_counts)
    {
        auto mb = [](const size_t doubles) { return doubles * sizeof(double) / 1048576.0; };
        size_t pes = 0, vm = 0, wg = 0;
        std::map<int, size_t> per_lam;
        for (const std::string &spe : present)
            for (int lam = 0; lam < lmax[spe] + 1; lam++)
            {
                const std::string k = spe + std::to_string(lam);
                const auto f = feat_index.find(k), pr = proj_index.find(k);
                const size_t p = (f == feat_index.end()) ? 0 : f->second.rows * f->second.cols;
                const size_t v = (pr == proj_index.end()) ? 0 : pr->second.rows * pr->second.cols;
                pes += p; vm += v; per_lam[lam] += p + v;
            }
        for (const auto &[lam, w] : wigner3j) wg += w.size();
        size_t worst = 0;
        for (const auto &[lam, sz] : per_lam) worst = std::max(worst, sz);
        std::cout << "[model] features " << mb(pes) << " MB + projectors " << mb(vm)
                  << " MB + wigner " << mb(wg) << " MB + weights " << mb(weights.size())
                  << " MB = " << mb(pes + vm + wg + weights.size()) << " MB if held whole;"
                  << " lazily, the largest lambda is " << mb(worst) << " MB" << std::endl;
        std::cout << "[model] indexed in "
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - _t_model).count()
                  << " s" << std::endl;
        std::cout << "[model] per lambda:";
        for (const auto &[lam, sz] : per_lam) std::cout << " l" << lam << "=" << mb(sz);
        std::cout << " MB" << std::endl;
    }
}


// Fetch the model matrices for one lambda, use them, drop them. Each lambda is
// visited once and nothing above it reads them again: the weight accounting
// downstream works off psi_nm and the projector shapes, which outlive the matrices
void SALTEDPredictor::load_model_lambda(const int lam)
{
    if (!model_file) return;
    for (const std::string &spe : model_species)
    {
        if (lam > lmax[spe]) continue;
        const std::string key = spe + std::to_string(lam);
        if (power_env_sparse.find(key) != power_env_sparse.end()) continue;
        const auto pr = proj_index.find(key);
        const auto ft = feat_index.find(key);
        if (pr == proj_index.end() || ft == feat_index.end()) continue;
        Vmat[key] = model_file->load_block(pr->second);
        dMatrix2 feats = model_file->load_block(ft->second);
        if (config.zeta == 1.0)
            power_env_sparse[key] = dot(Vmat[key], feats, true, false);
        else
            power_env_sparse[key] = std::move(feats);
    }
}

void SALTEDPredictor::free_model_lambda(const int lam)
{
    if (!model_file) return;
    for (const std::string &spe : model_species)
    {
        const std::string key = spe + std::to_string(lam);
        power_env_sparse.erase(key);
        Vmat.erase(key);
    }
}

vec SALTEDPredictor::predict()
{
    using namespace std;
    const auto _t_predict_start = std::chrono::steady_clock::now();
    auto _elapsed = [](const std::chrono::steady_clock::time_point &from)
    { return std::chrono::duration<double>(std::chrono::steady_clock::now() - from).count(); };
    double _t_equicomb = 0.0, _t_kernels = 0.0, _t_model_io = 0.0;
    // Compute equivariant descriptors for each lambda value entering the SPH expansion of the electron density
    // How many lambda blocks are alive at once. A block is natoms * (2*lam+1) *
    // featsize doubles and holding all of them sums to (nang+1)^2 times a single
    // block; fewer bounds that, at the cost of revisiting each species per group
    const int lmax_max = SALTED_Utils::get_lmax_max(lmax);
    const int lam_group = (lam_group_limit > 0) ? lam_group_limit : (lmax_max + 1);
    ivec featsize(lmax_max + 1);
    std::vector<std::vector<dMatrix2>> psi_nm(config.species.size());
    for (int spe_idx = 0; spe_idx < (int)config.species.size(); spe_idx++)
        psi_nm[spe_idx].resize(lmax[config.species[spe_idx]] + 1);
    // The only quantity that crosses lambda: set at lam = 0, reused above it when zeta != 1
    std::vector<dMatrix2> kernell0(config.species.size());
    for (int lam0 = 0; lam0 <= lmax_max; lam0 += lam_group)
    {
        const int lam1 = std::min(lam0 + lam_group, lmax_max + 1);
        vec2 pg(lam1 - lam0);
        const auto _t_eq = std::chrono::steady_clock::now();
        for (int lam = lam0; lam < lam1; lam++)
        {
        int llmax = 0;
        unordered_map<int, ivec> lvalues{};
        for (int l1 = 0; l1 < config.nang1 + 1; l1++)
        {
            for (int l2 = 0; l2 < config.nang2 + 1; l2++)
            {
                // keep only even combination to enforce inversion symmetry
                if ((lam + l1 + l2) % 2 == 0)
                {
                    if (abs(l2 - lam) <= l1 && l1 <= (l2 + lam))
                    {
                        lvalues[llmax] = {l1, l2};
                        llmax += 1;
                    }
                }
            }
        }
        // Fill dense array from dictionary
        ivec2 llvec(llmax, ivec(2));
        for (int i = 0; i < llmax; i++)
        {
            llvec[i] = lvalues[i];
        }

        cvec2 c2r = SALTED_Utils::complex_to_real_transformation({2 * lam + 1})[0];

        featsize[lam] = config.nspe1 * config.nspe2 * config.nrad1 * config.nrad2 * llmax;
        vec &p = pg[lam - lam0];
        ivec2 llvec_t = transpose<int>(llvec);
        if (config.sparsify)
        {
            int nfps = static_cast<int>(vfps[lam].size());
            p.assign((size_t)natoms * ((size_t)2 * lam + 1) * nfps, 0.0);
            equicomb(natoms, (config.nspe1 * config.nrad1), (config.nspe2 * config.nrad2), v1, v2, wigner3j[lam], llvec_t, lam, c2r, featsize[lam], nfps, vfps[lam], p, v2_is_conj_of_v1);
            featsize[lam] = nfps;
        }
        else
        {
            p.assign((size_t)natoms * ((size_t)2 * lam + 1) * featsize[lam], 0.0);
            equicomb(natoms, (config.nspe1 * config.nrad1), (config.nspe2 * config.nrad2), v1, v2, wigner3j[lam], llmax, llvec_t, lam, c2r, featsize[lam], p, v2_is_conj_of_v1);
        }
        }
        _t_equicomb += _elapsed(_t_eq);
        // After the descriptors, not before: holding the model across equicomb, where
        // memory is already highest, would raise the peak for nothing
        const auto _t_ld = std::chrono::steady_clock::now();
        for (int lam = lam0; lam < lam1; lam++) load_model_lambda(lam);
        _t_model_io += _elapsed(_t_ld);
        const auto _t_kn = std::chrono::steady_clock::now();
        // Species-outer within the group, so a species keeps its sparse matrices hot
        for (int spe_idx = 0; spe_idx < (int)config.species.size(); spe_idx++)
        {
            const string spe = config.species[spe_idx];
            if (atom_idx.find(spe) == atom_idx.end()) continue;
            for (int lam = lam0; lam < lam1; lam++)
            {
                if (lam > lmax[spe]) continue;
        {
            int lam2_1 = 2 * lam + 1;
            int row_size = featsize[lam] * lam2_1; // Size of a block of rows

            dMatrix2 pvec_lam(atom_idx[spe].size() * lam2_1, featsize[lam]);
            dMatrixRef2 _pvec(pg[lam - lam0].data(), natoms, featsize[lam] * lam2_1);
            double* pvec_ptr = pvec_lam.data();
            for (const int idx : atom_idx[spe])
            {
                auto _temp = Kokkos::submdspan(_pvec, idx, Kokkos::full_extent);
                std::copy(_temp.data_handle(), _temp.data_handle() + row_size, pvec_ptr);
                pvec_ptr += row_size;
            }
            dMatrix2 kernel_nm = dot(pvec_lam, power_env_sparse[spe + to_string(lam)], false, true);

            if (config.zeta == 1)
            {
                psi_nm[spe_idx][lam] = kernel_nm;
            }
            else {

                if (lam == 0)
                {
                    kernell0[spe_idx] = kernel_nm;
                    kernel_nm = elementWiseExponentiation(kernel_nm, config.zeta);
                }
                else
                {
                    for (size_t i1 = 0; i1 < natom_dict[spe]; ++i1)
                    {
                        for (size_t i2 = 0; i2 < Mspe[spe]; ++i2)
                        {
                            double scale_factor = pow(kernell0[spe_idx](i1, i2), config.zeta - 1);
                            size_t base_i = i1 * lam2_1;
                            size_t base_j = i2 * lam2_1;
                            for (size_t i = 0; i < lam2_1; ++i)
                            {
                                for (size_t j = 0; j < lam2_1; ++j)
                                {
                                    kernel_nm(base_i + i, base_j + j) *= scale_factor;
                                }
                            }
                        }
                    }
                }
                psi_nm[spe_idx][lam] = dot(kernel_nm, Vmat[spe + to_string(lam)], false, false);
            }
        }
            }
        }
        _t_kernels += _elapsed(_t_kn);
        for (int lam = lam0; lam < lam1; lam++) free_model_lambda(lam);
    }

    unordered_map<string, dMatrix1> C{};
    unordered_map<string, int> ispe{};
    int isize = 0;
    for (int spe_idx = 0; spe_idx < config.species.size(); spe_idx++)
    {
        const string spe = config.species[spe_idx];
        if (atom_idx.find(spe) == atom_idx.end())
        {
            for (int l = 0; l < lmax[spe] + 1; ++l)
            {
                // Never loaded for an absent species; its shape was, and that is all this needs
                const auto dim_it = proj_dims.find(spe + to_string(l));
                if (dim_it == proj_dims.end() || dim_it->second[1] == 0)
                {
                   std::cout << "The projector for species " << spe << " and l = " << l << " does not exist. This is a problem with the model, not NoSpherA2." << endl;
                   std::cout << "Continuing with the next species..., make sure there is no: " << spe << " in the structure you are trying to predict!!!!" << endl;
                    break;
                }

                // for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
                //{
                //     isize += static_cast<int>(Vmat[spe + to_string(l)][0].size());
                // }
                isize += static_cast<int>(dim_it->second[1]) * nmax[spe + to_string(l)];
            }
            continue;
        }
        ispe[spe] = 0;
        for (int l = 0; l < lmax[spe] + 1; ++l)
        {
            for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
            {
                // int Mcut = static_cast<int>(psi_nm[spe + to_string(l)][0].size());
                int Mcut = static_cast<int>(psi_nm[spe_idx][l].extent(1));
                // Check if isize + Mcut > weights.size()
                err_chekf(isize + Mcut <= weights.size(), "isize + Mcut > weights.size()", std::cout);

                dMatrix1 weights_subset(Mcut);
                std::copy(weights.data() + isize, weights.data() + isize + Mcut, weights_subset.data());

                C[spe + to_string(l) + to_string(n)] = dot(psi_nm[spe_idx][l], weights_subset, false);

                isize += Mcut;
            }
        }
    }
    psi_nm.clear();
    psi_nm.shrink_to_fit();


    int Tsize = 0;
    for (int iat = 0; iat < natoms; iat++)
    {
        string spe = atomic_symbols[iat];
        for (int l = 0; l < lmax[spe] + 1; l++)
        {
            for (int n = 0; n < nmax[spe + to_string(l)]; n++)
            {
                Tsize += 2 * l + 1;
            }
        }
    }

    vec Av_coeffs(Tsize, 0.0);

    // fill vector of predictions
    int i = 0;
    vec pred_coefs(Tsize, 0.0);
    for (int iat = 0; iat < natoms; ++iat)
    {
        string spe = atomic_symbols[iat];
        for (int l = 0; l < lmax[spe] + 1; ++l)
        {
            for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
            {
                // for (int ind = 0; ind < 2 * l + 1; ++ind)
                //{
                //     pred_coefs[i + ind] = C[spe + to_string(l) + to_string(n)][ispe[spe] * (2 * l + 1) + ind];
                // }
                std::copy_n(C[spe + to_string(l) + to_string(n)].data() + ispe[spe] * (2 * l + 1), 2 * l + 1, pred_coefs.begin() + i);

                if (config.average && l == 0)
                {
                    Av_coeffs[i] = av_coefs[spe][n];
                }
                i += 2 * l + 1;
            }
        }
        ispe[spe] += 1;
    }

    if (config.average)
    {
        for (i = 0; i < Tsize; i++)
        {
            pred_coefs[i] += Av_coeffs[i];
        }
    }

    //std::cout << "          ... done!\nNumber of predicted coefficients: " << pred_coefs.size() << endl;
    // npy::npy_data<double> coeffs;
    // coeffs.data = pred_coefs;
    // coeffs.fortran_order = false;
    // coeffs.shape = { unsigned long(pred_coefs.size()) };
    // npy::write_npy("folder_model.npy", coeffs);
    if (ProgressBar::report_counts)
    {
        const double total = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - _t_predict_start).count();
        std::cout << "[stages] predict " << total << " s = equicomb " << _t_equicomb
                  << " s (" << (100.0 * _t_equicomb / total) << "%) + kernels "
                  << _t_kernels << " s (" << (100.0 * _t_kernels / total)
                  << "%) + model I/O " << _t_model_io << " s ("
                  << (100.0 * _t_model_io / total) << "%) + rest "
                  << (total - _t_equicomb - _t_kernels - _t_model_io) << " s" << std::endl;
    }
    return pred_coefs;
}

vec SALTEDPredictor::gen_SALTED_densities()
{
    using namespace std;
    if (coef_file != "")
    {
        std::vector<float> coefs{};
        std::cout << "Reading coefficients from file: " << coef_file << endl;
        read_npy<float>(coef_file, coefs);
        vec double_coefs(coefs.size());
        for (int i = 0; i < coefs.size(); i++)
        {
            double_coefs[i] = static_cast<double>(coefs[i]);
        }
        return double_coefs;
    }

    // Run generation of tsc file
    _time_point start;
    if (debug)
        start = get_time();

    setup_atomic_environment();


    read_model_data();


    vec coefs = predict();

    // File VERSION 3 models carry an optional NORMC block asking for the
    // electron count to be constrained. Applied here rather than at each call
    // site so the tsc, the charge table and the cubes all see the same density.
    // V2 models have no such block, so they are untouched.
    if (model_file && model_file->charge_constraint_defined())
    {
        const auto entries = model_file->read_charge_constraint();
        const auto mode_it = entries.find("MODE");
        const int mode = (mode_it != entries.end() && !mode_it->second.empty())
                             ? static_cast<int>(std::lround(mode_it->second[0]))
                             : 0;
        if (mode == 1)
            apply_charge_constraint(wavy.get_atoms(), coefs, std::cout);
        else if (mode != 0)
            std::cout << "Unknown charge-constraint mode " << mode
                      << " in the model file; leaving the density alone." << std::endl;
    }

    shrink_intermediate_vectors();
    return coefs;
}

void SALTEDPredictor::shrink_intermediate_vectors()
{
    v1.clear();
    v2.clear();
    weights.clear();
    Vmat.clear();
    natom_dict.clear();
    lmax.clear();
    nmax.clear();
    Mspe.clear();
    vfps.clear();
    wigner3j.clear();
    av_coefs.clear();
    power_env_sparse.clear();
    featsize.clear();
    v1.shrink_to_fit();
    v2.shrink_to_fit();
    weights.shrink_to_fit();
    std::unordered_map<std::string, dMatrix2> umap;
    std::unordered_map<std::string, int> umap2;
    std::unordered_map<std::string, dMatrix1> umap3;
    Vmat.swap(umap);
    natom_dict.swap(umap2);
    lmax.swap(umap2);
    nmax.swap(umap2);
    power_env_sparse.swap(umap);
};
