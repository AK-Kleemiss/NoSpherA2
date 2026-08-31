#include "pch.h"
#include "SALTED_predictor.h"
#include "SALTED_utilities.h"
#include "SALTED_equicomb.h"
#include "nos_math.h"
#include "constants.h"
#include "wfn_class.h"
#include "basis_set.h"
#include <filesystem>


SALTEDPredictor::SALTEDPredictor(WFN wavy_in, options& opt_in) : wavy(wavy_in) //We copy the input WFN to modify it later
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
            return;
        }
        std::cout << "No SALTED binary file found in directory: " << opt_in.salted_model_dir << std::endl;
        exit(1);
    }

    if (opt_in.debug) std::cout << "Using SALTED Binary file: " << config.salted_filename << std::endl;
    _path = _path / config.salted_filename;
    SALTED_BINARY_FILE file = SALTED_BINARY_FILE(_path);
    file.populate_config(config);

    SALTED_Utils::filter_input(wavy, opt_in, config);
    
    //wavy.write_xyz("temp_rascaline.xyz"); //Also this
    //config.predict_filename = "temp_rascaline.xyz";

    natoms = wavy.get_ncen();
    
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
        if (label == "D" || label == "d")
        {
            label = "H";
        }
        atomic_symbols.emplace_back(label);
    }

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

    featomic::SimpleSystem featomic_system = SALTED_Utils::gen_featomic_system(wavy);
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
        // Same hyperparameters, so the second descriptor set would be an exact
        // duplicate of the first. Record that instead of copying it: equicomb
        // then reads conj(v1) where it would have read v2. Conjugation only
        // flips the sign of the imaginary part, which is exact in IEEE
        // arithmetic, so the result is bit-identical - and v2 is a copy of an
        // array that runs to gigabytes on a protein.
        v2_is_conj_of_v1 = true;
        v2.clear();
    }

    // Conjugate v2 once here rather than per use in equicomb. Skipped when v2
    // is just conj(v1): there is nothing to conjugate, equicomb applies it.
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
    // How long the model costs to get off disk decides whether a lazy loader
    // needs to be asynchronous or merely lazy: the total bytes are the same
    // either way, only the moment changes.
    const auto _t_model = std::chrono::steady_clock::now();
    const std::filesystem::path _SALTEDpath = SALTED_DIR / config.salted_filename;
    // Kept open for the whole prediction: the feature and projector matrices are
    // fetched lambda by lambda as the loop reaches them, not up front.
    model_file = std::make_unique<SALTED_BINARY_FILE>(_SALTEDpath);
    SALTED_BINARY_FILE &file = *model_file;
    if (config.field) {
        err_not_impl_f("Calculations using 'Field = True' are not yet supported", std::cout);
    }
    weights = file.read_weights();
    wigner3j = file.read_wigners();

    if (config.average) av_coefs = file.read_averages();
    if (config.sparsify) vfps = file.read_fps();


    // Only the species this structure actually contains need their model data.
    // The model knows every species it was trained on; a protein uses five of
    // them, alanine four. The rest are read for their SHAPE alone, because
    // `weights` is one flat vector laid out over every species the model knows
    // and the offset of a present species depends on the sizes of the absent ones
    // in front of it. Shapes are two integers; the payload is hundreds of MB.
    std::unordered_set<std::string> present;
    for (const std::string &spe : config.species)
        if (atom_idx.find(spe) != atom_idx.end()) present.insert(spe);

    // Nothing large is read here. Both blocks are INDEXED - offset and shape per
    // (species, lambda) - and the matrices are fetched in load_model_lambda()
    // when the prediction reaches that lambda, then dropped again. Each block is
    // used exactly once per run, so the bytes read are the same as before.
    feat_index = file.index_lambda_based_data("FEATS");
    proj_index = file.index_lambda_based_data("PROJ");
    model_species = present;

    // The two things the rest of the run needs from the shapes alone:
    //  - Mspe, the number of sparse environments of a present species;
    //  - the projector width of EVERY species the model knows, present or not,
    //    because `weights` is one flat vector laid out over all of them.
    for (const auto &[k, ref] : proj_index)
        proj_dims[k] = { ref.rows, ref.cols };
    for (const std::string &spe : present)
    {
        const auto it = feat_index.find(spe + "0");
        if (it != feat_index.end()) Mspe[spe] = static_cast<int>(it->second.rows);
    }

    // What the model actually costs, and how it splits by lambda. The prediction
    // walks one lambda at a time, so the per-lambda figure is the floor a lazy
    // loader could reach - worth knowing before writing one.
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


// Fetch the model matrices for one lambda, use them, drop them.
//
// The prediction visits each lambda exactly once, and within a lambda it needs
// power_env_sparse and Vmat for every species present. Nothing above that
// lambda is touched again afterwards - the weight accounting downstream works
// off psi_nm and the projector SHAPES, both of which outlive the matrices - so
// the blocks can be let go as soon as the kernels for that lambda are built.
//
// The bytes read over a run are exactly what they were when everything was
// loaded up front; only the moment changes. What changes with it is the peak:
// the largest single lambda instead of all of them at once.
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
    // How many lambda blocks are alive at once.
    //
    // A block is natoms * (2*lam+1) * featsize doubles; holding all of them sums
    // to (nang+1)^2 times a single block - 13.9 GB for 8,566 atoms, most of what
    // the prediction uses. Holding fewer bounds that, at the cost of revisiting
    // each species once per group instead of once in total.
    //
    // Streaming the tsc is the signal that memory is the binding constraint, so
    // -tsc_block selects the lean setting; otherwise everything is held, which is
    // the original loop order exactly.
    const int lmax_max = SALTED_Utils::get_lmax_max(lmax);
    ivec featsize(lmax_max + 1);
    std::vector<std::vector<dMatrix2>> psi_nm(config.species.size());
    for (int spe_idx = 0; spe_idx < (int)config.species.size(); spe_idx++)
        psi_nm[spe_idx].resize(lmax[config.species[spe_idx]] + 1);
    // The only quantity that crosses lambda: set at lam = 0, reused above it when
    // zeta != 1. One per species, and small.
    std::vector<dMatrix2> kernell0(config.species.size());
    for (int lam = 0; lam <= lmax_max; lam++)
    {
        vec p;
        const auto _t_eq = std::chrono::steady_clock::now();
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
                        lvalues[llmax] = { l1, l2 };
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

        cvec2 c2r = SALTED_Utils::complex_to_real_transformation({ 2 * lam + 1 })[0];

        featsize[lam] = config.nspe1 * config.nspe2 * config.nrad1 * config.nrad2 * llmax;
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
        _t_equicomb += _elapsed(_t_eq);
        // Deliberately after the descriptors, not before: the model is only needed
        // by the kernels, and the descriptor stage is where memory is already
        // highest. Holding it across equicomb would raise the peak for nothing.
        const auto _t_ld = std::chrono::steady_clock::now();
        load_model_lambda(lam);
        _t_model_io += _elapsed(_t_ld);
        const auto _t_kn = std::chrono::steady_clock::now();
        // Species-outer within the group, so a species walks several consecutive
        // lambdas and keeps its sparse matrices hot.
        for (int spe_idx = 0; spe_idx < (int)config.species.size(); spe_idx++)
        {
            const string spe = config.species[spe_idx];
            if (atom_idx.find(spe) == atom_idx.end()) continue;
            if (lam > lmax[spe]) continue;

            int lam2_1 = 2 * lam + 1;
            int row_size = featsize[lam] * lam2_1; // Size of a block of rows

            dMatrix2 pvec_lam(atom_idx[spe].size() * lam2_1, featsize[lam]);
            dMatrixRef2 _pvec(p.data(), natoms, featsize[lam] * lam2_1);
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
        _t_kernels += _elapsed(_t_kn);
        // Nothing above this lambda reads these matrices again.
        free_model_lambda(lam);
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
                // The projector data for an absent species was never loaded; its
                // shape was, and the shape is all this accounting needs.
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
