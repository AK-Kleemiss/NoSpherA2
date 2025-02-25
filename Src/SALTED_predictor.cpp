#include "pch.h"
#include "SALTED_predictor.h"
#include "SALTED_utilities.h"
#include "SALTED_equicomb.h"
#include "nos_math.h"
#include "constants.h"
#include "wfn_class.h"
#include <filesystem>

SALTEDPredictor::SALTEDPredictor(const WFN &wavy_in, options &opt_in)
{
    std::filesystem::path _path = opt_in.SALTED_DIR;
    SALTED_DIR = opt_in.SALTED_DIR;
    debug = opt_in.debug;

    config.salted_filename = find_first_salted_file(opt_in.SALTED_DIR);

	if (config.salted_filename == "") {
        std::cout << "No SALTED binary file found in directory: " << opt_in.SALTED_DIR << std::endl;
        exit(1);
    }

    if (opt_in.debug) std::cout << "Using SALTED Binary file: " << config.salted_filename << std::endl;
    _path = _path / config.salted_filename;
	SALTED_BINARY_FILE file = SALTED_BINARY_FILE(_path);
	file.populate_config(config);

    bool i_know_all = true;
#pragma omp parallel for reduction(&& : i_know_all)
    for (int a = 0; a < wavy_in.get_ncen(); a++)
        if (find(config.neighspe1.begin(), config.neighspe1.end(), std::string(constants::atnr2letter(wavy_in.get_atom_charge(a)))) == config.neighspe1.end())
            i_know_all = false;
    if (!i_know_all)
    {
        std::cout << "WARNING: Not all species in the structure are known to the model. The following species are not known: ";
        for (int a = 0; a < wavy_in.get_ncen(); a++)
        {
            if (find(config.neighspe1.begin(), config.neighspe1.end(), std::string(constants::atnr2letter(wavy_in.get_atom_charge(a)))) == config.neighspe1.end())
            {
                std::cout << constants::atnr2letter(wavy_in.get_atom_charge(a)) << " ";
            }
        }
        std::cout << "\nI will fill out these atoms using spherical Thakkar densities!\n";
        wavy = wavy_in; // make a copy of initial wavefunction, to leave the initial one untouched!
        for (int a = wavy_in.get_ncen() - 1; a >= 0; a--)
        {
            if (find(config.neighspe1.begin(), config.neighspe1.end(), std::string(constants::atnr2letter(wavy_in.get_atom_charge(a)))) == config.neighspe1.end())
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
}

const std::string SALTEDPredictor::get_dfbasis_name() const
{
    return config.dfbasis;
}

void calculateConjugate(std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> &v2)
{
#pragma omp parallel for
    for (int i = 0; i < v2.size(); ++i)
    {
        auto &vec3d = v2[i];
        for (auto &vec2d : vec3d)
        {
            for (auto &vec1d : vec2d)
            {
                std::transform(vec1d.begin(), vec1d.end(), vec1d.begin(), [](const std::complex<double> &val)
                               { return std::conj(val); });
            }
        }
    }
}

void SALTEDPredictor::setup_atomic_environment()
{
    const std::array<std::vector<primitive>, 118> *bs = wavy.get_basis_set_ptr();
    SALTED_Utils::set_lmax_nmax(lmax, nmax, *bs, config.species);

    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        atomic_symbols.push_back(wavy.get_atom_label(i));
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

    // RASCALINE (Generate descriptors)
    v1 = Rascaline_Descriptors(
             config.predict_filename,
             config.nrad1,
             config.nang1,
             config.sig1,
             config.rcut1,
             natoms,
             config.neighspe1,
             config.species)
             .calculate_expansion_coeffs();
    if ((config.nrad2 != config.nrad1) || (config.nang2 != config.nang1) || (config.sig2 != config.sig1) || (config.rcut2 != config.rcut1) || (config.neighspe2 != config.neighspe1))
    {
        v2 = Rascaline_Descriptors(
                 config.predict_filename,
                 config.nrad2,
                 config.nang2,
                 config.sig2,
                 config.rcut2,
                 natoms,
                 config.neighspe2,
                 config.species)
                 .calculate_expansion_coeffs();
    }
    else
    {
        v2 = v1;
    }

    // Calculate the conjugate of v2 and store it back in v2, to avoid recalculating it in the equicomb function
    calculateConjugate(v2);
    // END RASCALINE
}


void SALTEDPredictor::read_model_data() {
    const std::filesystem::path _SALTEDpath = SALTED_DIR / config.salted_filename;
    SALTED_BINARY_FILE file = SALTED_BINARY_FILE(_SALTEDpath);
    if (config.field) {
        err_not_impl_f("Calculations using 'Field = True' are not yet supported", std::cout);
    }
    weights = file.read_weights();
	wigner3j = file.read_wigners();

    if (config.average) av_coefs = file.read_averages();
	if (config.sparsify) vfps = file.read_fps();


	std::unordered_map<std::string, dMatrix2> features = file.read_features();
    Vmat = file.read_projectors();
    std::string key;
	for (std::string spe : config.species) {
		for (int lam = 0; lam < lmax[spe] + 1; lam++) {
			key = spe + std::to_string(lam);
			if (lam == 0) Mspe[spe] = features[key].extent(0);

            if (config.zeta == 1.0) {
				power_env_sparse[key] = dot(Vmat[key], features[key], true, false); //Transpose the first matrix
			}
            else {
                power_env_sparse[key] = features[key];
            }
		}
	}
}


vec SALTEDPredictor::predict()
{
    using namespace std;
    // Compute equivariant descriptors for each lambda value entering the SPH expansion of the electron density
    vec2 pvec(SALTED_Utils::get_lmax_max(lmax) + 1);
    for (int lam = 0; lam < SALTED_Utils::get_lmax_max(lmax) + 1; lam++)
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
        vec p;
        ivec2 llvec_t = transpose<int>(llvec);
        if (config.sparsify)
        {
            int nfps = static_cast<int>(vfps[lam].size());
            equicomb(natoms, (config.nspe1 * config.nrad1), (config.nspe2 * config.nrad2), v1, v2, wigner3j[lam], llvec_t, lam, c2r, featsize[lam], nfps, vfps[lam], p);
            featsize[lam] = nfps;
        }
        else
        {
            equicomb(natoms, (config.nspe1 * config.nrad1), (config.nspe2 * config.nrad2), v1, v2, wigner3j[lam], llmax, llvec_t, lam, c2r, featsize[lam], p);
        }
        pvec[lam] = p;
    }

    std::vector<std::vector<dMatrix2>> psi_nm(config.species.size());
    for (int spe_idx = 0; spe_idx < config.species.size(); spe_idx++)
    {
        const string spe = config.species[spe_idx];
        psi_nm[spe_idx].resize(lmax[spe] + 1);

        if (atom_idx.find(spe) == atom_idx.end())
        {
            continue;
        }
        dMatrix2 kernell0_nm;
        for (int lam = 0; lam < lmax[spe] + 1; ++lam)
        {
            int lam2_1 = 2 * lam + 1;
            int row_size = featsize[lam] * lam2_1; // Size of a block of rows

            dMatrix2 pvec_lam(atom_idx[spe].size() * lam2_1, featsize[lam]);
            for (int idx = 0; idx < atom_idx[spe].size(); idx++)
            {
                int start_idx = atom_idx[spe][idx] * row_size;     // Start index in pvec_flat
                int end_idx = start_idx + row_size; // End index in pvec_flat
                // Copy the block directly into flatVec2
                std::copy(pvec[lam].begin() + start_idx, pvec[lam].begin() + end_idx, pvec_lam.data() + idx*row_size);
            }
            dMatrix2 kernel_nm = dot(pvec_lam, power_env_sparse[spe + to_string(lam)], false, true); // I AM NOT SURE THIS WILL USE THE RIGHT SIZES....

            if (config.zeta == 1)
            {
                psi_nm[spe_idx][lam] = kernel_nm;
            }
            else {

                if (lam == 0)
                {
                    kernell0_nm = kernel_nm;
                    kernel_nm = elementWiseExponentiation(kernel_nm, config.zeta);
                }
                else
                {
                    for (size_t i1 = 0; i1 < natom_dict[spe]; ++i1)
                    {
                        for (size_t i2 = 0; i2 < Mspe[spe]; ++i2)
                        {
                            for (size_t i = 0; i < lam2_1; ++i)
                            {
                                for (size_t j = 0; j < lam2_1; ++j)
                                {
                                    kernel_nm( i1* lam2_1 + i,i2* lam2_1 + j ) *= pow(kernell0_nm( i1,i2 ), config.zeta - 1);
                                }
                            }
                        }
                    }
                }
                psi_nm[spe_idx][lam] = dot(kernel_nm, Vmat[spe + to_string(lam)], false, false);
            }
        }
    }
    pvec.clear();
    pvec.shrink_to_fit();

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
                // Check if Vmat[spe + to_string(l)][0] exists
                if (Vmat[spe + to_string(l)].size() == 0)
                {
                    cout << "The projector for species " << spe << " and l = " << l << " does not exist. This is a problem with the model, not NoSpherA2." << endl;
                    cout << "Continuing with the next species..., make sure there is no: " << spe << " in the structure you are trying to predict!!!!" << endl;
                    l = lmax[spe] + 1;
                    continue;
                }

                // for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
                //{
                //     isize += static_cast<int>(Vmat[spe + to_string(l)][0].size());
                // }
                isize += static_cast<int>(Vmat[spe + to_string(l)].extent(1)) * nmax[spe + to_string(l)];
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

    // cout << "          ... done!\nNumber of predicted coefficients: " << pred_coefs.size() << endl;
    // npy::npy_data<double> coeffs;
    // coeffs.data = pred_coefs;
    // coeffs.fortran_order = false;
    // coeffs.shape = { unsigned long(pred_coefs.size()) };
    // npy::write_npy("folder_model.npy", coeffs);
    return pred_coefs;
}

vec SALTEDPredictor::gen_SALTED_densities()
{
    using namespace std;
    if (coef_file != "")
    {
        vec coefs{};
        cout << "Reading coefficients from file: " << coef_file << endl;
        read_npy<double>(coef_file, coefs);
        return coefs;
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
    ;
    vfps.clear();
    ;
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