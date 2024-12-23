#include "SALTED_predictor.h"
#include "constants.h"
#include "pch.h"
#include <filesystem>

#ifdef _WIN32
#include "DLL_Helper.h"
#endif


//-SALTED D:\Models\Iron_Complex -cif mohrs_salt_IAM.cif -wfn mohrs_salt_IAM.xyz  -cpus 8 -hkl_min_max -14 14 -12 27 -19 20
// std::string find_first_h5_file(const std::string& directory_path)
SALTEDPredictor::SALTEDPredictor(const WFN &wavy_in, options &opt_in) : _opt(opt_in)
{
    std::filesystem::path _path = _opt.SALTED_DIR;

    config.h5_filename = find_first_h5_file(_opt.SALTED_DIR);

    if (config.h5_filename == "")
    {
        if (_opt.debug)
        {
            std::cout << "No HDF5 file found in the SALTED directory. Using inputs.txt instead." << std::endl;
        }
        _path = _path / "inputs.txt";
        if (_opt.debug)
        {
            std::cout << "Using inputs file: " << _path << std::endl;
        }
        config.populateFromFile(_path);
    }
    else
    {
        if (_opt.debug)
        {
            std::cout << "Using HDF5 file: " << config.h5_filename << std::endl;
        }

        // If RAS is enabled (i.e. hdf5 is enabled) read the contents of the hdf5 file
#if has_RAS
        _path = _path / config.h5_filename;
        H5::H5File config_file(_path, H5F_ACC_RDONLY);
        config.populateFromFile(config_file);
#else
        err_not_impl_f("HDF5 files are not supported by this build", std::cout);
#endif
    }
    bool i_know_all = true;
#pragma omp parallel for reduction(&& : i_know_all)
    for (int a = 0; a < wavy_in.get_ncen(); a++)
        if (find(config.neighspe1.begin(), config.neighspe1.end(), std::string(constants::atnr2letter(wavy_in.atoms[a].charge))) == config.neighspe1.end())
            i_know_all = false;
    if (!i_know_all)
    {
        std::cout << "WARNING: Not all species in the structure are known to the model. The following species are not known: ";
        for (int a = 0; a < wavy_in.get_ncen(); a++)
        {
            if (find(config.neighspe1.begin(), config.neighspe1.end(), std::string(constants::atnr2letter(wavy_in.atoms[a].charge))) == config.neighspe1.end())
            {
                std::cout << constants::atnr2letter(wavy_in.atoms[a].charge) << " ";
            }
        }
        std::cout << "\nI will fill out these atoms using spherical Thakkar densities!\n";
        wavy = wavy_in; // make a copy of initial wavefunction, to leave the initial one untouched!
        for (int a = wavy_in.get_ncen() - 1; a >= 0; a--)
        {
            if (find(config.neighspe1.begin(), config.neighspe1.end(), std::string(constants::atnr2letter(wavy_in.atoms[a].charge))) == config.neighspe1.end())
            {
                wavy.erase_atom(a);
            }
        }
        // remove all known basis sets, to not get problems with newly loaded ones
        for (int a = 0; a < wavy.get_ncen(); a++)
        {
            wavy.atoms[a].basis_set.clear();
        }
        std::filesystem::path new_fn = wavy.get_path().parent_path() / "SALTED_temp.xyz";
        wavy.write_xyz(new_fn);
        wavy.set_path(new_fn);
        _opt.needs_Thakkar_fill = true;
    }
    else
    {
        wavy = wavy_in;
    }
    wavy.write_xyz("temp_rascaline.xyz");
	config.predict_filename = "temp_rascaline.xyz";
    if (wavy.get_nmo() != 0)
        wavy.clear_MOs(); //Delete unneccesarry MOs, since we are predicting anyway.
}

SALTEDPredictor::SALTEDPredictor() : _opt(*(new options())) {}


SALTEDPredictor::~SALTEDPredictor()
{
    unload_BLAS();
}

void SALTEDPredictor::load_BLAS()
{
#if has_RAS
#ifdef _WIN32
    _putenv_s("OPENBLAS_NUM_THREADS", std::to_string(_opt.threads).c_str());
    typedef void (*ExampleFunctionType)(void);
    BLAS_pointer = math_load_BLAS(_opt.threads);
    has_BLAS = BLAS_pointer != NULL;
#else
    std::string nums = "OPENBLAS_NUM_THREADS=" + std::to_string(_opt.threads);
    char* env = strdup(nums.c_str());
    putenv(env);
    has_BLAS = true;
#endif
#endif
}

void SALTEDPredictor::unload_BLAS()
{
#ifdef _WIN32
    math_unload_BLAS(BLAS_pointer);
#endif
}

const std::string SALTEDPredictor::get_dfbasis_name()
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

    for (int i = 0; i < wavy.atoms.size(); i++)
    {
        atomic_symbols.push_back(wavy.atoms[i].label);
    }
    // # Define system excluding atoms that belong to species not listed in SALTED input
    atomic_symbols = SALTED_Utils::filter_species(atomic_symbols, config.species);

    // Print all Atomic symbols
    if (_opt.debug)
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
#if has_RAS
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
#else
    err_not_impl_f("RASCALINE is not supported by this build", std::cout);
#endif
    // END RASCALINE
}

void SALTEDPredictor::read_model_data()
{
    using namespace std;
    // Define zeta as a string with one decimal
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(1) << config.zeta;
    std::string zeta_str = stream.str();
#if has_RAS == 1
    H5::H5File features(_opt.SALTED_DIR / "GPR_data" / ("FEAT_M - " + std::to_string(config.Menv) + ".h5"), H5F_ACC_RDONLY);
    H5::H5File projectors(_opt.SALTED_DIR / "GPR_data" / ("projector_M" + std::to_string(config.Menv) + "_zeta" + zeta_str + ".h5"), H5F_ACC_RDONLY);
    std::vector<hsize_t> dims_out_descrip;
    std::vector<hsize_t> dims_out_proj;
    for (string spe : config.species)
    {
        for (int lam = 0; lam < lmax[spe] + 1; lam++)
        {

            vec temp_power = readHDF5<double>(features, "sparse_descriptors/" + spe + "/" + to_string(lam), dims_out_descrip);
            power_env_sparse[spe + std::to_string(lam)] = temp_power;
            vec temp_proj = readHDF5<double>(projectors, "projectors/" + spe + "/" + to_string(lam), dims_out_proj);
            Vmat[spe + std::to_string(lam)] = reshape(temp_proj, Shape2D{(int)dims_out_proj[0], (int)dims_out_proj[1]});

            if (lam == 0)
            {
                Mspe[spe] = (int)dims_out_descrip[0];
            }
            if (config.zeta == 1)
            {
                load_BLAS();
                power_env_sparse[spe + std::to_string(lam)] = flatten(dot<double>(temp_proj, temp_power, (int)dims_out_proj[0], (int)dims_out_proj[1], (int)dims_out_descrip[0], (int)dims_out_descrip[1], true, false));
                unload_BLAS();
            }
        }
    }
    features.close();
    projectors.close();
#else
    err_not_impl_f("Not possible wihtout Rascaline and BLAS", std::cout);
#endif

    for (int lam = 0; lam < SALTED_Utils::get_lmax_max(lmax) + 1; lam++)
    {
        wigner3j[lam] = readVectorFromFile<double>(_opt.SALTED_DIR / "wigners" / ("wigner_lam - " + to_string(lam) + "_lmax1 - " + to_string(config.nang1) + "_lmax2 - " + to_string(config.nang2) + ".dat"));
    }

    if (config.sparsify)
    {
        filesystem::path path = _opt.SALTED_DIR / "GPR_data" / ("fps" + to_string(config.ncut) + "-");
        vfps = read_fps<int64_t>(path, SALTED_Utils::get_lmax_max(lmax));
    };
    if (config.average)
    {
        for (string spe : config.species)
        {
            filesystem::path path = _opt.SALTED_DIR / "averages" / ("averages_" + spe + ".npy");
            read_npy(path, av_coefs[spe]);
        }
    }

    int ntrain = static_cast<int>(config.Ntrain * config.trainfrac);

    if (config.field)
    {
        cout << "Field" << endl;
        err_not_impl_f("Calculations using 'Field = True' are not yet supported", std::cout);
    }
    else
    {
        filesystem::path path = _opt.SALTED_DIR / "GPR_data"/ ("weights_N" + to_string(ntrain) + "_reg-6.npy");
        read_npy(path, weights);
    }
}

#if has_RAS
void SALTEDPredictor::read_model_data_h5()
{
    using namespace std;
    const filesystem::path _H5path = _opt.SALTED_DIR / config.h5_filename;
    H5::H5File input(_H5path, H5F_ACC_RDONLY);
    vector<hsize_t> dims_out_descrip;
    vector<hsize_t> dims_out_proj;
    for (string spe : config.species)
    {
        for (int lam = 0; lam < lmax[spe] + 1; lam++)
        {
            filesystem::path spar_descrip = "sparse_descriptors";
            spar_descrip = spar_descrip/ spe / to_string(lam);
            filesystem::path proj = "projectors";
            proj = proj / spe / to_string(lam);
            vec temp_power = readHDF5<double>(input, spar_descrip, dims_out_descrip);
            power_env_sparse[spe + to_string(lam)] = temp_power;
            vec temp_proj = readHDF5<double>(input, proj, dims_out_proj);
            Vmat[spe + to_string(lam)] = reshape(temp_proj, Shape2D{(int)dims_out_proj[0], (int)dims_out_proj[1]});

            if (lam == 0)
            {
                Mspe[spe] = (int)dims_out_descrip[0];
            }
            if (config.zeta == 1)
            {
                load_BLAS();
                power_env_sparse[spe + to_string(lam)] = flatten(dot<double>(temp_proj, temp_power, (int)dims_out_proj[0], (int)dims_out_proj[1], (int)dims_out_descrip[0], (int)dims_out_descrip[1], true, false));
                unload_BLAS();
            }
        }
    }

    vector<hsize_t> dims_out_temp;
    filesystem::path wigner = "wigners";
    for (int lam = 0; lam < SALTED_Utils::get_lmax_max(lmax) + 1; lam++)
    {
        wigner3j[lam] = readHDF5<double>(input, wigner / ("lam-" + to_string(lam)), dims_out_temp);
    }

    if (config.sparsify)
    {
        filesystem::path path = "fps";
        for (int lam = 0; lam < SALTED_Utils::get_lmax_max(lmax) + 1; lam++)
        {
            vfps[lam] = readHDF5<int64_t>(input, path/("lam-" + to_string(lam)), dims_out_temp);
        }
    };
    if (config.average)
    {
        for (string spe : config.species)
        {
            av_coefs[spe] = readHDF5<double>(input, "averages/" + spe, dims_out_temp);
        }
    }

    if (config.field)
    {
        cout << "Field" << endl;
        err_not_impl_f("Calculations using 'Field = True' are not yet supported", std::cout);
    }
    else
    {
        weights = readHDF5<double>(input, "weights", dims_out_temp);
    }
}
#endif

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

    load_BLAS();

    std::vector<vec3> psi_nm(config.species.size());
    for (int spe_idx = 0; spe_idx < config.species.size(); spe_idx++)
    {
        const string spe = config.species[spe_idx];
        psi_nm[spe_idx].resize(lmax[spe] + 1);

        if (atom_idx.find(spe) == atom_idx.end())
        {
            continue;
        }
        vec2 kernell0_nm;
        for (int lam = 0; lam < lmax[spe] + 1; ++lam)
        {
            int lam2_1 = 2 * lam + 1;
            int row_size = featsize[lam] * lam2_1; // Size of a block of rows

            vec pvec_lam;
            pvec_lam.reserve(atom_idx[spe].size() * row_size); // Preallocate memory for pvec_lam to avoid reallocations

            for (int idx : atom_idx[spe])
            {
                int start_idx = idx * row_size;     // Start index in pvec_flat
                int end_idx = start_idx + row_size; // End index in pvec_flat

                // Copy the block directly into flatVec2
                pvec_lam.insert(pvec_lam.end(), pvec[lam].begin() + start_idx, pvec[lam].begin() + end_idx);
            }

            vec2 kernel_nm = dot<double>(pvec_lam, power_env_sparse[spe + to_string(lam)], natom_dict[spe] * lam2_1, featsize[lam], Mspe[spe] * lam2_1, featsize[lam], false, true);

            if (config.zeta == 1)
            {
                psi_nm[spe_idx][lam] = kernel_nm;
                continue;
            }

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
                                kernel_nm[i1 * lam2_1 + i][i2 * lam2_1 + j] *= pow(kernell0_nm[i1][i2], config.zeta - 1);
                            }
                        }
                    }
                }
            }
            psi_nm[spe_idx][lam] = dot<double>(kernel_nm, Vmat[spe + to_string(lam)], false, false);
        }
    }
    pvec.clear();
    pvec.shrink_to_fit();

    unordered_map<string, vec> C{};
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
                isize += static_cast<int>(Vmat[spe + to_string(l)][0].size()) * nmax[spe + to_string(l)];
            }
            continue;
        }
        ispe[spe] = 0;
        for (int l = 0; l < lmax[spe] + 1; ++l)
        {
            for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
            {
                // int Mcut = static_cast<int>(psi_nm[spe + to_string(l)][0].size());
                int Mcut = static_cast<int>(psi_nm[spe_idx][l][0].size());
                // Check if isize + Mcut > weights.size()
                err_chekf(isize + Mcut <= weights.size(), "isize + Mcut > weights.size()", std::cout);
                vec weights_subset(weights.begin() + isize, weights.begin() + isize + Mcut);
                C[spe + to_string(l) + to_string(n)] = dot(psi_nm[spe_idx][l], weights_subset, false);

                isize += Mcut;
            }
        }
    }
    psi_nm.clear();
    psi_nm.shrink_to_fit();
    // std::unordered_map<string, vec2>  temp;
    // psi_nm.swap(temp);

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
                std::copy_n(C[spe + to_string(l) + to_string(n)].begin() + ispe[spe] * (2 * l + 1), 2 * l + 1, pred_coefs.begin() + i);

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
    if (_opt.coef_file != "")
    {
        vec coefs{};
        cout << "Reading coefficients from file: " << _opt.coef_file << endl;
        read_npy<double>(_opt.coef_file, coefs);
        return coefs;
    }

    // Run generation of tsc file
    time_point start;
    if (_opt.debug)
        start = get_time();

    setup_atomic_environment();

    if (!config.from_h5)
    {
        read_model_data();
    }
    else
    {
#if has_RAS
        read_model_data_h5();
#else
        err_not_impl_f("HDF5 files are not supported by this build", std::cout);
#endif
    }

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
    std::unordered_map<std::string, vec2> umap;
    std::unordered_map<std::string, int> umap2;
    std::unordered_map<std::string, vec> umap3;
    Vmat.swap(umap);
    natom_dict.swap(umap2);
    lmax.swap(umap2);
    nmax.swap(umap2);
    power_env_sparse.swap(umap3);
};