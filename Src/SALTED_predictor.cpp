#include "SALTED_predictor.h"
#include <filesystem>

//std::string find_first_h5_file(const std::string& directory_path)
SALTEDPredictor::SALTEDPredictor(const WFN &wavy_in, const options &opt_in) : wavy(wavy_in), opt(opt_in)
{
    std::string _path = opt.SALTED_DIR;

    config.h5_filename = find_first_h5_file(opt.SALTED_DIR);
    
    if (config.h5_filename == "")
	{
        std::string _f_path("inputs.txt");
        join_path(_path, _f_path);
        config.populateFromFile(_path);
	}
	else
	{
#if has_RAS
    join_path(_path, config.h5_filename);
    H5::H5File config_file(_path, H5F_ACC_RDONLY);
    config.populateFromFile(config_file);
#else
    err_not_impl_f("HDF5 files are not supported by this build", std::cout);
#endif
	}
    config.predict_filename = wavy.get_path();


}

const std::string SALTEDPredictor::get_dfbasis_name()
{
	return config.dfbasis;
}

void calculateConjugate(std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& v2) {
    std::for_each(std::execution::par, v2.begin(), v2.end(), [](auto& vec3d) {
        std::for_each(vec3d.begin(), vec3d.end(), [](auto& vec2d) {
            std::for_each(vec2d.begin(), vec2d.end(), [](auto& vec1d) {
                std::transform(vec1d.begin(), vec1d.end(), vec1d.begin(), [](const std::complex<double>& val) {
                    return std::conj(val);
                    });
                });
            });
        });
}


void SALTEDPredictor::setup_atomic_environment()
{

    SALTED_Utils::set_lmax_nmax(lmax, nmax, *(wavy.get_basis_set_ptr()), config.species);

    for (int i = 0; i < wavy.atoms.size(); i++)
    {
        atomic_symbols.push_back(wavy.atoms[i].label);
    }
    // # Define system excluding atoms that belong to species not listed in SALTED input
    atomic_symbols = SALTED_Utils::filter_species(atomic_symbols, config.species);

    // Print all Atomic symbols
    if (opt.debug)
    {
        std::cout << "Atomic symbols: ";
        for (const auto& symbol : atomic_symbols)
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
        config.species
    ).calculate_expansion_coeffs();
    if ((config.nrad2 != config.nrad1) || (config.nang2 != config.nang1) || (config.sig2 != config.sig1) || (config.rcut2 != config.rcut1) || (config.neighspe2 != config.neighspe1)) {
        v2 = Rascaline_Descriptors(
            config.predict_filename, 
            config.nrad2, 
            config.nang2, 
            config.sig2, 
            config.rcut2, 
            natoms, 
            config.neighspe2, 
            config.species
        ).calculate_expansion_coeffs();
    }
    else{
		v2 = v1;
	}

    //Calculate the conjugate of v2 and store it back in v2, to avoid recalculating it in the equicomb function
    calculateConjugate(v2);
//#pragma omp parallel for
//    for (int i = 0; i < v2.size(); ++i)
//	{
//		for (int j = 0; j < v2[0].size(); ++j)
//		{
//			for (int k = 0; k < v2[0][0].size(); ++k)
//			{
//				for (int l = 0; l < v2[0][0][0].size(); ++l)
//				{
//					v2[i][j][k][l] = conj(v2[i][j][k][l]);
//				}
//			}
//		}
//	}
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

    H5::H5File features(opt.SALTED_DIR + "/GPR_data/FEAT_M-" + std::to_string(config.Menv) + ".h5", H5F_ACC_RDONLY);
    H5::H5File projectors(opt.SALTED_DIR + "/GPR_data/projector_M" + std::to_string(config.Menv) + "_zeta" + zeta_str + ".h5", H5F_ACC_RDONLY);
    std::vector<hsize_t> dims_out_descrip;
    std::vector<hsize_t> dims_out_proj;
    for (string spe : config.species)
    {
        for (int lam = 0; lam < lmax[spe] + 1; lam++)
        {
            vec temp_power = readHDF5<double>(features, "sparse_descriptors/" + spe + "/" + to_string(lam), dims_out_descrip);
            power_env_sparse[spe + std::to_string(lam)] = temp_power;
            vec temp_proj = readHDF5<double>(projectors, "projectors/" + spe + "/" + to_string(lam), dims_out_proj);
            Vmat[spe + std::to_string(lam)] = reshape(temp_proj, Shape2D{ (int)dims_out_proj[0], (int)dims_out_proj[1] });

            if (lam == 0)
            {
                Mspe[spe] = (int)dims_out_descrip[0];
            }
            if (config.zeta == 1)
            {
                power_env_sparse[spe + std::to_string(lam)] = flatten(dot<double>(temp_proj, temp_power, (int)dims_out_proj[0], (int)dims_out_proj[1], (int)dims_out_descrip[0], (int)dims_out_descrip[1], true, false));
            }
        }
    }
    features.close();
    projectors.close();

    for (int lam = 0; lam < SALTED_Utils::get_lmax_max(lmax) + 1; lam++)
    {
        wigner3j[lam] = readVectorFromFile<double>(opt.SALTED_DIR + "/wigners/wigner_lam-" + to_string(lam) + "_lmax1-" + to_string(config.nang1) + "_lmax2-" + to_string(config.nang2) + ".dat");
    }
    


    if (config.sparsify)
    {
        std::string path = opt.SALTED_DIR;
        join_path(path, { "GPR_data", "fps" + to_string(config.ncut) + "-" });
        vfps = read_fps<int64_t>(path, SALTED_Utils::get_lmax_max(lmax));
    };
    if (config.average) {
        for (string spe : config.species)
        {
            string path = opt.SALTED_DIR;
            join_path(path, { "averages", "averages_" + spe + ".npy" });
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
        string path = opt.SALTED_DIR;
        join_path(path, { "GPR_data", "weights_N" + to_string(ntrain) + "_reg-6.npy"});
        read_npy(path, weights);
    }
}

#if has_RAS
void SALTEDPredictor::read_model_data_h5() {
    using namespace std;
    string _H5path = opt.SALTED_DIR;
    join_path(_H5path, config.h5_filename);
    H5::H5File input(_H5path, H5F_ACC_RDONLY);
    vector<hsize_t> dims_out_descrip;
    vector<hsize_t> dims_out_proj;
    for (string spe : config.species)
    {
        for (int lam = 0; lam < lmax[spe] + 1; lam++)
        {
            vec temp_power = readHDF5<double>(input, "sparse_descriptors/" + spe + "/" + to_string(lam), dims_out_descrip);
            power_env_sparse[spe + to_string(lam)] = temp_power;
            vec temp_proj = readHDF5<double>(input, "projectors/" + spe + "/" + to_string(lam), dims_out_proj);
            Vmat[spe + to_string(lam)] = reshape(temp_proj, Shape2D{(int)dims_out_proj[0], (int)dims_out_proj[1]});

            if (lam == 0)
            {
                Mspe[spe] = (int)dims_out_descrip[0];
            }
            if (config.zeta == 1)
            {
                power_env_sparse[spe + to_string(lam)] = flatten(dot<double>(temp_proj, temp_power, (int)dims_out_proj[0], (int)dims_out_proj[1], (int)dims_out_descrip[0], (int)dims_out_descrip[1], true, false));
            }
        }
    }

    vector<hsize_t> dims_out_temp;
    for (int lam = 0; lam < SALTED_Utils::get_lmax_max(lmax) + 1; lam++)
    {
        wigner3j[lam] = readHDF5<double>(input, "wigners/lam-" + to_string(lam), dims_out_temp);
    }
    
    if (config.sparsify)
    {
        for (int lam = 0; lam < SALTED_Utils::get_lmax_max(lmax) + 1; lam++)
        {
            vfps[lam] = readHDF5<int64_t>(input, "fps/lam-" + to_string(lam), dims_out_temp);
        }
    };
    if (config.average) {
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
    unordered_map<int, vec> pvec{};
    for (int lam = 0; lam < SALTED_Utils::get_lmax_max(lmax) + 1; lam++)
    {
        int llmax = 0;
        unordered_map<int, ivec> lvalues{};
        for (int l1 = 0; l1 < config.nang1 + 1; l1++)
        {
            for (int l2 = 0; l2 < config.nang2 + 1; l2++)
            {
                // keep only even combination to enforce inversion symmetryc
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
        vec p_test;
        ivec2 llvec_t = transpose<int>(llvec);
        if (config.sparsify)
        {
            int nfps = static_cast<int>(vfps[lam].size());
            equicomb(natoms, (config.nspe1 * config.nrad1), (config.nspe2 * config.nrad2), v1, v2, wigner3j[lam], llvec_t, lam, c2r, featsize[lam], nfps, vfps[lam], p);
            featsize[lam] = nfps;
        }
        else
        {
            equicomb(natoms, config.nang1, config.nang2, (config.nspe1 * config.nrad1), (config.nspe2 * config.nrad2), v1, v2, wigner3j[lam], llmax, llvec_t, lam, c2r, featsize[lam], p);
        }


        pvec[lam] = p;
    }



    unordered_map<string, vec2> psi_nm{};
    for (const string &spe : config.species)
    {
        if (atom_idx.find(spe) == atom_idx.end())
        {
            continue;
        }
        vec2 kernell0_nm;
        for (int lam = 0; lam < lmax[spe] + 1; ++lam) {
            int lam2_1  = 2 * lam + 1;
            
            //This block collects the rows of pvec corresponding to the atoms of species spe
            vec pvec_lam;
            pvec_lam.reserve(atom_idx[spe].size() * lam2_1 * featsize[lam]);  // Preallocate memory for pvec_lam to avoid reallocations
            int row_size = featsize[lam] * lam2_1;  // Size of a block of rows

            for (int idx : atom_idx[spe]) {
                int start_idx = idx * row_size;  // Start index in pvec_flat
                int end_idx = start_idx + row_size;  // End index in pvec_flat

                // Copy the block directly into flatVec2
                pvec_lam.insert(pvec_lam.end(), pvec[lam].begin() + start_idx, pvec[lam].begin() + end_idx);
            }

            vec2 kernel_nm = dot<double>(pvec_lam, power_env_sparse[spe + to_string(lam)], natom_dict[spe] * lam2_1, featsize[lam], Mspe[spe] * lam2_1, featsize[lam], false, true);

            if (config.zeta == 1) {
                psi_nm[spe + to_string(lam)] = kernel_nm;
                continue;
            }

            if (lam == 0) {
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
            psi_nm[spe + to_string(lam)] = dot<double>(kernel_nm, Vmat[spe + to_string(lam)], false, false);
        }
    }
    pvec.clear();

    unordered_map<string, vec> C{};
    unordered_map<string, int> ispe{};
    int isize = 0;
    for (const string &spe : config.species)
    {
        if (atom_idx.find(spe) == atom_idx.end())
        {
            for (int l = 0; l < lmax[spe] + 1; ++l)
            {
                //Check if Vmat[spe + to_string(l)][0] exists
                if (Vmat[spe + to_string(l)].size() == 0) {
                    cout << "The projector for species " << spe << " and l = " << l << " does not exist. This is a problem with the model, not NoSpherA2." << endl;
                    cout << "Continuing with the next species..., make sure there is no: " << spe << " in the structure you are trying to predict!!!!" << endl;
                    l = lmax[spe] + 1;
                    continue;
                }
               
                //for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
                //{
                //    isize += static_cast<int>(Vmat[spe + to_string(l)][0].size());
                //}
                isize += static_cast<int>(Vmat[spe + to_string(l)][0].size()) * nmax[spe + to_string(l)];
            }
            continue;
        }
        ispe[spe] = 0;
        for (int l = 0; l < lmax[spe] + 1; ++l)
        {
            for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
            {
                int Mcut = static_cast<int>(psi_nm[spe + to_string(l)][0].size());
                //Check if isize + Mcut > weights.size()
                err_chekf(isize + Mcut <= weights.size(), "isize + Mcut > weights.size()", std::cout);
                vec weights_subset(weights.begin() + isize, weights.begin() + isize + Mcut);
                C[spe + to_string(l) + to_string(n)] = dot(psi_nm[spe + to_string(l)], weights_subset, false);

                isize += Mcut;
            }
        }
    }

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
                //for (int ind = 0; ind < 2 * l + 1; ++ind)
                //{
                //    pred_coefs[i + ind] = C[spe + to_string(l) + to_string(n)][ispe[spe] * (2 * l + 1) + ind];
                //}
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

    //cout << "          ... done!\nNumber of predicted coefficients: " << pred_coefs.size() << endl;
    //npy::npy_data<double> coeffs;
    //coeffs.data = pred_coefs;
    //coeffs.fortran_order = false;
    //coeffs.shape = { unsigned long(pred_coefs.size()) };
    //npy::write_npy("folder_model.npy", coeffs);
    return pred_coefs;
}

vec SALTEDPredictor::gen_SALTED_densities()
{
    using namespace std;
    if (opt.coef_file != "")
	{
        vec coefs{};
		cout << "Reading coefficients from file: " << opt.coef_file << endl;
		read_npy<double>(opt.coef_file, coefs);
		return coefs;
	}

    // Run generation of tsc file
    time_point start;
    if (opt.debug) start = get_time();

    setup_atomic_environment();

    if (!config.from_h5) {
        read_model_data();
    }
    else {
#if has_RAS
        read_model_data_h5();
#else
        err_not_impl_f("HDF5 files are not supported by this build", std::cout);
#endif
        
    }


    vec coefs = predict();


    if (opt.debug)
    {
        time_point end = get_time();
        long long int dur = get_sec(start, end);
        cout << "Finished ML Density Prediction!" << endl;
        if (dur < 1)
            cout << "Time for SALTED prediction: " << fixed << setprecision(0) << get_msec(start, end) << " ms" << endl;
        else
            cout << "Time for SALTED prediction: " << fixed << setprecision(0) << dur << " s" << endl;

        cout << "Number of coefficients: " << coefs.size() << endl;

        if (opt.wfn == string("test_cysteine2.xyz"))
        {
            vector<unsigned long> shape{};
            bool fortran_order;
            vec ref_coefs{};
            npy::LoadArrayFromNumpy("test_cysteine.npy", shape, fortran_order, ref_coefs);
            // Compare coefs with the reference
            vector<double> diff_vec;
            double diff = 0.0;
            for (int i = 0; i < coefs.size(); i++)
            {
                diff += abs(coefs[i] - ref_coefs[i]);
                if (abs(coefs[i] - ref_coefs[i]) > 1e-4 || isnan(coefs[i] - ref_coefs[i]))
                {
                    cout << "Difference in coef " << fixed << setprecision(3) << i << " : " << coefs[i] << " - " << ref_coefs[i] << " = " << abs(coefs[i] - ref_coefs[i]) << endl;
                }
                diff_vec.push_back((coefs[i] / ref_coefs[i]) - 1);
            }
            cout << "Difference between calculated and reference coefs: " << fixed << setprecision(3) << diff << endl;
            cout << "Maximum ((pred / ref) -1): " << fixed << setprecision(3) << *max_element(diff_vec.begin(), diff_vec.end()) << endl;
            if (diff > 0.1 || isnan(diff))
            {
                cout << "WARNING: The difference between the calculated and reference coefficients is too large!" << endl;
                exit(1);
            }
        }
    }

    return coefs;
}