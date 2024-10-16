#include "SALTED_predictor.h"
#include <filesystem>

// std::string find_first_h5_file(const std::string& directory_path)
SALTEDPredictor::SALTEDPredictor(const WFN &wavy_in, const options &opt_in) : wavy(wavy_in), opt(opt_in)
{
    std::string _path = opt.SALTED_DIR;

    config.h5_filename = find_first_h5_file(opt.SALTED_DIR);
    config.predict_filename = wavy.get_path();

    if (config.h5_filename == "")
    {
        if (opt.debug)
        {
            std::cout << "No HDF5 file found in the SALTED directory. Using inputs.txt instead." << std::endl;
        }
        std::string _f_path("inputs.txt");
        join_path(_path, _f_path);
        config.populateFromFile(_path);
    }
    else
    {
        if (opt.debug)
        {
            std::cout << "Using HDF5 file: " << config.h5_filename << std::endl;
        }
#if has_RAS
        join_path(_path, config.h5_filename);
        H5::H5File config_file(_path, H5F_ACC_RDONLY);
        config.populateFromFile(config_file);
#else
        err_not_impl_f("HDF5 files are not supported by this build", std::cout);
#endif
    }
    config.neighspe1 = { "H", "O" , "C", "N", "S" };
    config.neighspe2 = { "H", "O" , "C", "N", "S" };
    config.species = { "H", "O" , "C", "N", "S" };
    config.nspe1 = 5;
    config.nspe2 = 5;
    config.nang1 = 6;
    config.nang2 = 6;
    config.nrad1 = 5;
    config.nrad2 = 5;
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
            Vmat[spe + std::to_string(lam)] = reshape(temp_proj, Shape2D{(int)dims_out_proj[0], (int)dims_out_proj[1]});

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
        wigner3j_old[lam] = readVectorFromFile<double>(opt.SALTED_DIR + "/wigners/wigner_lam-" + to_string(lam) + "_lmax1-" + to_string(config.nang1) + "_lmax2-" + to_string(config.nang2) + ".dat");
    }
    


    if (config.sparsify)
    {
        std::string path = opt.SALTED_DIR;
        join_path(path, {"GPR_data", "fps" + to_string(config.ncut) + "-"});
        vfps = read_fps<int64_t>(path, SALTED_Utils::get_lmax_max(lmax));
    };
    if (config.average)
    {
        for (string spe : config.species)
        {
            string path = opt.SALTED_DIR;
            join_path(path, {"averages", "averages_" + spe + ".npy"});
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
        join_path(path, {"GPR_data", "weights_N" + to_string(ntrain) + "_reg-6.npy"});
        read_npy(path, weights);
    }
}

#if has_RAS
void SALTEDPredictor::read_model_data_h5()
{
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
        wigner3j_old[lam] = readHDF5<double>(input, "wigners_old/lam-" + to_string(lam), dims_out_temp);
    }

    for (int lam = 0; lam < SALTED_Utils::get_lmax_max(lmax) + 1; lam++)
    {
        wigner3j_new[lam] = readHDF5<double>(input, "wigners_new/lam-" + to_string(lam), dims_out_temp);
    }

    if (config.sparsify)
    {
        for (int lam = 0; lam < SALTED_Utils::get_lmax_max(lmax) + 1; lam++)
        {
            vfps[lam] = readHDF5<int64_t>(input, "fps/lam-" + to_string(lam), dims_out_temp);
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

    //Get the maximum neccecary lambda value
    int lam_max = SALTED_Utils::get_lmax_max(lmax);

    // Compute equivariant descriptors for each lambda value entering the SPH expansion of the electron density
    unordered_map<int, vec> pvec{};

    //Old Version of equicomb
    //if (!config.fast_equicomb)
    {
        for (int lam = 0; lam < lam_max + 1; lam++)
        {
            std::cout << "Calculating descriptors for l = " << lam << std::endl;
            ivec2 llvec = calc_llvec(config.nang1, config.nang2, lam);
            int llmax = llvec.size();

            cvec2 c2r = SALTED_Utils::complex_to_real_transformation({ 2 * lam + 1 })[0];

            featsize[lam] = config.nspe1 * config.nspe2 * config.nrad1 * config.nrad2 * llmax;
            vec p;
            vec p_test;
            ivec2 llvec_t = transpose<int>(llvec);
            if (config.sparsify)
            {
                int nfps = static_cast<int>(vfps[lam].size());
                equicomb(natoms, (config.nspe1 * config.nrad1), (config.nspe2 * config.nrad2), v1, v2, wigner3j_old[lam], llvec_t, lam, c2r, featsize[lam], nfps, vfps[lam], p);
                featsize[lam] = nfps;
            }
            else
            {
                equicomb(natoms, config.nang1, config.nang2, (config.nspe1 * config.nrad1), (config.nspe2 * config.nrad2), v1, v2, wigner3j_old[lam], llvec_t, lam, c2r, featsize[lam], p);
            }

            pvec[lam] = p;
        }
	}
    //else
    unordered_map<int, vec> old_pvec = pvec;
	pvec.clear();
    {
        //Matrix of size l1 x l2 with the highest lam value for each combination (initialized as -1 )
        ivec2 max_lam_per_l1l2(config.nang1 + 1, ivec(config.nang2 + 1, -1));

        //For every lam the corresponding l1 and l2 values (basically llvec from the old version)
        std::unordered_map<int, ivec2> l1l2_per_lam;
        get_angular_indexes_symmetric_per_lambda(lam_max, config.nang1, config.nang2, max_lam_per_l1l2, l1l2_per_lam);

        std::cout << "Calculating combinations of v1 and v2..." << std::endl;
        std::unordered_map<std::string, ivec> feats_per_l1l2;
        std::unordered_map<int, cvec> prod_vec_map;
        equicomb_vec_multiply(natoms, (config.nspe1 * config.nrad1), (config.nspe2 * config.nrad2), config.nang1, config.nang2, v1, v2, max_lam_per_l1l2, lam_max, prod_vec_map, feats_per_l1l2);

        //int featsize_new = (config.nang1 + 1) * (config.nang2 + 1) * config.nspe1 * config.nrad1 * config.nspe2 * config.nrad2;
        const int n_feats_per_atom = config.nrad1 * config.nspe1 * config.nrad2 * config.nspe1;
        const int n_feats = n_feats_per_atom * natoms;

        //#pragma omp parallel
        for (int lam = 0; lam < lam_max + 1; lam++) {
            const int l21 = 2 * lam + 1;
            const ivec2 l1l2 = l1l2_per_lam[lam];
            const int l1l2_size = l1l2.size();

            // Wigner 3j coefficients
            vec wigner3j = wigner3j_new[lam];
            cvec2 c2r_lam = SALTED_Utils::complex_to_real_transformation({ 2 * lam + 1 })[0];
            cvec c2r_lam_flat = flatten(c2r_lam);

            int iwig = 0;
            vec normfacts(natoms, 0.0);
            int featsize_orig = config.nspe1 * config.nspe2 * config.nrad1 * config.nrad2 * l1l2_size;

            std::unordered_map<int,cvec2> p_temp;

            for (int l1l2_index = 0; l1l2_index < l1l2_size; l1l2_index++) {
                ivec vec = l1l2[l1l2_index];

                int size_inner = 2 * vec[0] + 1;
                int max_lam_inner = max_lam_per_l1l2[vec[0]][vec[1]];
                int start = (max_lam_inner - lam) * size_inner;
                int end = (max_lam_inner + lam + 1) * size_inner;

                ivec feat_indices = feats_per_l1l2[std::to_string(vec[0]) + std::to_string(vec[1])];

                cvec calc_vec(n_feats * l21, 0.0);

				//calculate the product of v1,v1 with the corresponding wigner3j-symbols and sum
    //#pragma omp for
                for (int i = 0; i < n_feats; i++) {
                    cvec prod(size_inner * l21, 0.0);
                    cvec::iterator temp = prod_vec_map[feat_indices[i]].begin();
                    std::transform(temp + start, temp + end, wigner3j.begin() + iwig, prod.begin(), std::multiplies<cdouble>());
                    //Calc the sum over every size_inner for every l21
                    for (int l = 0; l < l21; l++) {
                        for (int inner = 0; inner < size_inner; inner++) {
                            calc_vec[i * l21 + l] += prod[l * size_inner + inner];
                        }
                    }
                }

                cvec2 res = dot(calc_vec, c2r_lam_flat, n_feats, l21, l21, l21, false, true);
                
				//Calculate the norm factor for each atom
                //#pragma omp for
                for (int iat = 0; iat < natoms; ++iat) {
                    double normfact_sum = 0.0;
                    int base_idx = iat * n_feats_per_atom;

                    for (int i = 0; i < n_feats_per_atom; ++i) {
                        const cvec& res_row = res[base_idx + i];
                        for (int l = 0; l < l21; ++l) {
                            normfact_sum += std::real(res_row[l]) * std::real(res_row[l]);
                        }
                    }
                    normfacts[iat] += normfact_sum;
                }
                iwig += l21 * size_inner;
                p_temp[l1l2_index] = res;
            }
            
			//p of size (n1 * n2 * llmax,l21, natoms)
            vec3 p;
            if (config.sparsify) {
                p.assign(natoms, vec2(l21, vec(vfps[lam].size(), 0.0)));
            }
            else {
				p.assign(natoms, vec2(l21, vec(featsize_orig, 0.0)));
            }
            //ptemp of size(llmax, natoms*n1*n2, l21)

            for (int iat = 0; iat < natoms; ++iat) {
                normfacts[iat] = std::sqrt(normfacts[iat]);
            }
			int rad1spe1 = config.nrad1 * config.nspe1;
			int rad2spe2 = config.nrad2 * config.nspe2;
            //for (int atom = 0; atom < natoms; ++atom) {
            //    for (int i1 = 0; i1 < rad1spe1; ++i1) {
            //        for (int i2 = 0; i2 < rad2spe2; ++i2) {
            //            for (int ll = 0; ll < l1l2_size; ++ll) {
            //                int indxp = i1 * rad2spe2 * l1l2_size + i2 * l1l2_size + ll;
     
            //                int indxp_temp = atom * n_feats_per_atom + i1 * rad2spe2 + i2;
            //                for (int l = 0; l < l21; ++l) {
            //                    p[atom][l][indxp] = std::real(p_temp[ll][indxp_temp][l]) / normfacts[atom];
            //                }
            //            }
            //        }
            //    }
            //}
            for (int atom = 0; atom < natoms; ++atom) {
                for (int fps : vfps[lam]) {
					int i1 = fps / rad2spe2;
					int i2 = fps % rad2spe2;
					int ll = (fps / rad2spe2) % l1l2_size;


                    int indxp_temp = atom * n_feats_per_atom + i1 * rad2spe2 + i2;
                    for (int l = 0; l < l21; ++l) {
                        p[atom][l][fps] = std::real(p_temp[ll][indxp_temp][l]) / normfacts[atom];
                    }
                }
            }


            pvec[lam] = flatten(p);
        }
    }

	for (int lam = 0; lam < lam_max + 1; lam++)
	{
		std::cout << "Old: " << old_pvec[lam].size() << " New: " << pvec[lam].size() << std::endl;
		for (int i = 0; i < pvec[lam].size(); i++)
		{
			if (std::abs(old_pvec[lam][i] - pvec[lam][i]) > 1e-6)
			{
                std::cout << "Error at index " << i << " for lam " << lam  << "  " << std::to_string(std::abs(old_pvec[lam][i] - pvec[lam][i])) << std::endl;
                //find if the old_pvec[lam][i] is anywhere in pvec[lam] and if so print the index
                bool found = false;
				for (int j = i; j < pvec[lam].size(); j++)
				{
					if (std::abs(old_pvec[lam][i] - pvec[lam][j]) < 1e-6)
					{
                        found = true;
                        break;
                    }
				}
				if (!found)
				{
					std::cout << "Not found" << std::endl;
				}
			}
		}
	}


    unordered_map<string, vec2> psi_nm{};
    for (const string &spe : config.species)
    {
        if (atom_idx.find(spe) == atom_idx.end())
        {
            continue;
        }
        vec2 kernell0_nm;
        for (int lam = 0; lam < lmax[spe] + 1; ++lam)
        {
            int lam2_1 = 2 * lam + 1;

            // This block collects the rows of pvec corresponding to the atoms of species spe
            vec pvec_lam;
            pvec_lam.reserve(atom_idx[spe].size() * lam2_1 * featsize[lam]); // Preallocate memory for pvec_lam to avoid reallocations
            int row_size = featsize[lam] * lam2_1;                           // Size of a block of rows

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
                psi_nm[spe + to_string(lam)] = kernel_nm;
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
                int Mcut = static_cast<int>(psi_nm[spe + to_string(l)][0].size());
                // Check if isize + Mcut > weights.size()
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
    if (opt.coef_file != "")
    {
        vec coefs{};
        cout << "Reading coefficients from file: " << opt.coef_file << endl;
        read_npy<double>(opt.coef_file, coefs);
        return coefs;
    }

    // Run generation of tsc file
    time_point start;
    if (opt.debug)
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