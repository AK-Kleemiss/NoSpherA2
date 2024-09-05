#include "SALTED_predictor.h"
#include <filesystem>

using namespace std;
//std::string find_first_h5_file(const std::string& directory_path)
SALTEDPredictor::SALTEDPredictor(const WFN &wavy_in, const options &opt_in) : wavy(wavy_in), opt(opt_in)
{
    string _path = this->opt.SALTED_DIR;

    this->config.h5_filename = find_first_h5_file(this->opt.SALTED_DIR);
    
    if (this->config.h5_filename == "")
	{
        string _f_path("inputs.txt");
        join_path(_path, _f_path);
        this->config.populateFromFile(_path);
	}
	else
	{
#if has_RAS
    join_path(_path, this->config.h5_filename);
    H5::H5File config_file(_path, H5F_ACC_RDONLY);
    this->config.populateFromFile(config_file);
#else
    err_not_impl_f("HDF5 files are not supported by this build", std::cout);
#endif
	}
    this->config.predict_filename = this->wavy.get_path();


}

const string SALTEDPredictor::get_dfbasis_name()
{
	return this->config.dfbasis;
}



void SALTEDPredictor::setup_atomic_environment()
{

    SALTED_Utils::set_lmax_nmax(this->lmax, this->nmax, *(this->wavy.get_basis_set_ptr()), this->config.species);

    for (int i = 0; i < this->wavy.atoms.size(); i++)
    {
        this->atomic_symbols.push_back(this->wavy.atoms[i].label);
    }
    // # Define system excluding atoms that belong to species not listed in SALTED input
    this->atomic_symbols = SALTED_Utils::filter_species(this->atomic_symbols, this->config.species);

    // Print all Atomic symbols
    if (this->opt.debug)
    {
        cout << "Atomic symbols: ";
        for (const auto& symbol : atomic_symbols)
        {
            cout << symbol << " ";
        }
        cout << endl;
    }

    this->natoms = static_cast<int>(this->atomic_symbols.size());
    for (int i = 0; i < this->atomic_symbols.size(); i++)
    {
        this->atom_idx[this->atomic_symbols[i]].push_back(i);
        this->natom_dict[this->atomic_symbols[i]] += 1;
    }

    // RASCALINE (Generate descriptors)
#if has_RAS
    v1 = Rascaline_Descriptors(this->config.predict_filename, this->config.nrad1, this->config.nang1, this->config.sig1, this->config.rcut1, this->natoms, this->config.neighspe1, this->config.species).calculate_expansion_coeffs();
    if ((this->config.nrad2 != this->config.nrad1) || (this->config.nang2 != this->config.nang1) || (this->config.sig2 != this->config.sig1) || (this->config.rcut2 != this->config.rcut1) || (this->config.neighspe2 != this->config.neighspe1)) {
        v2 = Rascaline_Descriptors(this->config.predict_filename, this->config.nrad2, this->config.nang2, this->config.sig2, this->config.rcut2, this->natoms, this->config.neighspe2, this->config.species).calculate_expansion_coeffs();
    }else{
		v2 = v1;
	}
#else
    err_not_impl_f("RASCALINE is not supported by this build", std::cout);
#endif
    // END RASCALINE
}

void SALTEDPredictor::read_model_data()
{
    // Define zeta as a string with one decimal
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(1) << this->config.zeta;
    std::string zeta_str = stream.str();

    H5::H5File features(this->opt.SALTED_DIR + "/GPR_data/FEAT_M-" + to_string(this->config.Menv) + ".h5", H5F_ACC_RDONLY);
    H5::H5File projectors(this->opt.SALTED_DIR + "/GPR_data/projector_M" + to_string(this->config.Menv) + "_zeta" + zeta_str + ".h5", H5F_ACC_RDONLY);
    vector<hsize_t> dims_out_descrip;
    vector<hsize_t> dims_out_proj;
    for (string spe : config.species)
    {
        for (int lam = 0; lam < lmax[spe] + 1; lam++)
        {
            vec temp_power = readHDF5<double>(features, "sparse_descriptors/" + spe + "/" + to_string(lam), dims_out_descrip);
            this->power_env_sparse[spe + to_string(lam)] = temp_power;
            vec temp_proj = readHDF5<double>(projectors, "projectors/" + spe + "/" + to_string(lam), dims_out_proj);
            this->Vmat[spe + to_string(lam)] = reshape(temp_proj, Shape2D{ dims_out_proj[0], dims_out_proj[1] });

            if (lam == 0)
            {
                this->Mspe[spe] = dims_out_descrip[0];
            }
            if (this->config.zeta == 1)
            {
                this->power_env_sparse[spe + to_string(lam)] = flatten(dot<double>(temp_proj, temp_power, dims_out_proj[0], dims_out_proj[1], dims_out_descrip[0], dims_out_descrip[1], true, false));
            }
        }
    }
    features.close();
    projectors.close();

    for (int lam = 0; lam < SALTED_Utils::get_lmax_max(lmax) + 1; lam++)
    {
        this->wigner3j[lam] = readVectorFromFile<double>(this->opt.SALTED_DIR + "/wigners/wigner_lam-" + to_string(lam) + "_lmax1-" + to_string(this->config.nang1) + "_lmax2-" + to_string(this->config.nang2) + ".dat");
    }
    


    if (config.sparsify)
    {
        this->vfps = read_fps<int64_t>(this->opt.SALTED_DIR + "/GPR_data/fps" + to_string(this->config.ncut) + "-", SALTED_Utils::get_lmax_max(this->lmax));
    };
    if (config.average) {
        for (string spe : this->config.species)
        {
            read_npy(this->opt.SALTED_DIR + "/averages/averages_" + spe + ".npy", this->av_coefs[spe]);
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
        read_npy(this->opt.SALTED_DIR + "/GPR_data/weights_N" + to_string(ntrain) + "_reg-6.npy", this->weights);
    }
}

#if has_RAS
void SALTEDPredictor::read_model_data_h5() {
    string _H5path = this->opt.SALTED_DIR;
    join_path(_H5path, this->config.h5_filename);
    H5::H5File input(_H5path, H5F_ACC_RDONLY);
    vector<hsize_t> dims_out_descrip;
    vector<hsize_t> dims_out_proj;
    for (string spe : config.species)
    {
        for (int lam = 0; lam < lmax[spe] + 1; lam++)
        {
            vec temp_power = readHDF5<double>(input, "sparse_descriptors/" + spe + "/" + to_string(lam), dims_out_descrip);
            this->power_env_sparse[spe + to_string(lam)] = temp_power;
            vec temp_proj = readHDF5<double>(input, "projectors/" + spe + "/" + to_string(lam), dims_out_proj);
            this->Vmat[spe + to_string(lam)] = reshape(temp_proj, Shape2D{dims_out_proj[0], dims_out_proj[1]});

            if (lam == 0)
            {
                this->Mspe[spe] = dims_out_descrip[0];
            }
            if (this->config.zeta == 1)
            {
                this->power_env_sparse[spe + to_string(lam)] = flatten(dot<double>(temp_proj, temp_power,dims_out_proj[0], dims_out_proj[1], dims_out_descrip[0], dims_out_descrip[1], true, false));
            }
        }
    }

    vector<hsize_t> dims_out_temp;
    for (int lam = 0; lam < SALTED_Utils::get_lmax_max(lmax) + 1; lam++)
    {
        this->wigner3j[lam] = readHDF5<double>(input, "wigners/lam-" + to_string(lam), dims_out_temp);
    }
    
    if (config.sparsify)
    {
        for (int lam = 0; lam < SALTED_Utils::get_lmax_max(this->lmax) + 1; lam++)
        {
            this->vfps[lam] = readHDF5<int64_t>(input, "fps/lam-" + to_string(lam), dims_out_temp);
        }
    };
    if (config.average) {
        for (string spe : this->config.species)
        {
            this->av_coefs[spe] = readHDF5<double>(input, "averages/" + spe, dims_out_temp);
        }
    }

    if (config.field)
    {
        cout << "Field" << endl;
        err_not_impl_f("Calculations using 'Field = True' are not yet supported", std::cout);
    }
    else
    {
        this->weights = readHDF5<double>(input, "weights", dims_out_temp);
    }
}
#endif

vec SALTEDPredictor::predict()
{
    // Compute equivariant descriptors for each lambda value entering the SPH expansion of the electron density
    unordered_map<int, vec> pvec{};
    for (int lam = 0; lam < SALTED_Utils::get_lmax_max(lmax) + 1; lam++)
    {
        int llmax = 0;
        unordered_map<int, ivec> lvalues{};
        for (int l1 = 0; l1 < this->config.nang1 + 1; l1++)
        {
            for (int l2 = 0; l2 < this->config.nang2 + 1; l2++)
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

        this->featsize[lam] = this->config.nspe1 * this->config.nspe2 * this->config.nrad1 * this->config.nrad2 * llmax;
        vec p;
        ivec2 llvec_t = transpose<int>(llvec);
        if (this->config.sparsify)
        {
            int nfps = static_cast<int>(vfps[lam].size());
            equicomb(this->natoms, this->config.nang1, this->config.nang2, (this->config.nspe1 * this->config.nrad1), (this->config.nspe2 * this->config.nrad2), this->v1, this->v2, this->wigner3j[lam], llmax, llvec_t, lam, c2r, this->featsize[lam], nfps, this->vfps[lam], p);
            this->featsize[lam] = nfps;
        }
        else
        {
            equicomb(this->natoms, this->config.nang1, this->config.nang2, (this->config.nspe1 * this->config.nrad1), (this->config.nspe2 * this->config.nrad2), this->v1, this->v2, this->wigner3j[lam], llmax, llvec_t, lam, c2r, this->featsize[lam], p);
        }
        pvec[lam] = p;
    }

    unordered_map<string, vec2> psi_nm{};
    for (const string &spe : this->config.species)
    {
        if (this->atom_idx.find(spe) == this->atom_idx.end())
        {
            continue;
        }
        vec2 kernell0_nm;
        for (int lam = 0; lam < lmax[spe] + 1; ++lam) {
            int lam2_1  = 2 * lam + 1;
            
            //This block collects the rows of pvec corresponding to the atoms of species spe
            vec pvec_lam;
            pvec_lam.reserve(this->atom_idx[spe].size() * lam2_1 * this->featsize[lam]);  // Preallocate memory for pvec_lam to avoid reallocations
            int row_size = this->featsize[lam] * lam2_1;  // Size of a block of rows

            for (int idx : this->atom_idx[spe]) {
                int start_idx = idx * row_size;  // Start index in pvec_flat
                int end_idx = start_idx + row_size;  // End index in pvec_flat

                // Copy the block directly into flatVec2
                pvec_lam.insert(pvec_lam.end(), pvec[lam].begin() + start_idx, pvec[lam].begin() + end_idx);
            }

            vec2 kernel_nm = dot<double>(pvec_lam, this->power_env_sparse[spe + to_string(lam)], this->natom_dict[spe] * lam2_1, this->featsize[lam], this->Mspe[spe] * lam2_1, this->featsize[lam], false, true);

            if (config.zeta == 1) {
                psi_nm[spe + to_string(lam)] = kernel_nm;
                continue;
            }

            if (lam == 0) {
                kernell0_nm = kernel_nm;
                kernel_nm = elementWiseExponentiation(kernel_nm, this->config.zeta);
            }
            else
            {
                for (size_t i1 = 0; i1 < this->natom_dict[spe]; ++i1)
                {
                    for (size_t i2 = 0; i2 < this->Mspe[spe]; ++i2)
                    {
                        for (size_t i = 0; i < lam2_1; ++i)
                        {
                            for (size_t j = 0; j < lam2_1; ++j)
                            {
                                kernel_nm[i1 * lam2_1 + i][i2 * lam2_1 + j] *= pow(kernell0_nm[i1][i2], this->config.zeta - 1);
							}
						}
					}
				}
            }
            psi_nm[spe + to_string(lam)] = dot<double>(kernel_nm, this->Vmat[spe + to_string(lam)], false, false);
        }
    }
    pvec.clear();

    unordered_map<string, vec> C{};
    unordered_map<string, int> ispe{};
    int isize = 0;
    for (const string &spe : this->config.species)
    {
        if (this->atom_idx.find(spe) == this->atom_idx.end())
        {
            for (int l = 0; l < this->lmax[spe] + 1; ++l)
            {
                //Check if Vmat[spe + to_string(l)][0] exists
                if (this->Vmat[spe + to_string(l)].size() == 0) {
                    cout << "The projector for species " << spe << " and l = " << l << " does not exist. This is a problem with the model, not NoSpherA2." << endl;
                    cout << "Continuing with the next species..., make sure there is no: " << spe << " in the structure you are trying to predict!!!!" << endl;
                    l = this->lmax[spe] + 1;
                    continue;
                }
               
                //for (int n = 0; n < this->nmax[spe + to_string(l)]; ++n)
                //{
                //    isize += static_cast<int>(this->Vmat[spe + to_string(l)][0].size());
                //}
                isize += static_cast<int>(this->Vmat[spe + to_string(l)][0].size()) * this->nmax[spe + to_string(l)];
            }
            continue;
        }
        ispe[spe] = 0;
        for (int l = 0; l < this->lmax[spe] + 1; ++l)
        {
            for (int n = 0; n < this->nmax[spe + to_string(l)]; ++n)
            {
                int Mcut = static_cast<int>(psi_nm[spe + to_string(l)][0].size());
                //Check if isize + Mcut > weights.size()
                err_chekf(isize + Mcut <= this->weights.size(), "isize + Mcut > weights.size()", std::cout);
                vec weights_subset(this->weights.begin() + isize, this->weights.begin() + isize + Mcut);
                C[spe + to_string(l) + to_string(n)] = dot(psi_nm[spe + to_string(l)], weights_subset);

                isize += Mcut;
            }
        }
    }

    int Tsize = 0;
    for (int iat = 0; iat < this->natoms; iat++)
    {
        string spe = this->atomic_symbols[iat];
        for (int l = 0; l < this->lmax[spe] + 1; l++)
        {
            for (int n = 0; n < this->nmax[spe + to_string(l)]; n++)
            {
                Tsize += 2 * l + 1;
            }
        }
    }

    vec Av_coeffs(Tsize, 0.0);

    // fill vector of predictions
    int i = 0;
    vec pred_coefs(Tsize, 0.0);
    for (int iat = 0; iat < this->natoms; ++iat)
    {
        string spe = this->atomic_symbols[iat];
        for (int l = 0; l < this->lmax[spe] + 1; ++l)
        {
            for (int n = 0; n < this->nmax[spe + to_string(l)]; ++n)
            {
                //for (int ind = 0; ind < 2 * l + 1; ++ind)
                //{
                //    pred_coefs[i + ind] = C[spe + to_string(l) + to_string(n)][ispe[spe] * (2 * l + 1) + ind];
                //}
                std::copy_n(C[spe + to_string(l) + to_string(n)].begin() + ispe[spe] * (2 * l + 1), 2 * l + 1, pred_coefs.begin() + i);

                if (this->config.average && l == 0)
                {
                    Av_coeffs[i] = this->av_coefs[spe][n];
                }
                i += 2 * l + 1;
            }
        }
        ispe[spe] += 1;
    }

    if (this->config.average)
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
    if (this->opt.coef_file != "")
	{
        vec coefs{};
		cout << "Reading coefficients from file: " << this->opt.coef_file << endl;
		read_npy<double>(this->opt.coef_file, coefs);
		return coefs;
	}

    // Run generation of tsc file
    time_point start;
    if (this->opt.debug) start = get_time();

    this->setup_atomic_environment();

    if (!this->config.from_h5) {
        this->read_model_data();
    }
    else {
#if has_RAS
        this->read_model_data_h5();
#else
        err_not_impl_f("HDF5 files are not supported by this build", std::cout);
#endif
        
    }


    vec coefs = this->predict();


    if (this->opt.debug)
    {
        time_point end = get_time();
        long long int dur = get_sec(start, end);
        cout << "Finished ML Density Prediction!" << endl;
        if (dur < 1)
            cout << "Time for SALTED prediction: " << fixed << setprecision(0) << get_msec(start, end) << " ms" << endl;
        else
            cout << "Time for SALTED prediction: " << fixed << setprecision(0) << dur << " s" << endl;

        cout << "Number of coefficients: " << coefs.size() << endl;

        if (this->opt.wfn == string("test_cysteine.xyz"))
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