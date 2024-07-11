/*
Original Code by:
    Authors:    Andrea Grisafi, Alberto Fabrizio, David M. Wilkins, Benjamin A. R. Meyer, Clemence Corminboeuf, Michele Ceriotti, Zekun Lou, Alan Lewis
    Date:       2024
    License:    GPL-3.0
    GitHub:     https://github.com/andreagrisafi/SALTED
    Article:    https://pubs.acs.org/doi/10.1021/acscentsci.8b00551
*/

#include "ML_predict.h"

using namespace std;
using namespace SALTED_Utils;

vec predict(const WFN& wavy, const string model_folder)
{
    // Read the configuration file
    cout << "Starting prediction... " << flush;
    Config config;
    string _f_path("inputs.txt");
    string _path = model_folder;
    join_path(_path, _f_path);
    config.populateFromFile(_path);

    config.predict_filename = wavy.get_path();
    int nspe1 = static_cast<int>(config.neighspe1.size());
    int nspe2 = static_cast<int>(config.neighspe2.size());
    /* bool average = true;
     bool sparsify = true;
     int ncut = 1000;
     bool field = false;
     double rcut1 = 4.0;
     int nang1 = 7;
     int nang2 = 7;
     int nrad1 = 7;
     int nrad2 = 7;
     double sig1 = 0.3;

     float zeta = 2.0;
     int Ntrain = 600;
     float trainfrac = 0.8;
     std::vector<std::string> neighspe1 = { "H", "C", "N", "O", "S"};
     std::vector<std::string> species = { "H", "C", "N", "O", "S" };
     int nspe1 = neighspe1.size();
     int nspe2 = nspe2;*/

    vector<string> atomic_symbols{};
    for (int i = 0; i < wavy.atoms.size(); i++)
    {
        atomic_symbols.push_back(wavy.atoms[i].label);
    }
    //# Define system excluding atoms that belong to species not listed in SALTED input
    atomic_symbols = filter_species(atomic_symbols, config.species);
    //// Output the result
    // for (const auto& symbol : atomic_symbols) {
    //     std::cout << symbol << " ";
    // }
    // std::cout << std::endl;

    int natoms = static_cast<int>(atomic_symbols.size());
    unordered_map<string, vector<int>> atom_idx{};
    unordered_map<string, int> natom_dict{};
    for (int i = 0; i < atomic_symbols.size(); i++)
    {
        atom_idx[atomic_symbols[i]].push_back(i);
        natom_dict[atomic_symbols[i]] += 1;
    }

    unordered_map<string, int> lmax{};
    unordered_map<string, int> nmax{};
    set_lmax_nmax(lmax, nmax, QZVP_JKfit, config.species);
    


    // RASCALINE (Generate descriptors)
#if defined(_WIN32) || defined(__RASCALINE__)
    cvec4 v1 = Rascaline_Descriptors(config.predict_filename, config.nrad1, config.nang1, config.sig1, config.rcut1, natoms, config.neighspe1, config.species).calculate_expansion_coeffs();
    cvec4 v2 = Rascaline_Descriptors(config.predict_filename, config.nrad2, config.nang2, config.sig2, config.rcut2, natoms, config.neighspe2, config.species).calculate_expansion_coeffs();
#else
    err_not_impl_f("RASCALINE is not supported by this build", std::cout);
    cvec4 v2, v1;
#endif
    // END RASCALINE
    // 
    // Read Model variables
    unordered_map<string, vec2> Vmat{};
    unordered_map<string, int> Mspe{};
    unordered_map<string, vec2> power_env_sparse{};
    // Define zeta as a string with one decimal
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(1) << config.zeta;
    std::string zeta_str = stream.str();
#if defined(_WIN32) || defined(__RASCALINE__)
    H5::H5File features(model_folder + "/GPR_data/FEAT_M-" + to_string(config.Menv) + ".h5", H5F_ACC_RDONLY);
    H5::H5File projectors(model_folder + "/GPR_data/projector_M" + to_string(config.Menv) + "_zeta" + zeta_str + ".h5", H5F_ACC_RDONLY);
    for (string spe : config.species)
    {
        for (int lam = 0; lam < lmax[spe] + 1; lam++)
        {
            power_env_sparse[spe + to_string(lam)] = readHDF5(features, "sparse_descriptors/" + spe + "/" + to_string(lam));
            Vmat[spe + to_string(lam)] = readHDF5(projectors, "projectors/" + spe + "/" + to_string(lam));

            if (lam == 0)
            {
                Mspe[spe] = static_cast<int>(power_env_sparse[spe + to_string(lam)].size());
            }
            if (config.zeta == 1)
            {
                vec2 Vmat_t = transpose<double>(Vmat[spe + to_string(lam)]);
                power_env_sparse[spe + to_string(lam)] = dot<double>(Vmat_t, power_env_sparse[spe + to_string(lam)]);
            }
        }
    }
    features.close();
    projectors.close();
#endif

    unordered_map<int, vector<int64_t>> vfps{};
    if (config.sparsify)
    {
        vfps = read_fps<int64_t>(model_folder + "/GPR_data/fps" + to_string(config.ncut) + "-", get_lmax_max(lmax));
    };

    int ntrain = static_cast<int>(config.Ntrain * config.trainfrac);
    vec weights{};
    if (config.field)
    {
        cout << "Field" << endl;
        err_not_impl_f("Calculations using 'Field = True' are not yet supported", std::cout);
    }
    else
    {
        read_npy(model_folder + "/GPR_data/weights_N" + to_string(ntrain) + "_reg-6.npy", weights);
    }
    // END READ MODEL VARIABLES

    // Compute equivariant descriptors for each lambda value entering the SPH expansion of the electron density
    vec2 pvec_l0{};
    unordered_map<int, vec3> pvec{};
    for (int lam = 0; lam < get_lmax_max(lmax) + 1; lam++)
    {
        int llmax = 0;
        unordered_map<int, vector<int>> lvalues{};
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

        vector<double> wigner3j = readVectorFromFile<double>(model_folder + "/wigners/wigner_lam-" + to_string(lam) + "_lmax1-" + to_string(config.nang1) + "_lmax2-" + to_string(config.nang2) + ".dat");

        cvec2 c2r = SALTED_Utils::complex_to_real_transformation({2 * lam + 1})[0];

        int featsize = nspe1 * nspe2 * config.nrad1 * config.nrad2 * llmax;
        vec3 p;
        ivec2 llvec_t = transpose<int>(llvec);
        if (config.sparsify)
        {
            int nfps = static_cast<int>(vfps[lam].size());
            equicomb(natoms, config.nang1, config.nang2, (nspe1 * config.nrad1), (nspe2 * config.nrad2), v1, v2, wigner3j, llmax, llvec_t, lam, c2r, featsize, nfps, vfps[lam], p);
            featsize = config.ncut;
        }
        else
        {
            equicomb(natoms, config.nang1, config.nang2, (nspe1 * config.nrad1), (nspe2 * config.nrad2), v1, v2, wigner3j, llmax, llvec_t, lam, c2r, featsize, p);
        }
        wigner3j.clear();

        vec3 p_t = reorder3D<double>(p);
        vec flat_p = flatten<double>(p_t);

        if (lam == 0)
        {
            pvec_l0 = reshape<double>(flat_p, Shape2D{natoms, featsize});
        }
        else
        {
            pvec[lam] = reshape<double>(flat_p, Shape3D{natoms, 2 * lam + 1, featsize});
        }
        flat_p.clear();
    }

    unordered_map<string, vec2> psi_nm{};

    for (const string &spe : config.species)
    {
        if (atom_idx.find(spe) == atom_idx.end())
        {
            continue;
        }
        vec2 power_env_sparse_t = transpose<double>(power_env_sparse[spe + "0"]);
        vec2 kernell0_nm;
        if (config.zeta == 1)
        {
            psi_nm[spe + "0"] = dot<double>(collectRows(pvec_l0, atom_idx[spe]), power_env_sparse_t);
        }
        else
        {
            kernell0_nm = dot<double>(collectRows(pvec_l0, atom_idx[spe]), power_env_sparse_t);
            vec2 kernel_nm = elementWiseExponentiation(kernell0_nm, config.zeta);
            psi_nm[spe + "0"] = dot<double>(kernel_nm, Vmat[spe + "0"]);
        }

        for (int lam = 1; lam <= lmax[spe]; ++lam)
        {
            int featsize = static_cast<int>(pvec[lam][0][0].size());
            power_env_sparse_t = transpose<double>(power_env_sparse[spe + to_string(lam)]);
            vec3 pVec_Rows = collectRows(pvec[lam], atom_idx[spe]);
            vec2 pvec_lam = reshape<double>(flatten<double>(pVec_Rows), Shape2D{natom_dict[spe] * (2 * lam + 1), featsize});

            if (config.zeta == 1)
            {
                psi_nm[spe + to_string(lam)] = dot<double>(pvec_lam, power_env_sparse_t);
            }
            else
            {
                vec2 kernel_nm = dot<double>(pvec_lam, power_env_sparse_t);
                for (size_t i1 = 0; i1 < natom_dict[spe]; ++i1)
                {
                    for (size_t i2 = 0; i2 < Mspe[spe]; ++i2)
                    {
                        for (size_t i = 0; i < 2 * lam + 1; ++i)
                        {
                            for (size_t j = 0; j < 2 * lam + 1; ++j)
                            {
                                kernel_nm[i1 * (2 * lam + 1) + i][i2 * (2 * lam + 1) + j] *= pow(kernell0_nm[i1][i2], config.zeta - 1);
                            }
                        }
                    }
                }
                psi_nm[spe + to_string(lam)] = dot<double>(kernel_nm, Vmat[spe + to_string(lam)]);
            }
        }
    }
    
    unordered_map<string, vector<double>> C{};
    unordered_map<string, int> ispe{};
    int isize = 0;
    for (const string &spe : config.species)
    {
        if (atom_idx.find(spe) == atom_idx.end())
        {
            for (int l = 0; l <= lmax[spe]; ++l)
            {
                for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
                {
                    isize += static_cast<int>(Vmat[spe + to_string(l)].size());
                }
            }
            continue;
        }
        ispe[spe] = 0;
        for (int l = 0; l <= lmax[spe]; ++l)
        {
            for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
            {
                int Mcut = static_cast<int>(psi_nm[spe + to_string(l)][0].size());
                vector<double> weights_subset(weights.begin() + isize, weights.begin() + isize + Mcut);
                C[spe + to_string(l) + to_string(n)] = dot(psi_nm[spe + to_string(l)], weights_subset);
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

    vector<double> Av_coeffs(Tsize, 0.0);
    unordered_map<string, vector<double>> av_coefs{};
    if (config.average)
    {
        for (string spe : config.species)
        {
            read_npy(model_folder + "/averages/averages_" + spe + ".npy", av_coefs[spe]);
        }
    }

    // fill vector of predictions
    int i = 0;
    vector<double> pred_coefs(Tsize, 0.0);
    for (int iat = 0; iat < natoms; ++iat)
    {
        string spe = atomic_symbols[iat];
        for (int l = 0; l <= lmax[spe]; ++l)
        {
            for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
            {
                for (int ind = 0; ind < 2 * l + 1; ++ind)
                {
                    pred_coefs[i + ind] = C[spe + to_string(l) + to_string(n)][ispe[spe] * (2 * l + 1) + ind];
                }
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

    cout << "          ... done!\nNumber of predicted coefficients: " << pred_coefs.size() << endl;
    // npy::npy_data<double> coeffs;
    // coeffs.data = pred_coefs;
    // coeffs.fortran_order = false;
    // coeffs.shape = { unsigned long(pred_coefs.size()) };
    // npy::write_npy("coeffs_by_black_magic.npy", coeffs);
    return pred_coefs;
}

//Wrapper function to generate the ML density prediction
vec gen_SALTED_densities(const WFN& wave, options opt, time_point& start, time_point& end_SALTED) {
    // Run generation of tsc file
    if (opt.debug) cout << "Finished ML Density Prediction!" << endl;
    vec coefs = predict(wave, opt.SALTED_DIR);
    end_SALTED = get_time();
    if (opt.debug) {
        long long int dur = get_sec(start, end_SALTED);
        cout << "Finished ML Density Prediction!" << endl;
        if (dur < 1)
            cout << "Time for SALTED prediction: " << fixed << setprecision(0) << get_msec(start, end_SALTED) << " ms" << endl;
        else
            cout << "Time for SALTED prediction: " << fixed << setprecision(0) << dur << " s" << endl;
    }
    if (opt.debug && (opt.wfn == string("test_cysteine.xyz") || opt.wfn == string("test_sucrose.xyz"))) {
        vector<unsigned long> shape{};
        bool fortran_order;
        vec ref_coefs{};
        //Depending on the value of opt.wfn read either the cysteine or the sucrose reference coefs
        if (opt.wfn == string("test_sucrose.xyz"))
            npy::LoadArrayFromNumpy("sucrose_ref_Combined_v1.npy", shape, fortran_order, ref_coefs);
        else
            npy::LoadArrayFromNumpy("cysteine_ref_Cysteine.npy", shape, fortran_order, ref_coefs);
        // Compare coefs with the reference
        vector<double> diff_vec;
        double diff = 0.0;
        for (int i = 0; i < coefs.size(); i++)
        {
            diff += abs(coefs[i] - ref_coefs[i]);
            if (abs(coefs[i] - ref_coefs[i]) > 1e-4) {
                cout << "Difference in coef " << fixed << setprecision(3) << i << " : " << coefs[i] << " - " << ref_coefs[i] << " = " << abs(coefs[i] - ref_coefs[i]) << endl;
            }
            diff_vec.push_back((coefs[i] / ref_coefs[i]) - 1);
        }
        cout << "Difference between calculated and reference coefs: " << fixed << setprecision(3) << diff << endl;
        cout << "Maximum ((pred / ref) -1): " << fixed << setprecision(3) << *max_element(diff_vec.begin(), diff_vec.end()) << endl;
        if (diff > 0.1) {
			cout << "WARNING: The difference between the calculated and reference coefficients is too large!" << endl;
            exit(1);
        }
    }

    return coefs;
}
