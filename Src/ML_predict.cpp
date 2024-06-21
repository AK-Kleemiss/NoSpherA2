#include "ML_predict.h"
using namespace std;

struct Config {
    string saltedname;
    string saltedpath;
    vector<string> species;
    int Menv;
    double zeta;
    double regul;
    int ncut;
    bool sparsify;
    string qmcode;
    string qmbasis;
    string dfbasis;
    bool field;
    int Ntrain;
    double trainfrac;
};

struct RadialBasis {
    std::string type;
    double spline_accuracy;
};

struct CutoffFunction {
    std::string type;
    double width;
};

struct HyperParametersDensity {
    double cutoff;
    int max_radial;
    int max_angular;
    double atomic_gaussian_width;
    double center_atom_weight;
    RadialBasis radial_basis;
    CutoffFunction cutoff_function;
};

// Define types for simplicity
using ComplexMatrix = vector<vector<complex<double>>>;

void equicombsparse(int natoms, int nang1, int nang2, int nrad1, int nrad2,
    vector<vector<vector<vector<complex<double>>>>>& v1,
    vector<vector<vector<vector<complex<double>>>>>& v2,
    int wigdim, const vector<double>& w3j, int llmax,
    const vector<vector<int>>& llvec, int lam,
    const vector<vector<complex<double>>>& c2r,
    int featsize, int nfps, const vector<long>& vfps,
    vector<vector<vector<double>>>& p) {
    p = vector<vector<vector<double>>>(2 * lam + 1, vector<vector<double>>(nfps, vector<double>(natoms, 0.0)));

#pragma omp parallel for default(none) shared(natoms, nang1, nang2, nrad1, nrad2, v1, v2, wigdim, w3j, llmax, llvec, lam, c2r, nfps, vfps, featsize, p)
    for (int iat = 0; iat < natoms; ++iat) {
        double inner = 0.0;
        vector<vector<double>> ptemp(2 * lam + 1, vector<double>(featsize, 0.0));
        int ifeat = 0;
        for (int n1 = 0; n1 < nrad1; ++n1) {
            for (int n2 = 0; n2 < nrad2; ++n2) {
                int iwig = 0;
                for (int il = 0; il < llmax; ++il) {
                    int l1 = llvec[il][0];
                    int l2 = llvec[il][1];
                    vector<complex<double>> pcmplx(2 * lam + 1, complex<double>(0.0, 0.0));
                    for (int imu = 0; imu < 2 * lam + 1; ++imu) {
                        int mu = imu - lam;
                        for (int im1 = 0; im1 < 2 * l1 + 1; ++im1) {
                            int m1 = im1 - l1;
                            int m2 = m1 - mu;
                            if (abs(m2) <= l2) {
                                int im2 = m2 + l2;
                                pcmplx[imu] += w3j[iwig] * v1[im1][l1][n1][iat] * conj(v2[im2][l2][n2][iat]);
                                iwig++;
                            }
                        }
                    }
                    vector<double> preal(2 * lam + 1);
                    for (int imu = 0; imu < 2 * lam + 1; ++imu) {
                        for (int j = 0; j < 2 * lam + 1; ++j) {
                            preal[imu] += real(c2r[imu][j] * pcmplx[j]);
                        }
                        inner += pow(preal[imu], 2);
                        ptemp[imu][ifeat] = preal[imu];
                    }
                    ifeat++;
                }
            }
        }
        double normfact = sqrt(inner);
        for (int n = 0; n < nfps; ++n) {
            int ifps = vfps[n];
            for (int imu = 0; imu < 2 * lam + 1; ++imu) {
                p[imu][n][iat] = ptemp[imu][ifps] / normfact;
            }
        }
    }
}


void equicomb(int natoms, int nang1, int nang2, int nrad1, int nrad2,
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& v1,
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>& v2,
    int wigdim, std::vector<double>& w3j, int llmax,
    std::vector<std::vector<int>>& llvec, int lam,
    std::vector<std::vector<std::complex<double>>>& c2r, int featsize,
    std::vector<std::vector<std::vector<double>>>& p) {
    // Initialize p with zeros
    p = std::vector<std::vector<std::vector<double>>>(2 * lam + 1, std::vector<std::vector<double>>(featsize, std::vector<double>(natoms, 0.0)));

    // Parallel region
#pragma omp parallel for default(none) shared(natoms, nang1, nang2, nrad1, nrad2, v1, v2, wigdim, w3j, llmax, llvec, lam, c2r, featsize, p)
    for (int iat = 0; iat < natoms; ++iat) {
        double inner = 0.0;
        std::vector<std::vector<double>> ptemp(2 * lam + 1, std::vector<double>(featsize, 0.0));
        int ifeat = 0;
        for (int n1 = 0; n1 < nrad1; ++n1) {
            for (int n2 = 0; n2 < nrad2; ++n2) {
                int iwig = 0;
                for (int il = 0; il < llmax; ++il) {
                    int l1 = llvec[0][il];
                    int l2 = llvec[1][il];
                    std::vector<std::complex<double>> pcmplx(2 * lam + 1, std::complex<double>(0.0, 0.0));
                    for (int imu = 0; imu < 2 * lam + 1; ++imu) {
                        int mu = imu - lam;
                        for (int im1 = 0; im1 < 2 * l1 + 1; ++im1) {
                            int m1 = im1 - l1;
                            int m2 = m1 - mu;
                            if (std::abs(m2) <= l2) {
                                int im2 = m2 + l2;
                                pcmplx[imu] += w3j[iwig] * v1[im1][l1][n1][iat] * std::conj(v2[im2][l2][n2][iat]);
                                iwig++;
                            }
                        }
                    }
                    std::vector<double> preal(2 * lam + 1);
                    // Matrix multiplication c2r * pcmplx
                    for (int i = 0; i < 2 * lam + 1; ++i) {
                        std::complex<double> sum = 0.0;
                        for (int j = 0; j < 2 * lam + 1; ++j) {
                            sum += c2r[i][j] * pcmplx[j];
                        }
                        preal[i] = std::real(sum);
                        inner += std::norm(sum);
                        ptemp[i][ifeat] = preal[i];
                    }
                    ifeat++;
                }
            }
        }
        double normfact = std::sqrt(inner);
        for (int ifeat = 0; ifeat < featsize; ++ifeat) {
            for (int imu = 0; imu < 2 * lam + 1; ++imu) {
                p[imu][ifeat][iat] = ptemp[imu][ifeat] / normfact;
            }
        }
    }
}


vector<ComplexMatrix> complex_to_real_transformation(vector<int> sizes) {
    vector<ComplexMatrix> matrices{};
    for (int i = 0; i < sizes.size(); i++) {
        int lval = (sizes[i] - 1) / 2;
        int st = pow(-1.0, (lval + 1));
        ComplexMatrix transformed_matrix(sizes[i], vector<complex<double>>(sizes[i], 0.0));
        for (int j = 0; j < lval; j++) {
            transformed_matrix[j][j] = complex<double>(0.0, 1.0);
            transformed_matrix[j][sizes[i] - j - 1] = complex<double>(0.0, st);
            transformed_matrix[sizes[i] - j - 1][j] = complex<double>(1.0, 0.0);
            transformed_matrix[sizes[i] - j - 1][sizes[i] - j - 1] = complex<double>(st * -1.0, 0.0);
            st = st * -1.0;
        }
        transformed_matrix[lval][lval] = sqrt(2.0);
        // Divide each element by sqrt(2.0)
        for (auto& row : transformed_matrix) {
            for (auto& elem : row) {
                elem /= sqrt(2.0);
            }
        }
        matrices.push_back(transformed_matrix);
    }
    return matrices;

}

Config parse_input() {
    // Implement this function to parse the input configuration
    Config inp;
    // Dummy implementation
    return inp;
}

std::string to_json(const HyperParametersDensity& params) {
    std::ostringstream oss;
    oss << "{\n"
        << "  \"cutoff\": " << params.cutoff << ",\n"
        << "  \"max_radial\": " << params.max_radial << ",\n"
        << "  \"max_angular\": " << params.max_angular << ",\n"
        << "  \"atomic_gaussian_width\": " << params.atomic_gaussian_width << ",\n"
        << "  \"center_atom_weight\": " << params.center_atom_weight << ",\n"
        << "  \"radial_basis\": {\n"
        << "    \"Gto\": {\n"
        << "      \"spline_accuracy\": " << params.radial_basis.spline_accuracy << "\n"
        << "    }\n"
        << "  },\n"
        << "  \"cutoff_function\": {\n"
        << "    \"ShiftedCosine\": {\n"
        << "      \"width\": " << params.cutoff_function.width << "\n"
        << "    }\n"
        << "  }\n"
        << "}";
    return oss.str();
}

string gen_parameters(double rcut, int nrad, int nang, double atomic_gaussian_width = 0.3, double center_atom_weight = 1.0, double spline_accuracy = 1e-6, double cutoff_width = 0.1) {
    std::ostringstream oss;
    HyperParametersDensity hyper_parameters_density = {
            rcut,                        // cutoff
            nrad,                       // max_radial
            nang,                       // max_angular
            atomic_gaussian_width,      // atomic_gaussian_width
            center_atom_weight,         // center_atom_weight
            {"Gto", spline_accuracy},              // radial_basis
            {"ShiftedCosine", cutoff_width}      // cutoff_function
    };
    string json_string = to_json(hyper_parameters_density);
    std::cout << json_string << std::endl;

    return json_string;
}

metatensor::TensorMap get_feats_projs(string filepath, double rcut1, int nrad1, int nang1, double atomic_gaussian_width, std::vector<std::string> neighspe, std::vector<std::string> species) {
    rascaline::BasicSystems system = rascaline::BasicSystems(filepath);
    //Construct the parameters for the calculator from the inputs given
    string temp_p= gen_parameters(rcut1, nrad1, nang1, atomic_gaussian_width);
    const char* parameters = temp_p.c_str();

    int nspe1 = neighspe.size();
	std::vector<std::vector<int32_t>> keys_array; 

	for (int l = 0; l <= nang1 + 1; ++l) {
		for (const std::string& specen : species) {
			for (const std::string& speneigh : neighspe) {
				// Directly emplace back initializer_lists into keys_array
				keys_array.emplace_back(std::vector<int32_t>{l, 1, get_Z_from_label(specen.c_str()) + 1, get_Z_from_label(speneigh.c_str()) + 1});
			}
		}
	}

	// Assuming metatensor::Labels expects a flat sequence of integers for each label
	std::vector<int32_t> flattened_keys;
	for (const auto& subVector : keys_array) {
		flattened_keys.insert(flattened_keys.end(), subVector.begin(), subVector.end());
	}

	// Convert keys_array to rascaline::Labels
	std::vector<std::string> names = { "o3_lambda", "o3_sigma", "center_type", "neighbor_type" };
	metatensor::Labels keys_selection(names, flattened_keys.data(), flattened_keys.size() / names.size());

	// create the calculator with its name and parameters
	rascaline::Calculator calculator = rascaline::Calculator("spherical_expansion", parameters);

	rascaline::CalculationOptions calc_opts;
	calc_opts.selected_keys = keys_selection;
	// run the calculation
	metatensor::TensorMap descriptor = calculator.compute(system, calc_opts);

	// The descriptor is a metatensor `TensorMap`, containing multiple blocks.
	// We can transform it to a single block containing a dense representation,
	// with one sample for each atom-centered environment.
	descriptor = descriptor.keys_to_samples("center_type");
	descriptor = descriptor.keys_to_properties("neighbor_type");

	return descriptor;
}




void predict() {

    WFN wavy(0);
    wavy.read_xyz("sucrose.xyz", cout, false);

    //TODO: sanitize input to only inclue atoms that are included in the model
    vector<string> atomic_symbols{};
    for (int i = 0; i < wavy.atoms.size(); i++) {
        atomic_symbols.push_back(wavy.atoms[i].label);
    }
    int n_atoms = atomic_symbols.size();


    double rcut1 = 4.0;
    int nrad1 = 5;
    int nang1 = 7;
    int nspe1 = 3;
    double sig1 = 0.3;
    std::vector<std::string> neighspe1 = { "H", "C", "O" };
    std::vector<std::string> species = { "H", "C", "O" };
    metatensor::TensorMap spx = get_feats_projs("sucrose.xyz", rcut1, nrad1, nang1, sig1, neighspe1, species);
    vector<vector<vector<vector<complex<double>>>>> omega2(nang1 + 1, vector<vector<vector<complex<double>>>>(n_atoms, vector<vector<complex<double>>>(2 * nang1 + 1, vector<complex<double>>(nspe1 * nrad1, { 0.0, 0.0 }))));
    for (int l = 0; l <= nang1 + 1; ++l) {
        ComplexMatrix c2r = complex_to_real_transformation({ (2 * l) + 1 })[0];
        metatensor::TensorBlock spx_block = spx.block_by_id(l);
        auto test = spx.blocks_matching(metatensor::Labels({"o3_lmabda"}, static_cast<const int32_t*>(nullptr), 0));
        metatensor::NDArray<double> spx_values = spx_block.values();
        // Perform the einsum operation
        for (int a = 0; a < n_atoms; ++a) {
            for (int r = 0; r < nspe1 * nrad1; ++r) {
                for (int d = 0; d < 2 * l + 1; ++d) {
                    for (int c = 0; c < 2 * l + 1; ++c) {
                        for (int k = 0; k < 2 * l + 1; ++k) {
                            complex<double> g = omega2[l][a][d][r];
                            complex<double> t1 = conj(c2r[c][k]);
                            auto t2 = spx_values(a, r, k);
                            //omega2[l][a][d][r] += conj(c2r[c][k]) * spx_values(a, r, k, d);
                        }
                    }
                }
            }
        }
        c2r.clear();
	}
    
}