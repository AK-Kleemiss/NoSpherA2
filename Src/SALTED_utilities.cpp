#include "SALTED_utilities.h"

using namespace std;


vector<cvec2> SALTED_Utils::complex_to_real_transformation(vector<int> sizes)
{
    vector<cvec2> matrices{};
    for (int i = 0; i < sizes.size(); i++)
    {
        int lval = (sizes[i] - 1) / 2;
        int st = (lval & 1) ? 1 : -1;;
        ;
        cvec2 transformed_matrix(sizes[i], vector<complex<double>>(sizes[i], 0.0));
        for (int j = 0; j < lval; j++)
        {
            transformed_matrix[j][j] = complex<double>(0.0, 1.0);
            transformed_matrix[j][sizes[i] - j - 1] = complex<double>(0.0, st);
            transformed_matrix[sizes[i] - j - 1][j] = complex<double>(1.0, 0.0);
            transformed_matrix[sizes[i] - j - 1][sizes[i] - j - 1] = complex<double>(-st, 0.0);
            st = -st;
        }
        transformed_matrix[lval][lval] = sqrt(2.0);
        // Divide each element by sqrt(2.0)
        for (auto& row : transformed_matrix)
        {
            for (auto& elem : row)
            {
                elem /= sqrt(2.0);
            }
        }
        matrices.push_back(transformed_matrix);
    }
    return matrices;
}


int SALTED_Utils::get_lmax_max(unordered_map<string, int>& lmax)
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

void SALTED_Utils::set_lmax_nmax(unordered_map<string, int>& lmax, unordered_map<string, int>& nmax, array<vector<primitive>, 35> basis_set, vector<string> species)
{
    // lmax = {"C": 5, "H":2,...} with the numbers beeing the maximum angular momentum (type) for the given atom
    // nmax = {C0: 10, C1: 7, ...} with the numbers beeing the maximum number of primitives for the given atom and type

    for (auto& spe : species)
    {
        int atom_index = get_Z_from_label(spe.c_str());
        // get the last element of the basis set for the given atom
        lmax[spe] = basis_set[atom_index].back().type;
        // initialize nmax with symbol + type
        for (int i = 0; i < basis_set[atom_index].back().type + 1; i++)
        {
            nmax[spe + to_string(i)] = 0;
        }
        // count the number of primitives for the given atom and type
        for (int i = 0; i < basis_set[atom_index].size(); i++)
        {
            nmax[spe + to_string(basis_set[atom_index][i].type)] += 1;
        }
    }
}



// Function to filter out atoms that belong to species not available for the model selected
std::vector<std::string> SALTED_Utils::filter_species(const std::vector<std::string>& atomic_symbols, const std::vector<std::string>& species)
{
    std::vector<std::string> filtered_symbols;
    std::set<std::string> excluded_species;

    // Convert species vector to a set for efficient lookup
    std::set<std::string> species_set(species.begin(), species.end());

    // Find all species that are not in the input species set
    for (const auto& symbol : atomic_symbols)
    {
        if (species_set.find(symbol) == species_set.end())
        {
            excluded_species.insert(symbol);
        }
    }

    // Filter out excluded species from atomic_symbols
    for (const auto& symbol : atomic_symbols)
    {
        if (excluded_species.find(symbol) == excluded_species.end())
        {
            filtered_symbols.push_back(symbol);
        }
    }

    return filtered_symbols;
}


#if defined(_WIN32) || defined(__RASCALINE__)
std::string Rascaline_Descriptors::to_json(const HyperParametersDensity& params)
{
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

string Rascaline_Descriptors::gen_parameters()
{
    std::ostringstream oss;
    HyperParametersDensity hyper_parameters_density = {
        this->rcut,                           // cutoff
        this->nrad,                           // max_radial
        this->nang,                           // max_angular
        this->atomic_gaussian_width,          // atomic_gaussian_width
        this->center_atom_weight,             // center_atom_weight
        {"Gto", this->spline_accuracy},       // radial_basis
        {"ShiftedCosine", this->cutoff_width} // cutoff_function
    };
    string json_string = to_json(hyper_parameters_density);
    return json_string;
}


Rascaline_Descriptors::Rascaline_Descriptors(std::string filepath, int nrad, int nang, double atomic_gaussian_width,
    double rcut, int n_atoms, std::vector<std::string> neighspe, std::vector<std::string> species,
    double center_atom_weight, double spline_accuracy, double cutoff_width) {
    this->filepath = filepath;
    this->nrad = nrad;
    this->nang = nang;
    this->atomic_gaussian_width = atomic_gaussian_width;
    this->rcut = rcut;
    this->n_atoms = n_atoms;
    this->neighspe = neighspe;
    this->species = species;
    this->center_atom_weight = center_atom_weight;
    this->spline_accuracy = spline_accuracy;
    this->cutoff_width = cutoff_width;
    this->nspe = neighspe.size();
}

// RASCALINE1
metatensor::TensorMap Rascaline_Descriptors::get_feats_projs()
{
    rascaline::BasicSystems system = rascaline::BasicSystems(this->filepath);
    // Construct the parameters for the calculator from the inputs given
    string temp_p = gen_parameters();
    const char* parameters = temp_p.c_str();

    // size_t nspe1 = neighspe.size();
    std::vector<std::vector<int32_t>> keys_array;
    for (int l = 0; l < this->nang + 1; ++l)
    {
        for (const std::string& specen : this->species)
        {
            for (const std::string& speneigh : this->neighspe)
            {
                // Directly emplace back initializer_lists into keys_array
                keys_array.emplace_back(std::vector<int32_t>{l, 1, get_Z_from_label(specen.c_str()) + 1, get_Z_from_label(speneigh.c_str()) + 1});
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
    // descriptor.save("spx_pred.npy");

    return descriptor;
}

// RASCALINE2
cvec4 Rascaline_Descriptors::get_expansion_coeffs(vector<uint8_t> descriptor_buffer)
{
    metatensor::TensorMap descriptor = metatensor::TensorMap::load_buffer(descriptor_buffer);
    vector<vector<vector<vector<complex<double>>>>> omega(this->nang + 1, vector<vector<vector<complex<double>>>>(this->n_atoms, vector<vector<complex<double>>>(2 * this->nang + 1, vector<complex<double>>(this->nspe * this->nrad, { 0.0, 0.0 }))));
    for (int l = 0; l < nang + 1; ++l)
    {
        cvec2 c2r = SALTED_Utils::complex_to_real_transformation({ (2 * l) + 1 })[0];
        metatensor::TensorBlock descriptor_block = descriptor.block_by_id(l);
        metatensor::NDArray<double> descriptor_values = descriptor_block.values();

        // Perform the matrix multiplication and assignment
        for (int a = 0; a < this->n_atoms; ++a)
        {
            for (int c = 0; c < 2 * l + 1; ++c)
            {
                for (int d = 0; d < this->nspe * this->nrad; ++d)
                {
                    omega[l][a][c][d] = 0.0;
                    for (int r = 0; r < 2 * l + 1; ++r)
                    {
                        omega[l][a][c][d] += conj(c2r[r][c]) * descriptor_values(a, r, d);
                    }
                }
            }
        }
        c2r.clear();
    }

    cvec4 expansion_coeffs(2 * this->nang + 1, std::vector<cvec2>(this->nang + 1, cvec2(this->nspe * this->nrad, cvec(this->n_atoms))));

    for (int a = 0; a < omega.size(); ++a)
    {
        for (int b = 0; b < omega[0].size(); ++b)
        {
            for (int c = 0; c < omega[0][0].size(); ++c)
            {
                for (int d = 0; d < omega[0][0][0].size(); ++d)
                {
                    expansion_coeffs[c][a][d][b] = omega[a][b][c][d];
                }
            }
        }
    }
    return expansion_coeffs;
}

cvec4 Rascaline_Descriptors::calculate_expansion_coeffs()
{
	
    metatensor::TensorMap descriptor = get_feats_projs();
    std::vector<uint8_t> descriptor_buffer = descriptor.save_buffer();
    return get_expansion_coeffs(descriptor_buffer);
}
#endif