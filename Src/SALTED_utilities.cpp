#include "pch.h"
#include "SALTED_utilities.h"
#include "constants.h"
#include "atoms.h"
#include "cube.h"
#include "wfn_class.h"
#include "metatensor.hpp"
#include "featomic.hpp"

std::vector<cvec2> SALTED_Utils::complex_to_real_transformation(std::vector<int> sizes)
{
	const double sqrt_2 = sqrt(2.0);
    using namespace std;
    vector<cvec2> matrices{};
    for (int i = 0; i < sizes.size(); i++)
    {
        int lval = (sizes[i] - 1) / 2;
        int st = (lval & 1) ? 1 : -1;

        cvec2 transformed_matrix(sizes[i], cvec(sizes[i], 0.0));
        for (int j = 0; j < lval; j++)
        {
            transformed_matrix[j][j] = complex<double>(0.0, 1.0);
            transformed_matrix[j][sizes[i] - j - 1] = complex<double>(0.0, st);
            transformed_matrix[sizes[i] - j - 1][j] = complex<double>(1.0, 0.0);
            transformed_matrix[sizes[i] - j - 1][sizes[i] - j - 1] = complex<double>(-st, 0.0);
            st = -st;
        }
        transformed_matrix[lval][lval] = sqrt_2;
        // Divide each element by sqrt(2.0)
        for (auto &row : transformed_matrix)
        {
            for (auto &elem : row)
            {
                elem /= sqrt_2;
            }
        }
        matrices.push_back(transformed_matrix);
    }
    return matrices;
}

int SALTED_Utils::get_lmax_max(std::unordered_map<std::string, int> &lmax)
{
    int lmax_max = 0;
    for (auto &[key, value] : lmax)
    {
        if (value > lmax_max)
        {
            lmax_max = value;
        }
    }
    return lmax_max;
}

void SALTED_Utils::set_lmax_nmax(std::unordered_map<std::string, int> &lmax, std::unordered_map<std::string, int> &nmax, const std::array<std::vector<primitive>, 118> &basis_set, std::vector<std::string> species)
{
    // lmax = {"C": 5, "H":2,...} with the numbers beeing the maximum angular momentum (type) for the given atom
    // nmax = {C0: 10, C1: 7, ...} with the numbers beeing the maximum number of primitives for the given atom and type

    for (auto &spe : species)
    {
        int atom_index = constants::get_Z_from_label(spe.c_str());
        // get the last element of the basis set for the given atom
        lmax[spe] = basis_set[atom_index].back().get_type();
        // initialize nmax with symbol + type
        for (int i = 0; i < basis_set[atom_index].back().get_type() + 1; i++)
        {
            nmax[spe + std::to_string(i)] = 0;
        }
        // count the number of primitives for the given atom and type
        for (int i = 0; i < basis_set[atom_index].size(); i++)
        {
            nmax[spe + std::to_string(basis_set[atom_index][i].get_type())] += 1;
        }
    }
}

// Function to filter out atoms that belong to species not available for the model selected
std::vector<std::string> SALTED_Utils::filter_species(const std::vector<std::string> &atomic_symbols, const std::vector<std::string> &species)
{
    std::vector<std::string> filtered_symbols;
    std::set<std::string> excluded_species;

    // Convert species vector to a set for efficient lookup
    std::set<std::string> species_set(species.begin(), species.end());

    // Find all species that are not in the input species set
    for (const auto &symbol : atomic_symbols)
    {
        if (species_set.find(symbol) == species_set.end())
        {
            excluded_species.insert(symbol);
        }
    }

    // Print out the excluded species
    if (!excluded_species.empty())
    {
        std::cout << "Excluded species: ";
        for (const auto &_species : excluded_species)
        {
            std::cout << _species << " ";
        }
        std::cout << std::endl;
        err_not_impl_f("This Model does not contain all neccecary molecules to predict this structure\n", std::cout);
    }

    // Filter out excluded species from atomic_symbols
    for (const auto &symbol : atomic_symbols)
    {
        if (excluded_species.find(symbol) == excluded_species.end())
        {
            filtered_symbols.push_back(symbol);
        }
    }

    return filtered_symbols;
}

std::string Rascaline_Descriptors::to_json(const HyperParametersDensity &params)
{
    std::ostringstream oss;
    oss << "{\n"
        << "  \"cutoff\": {\n    \"radius\": " << params.cutoff << ",\n    \"smoothing\": {\"type\": \"ShiftedCosine\", \"width\": " << params.cutoff_function.width << " }\n  }, \n"
        << "  \"density\": {\n    \"type\": \"Gaussian\",\n    \"width\": "<< params.atomic_gaussian_width << ",\n    \"center_atom_weight\": "<< params.center_atom_weight << "\n  },\n"
        << "  \"basis\": {\n    \"type\": \"TensorProduct\",\n    \"max_angular\": " << params.max_angular << ",\n    \"radial\": {\"type\": \"Gto\", \"max_radial\": " << params.max_radial << "}\n  }\n"
        << "}";
    return oss.str();
}

std::string Rascaline_Descriptors::gen_parameters()
{
    HyperParametersDensity hyper_parameters_density = {
        this->rcut,                           // cutoff
        this->nrad,                           // max_radial
        this->nang,                           // max_angular
        this->atomic_gaussian_width,          // atomic_gaussian_width
        this->center_atom_weight,             // center_atom_weight
        {"Gto", this->spline_accuracy},       // radial_basis
        {"ShiftedCosine", this->cutoff_width} // cutoff_function
    };
    std::string json_string = to_json(hyper_parameters_density);
    return json_string;
}

// How about this?
Rascaline_Descriptors::Rascaline_Descriptors(const std::filesystem::path &filepath, const int &nrad, const int &nang, const double &atomic_gaussian_width,
                                             const double &rcut, const int &n_atoms, const std::vector<std::string> &neighspe, const std::vector<std::string> &species,
                                             const double &center_atom_weight, const double &spline_accuracy, const double &cutoff_width)
{
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
    this->nspe = (int)neighspe.size();
}

// RASCALINE1
metatensor::TensorMap Rascaline_Descriptors::get_feats_projs()
{
    featomic::SimpleSystem system;
    WFN wfn = WFN(this->filepath.c_str());
    for (const atom& a : *wfn.get_atoms_ptr())
    {
        system.add_atom(a.get_charge(), a.get_coords());
    }
    // Construct the parameters for the calculator from the inputs given
    std::string temp_p = gen_parameters();
    const char *parameters = temp_p.c_str();

    // size_t nspe1 = neighspe.size();
    std::vector<std::vector<int32_t>> keys_array;
    for (int l = 0; l < this->nang + 1; ++l)
    {
        for (const std::string &specen : this->species)
        {
            for (const std::string &speneigh : this->neighspe)
            {
                // Directly emplace back initializer_lists into keys_array
                keys_array.emplace_back(std::vector<int32_t>{l, 1, constants::get_Z_from_label(specen.c_str()) + 1, constants::get_Z_from_label(speneigh.c_str()) + 1});
            }
        }
    }

    // Assuming metatensor::Labels expects a flat sequence of integers for each label
    std::vector<int32_t> flattened_keys;
    for (const auto &subVector : keys_array)
    {
        flattened_keys.insert(flattened_keys.end(), subVector.begin(), subVector.end());
    }

    // Convert keys_array to rascaline::Labels
    std::vector<std::string> names = {"o3_lambda", "o3_sigma", "center_type", "neighbor_type"};
    metatensor::Labels keys_selection(names, flattened_keys.data(), flattened_keys.size() / names.size());

    // create the calculator with its name and parameters
    auto calculator = featomic::Calculator("spherical_expansion", parameters);

    featomic::CalculationOptions calc_opts;
    calc_opts.selected_keys = keys_selection;
    calc_opts.use_native_system = true;
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
cvec4 Rascaline_Descriptors::get_expansion_coeffs(std::vector<uint8_t> descriptor_buffer)
{

    metatensor::TensorMap descriptor = metatensor::TensorMap::load_buffer(descriptor_buffer);
    // cvec4 omega(this->nang + 1, std::vector<cvec2>(this->n_atoms, cvec2(2 * this->nang + 1, cvec(this->nspe * this->nrad, {0.0, 0.0}))));
    cvec4 omega(this->n_atoms, std::vector<cvec2>(this->nspe * this->nrad, cvec2(this->nang + 1, cvec(2 * this->nang + 1, {0.0, 0.0}))));
    for (int l = 0; l < nang + 1; ++l)
    {
        cvec2 c2r = SALTED_Utils::complex_to_real_transformation({(2 * l) + 1})[0];
        metatensor::TensorBlock descriptor_block = descriptor.block_by_id(l);
        metatensor::NDArray<double> descriptor_values = descriptor_block.values();

        // Perform the matrix multiplication and assignment
        for (int a = 0; a < this->n_atoms; ++a)
        {
            for (int c = 0; c < 2 * l + 1; ++c)
            {
                for (int d = 0; d < this->nspe * this->nrad; ++d)
                {
                    // omega[l][a][c][d] = 0.0;
                    omega[a][d][l][c] = 0.0;
                    for (int r = 0; r < 2 * l + 1; ++r)
                    {
                        // omega[l][a][c][d] += conj(c2r[r][c]) * descriptor_values(a, r, d);
                        omega[a][d][l][c] += conj(c2r[r][c]) * descriptor_values(a, r, d);
                    }
                } /*cvec* v2_ptr = (cvec*)&v2[iat][l2][n2];*/
            }
        }
        c2r.clear();
    }

    return omega;
}

cvec4 Rascaline_Descriptors::calculate_expansion_coeffs()
{

    metatensor::TensorMap descriptor = get_feats_projs();
    std::vector<uint8_t> descriptor_buffer = descriptor.save_buffer();
    return get_expansion_coeffs(descriptor_buffer);
}


const double calc_density_ML(const double& x,
    const double& y,
    const double& z,
    const vec& coefficients,
    const std::vector<atom>& atoms)
{
    double dens = 0, radial;
    int coef_counter = 0;
    int e = 0, size = 0;
    basis_set_entry bf;
    primitive p;

    for (int a = 0; a < atoms.size(); a++)
    {
        size = (int)atoms[a].get_basis_set_size();
        double d[4]{
            x - atoms[a].get_coordinate(0),
            y - atoms[a].get_coordinate(1),
            z - atoms[a].get_coordinate(2), 0.0 };
        // store r in last element
        d[3] = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (d[3] < -46.0517)
        { // corresponds to cutoff of ex ~< 1E-20
            for (e = 0; e < size; e++)
            {
                bf = atoms[a].get_basis_set_entry(e);
                coef_counter += (2 * bf.get_type() + 1);
            }
            continue;
        }
        // normalize distances for spherical harmonic
        for (e = 0; e < 3; e++)
            d[e] /= d[3];
        int type = -1, prim = 0;
        for (int shell = 0; shell < atoms[a].get_shellcount().size(); shell++) {
            radial = 0;
            int type = atoms[a].get_basis_set_entry(prim).get_type();

            for (int e = 0; e < atoms[a].get_shellcount()[shell]; e++, prim++) {
                bf = atoms[a].get_basis_set_entry(prim);
                radial += gaussian_radial(bf.get_primitive(), d[3]) * bf.get_coefficient();
            }

            if (radial < 1E-10)
            {
                coef_counter += (2 * type + 1);
                continue;
            }

            dens += radial * constants::spherical_harmonic(type, d, &coefficients[coef_counter]);
            coef_counter += (2 * type + 1);
        }
    }
    // err_checkf(coef_counter == exp_coefs, "WRONG NUMBER OF COEFFICIENTS! " + std::to_string(coef_counter) + " vs. " + std::to_string(exp_coefs), std::cout);
    return dens;
}

const double calc_density_ML(const double& x,
    const double& y,
    const double& z,
    const vec& coefficients,
    const std::vector<atom>& atoms,
    const int& atom_nr)
{
    double dens = 0, radial = 0;
    int coef_counter = 0;
    int e = 0, size = 0;

    for (int a = 0; a < atoms.size(); a++)
    {
        if (a != atom_nr) {
            for (e = 0; e < size; e++)
            {
                coef_counter += (2 * atoms[a].get_basis_set_type(e) + 1);
            }
            continue;
        }
        size = (int)atoms[a].get_basis_set_size();

        basis_set_entry bf;
        double d[4]{
            x - atoms[a].get_coordinate(0),
            y - atoms[a].get_coordinate(1),
            z - atoms[a].get_coordinate(2), 0.0 };
        // store r in last element
        d[3] = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        // normalize distances for spherical harmonic
        for (e = 0; e < 3; e++)
            d[e] /= d[3];
        int type = -1, prim = 0;
        for (int shell = 0; shell < atoms[a].get_shellcount().size(); shell++) {
            radial = 0;
            int type = atoms[a].get_basis_set_entry(prim).get_type();

            for (int e = 0; e < atoms[a].get_shellcount()[shell]; e++, prim++) {
                bf = atoms[a].get_basis_set_entry(prim);
                radial += gaussian_radial(bf.get_primitive(), d[3]) * bf.get_coefficient();
            }

            if (radial < 1E-10)
            {
                coef_counter += (2 * type + 1);
                continue;
            }

            dens += radial * constants::spherical_harmonic(type, d, &coefficients[coef_counter]);
            coef_counter += (2 * type + 1);
        }
        return dens;
    }
}



double helper_even(const int l) {
    if (l == 0) {
        return 1.0;
    }
    double res = 1.0;
    for (int i = 0; i <= l; i += 2) {
        res *= i;
    }
    return res;
}

/**
 * Calculates the atomic density for a given list of atoms and coefficients.
 *
 * @param atoms The list of atoms.
 * @param coefs The coefficients used in the calculation.
 * @return The atomic density for each atom.
 */
vec calc_atomic_density(const std::vector<atom> &atoms, const vec &coefs)
{
    int e = 0, size;
    double radial;
    basis_set_entry bf;

    vec atom_elecs(atoms.size(), 0.0);

    int coef_counter = 0;
    for (int a = 0; a < atoms.size(); a++)
    {

        int type = -1, prim = 0;
        for (int shell = 0; shell < atoms[a].get_shellcount().size(); shell++) {
            radial = 0;
            type = atoms[a].get_basis_set_entry(prim).get_type();
            if (type != 0)
            {
                coef_counter += (2 * type + 1); prim += atoms[a].get_shellcount()[shell]; //Skip functions and coefficients
                continue;
            }

            for (int e = 0; e < atoms[a].get_shellcount()[shell]; e++, prim++) {
                bf = atoms[a].get_basis_set_entry(prim);
                primitive p(a, bf.get_type(), bf.get_exponent(), bf.get_coefficient());
                radial += constants::PI / (2.0 * std::pow(p.get_exp(), 1.5)) * p.normalization_constant() * p.get_coef();
            }

            atom_elecs[a] += radial * coefs[coef_counter];
            coef_counter++;
        }
    }
    return atom_elecs;
}

cube calc_cube_ML(const vec data, WFN &dummy, const int atom_nr)
{
    double MinMax[6]{0, 0, 0, 0, 0, 0};
    int steps[3]{0, 0, 0};
    readxyzMinMax_fromWFN(dummy, MinMax, steps, 2.5, 0.1, true);
    cube CubeRho(steps[0], steps[1], steps[2], dummy.get_ncen(), true);
    CubeRho.give_parent_wfn(dummy);

    for (int i = 0; i < 3; i++)
    {
        CubeRho.set_origin(i, MinMax[i]);
        CubeRho.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
    }
    CubeRho.set_comment1("Calculated density using NoSpherA2 from ML Data");
    CubeRho.set_comment2("from " + dummy.get_path().string());
    CubeRho.set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_RI_rho.cube");

    _time_point start = get_time();

    ProgressBar *progress = new ProgressBar(CubeRho.get_size(0), 60, "=", " ", "Calculating Values");
    if (atom_nr != -1)
        std::cout << "Calculation for atom " << atom_nr << std::endl;

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < CubeRho.get_size(0); i++)
    {
        for (int j = 0; j < CubeRho.get_size(1); j++)
            for (int k = 0; k < CubeRho.get_size(2); k++)
            {

                vec PosGrid{
                    i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0),
                    i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1),
                    i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2)};

                if (atom_nr == -1)
                    CubeRho.set_value(i, j, k, calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.get_atoms()));
                else
                    CubeRho.set_value(i, j, k, calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.get_atoms(), atom_nr));
            }
        progress->update();
    }
    delete (progress);

    using namespace std;
    _time_point end = get_time();
    if (get_sec(start, end) < 60)
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
    CubeRho.calc_dv();
    std::cout << "Number of electrons: " << std::fixed << std::setprecision(4) << CubeRho.sum() << std::endl;
    return CubeRho;
};

void calc_cube_ML(const vec data, WFN& dummy, const int atom_nr, cube &CubeRho)
{
    CubeRho.set_comment1("Calculated density using NoSpherA2 from ML Data");
    CubeRho.set_comment2("from " + dummy.get_path().string());
    CubeRho.set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_RI_rho.cube");
    _time_point start = get_time();

    ProgressBar* progress = new ProgressBar(CubeRho.get_size(0), 60, "=", " ", "Calculating Values");
    if (atom_nr != -1)
        std::cout << "Calculation for atom " << atom_nr << std::endl;
    double dens = 0.0;
#pragma omp parallel for schedule(dynamic) private(dens)
    for (int i = 0; i < CubeRho.get_size(0); i++)
    {
        for (int j = 0; j < CubeRho.get_size(1); j++)
            for (int k = 0; k < CubeRho.get_size(2); k++)
            {

                vec PosGrid{
                    i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0),
                    i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1),
                    i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2) };

                dens = (atom_nr == -1)
                    ? calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.get_atoms())
                    : calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.get_atoms(), atom_nr);

                CubeRho.set_value(i, j, k, dens);

                //if (atom_nr == -1)
                //    CubeRho.set_value(i, j, k, calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.get_atoms()));
                //else
                //    CubeRho.set_value(i, j, k, calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.get_atoms(), atom_nr));
            }
        progress->update();
    }
    delete (progress);

    using namespace std;
    _time_point end = get_time();
    if (get_sec(start, end) < 60)
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
    CubeRho.calc_dv();
    std::cout << "Number of electrons: " << std::fixed << std::setprecision(4) << CubeRho.sum() << std::endl;
};