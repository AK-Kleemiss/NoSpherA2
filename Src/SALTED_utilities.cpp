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
        for (auto& row : transformed_matrix)
        {
            for (auto& elem : row)
            {
                elem /= sqrt_2;
            }
        }
        matrices.push_back(transformed_matrix);
    }
    return matrices;
}

int SALTED_Utils::get_lmax_max(std::unordered_map<std::string, int>& lmax)
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

void SALTED_Utils::set_lmax_nmax(std::unordered_map<std::string, int>& lmax, std::unordered_map<std::string, int>& nmax, const std::array<std::vector<primitive>, 118>& basis_set, std::vector<std::string> species)
{
    // lmax = {"C": 5, "H":2,...} with the numbers beeing the maximum angular momentum (type) for the given atom
    // nmax = {C0: 10, C1: 7, ...} with the numbers beeing the maximum number of primitives for the given atom and type

    for (auto& spe : species)
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

    // Print out the excluded species
    if (!excluded_species.empty())
    {
        std::cout << "Excluded species: ";
        for (const auto& _species : excluded_species)
        {
            std::cout << _species << " ";
        }
        std::cout << std::endl;
        err_not_impl_f("This Model does not contain all neccecary molecules to predict this structure\n", std::cout);
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

std::string SALTED_Utils::FeatomicHyperParameters::to_json() const
{
    std::ostringstream oss;
    oss << "{\n"
        << "  \"cutoff\": {  \
                    \"radius\": " << this->cutoff_radius << " , \"smoothing\": \
                                  {\"type\": \"" << this->cutoff_function.type << "\", \"width\": " << this->cutoff_function.width << "} }, \n"
        << "  \"density\": { \
                    \"type\": \"Gaussian\", \"width\": " << this->atomic_gaussian_width << ", \"center_atom_weight\": " << this->center_atom_weight << "},\n"
        << "  \"basis\": { \
                    \"type\": \"TensorProduct\", \"max_angular\": " << this->max_angular << ", \"radial\": {\"type\":  \"" << this->radial_basis.type << "\", \"max_radial\": " << this->max_radial << "} , \"spline_accuracy\": " << this->radial_basis.spline_accuracy << "}\n"
        << "}";
    return oss.str();
}


// Used to generate metatensor::TensorMap and save the buffer location into the descriptor_buffer
static metatensor::TensorMap get_feats_projs(featomic::SimpleSystem featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters)
{
    // size_t nspe1 = neighspe.size();
    std::vector<std::array<int32_t, 4>> keys_array;
    keys_array.reserve((parameters.max_angular + 1) * parameters.species.size() * parameters.neighspe.size());

    for (int l = 0; l < parameters.max_angular + 1; ++l)
    {
        for (const std::string& center_spe : parameters.species)
        {
            int32_t center_z = constants::get_Z_from_label(center_spe.c_str()) + 1;
            for (const std::string& neigh_spe : parameters.neighspe)
            {
                int32_t neigh_z = constants::get_Z_from_label(neigh_spe.c_str()) + 1;
                // Directly emplace back initializer_lists into keys_array
                keys_array.push_back({ l, 1, center_z, neigh_z});
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


    //create the calculator with its name and parameters
    //Do not ask me, why Featomic expects the max_radial to be one less than the actual number of radial basis functions, but it does, so here we are
    SALTED_Utils::FeatomicHyperParameters modif_param = parameters;
    modif_param.max_radial -= 1;
    auto calculator = featomic::Calculator("spherical_expansion", modif_param.to_json().c_str());

    featomic::CalculationOptions calc_opts;
    calc_opts.selected_keys = keys_selection;
    calc_opts.use_native_system = true;
    // run the calculation
    metatensor::TensorMap descriptor = calculator.compute(featomic_system, calc_opts);

    // The descriptor is a metatensor `TensorMap`, containing multiple blocks.
    // We can transform it to a single block containing a dense representation,
    // with one sample for each atom-centered environment.
    descriptor = descriptor.keys_to_samples("center_type");
    descriptor = descriptor.keys_to_properties("neighbor_type");
    // descriptor.save("spx_pred.npy");

    return descriptor;
}

// Reads the descriptor buffer and fills the expansion coefficients vector
static cvec4 get_expansion_coeffs(std::vector<uint8_t> descriptor_buffer, const featomic::SimpleSystem& featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters)
{
    int n_atoms = (int)featomic_system.size();
    int nspe = (int)parameters.species.size();
    metatensor::TensorMap descriptor = metatensor::TensorMap::load_buffer(descriptor_buffer);
    // cvec4 omega(this->nang + 1, std::vector<cvec2>(this->n_atoms, cvec2(2 * this->nang + 1, cvec(this->nspe * this->nrad, {0.0, 0.0}))));
    cvec4 omega(n_atoms, std::vector<cvec2>((size_t)nspe * parameters.max_radial, cvec2((size_t)parameters.max_angular + 1, cvec((size_t)2 * parameters.max_angular + 1, { 0.0, 0.0 }))));
    for (int l = 0; l < parameters.max_angular + 1; ++l)
    {
        cvec2 c2r = SALTED_Utils::complex_to_real_transformation({ (2 * l) + 1 })[0];
        metatensor::TensorBlock descriptor_block = descriptor.block_by_id(l);
        metatensor::NDArray<double> descriptor_values = descriptor_block.values();

        // Perform the matrix multiplication and assignment
        for (int a = 0; a < n_atoms; ++a)
        {
            for (int c = 0; c < 2 * l + 1; ++c)
            {
                for (int d = 0; d < nspe * parameters.max_radial; ++d)
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


cvec4 SALTED_Utils::calculate_SALTED_descriptors(const featomic::SimpleSystem& featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters)
{
    metatensor::TensorMap descriptor = get_feats_projs(featomic_system, parameters);
    std::vector<uint8_t> descriptor_buffer = descriptor.save_buffer();
    return get_expansion_coeffs(descriptor_buffer, featomic_system, parameters);
}


//FEATOMIC POWER Spectrum
metatensor::TensorMap SALTED_Utils::calculate_SOAP_Powerspectrum(featomic::SimpleSystem featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters) {
    // create the calculator with its name and parameters
    auto calculator = featomic::Calculator("soap_power_spectrum", parameters.to_json().c_str());

    std::vector<std::array<int32_t,3>> keys_array;
    for (const std::string& center_type : parameters.species)
    {
        int32_t z_center = constants::get_Z_from_label(center_type.c_str()) + 1;

        for (size_t i = 0; i < parameters.neighspe.size(); ++i)
        {
            int32_t z1 = constants::get_Z_from_label(parameters.neighspe[i].c_str()) + 1;

            for (size_t j = i; j < parameters.neighspe.size(); ++j)
            {
                int32_t z2 = constants::get_Z_from_label(parameters.neighspe[j].c_str()) + 1;

                keys_array.push_back({ z_center, z1, z2 });
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
    std::vector<std::string> names = {"center_type", "neighbor_1_type", "neighbor_2_type"};
    metatensor::Labels keys_selection(names, flattened_keys.data(), flattened_keys.size() / names.size());

    featomic::CalculationOptions calc_opts;
    calc_opts.use_native_system = true;
    calc_opts.selected_keys = keys_selection;
    // run the calculation
    // Initialize descriptor directly from computation result
    metatensor::TensorMap descriptor = calculator.compute(featomic_system, calc_opts);

    // The descriptor is a metatensor `TensorMap`, containing multiple blocks.
    // We can transform it to a single block containing a dense representation,
    // with one sample for each atom-centered environment.
    descriptor = descriptor.keys_to_samples("center_type");
    descriptor = descriptor.keys_to_properties(svec{ "neighbor_1_type" , "neighbor_2_type" });

    return descriptor;
}


const double calc_density_ML(const double& x,
    const double& y,
    const double& z,
    const vec& coefficients,
    const std::vector<atom>& atoms)
{
    double dens = 0, radial;
    int coef_counter = 0;
    unsigned int shell = 0, n_shells = 0, prim = 0;
    basis_set_entry bf;
    primitive p;

    for (int a = 0; a < atoms.size(); a++)
    {
        prim = 0;
        n_shells = static_cast<unsigned int>(atoms[a].get_shellcount().size());
        double d[4]{
            x - atoms[a].get_coordinate(0),
            y - atoms[a].get_coordinate(1),
            z - atoms[a].get_coordinate(2), 0.0 };
        // store r in last element
        d[3] = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (d[3] < -46.0517)
        { // corresponds to cutoff of ex ~< 1E-20
            for (shell = 0; shell < n_shells; shell++)
            {
                coef_counter += (2 * atoms[a].get_basis_set_type(prim) + 1);
                prim += atoms[a].get_shellcount()[shell];
            }
            continue;
        }
        // normalize distances for spherical harmonic
        for (int i = 0; i < 3; i++)
            d[i] /= d[3];
        
        for (int shell = 0; shell < n_shells; shell++) {
            radial = 0;
            int type = atoms[a].get_basis_set_entry(prim).get_type();

            for (unsigned int e = 0; e < atoms[a].get_shellcount()[shell]; e++, prim++) {
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
    unsigned int shell = 0, n_shells = 0;

    for (int a = 0; a < atoms.size(); a++)
    {
        n_shells = static_cast<unsigned int>(atoms[a].get_shellcount().size());
        unsigned int prim = 0;
        if (a != atom_nr) {
            for (shell = 0; shell < n_shells; shell++)
            {
                coef_counter += (2 * atoms[a].get_basis_set_type(prim) + 1);
                prim += atoms[a].get_shellcount()[shell];
            }
            continue;
        }
        
        basis_set_entry bf;
        double d[4]{
            x - atoms[a].get_coordinate(0),
            y - atoms[a].get_coordinate(1),
            z - atoms[a].get_coordinate(2), 0.0 };
        // store r in last element
        d[3] = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        // normalize distances for spherical harmonic
        for (int i = 0; i < 3; i++)
            d[i] /= d[3];
        for (shell = 0; shell < n_shells; shell++) {
            radial = 0;
            unsigned int type = atoms[a].get_basis_set_entry(prim).get_type();

            for (unsigned int e = 0; e < atoms[a].get_shellcount()[shell]; e++, prim++) {
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
    //This should never happen, but just in case
    std::cout << "Atom number " << atom_nr << " not found in the list of atoms." << std::endl;
    return -1;
}


/**
 * Calculates the atomic density for a given list of atoms and coefficients.
 *
 * @param atoms The list of atoms.
 * @param coefs The coefficients used in the calculation.
 * @return The atomic density for each atom.
 */
vec calc_atomic_density(const std::vector<atom>& atoms, const vec& coefs)
{
    double radial;
    basis_set_entry bf;

    vec atom_elecs(atoms.size(), 0.0);

    int coef_counter = 0;
    for (int a = 0; a < atoms.size(); a++)
    {

        int type = -1, prim = 0;
        for (unsigned int shell = 0; shell < atoms[a].get_shellcount().size(); shell++) {
            radial = 0;
            type = atoms[a].get_basis_set_entry(prim).get_type();
            if (type != 0)
            {
                coef_counter += (2 * type + 1); prim += atoms[a].get_shellcount()[shell]; //Skip functions and coefficients
                continue;
            }

            for (unsigned int e = 0; e < atoms[a].get_shellcount()[shell]; e++, prim++) {
                bf = atoms[a].get_basis_set_entry(prim);
                primitive p(a, bf.get_type(), bf.get_exponent(), bf.get_coefficient());
                radial += constants::PI / (2.0 * std::pow(p.get_exp(), 1.5)) * p.normalization_constant() * p.get_coef();
            }

            atom_elecs[a] += radial * coefs[coef_counter];
            coef_counter++;
        }
        atom_elecs[a] += atoms[a].get_ECP_electrons();
    }
    return atom_elecs;
}

void calc_cube_ML(const vec& data, WFN& dummy, cube& cube_data, const int& atom_nr)
{
    _time_point start = get_time();

    const int s1 = cube_data.get_size(0), s2 = cube_data.get_size(1), s3 = cube_data.get_size(2), total_size = s1 * s2 * s3;
    std::cout << "Lets go into the loop! There is " << total_size << " points" << std::endl;

    ProgressBar* progress = new ProgressBar(total_size, 60, "=", " ", "Calculating Values");
    vec v1{
        cube_data.get_vector(0, 0),
        cube_data.get_vector(1, 0),
        cube_data.get_vector(2, 0) },
        v2{
            cube_data.get_vector(0, 1),
            cube_data.get_vector(1, 1),
            cube_data.get_vector(2, 1) },
            v3{
                cube_data.get_vector(0, 2),
                cube_data.get_vector(1, 2),
                cube_data.get_vector(2, 2) },
                orig{
                    cube_data.get_origin(0),
                    cube_data.get_origin(1),
                    cube_data.get_origin(2) };

    if (atom_nr != -1)
        std::cout << "Calculation for atom " << atom_nr << std::endl;

    std::vector<atom> atoms = dummy.get_atoms();
#pragma omp parallel for schedule(dynamic)
    for (int index = 0; index < total_size; index++)
    {
        int i = index / (s2 * s3);
        int j = (index / s3) % s2;
        int k = index % s3;

        vec PosGrid{
            i * v1[0] + j * v2[0] + k * v3[0] + orig[0],
            i * v1[1] + j * v2[1] + k * v3[1] + orig[1],
            i * v1[2] + j * v2[2] + k * v3[2] + orig[2] };
        cube_data.set_value(i, j, k,
            (atom_nr == -1)
            ? calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, atoms)
            : calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, atoms, atom_nr));
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
    cube_data.calc_dv();
    std::cout << "Number of electrons: " << std::fixed << std::setprecision(4) << cube_data.sum() << std::endl;
};

cube calc_cube_ML(const vec& data, WFN& dummy, const int& atom_nr)
{
    properties_options opts;
    readxyzMinMax_fromWFN(dummy, opts, true);
    cube CubeRho(opts.NbSteps, dummy.get_ncen(), true);
    CubeRho.give_parent_wfn(dummy);

    for (int i = 0; i < 3; i++)
    {
        CubeRho.set_origin(i, opts.MinMax[i]);
        CubeRho.set_vector(i, i, (opts.MinMax[i + 3] - opts.MinMax[i]) / opts.NbSteps[i]);
    }
    CubeRho.set_comment1("Calculated density using NoSpherA2 from ML Data");
    CubeRho.set_comment2("from " + dummy.get_path().string());
    CubeRho.set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_RI_rho.cube");

    calc_cube_ML(data, dummy, CubeRho, atom_nr);

    return CubeRho;
};

#include "integrator.h"
#include "libCintMain.h"
#include "nos_math.h"
#include "npy.h"
void create_SALTED_training_data(const WFN& orbital, const WFN& aux) {
    std::cout << "Calculating density fitting coefficients..." << std::endl;
    DensityFitting::CONFIG config;
    config.analyze_quality = true;
    //config.restrain_type = DensityFitting::RESTRAINT_TYPE::SIMPLE_AND_TIK;
    //config.charge_scheme = DensityFitting::CHARGE_SCHEME::HIRSHFELD;
    //if (wavy->get_origin() == e_origin::ptb)
    //    config.restraint_strength = 1.0e-4;

    vec coefs = DensityFitting::density_fit(orbital, aux, config);

    vec overlap;
    Int_Params aux_basis(aux);
    compute2C<Overlap2C_SPH>(aux_basis, overlap);
    const int nao_max = aux_basis.get_nao();

    dMatrix1 coefs_vec(coefs.size());
    coefs_vec.container() = coefs;
    dMatrix2 overlap_mat(nao_max, nao_max);
    overlap_mat.container() = overlap;

    dMatrix1 proj = dot(overlap_mat, coefs_vec);

    npy::write_npy("coefficients.npy", 
        npy::npy_data<double>{
            coefs, 
            { static_cast<unsigned long>(coefs.size()) },
            false}
    );

    npy::write_npy("projections.npy",
        npy::npy_data<double>{
        proj.container(),
        { static_cast<unsigned long>(proj.size()) },
            false}
    );

    npy::write_npy("overlap.npy",
        npy::npy_data<double>{
        overlap,
        { static_cast<unsigned long>(nao_max), static_cast<unsigned long>(nao_max) },
            false}
    );

}
