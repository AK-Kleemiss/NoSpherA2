/**
 * @file scattering_factors.h
 * @brief This file contains the declarations of functions and classes related to scattering factors calculations.
 */

#pragma once

#include <string>
#include <vector>
#include <fstream>
#include "convenience.h"
#include "SALTED_predictor.h"

/**
 * @class WFN
 * @brief Class representing the wavefunction.
 */
class WFN;
template class tsc_block<int, cdouble>;
typedef tsc_block<int, cdouble> itsc_block;

/**
 * @class cell
 * @brief Class representing the unit cell.
 */
class cell;

/**
 * @struct options
 * @brief Struct representing the options for calculations.
 */
struct options;

/**
 * @brief Calculates the Thakkar scattering factors.
 * @param opt The options for calculations.
 * @param file The output stream to write the results to.
 * @param wave The wavefunction.
 * @return True if the calculation is successful, false otherwise.
 */
bool thakkar_sfac(const options &opt, std::ostream &file, const WFN &wave);

/**
 * @brief Calculates the MTC (Modified Thakkar) scattering factors.
 * @param opt The options for scattering factors calculations.
 * @param file The output stream to write the results to.
 * @param known_atoms The list of known atoms.
 * @param wave The list of wavefunctions.
 * @param nr The number of wavefunctions.
 * @return True if the calculation is successful, false otherwise.
 */
itsc_block MTC_thakkar_sfac(
    options &opt,
    std::ostream &file,
    svec &known_atoms,
    std::vector<WFN> &wave,
    const int &nr);

/**
 * @brief Calculates the scattering factors for HF (Hartree-Fock) method.
 * @param opt The options for scattering factors calculations.
 * @param wave The wavefunction.
 * @param file The output stream to write the results to.
 * @return True if the calculation is successful, false otherwise.
 */
bool calculate_scattering_factors_HF(const options &opt, const WFN &wave, std::ostream &file);

/**
 * @brief Calculates the scattering factors for ML method.
 * @param opt The options for scattering factors calculations.
 * @param wave The wavefunction.
 * @param file The output stream to write the results to.
 * @param exp_coefs The number of expansion coefficients.
 * @return True if the calculation is successful, false otherwise.
 */
bool calculate_scattering_factors_ML(
    options &opt,
    SALTEDPredictor &SP,
    std::ostream &file);

/**
 * @brief Calculates the MTC (Modified Thakkar) scattering factors.
 * @param opt The options for scattering factors calculations.
 * @param wave The list of wavefunctions.
 * @param file The output stream to write the results to.
 * @param known_atoms The list of known atoms.
 * @param nr The number of wavefunctions.
 * @param kpts The list of k-points.
 * @return The calculated scattering factors.
 */
itsc_block calculate_scattering_factors_MTC_SALTED(
    options& opt,
    SALTEDPredictor& SP,
    std::ostream& file,
    svec& known_atoms,
    const int& nr,
    vec2* kpts = NULL);

/**
 * @brief Calculates the MTC (Modified Thakkar) scattering factors.
 * @param opt The options for scattering factors calculations.
 * @param wave The list of wavefunctions.
 * @param file The output stream to write the results to.
 * @param known_atoms The list of known atoms.
 * @param nr The number of wavefunctions.
 * @param kpts The list of k-points.
 * @return The calculated scattering factors.
 */
itsc_block calculate_scattering_factors_MTC(
    options &opt,
    const std::vector<WFN> &wave,
    std::ostream &file,
    svec &known_atoms,
    const int &nr,
    vec2 *kpts = NULL);

/**
 * @brief Generates the hkl (Miller indices) list.
 * @param dmin The minimum d-spacing.
 * @param hkl The generated hkl list.
 * @param twin_law The list of twin laws.
 * @param unit_cell The unit cell.
 * @param file The output stream to write the results to.
 * @param debug Flag indicating whether to enable debug mode.
 */
void generate_hkl(
    const double &dmin,
    hkl_list &hkl,
    const vec2 &twin_law,
    cell &unit_cell,
    std::ostream &file,
    bool debug = false);

void generate_hkl(
    const ivec2 &hkl_min_max,
    hkl_list &hkl,
    const vec2 &twin_law,
    cell &unit_cell,
    std::ostream &file,
    bool debug = false,
    bool ED = false);

/**
 * @brief Generates the fractional hkl (Miller indices) list.
 * @param dmin The minimum d-spacing.
 * @param hkl The generated fractional hkl list.
 * @param twin_law The list of twin laws.
 * @param unit_cell The unit cell.
 * @param file The output stream to write the results to.
 * @param stepsize The step size for generating fractional hkl.
 * @param debug Flag indicating whether to enable debug mode.
 */
void generate_fractional_hkl(
    const double &dmin,
    hkl_list_d &hkl,
    const vec2 &twin_law,
    cell &unit_cell,
    std::ostream &file,
    double stepsize,
    bool debug);

/**
 * @brief Converts the scattering factor to the ED (Electron Density) single.
 * @tparam NumType The type of the charge.
 * @param neutralcharge The neutral charge.
 * @param sf The scattering factor.
 * @param k_vector The k-vector.
 * @param charge The charge.
 * @return The converted ED single.
 */
cdouble convert_to_ED_single(
    const int &neutralcharge,
    cdouble &sf,
    const double &k_vector,
    const int charge = 0);

/**
 * @brief Reads atoms from CIF (Crystallographic Information File).
 * @param cif_input The input stream for the CIF file.
 * @param input_groups The list of input groups.
 * @param unit_cell The unit cell.
 * @param wave The wavefunction.
 * @param known_atoms The list of known atoms.
 * @param atom_type_list The list of atom types.
 * @param asym_atom_to_type_list The list of asymmetric atoms to atom types.
 * @param asym_atom_list The list of asymmetric atoms.
 * @param needs_grid The list indicating whether an atom needs a grid.
 * @param file The output stream to write the results to.
 * @param debug Flag indicating whether to enable debug mode.
 */
svec read_atoms_from_CIF(
    std::ifstream &cif_input,
    const ivec &input_groups,
    const cell &unit_cell,
    const WFN &wave,
    const svec &known_atoms,
    ivec &atom_type_list,
    ivec &asym_atom_to_type_list,
    ivec &asym_atom_list,
    bvec &needs_grid,
    std::ostream &file,
    bvec &constant_atoms,
    const bool SALTED = false,
    const bool debug = false);

/**
 * @brief Makes Hirshfeld grids.
 * @param pbc The periodic boundary condition.
 * @param accuracy The accuracy.
 * @param unit_cell The unit cell.
 * @param wave The wavefunction.
 * @param atom_type_list The list of atom types.
 * @param asym_atom_list The list of asymmetric atoms.
 * @param needs_grid The list indicating whether an atom needs a grid.
 * @param d1 The list of d1 vectors.
 * @param d2 The list of d2 vectors.
 * @param d3 The list of d3 vectors.
 * @param dens The list of density vectors.
 * @param file The output stream to write the results to.
 * @param start The start time point.
 * @param end_becke The end time point for Becke grid generation.
 * @param end_prototypes The end time point for prototypes generation.
 * @param end_spherical The end time point for spherical grid generation.
 * @param end_prune The end time point for pruning grid points.
 * @param end_aspherical The end time point for aspherical grid generation.
 * @param debug Flag indicating whether to enable debug mode.
 * @param no_date Flag indicating whether to exclude the date from the output.
 * @return The number of Hirshfeld grids generated.
 */
int make_hirshfeld_grids(
    const int &pbc,
    const int &accuracy,
    cell &unit_cell,
    const WFN &wave,
    const ivec &atom_type_list,
    const ivec &cif2wfn_list,
    bvec &needs_grid,
    vec2 &d1,
    vec2 &d2,
    vec2 &d3,
    vec2 &dens,
    const svec &labels,
    std::ostream &file,
    std::vector<time_point> &time_points,
    std::vector<std::string> &time_descriptions,
    bool debug = false,
    bool no_date = false);

/**
 * @brief Adds ECP (Effective Core Potential) contribution to the scattering factors.
 * @param asym_atom_list The list of asymmetric atoms.
 * @param wave The wavefunction.
 * @param sf The scattering factors.
 * @param cell The unit cell.
 * @param hkl The hkl (Miller indices) list.
 * @param file The output stream to write the results to.
 * @param mode The mode for adding ECP contribution.
 * @param debug Flag indicating whether to enable debug mode.
 */
static void add_ECP_contribution(
    const ivec &asym_atom_list,
    const WFN &wave,
    cvec2 &sf,
    const cell &cell,
    hkl_list &hkl,
    std::ostream &file,
    const int &mode = 0,
    const bool debug = false);

/**
 * @brief Calculates the scattering factors.
 * @param points The number of points.
 * @param k_pt The list of k-points.
 * @param d1 The list of d1 vectors.
 * @param d2 The list of d2 vectors.
 * @param d3 The list of d3 vectors.
 * @param dens The list of density vectors.
 * @param sf The scattering factors.
 * @param file The output stream to write the results to.
 * @param start The start time point.
 * @param end1 The end time point.
 * @param debug Flag indicating whether to enable debug mode.
 */
void calc_SF(const int &points,
             vec2 &k_pt,
             vec2 &d1,
             vec2 &d2,
             vec2 &d3,
             vec2 &dens,
             cvec2 &sf,
             std::ostream &file,
             time_point &start,
             time_point &end1,
             bool debug,
             bool no_date = false);

double fourier_bessel_integral(
    const primitive &p,
    const double &H);

cdouble sfac_bessel(
    const primitive &p,
    const double *k_point,
    const vec &coefs);

/**
 * @brief Calculates the diffuse scattering factors.
 * @param opt The options for scattering factors calculations.
 * @param log_file The output stream to write the log to.
 */
void calc_sfac_diffuse(const options &opt, std::ostream &log_file);

#include "wfn_class.h"
#include "tsc_block.h"