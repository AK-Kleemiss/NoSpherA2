#pragma once

#include <string>
#include <fstream>
#include <vector>
#include <ostream>
#include <filesystem>

class WFN;
class cell;
struct options;
struct properties_options;
// NOTE: this enum is coupled BY POSITION to the emplace_back() sequence that
// builds the cube vector in properties_calculation() (properties.cpp). Nothing
// enforces the correspondence, so inserting an enumerator in the middle without
// inserting the matching emplace_back at the same index silently shifts every
// later property by one. Always append new entries at the end.
enum cube_type {
    Rho = 0,
    RDG = 1,
    Elf = 2,
    Eli = 3,
    Lap = 4,
    ESP = 5,
    MO_val = 6,
    HDEF = 7,
    DEF = 8,
    Hirsh = 9,
    Spin_Density = 10,
    spherical_density = 11,
    Fukui_plus = 12,
    Fukui_minus = 13,
    Fukui_zero = 14,
    Dual_Descriptor = 15
};

/**
 * Calculates the static deflection using the given parameters.
 *
 * @param CubeDEF The cube object representing the deflection.
 * @param CubeRho The cube object representing the density.
 * @param wavy The WFN object representing the wave.
 * @param radius The radius value for the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Static_Def(
    cube &CubeDEF,
    cube &CubeRho,
    const WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap = true);

/**
 * Calculates the static definition of a cube.
 *
 * This function calculates the static definition of a cube based on the given parameters.
 *
 * @param CubeDEF The cube object representing the static definition.
 * @param CubeRho The cube object representing the density.
 * @param CubeSpher The cube object representing the spherical coordinates.
 * @param wavy The WFN object representing the wave function.
 * @param radius The radius parameter for the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Static_Def(
    std::vector<cube> &Cubes,
    const WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap = true);
/**
 * Calculates the spherical density of a cube using the given WFN object.
 *
 * @param CubeSpher The cube object to store the spherical density.
 * @param wavy The WFN object containing the wavefunction data.
 * @param radius The radius parameter for the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Spherical_Dens(
    cube &CubeSpher,
    const WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap = true);
/**
 * Calculates the density (Rho) for a given cube and WFN object.
 *
 * @param CubeRho The cube object to store the calculated density.
 * @param wavy The WFN object containing the wavefunction information.
 * @param radius The radius parameter for the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Rho(
    cube &CubeRho,
    const WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap = true);

void Calc_Eli(
    cube &CubeRho,
    const WFN &wavy,
    double radius,
    std::ostream &file = std::cout,
    bool wrap = false);

/**
 * Calculates the density based on a wfn with spherical harmonicsand stores the result in the given cube.
 *
 * @param CubeRho The cube object to store the calculated spherical harmonics.
 * @param wavy The WFN object containing the input data for the calculation.
 * @param file The output stream to write the result to.
 */
void Calc_Rho_spherical_harmonics(
    cube &CubeRho,
    const WFN &wavy,
    std::ostream &file);
/**
 * Calculates the molecular orbital (MO) using the spherical harmonics.
 *
 * @param CubeMO The cube object to store the calculated spherical harmonics.
 * @param wavy The WFN object containing the wavefunction information.
 * @param MO The index of the molecular orbital to calculate the spherical harmonics for.
 * @param file The output stream to write the calculated spherical harmonics.
 */
void Calc_MO_spherical_harmonics(
    cube &CubeRho,
    const WFN &wavy,
    int MO,
    std::ostream &file,
    bool nodate = false);
/**
 * Calculates the properties of the given cubes and WFN object.
 *
 * @param CubeRho The cube object representing the electron density.
 * @param CubeRDG The cube object representing the reduced density gradient.
 * @param CubeElf The cube object representing the electron localization function.
 * @param CubeEli The cube object representing the electron localization index.
 * @param CubeLap The cube object representing the Laplacian of the electron density.
 * @param CubeESP The cube object representing the electrostatic potential.
 * @param wavy The WFN object representing the wavefunction.
 * @param radius The radius parameter for the calculation.
 * @param file The output stream to write the results to.
 * @param test A boolean flag indicating whether to run the function in test mode.
 */
void Calc_Prop(
    std::vector<cube> &Cubes,
    const WFN &wavy,
    double radius,
    std::ostream &file,
    bool test,
    bool wrap = true);
/**
 * Calculates the Electrostatic Potential (ESP) for a given cube and WFN object.
 *
 * @param CubeESP The cube object to store the calculated ESP.
 * @param wavy The WFN object containing the wavefunction information.
 * @param radius The radius parameter for the ESP calculation.
 * @param no_date A flag indicating whether to include the date in the output.
 * @param file The output stream to write the ESP results.
 */
void Calc_ESP(
    cube &CubeESP,
    const WFN &wavy,
    double radius,
    bool no_date,
    std::ostream &file,
    bool wrap = true);
/**
 * Calculates the molecular orbital (MO) for a given cube.
 *
 * @param CubeMO The cube object to store the calculated MO.
 * @param mo The index of the MO to calculate.
 * @param wavy The WFN object containing the wavefunction data.
 * @param radius The radius for the calculation.
 * @param file The output stream to write the calculated MO.
 */
void Calc_MO(
    cube &CubeMO,
    int mo,
    const WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap = true);
/**
 * Locates the frontier molecular orbitals (HOMO and LUMO) of a wavefunction.
 *
 * There is no stored HOMO/LUMO index in WFN, so the frontier pair is derived
 * from the occupation numbers. Among all occupied MOs the one with the highest
 * orbital energy is taken as the HOMO, and among all unoccupied MOs the one
 * with the lowest orbital energy as the LUMO. Some wavefunction formats carry
 * no meaningful orbital energies (every energy identically zero), in which case
 * the routine falls back to MO ordering: the last occupied and the first
 * unoccupied index.
 *
 * For unrestricted wavefunctions the frontier pair is taken across BOTH spin
 * manifolds, i.e. the global highest occupied and global lowest virtual. That
 * is a deliberate simplification - a spin-resolved Fukui function would need
 * one pair per operator - and is reported to the caller via @p unrestricted.
 *
 * @param wavy The wavefunction to inspect.
 * @param homo Receives the HOMO index, or -1 if no occupied MO exists.
 * @param lumo Receives the LUMO index, or -1 if no virtual MO exists.
 * @param unrestricted Receives whether the wavefunction is unrestricted.
 * @return true if BOTH a HOMO and a LUMO were found.
 */
bool find_frontier_orbitals(
    const WFN &wavy,
    int &homo,
    int &lumo,
    bool &unrestricted);
/**
 * Calculates the Fukui functions and the dual descriptor in the
 * frontier-orbital (frozen-orbital) approximation.
 *
 * f+(r) = |psi_LUMO(r)|^2   response to adding an electron; large where the
 *                           molecule is attacked by a NUCLEOPHILE.
 * f-(r) = |psi_HOMO(r)|^2   response to removing an electron; large where the
 *                           molecule is attacked by an ELECTROPHILE.
 * f0(r) = (f+ + f-) / 2     radical attack.
 * df(r) = f+ - f-           dual descriptor; > 0 electrophilic region,
 *                           < 0 nucleophilic region.
 *
 * All four fields derive from the same two orbital amplitudes, so they are
 * evaluated in a single pass over the grid. Only the cubes that are flagged as
 * loaded are filled. Values are accumulated (+=) rather than assigned, which is
 * what the periodic (wrapped) grid path requires.
 *
 * NOTE: Calc_MO returns the orbital AMPLITUDE psi, not the density |psi|^2.
 * The squaring is done here.
 *
 * @param Cubes The cube vector; Fukui_plus/minus/zero and Dual_Descriptor are filled.
 * @param wavy The WFN object containing the wavefunction data.
 * @param homo Index of the HOMO, from find_frontier_orbitals.
 * @param lumo Index of the LUMO, from find_frontier_orbitals.
 * @param radius The radius (Angstrom) around the atoms to evaluate.
 * @param file The output stream for timing information.
 * @param wrap Whether to wrap the grid periodically (set when a CIF is given).
 */
void Calc_Fukui(
    std::vector<cube> &Cubes,
    const WFN &wavy,
    int homo,
    int lumo,
    double radius,
    std::ostream &file,
    bool wrap = true);
/**
 * Condensed (atom-summed) Fukui functions under all five atomic partitions.
 *
 * f+_A and f-_A are indexed [partition][atom], with the partition index being
 * PartitionResults::CHARGE_ORDER (S_BECKE, S_TFVC, S_HIRSH, S_MBIS, S_EMBIS).
 */
// Spelled with explicit std types rather than the svec/vec2 aliases: those live
// in convenience.h, which properties.h deliberately does not include.
struct CondensedFukuiResults {
    bool valid = false;
    std::vector<std::string> labels;
    std::vector<std::vector<double>> f_plus;   // [partition][atom]
    std::vector<std::vector<double>> f_minus;  // [partition][atom]
};
/**
 * Integrates the frontier-orbital densities against every atomic partition
 * weight GridManager knows about, giving the condensed Fukui functions
 *
 *     f+_A = integral w_A(r) |psi_LUMO(r)|^2 dr
 *     f-_A = integral w_A(r) |psi_HOMO(r)|^2 dr
 *
 * for Becke, Hirshfeld, TFVC, MBIS and EMBIS at once.
 *
 * The point-wise Fukui function is dominated by near-nuclear density, so its
 * maxima track the heaviest atom rather than the reactive site. The condensed
 * form is what actually carries the chemistry, and it is the standard object in
 * the conceptual-DFT literature (De Proft et al., J. Comput. Chem. 2002, 23,
 * 1198 condense it with exactly the Hirshfeld weights used here).
 *
 * The partition weights are built from the GROUND-STATE density, as they must
 * be - MBIS and EMBIS are refined self-consistently against it. Only the
 * integrand is swapped for a frontier-orbital density.
 *
 * @param wavy The wavefunction, which must still carry its virtual orbitals.
 * @param homo HOMO index from find_frontier_orbitals.
 * @param lumo LUMO index from find_frontier_orbitals.
 * @param unit_cell Cell for the integration grid.
 * @param accuracy GridManager accuracy level (as used for atomic charges).
 * @param file Output stream for GridManager's own progress reporting.
 */
CondensedFukuiResults Calc_Condensed_Fukui(
    const WFN &wavy,
    int homo,
    int lumo,
    const cell &unit_cell,
    int accuracy,
    std::ostream &file);
/**
 * Prints the condensed Fukui table: f+, f- and the dual descriptor, one column
 * per partition, with the per-quantity sum over atoms as a normalisation check
 * (each should be 1).
 */
void print_condensed_fukui(
    const CondensedFukuiResults &r,
    std::ostream &file);
/**
 * Standalone conceptual-DFT reactivity analysis of one wavefunction.
 *
 * Reads opt.wfn, reports the frontier orbitals and the HOMO-LUMO gap, computes
 * the condensed Fukui functions under all five atomic partitions, names the most
 * electrophilic and most nucleophilic atom, and writes <stem>_fukui.dat.
 *
 * Deliberately does NOT touch the cube machinery: no grid, no radius, no
 * resolution and no CIF are needed, so this answers "where does this molecule
 * react" from a wavefunction alone. Driven by -fukui_analysis; the cube-producing
 * path is -fukui inside a properties run.
 *
 * Dispatched from run_app_impl, NOT from options::digest_options(), so that the
 * console streambuf can be restored first and the table actually reaches the
 * terminal rather than NoSpherA2.log.
 */
void fukui_analysis(options &opt, std::ostream &log2 = std::cout);
/**
 * Calculates the Spin density cube using the provided WFN object.
 *
 * @param Cube_S_Rho The cube object to store the calculated S_Rho cube.
 * @param wavy The WFN object used for the calculation.
 * @param file The output stream to write the results to.
 * @param nodate A boolean flag indicating whether to include the date in the output.
 */
void Calc_S_Rho(
    cube &Cube_S_Rho,
    const WFN &wavy,
    std::ostream &file,
    bool &nodate);
/**
 * Calculates the Hirshfeld Deformation Density for a given set of parameters.
 *
 * @param CubeHDEF The cube object representing the Hirshfeld Density.
 * @param CubeRho The cube object representing the electron density.
 * @param wavy The WFN object containing additional information.
 * @param radius The radius parameter for the calculation.
 * @param ignore_atom The index of the atom to ignore in the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Hirshfeld(
    std::vector<cube> &Cubes,
    const WFN &wavy,
    double radius,
    int ignore_atom,
    std::ostream &file,
    bool wrap = true);
/**
 * Calculates the Hirshfeld Deformation Density for a given set of parameters.
 *
 * @param CubeHDEF The cube object representing the Hirshfeld Density.
 * @param CubeRho The cube object representing the electron density.
 * @param CubeSpherical The cube object representing the spherical density.
 * @param wavy The WFN object containing the wavefunction information.
 * @param radius The radius parameter for the calculation.
 * @param ignore_atom The index of the atom to be ignored in the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Hirshfeld(
    cube &CubeHDEF,
    cube &CubeRho,
    cube &CubeSpherical,
    const WFN &wavy,
    double radius,
    int ignore,
    std::ostream &file,
    bool wrap = true);
/**
 * Calculates the Hirshfeld atom electron density.
 *
 * @param CubeHirsh The cube object to store the calculated Hirshfeld atom.
 * @param CubeRho The cube object containing the electron density.
 * @param CubeSpherical The cube object containing the spherical density.
 * @param wavy The WFN object containing the wavefunction information.
 * @param radius The radius parameter for the calculation.
 * @param ignore_atom The index of the atom to ignore in the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Hirshfeld_atom(
    std::vector<cube> &Cubes,
    const WFN &wavy,
    double radius,
    int ignore_atom,
    std::ostream &file,
    bool wrap = true);

void Calc_RhoEli(
    cube &CubeRho,
    cube &CubeEli,
    const WFN &wavy,
    double radius);

/**
 * Calculates the properties based on the given options.
 *
 * @param opt The options object containing the necessary data for the calculation.
 */
void properties_calculation(options &opt);

void promolecular_nci_analysis(
    const std::filesystem::path& xyz1,
    const std::filesystem::path& xyz2,
    const properties_options& opts,
    std::ostream& log);

/**
 * Combines the MO (Molecular Orbital) files based on the given options.
 *
 * @param opt The options specifying the details of the combination process.
 */
void do_combine_mo(options &opt);

void dipole_moments(options &opt, std::ostream &log2 = std::cout);
vec2 dipole_moments(WFN &wavy, cube &SPHER, double *MinMax, int *NbSteps, int threads, double radius, std::ostream &log2 = std::cout, bool debug = false);
void polarizabilities(options &opt, std::ostream &log2 = std::cout);

#include "wfn_class.h"
#include "cell.h"
