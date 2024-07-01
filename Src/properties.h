#pragma once

#include <string>
#include <fstream>
#include <vector>

class WFN;
class cell;
/**
 * Calculates the static deflection using the given parameters.
 *
 * @param CubeDEF The cube object representing the deflection.
 * @param CubeRho The cube object representing the density.
 * @param wavy The WFN object representing the wave.
 * @param cpus The number of CPUs to use for the calculation.
 * @param radius The radius value for the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Static_Def(
    cube &CubeDEF,
    cube &CubeRho,
    WFN &wavy,
    int cpus,
    double radius,
    std::ostream &file);

/**
 * Calculates the static definition of a cube.
 *
 * This function calculates the static definition of a cube based on the given parameters.
 *
 * @param CubeDEF The cube object representing the static definition.
 * @param CubeRho The cube object representing the density.
 * @param CubeSpher The cube object representing the spherical coordinates.
 * @param wavy The WFN object representing the wave function.
 * @param cpus The number of CPUs to be used for the calculation.
 * @param radius The radius parameter for the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Static_Def(
    cube &CubeDEF,
    cube &CubeRho,
    cube &CubeSpher,
    WFN &wavy,
    int cpus,
    double radius,
    std::ostream &file);
/**
 * Calculates the spherical density of a cube using the given WFN object.
 *
 * @param CubeSpher The cube object to store the spherical density.
 * @param wavy The WFN object containing the wavefunction data.
 * @param cpus The number of CPUs to use for the calculation.
 * @param radius The radius parameter for the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Spherical_Dens(
    cube &CubeSpher,
    WFN &wavy,
    int cpus,
    double radius,
    std::ostream &file);
/**
 * Calculates the density (Rho) for a given cube and WFN object.
 *
 * @param CubeRho The cube object to store the calculated density.
 * @param wavy The WFN object containing the wavefunction information.
 * @param cpus The number of CPUs to use for the calculation.
 * @param radius The radius parameter for the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Rho(
    cube &CubeRho,
    WFN &wavy,
    int cpus,
    double radius,
    std::ostream &file);
/**
 * Calculates the density matrix without any transformation.
 *
 * @param CubeRho The cube object to store the calculated density matrix.
 * @param wavy The WFN object containing the wavefunction data.
 * @param cpus The number of CPUs to use for the calculation.
 * @param radius The radius parameter for the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Rho_no_trans(
    cube &CubeRho,
    WFN &wavy,
    int cpus,
    double radius,
    std::ostream &file);
/**
 * Calculates the density based on a wfn with spherical harmonicsand stores the result in the given cube.
 *
 * @param CubeRho The cube object to store the calculated spherical harmonics.
 * @param wavy The WFN object containing the input data for the calculation.
 * @param cpus The number of CPUs to use for the calculation.
 * @param file The output stream to write the result to.
 */
void Calc_Rho_spherical_harmonics(
    cube &CubeRho,
    WFN &wavy,
    int cpus,
    std::ostream &file);
/**
 * Calculates the molecular orbital (MO) using the spherical harmonics.
 *
 * @param CubeMO The cube object to store the calculated spherical harmonics.
 * @param wavy The WFN object containing the wavefunction information.
 * @param cpus The number of CPUs to use for the calculation.
 * @param MO The index of the molecular orbital to calculate the spherical harmonics for.
 * @param file The output stream to write the calculated spherical harmonics.
 */
void Calc_MO_spherical_harmonics(
    cube &CubeRho,
    WFN &wavy,
    int cpus,
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
 * @param cpus The number of CPUs to be used for the calculation.
 * @param radius The radius parameter for the calculation.
 * @param file The output stream to write the results to.
 * @param test A boolean flag indicating whether to run the function in test mode.
 */
void Calc_Prop(
    cube &CubeRho,
    cube &CubeRDG,
    cube &CubeElf,
    cube &CubeEli,
    cube &CubeLap,
    cube &CubeESP,
    WFN &wavy,
    int cpus,
    double radius,
    std::ostream &file,
    bool test);
/**
 * Calculates the Electrostatic Potential (ESP) for a given cube and WFN object.
 *
 * @param CubeESP The cube object to store the calculated ESP.
 * @param wavy The WFN object containing the wavefunction information.
 * @param cpus The number of CPUs to use for the calculation.
 * @param radius The radius parameter for the ESP calculation.
 * @param no_date A flag indicating whether to include the date in the output.
 * @param file The output stream to write the ESP results.
 */
void Calc_ESP(
    cube &CubeESP,
    WFN &wavy,
    int cpus,
    double radius,
    bool no_date,
    std::ostream &file);
/**
 * Calculates the molecular orbital (MO) for a given cube.
 *
 * @param CubeMO The cube object to store the calculated MO.
 * @param mo The index of the MO to calculate.
 * @param wavy The WFN object containing the wavefunction data.
 * @param cpus The number of CPUs to use for the calculation.
 * @param radius The radius for the calculation.
 * @param file The output stream to write the calculated MO.
 */
void Calc_MO(
    cube &CubeMO,
    int mo,
    WFN &wavy,
    int cpus,
    double radius,
    std::ostream &file);
/**
 * Calculates the Spin density cube using the provided WFN object.
 *
 * @param Cube_S_Rho The cube object to store the calculated S_Rho cube.
 * @param wavy The WFN object used for the calculation.
 * @param cpus The number of CPUs to use for the calculation.
 * @param file The output stream to write the results to.
 * @param nodate A boolean flag indicating whether to include the date in the output.
 */
void Calc_S_Rho(
    cube &Cube_S_Rho,
    WFN &wavy,
    int cpus,
    std::ostream &file,
    bool nodate);
/**
 * Calculates the Hirshfeld Deformation Density for a given set of parameters.
 *
 * @param CubeHDEF The cube object representing the Hirshfeld Density.
 * @param CubeRho The cube object representing the electron density.
 * @param wavy The WFN object containing additional information.
 * @param cpus The number of CPUs to use for the calculation.
 * @param radius The radius parameter for the calculation.
 * @param ignore_atom The index of the atom to ignore in the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Hirshfeld(
    cube &CubeHDEF,
    cube &CubeRho,
    WFN &wavy,
    int cpus,
    double radius,
    int ignore,
    std::ostream &file);
/**
 * Calculates the Hirshfeld Deformation Density for a given set of parameters.
 *
 * @param CubeHDEF The cube object representing the Hirshfeld Density.
 * @param CubeRho The cube object representing the electron density.
 * @param CubeSpherical The cube object representing the spherical density.
 * @param wavy The WFN object containing the wavefunction information.
 * @param cpus The number of CPUs to be used for the calculation.
 * @param radius The radius parameter for the calculation.
 * @param ignore_atom The index of the atom to be ignored in the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Hirshfeld(
    cube &CubeHDEF,
    cube &CubeRho,
    cube &CubeSpherical,
    WFN &wavy,
    int cpus,
    double radius,
    int ignore,
    std::ostream &file);
/**
 * Calculates the Hirshfeld atom electron density.
 *
 * @param CubeHirsh The cube object to store the calculated Hirshfeld atom.
 * @param CubeRho The cube object containing the electron density.
 * @param CubeSpherical The cube object containing the spherical density.
 * @param wavy The WFN object containing the wavefunction information.
 * @param cpus The number of CPUs to use for the calculation.
 * @param radius The radius parameter for the calculation.
 * @param ignore_atom The index of the atom to ignore in the calculation.
 * @param file The output stream to write the results to.
 */
void Calc_Hirshfeld_atom(
    cube &CubeHirsh,
    cube &CubeRho,
    cube &CubeSpherical,
    WFN &wavy,
    int cpus,
    double radius,
    int ignore,
    std::ostream &file);

/**
 * Calculates the properties based on the given options.
 *
 * @param opt The options object containing the necessary data for the calculation.
 */
void properties_calculation(options &opt);

/**
 * Combines the MO (Molecular Orbital) files based on the given options.
 *
 * @param opt The options specifying the details of the combination process.
 */
void do_combine_mo(options &opt);

#include "wfn_class.h"
#include "cell.h"