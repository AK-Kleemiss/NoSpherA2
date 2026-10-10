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
// Coupled by position to the emplace_back sequence in properties_calculation(); append new entries at the end only
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
	bool wrap);

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
	bool wrap);
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
	bool wrap);
/**
 * Calculates the density (Rho) for a given cube and WFN object.
 *
 * @param CubeRho The cube object to store the calculated density.
 * @param wavy The WFN object containing the wavefunction information.
 * @param radius The radius parameter for the calculation.
 * @param file The output stream to write the results to.
 */
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
	bool wrap);
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
	bool wrap);
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
	bool wrap);
// HOMO = highest-energy occupied MO, LUMO = lowest-energy virtual (-1 if none); with all orbital energies zero
// it falls back to MO order. Unrestricted: one pair across both spins, a simplification flagged by unrestricted.
// Returns true if both exist
bool find_frontier_orbitals(
	const WFN &wavy,
	int &homo,
	int &lumo,
	bool &unrestricted);
// Frozen-orbital Fukui functions into the loaded cubes, one grid pass: f+ = |psi_LUMO|^2 (nucleophilic attack),
// f- = |psi_HOMO|^2 (electrophilic attack), f0 = (f+ + f-)/2, dual descriptor f+ - f- (> 0 electrophilic).
// Calc_MO gives the amplitude psi, squared here; values accumulate (+=) for the wrapped grid. radius in Angstrom
void Calc_Fukui(
	std::vector<cube> &Cubes,
	const WFN &wavy,
	int homo,
	int lumo,
	double radius,
	std::ostream &file,
	bool wrap);
// Partition index follows PartitionResults::CHARGE_ORDER (Becke, TFVC, Hirshfeld, MBIS, EMBIS)
struct CondensedFukuiResults {
	bool valid = false;
	std::vector<std::string> labels;
	vec2 f_plus;   // [partition][atom]
	vec2 f_minus;  // [partition][atom]
};
// f+_A = int w_A |psi_LUMO|^2, f-_A = int w_A |psi_HOMO|^2 for every GridManager partition (De Proft et al.,
// J. Comput. Chem. 2002, 23, 1198); the point-wise maxima sit on heavy-atom cores, the condensed form finds the
// reactive site. Weights stay those of the ground-state density (MBIS/EMBIS are refined against it); wavy needs its virtuals
CondensedFukuiResults Calc_Condensed_Fukui(
	const WFN &wavy,
	int homo,
	int lumo,
	const cell &unit_cell,
	int accuracy,
	std::ostream &file);
// One column per partition; the sums over atoms (each 1) check the normalisation
void print_condensed_fukui(
	const CondensedFukuiResults &r,
	std::ostream &file);
// -fukui_analysis: frontier orbitals, gap and condensed Fukui functions from opt.wfn alone, no cubes, into <stem>_fukui.dat.
// Dispatched from run_app_impl, not digest_options(), so the console streambuf is restored and the table reaches the terminal
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
	bool wrap);
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
	bool wrap);
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
	bool wrap);

//rho and ELI-D on two cubes from the orbitals. With `field` (a fitted density) rho comes
//from the fit; ELI-D still needs the orbitals and stays zero without them
struct density_field;
void Calc_RhoEli(
	cube &CubeRho,
	cube &CubeEli,
	const WFN &wavy,
	double radius,
	const density_field *field = nullptr);

/**
 * Calculates the properties based on the given options.
 *
 * @param opt The options object containing the necessary data for the calculation.
 */
void properties_calculation(options &opt);

void promolecular_nci_analysis(
	const pathvec& xyz_files,
	const properties_options& opts,
	std::ostream& log,
	const std::filesystem::path& cif = {}); // cif: grid = the unit cell instead of a box around the fragments

/**
 * Combines the MO (Molecular Orbital) files based on the given options.
 *
 * @param opt The options specifying the details of the combination process.
 */
void do_combine_mo(options &opt);

void dipole_moments(options &opt, std::ostream &log2 = std::cout);
void polarizabilities(options &opt, std::ostream &log2 = std::cout);

#include "wfn_class.h"
#include "cell.h"
#include "density_source.h"

void print_time(_time_point &start, _time_point &end, std::ostream &file);

//The same cube functions over any density source (density_source.h): the WFN, a Gaussian_Molecule, a Gaussian_Atom,
//a Centred atom model. Calc_Prop and Calc_ESP keep WFN overloads for the orbital kernels and the ESP pair table;
//here rho, its derivatives and the ESP come from the source alone, so ELF (orbitals) is refused with a message
inline bool near_any(const d3 &pos, const std::vector<d3> &centres, double radius_bohr)
{
	for (const d3 &c : centres)
		if (array_length(pos, c) < radius_bohr)
			return true;
	return false;
}
//Any point function on the grid within radius of the source; wrap sums the 27 periodic images (a CIF grid)
template <class S, typename F>
void Calc_Cube(cube &Cube, const S &src, F &&f, double radius, std::ostream &file, bool wrap, bool no_date = false)
{
	_time_point start = get_time();
	const double r = constants::ang2bohr(radius);
	const std::vector<d3> centres = source_positions(src);
	Cube.evaluate_on_grid([&](const d3 &pos) { return near_any(pos, centres, r) ? f(pos) : 0.0; }, wrap);
	if (!no_date)
	{
		_time_point end = get_time();
		print_time(start, end, file);
	}
}
template <class S>
void Calc_Rho(cube &CubeRho, const S &src, double radius, std::ostream &file, bool wrap)
{
	Calc_Cube(CubeRho, src, [&](const d3 &pos) { return calculate_density(src, pos); }, radius, file, wrap);
}
template <class S>
void Calc_Eli(cube &CubeEli, const S &src, double radius, std::ostream &file, bool wrap)
{
	Calc_Cube(CubeEli, src, [&](const d3 &pos) { return calculate_eli(src, pos); }, radius, file, wrap);
}
template <PointPotential S>
void Calc_ESP(cube &CubeESP, const S &src, double radius, bool no_date, std::ostream &file, bool wrap)
{
	Calc_Cube(CubeESP, src, [&](const d3 &pos) { return src.esp(pos); }, radius, file, wrap, no_date);
}
//RDG visualisations use 101.0 as the mask outside the calculation radius; with wrap every point is visited once per
//periodic image, so the mask is applied after the loop where rho stayed exactly zero. Rho becomes sign(lambda2) rho
inline void finish_signed_rho(std::vector<cube> &Cubes, const cube &rho_contrib)
{
	Cubes[cube_type::Rho] = rho_contrib;
	cube &rdg = Cubes[cube_type::RDG];
	for (int x = 0; x < rdg.get_size(0); x++)
		for (int y = 0; y < rdg.get_size(1); y++)
			for (int z = 0; z < rdg.get_size(2); z++)
				if (rho_contrib.get_value(x, y, z) == 0.0)
					rdg.set_value(x, y, z, 101.0);
}
//RDG, Laplacian and ELI-D (PC07) of the loaded cubes from the density of src; ELF needs orbitals and exits. The
//grid points within radius go to the batch forms of density_source.h in chunks (with the Hessian only when the RDG
//wants lambda2), so the fitted density runs its kernel on the GPU; with wrap every periodic image is a point
template <class S>
void Calc_Prop(std::vector<cube> &Cubes, const S &src, double radius, std::ostream &file, bool test, bool wrap)
{
	err_checkf(!Cubes[cube_type::Elf].get_loaded(), std::string("ELF needs orbitals, the ") + source_name(src) + " has none", file);
	const bool rdg = Cubes[cube_type::RDG].get_loaded(), lap_c = Cubes[cube_type::Lap].get_loaded(), eli_c = Cubes[cube_type::Eli].get_loaded();
	if (!rdg && !lap_c && !eli_c)
		return;
	_time_point start = get_time();
	const double r = constants::ang2bohr(radius);
	const std::vector<d3> centres = source_positions(src);
	cube rho_contrib(Cubes[cube_type::Rho]);
	rho_contrib.set_zero();
	const i3 n = rho_contrib.get_sizes();
	vec x, y, z, rho, gx, gy, gz, lap, H;
	std::vector<i3> idx;
	auto flush = [&]() {
		const int np = static_cast<int>(x.size());
		rho.resize(np), gx.resize(np), gy.resize(np), gz.resize(np), lap.resize(np);
		if (rdg)
		{
			H.resize(9 * static_cast<size_t>(np));
			calculate_hessian(src, np, x.data(), y.data(), z.data(), rho.data(), gx.data(), gy.data(), gz.data(), H.data());
			for (int p = 0; p < np; p++)
				lap[p] = H[9 * (size_t)p] + H[9 * (size_t)p + 4] + H[9 * (size_t)p + 8];
		}
		else
			calculate_density(src, np, x.data(), y.data(), z.data(), rho.data(), gx.data(), gy.data(), gz.data(), lap.data());
		for (int p = 0; p < np; p++)
		{
			const i3 &i = idx[p];
			const double g2 = gx[p] * gx[p] + gy[p] * gy[p] + gz[p] * gz[p];
			auto add = [&](cube &c, double v) { c.set_value(i[0], i[1], i[2], c.get_value(i[0], i[1], i[2]) + (std::isfinite(v) ? v : 0.0)); };
			if (lap_c)
				add(Cubes[cube_type::Lap], lap[p]);
			if (eli_c)
				add(Cubes[cube_type::Eli], aux_density::eli_from_density(rho[p], g2, lap[p]));
			double v = rho[p];
			if (rdg)
			{
				add(Cubes[cube_type::RDG], v > 0 ? constants::alpha_coef * std::sqrt(g2) / std::pow(v, constants::c_43) : 0.0);
				if (get_lambda_1(H.data() + 9 * (size_t)p) < 0)
					v = -v;
			}
			add(rho_contrib, v);
		}
		x.clear(), y.clear(), z.clear(), idx.clear();
	};
	//ponytail: 2^20 points per chunk = 136 MB of arrays on the host or the device; the gather is serial
	constexpr int chunk = 1 << 20;
	const int lo = wrap ? -1 : 0, hi = wrap ? 2 : 1;
	for (int i = lo * n[0]; i < hi * n[0]; i++)
		for (int j = lo * n[1]; j < hi * n[1]; j++)
			for (int k = lo * n[2]; k < hi * n[2]; k++)
			{
				const d3 pos = rho_contrib.get_pos(i, j, k);
				if (!near_any(pos, centres, r))
					continue;
				x.push_back(pos[0]), y.push_back(pos[1]), z.push_back(pos[2]);
				idx.push_back({ (i + n[0]) % n[0], (j + n[1]) % n[1], (k + n[2]) % n[2] });
				if (static_cast<int>(x.size()) == chunk)
					flush();
			}
	flush();
	if (rdg)
		finish_signed_rho(Cubes, rho_contrib);
	if (!test)
	{
		_time_point end = get_time();
		print_time(start, end, file);
	}
}
