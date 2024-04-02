#pragma once

#include <string>
#include <fstream>
#include <vector>

class WFN;
class cell;

void Calc_Static_Def(
	cube& CubeDEF,
	cube& CubeRho,
	WFN& wavy,
	int cpus,
	double radius,
	std::ostream& file
);
void Calc_Static_Def(
	cube& CubeDEF,
	cube& CubeRho,
	cube& CubeSpher,
	WFN& wavy,
	int cpus,
	double radius,
	std::ostream& file
);
void Calc_Spherical_Dens(
	cube& CubeSpher,
	WFN& wavy,
	int cpus,
	double radius,
	std::ostream& file
);
void Calc_Rho(
	cube& CubeRho,
	WFN& wavy,
	int cpus,
	double radius,
	std::ostream& file
);
void Calc_Rho_spherical_harmonics(
	cube& CubeRho,
	WFN& wavy,
	int cpus,
	std::ostream& file
);
void Calc_MO_spherical_harmonics(
	cube& CubeRho,
	WFN& wavy,
	int cpus,
	int MO,
	std::ostream& file
);
void Calc_Prop(
	cube& CubeRho,
	cube& CubeRDG,
	cube& CubeElf,
	cube& CubeEli,
	cube& CubeLap,
	cube& CubeESP,
	WFN& wavy,
	int cpus,
	double radius,
	std::ostream& file,
	bool test
);
void Calc_ESP(
	cube& CubeESP,
	WFN& wavy,
	int cpus,
	double radius,
	bool no_date,
	std::ostream& file
);
void Calc_MO(
	cube& CubeMO,
	int mo,
	WFN& wavy,
	int cpus,
	double radius,
	std::ostream& file
);
void Calc_S_Rho(
	cube& Cube_S_Rho,
	WFN& wavy,
	int cpus,
	std::ostream& file,
	bool nodate
);
void Calc_Hirshfeld(
	cube& CubeHDEF,
	cube& CubeRho,
	WFN& wavy,
	int cpus,
	double radius,
	int ignore,
	std::ostream& file
);
void Calc_Hirshfeld(
	cube& CubeHDEF,
	cube& CubeRho,
	cube& CubeSpherical,
	WFN& wavy,
	int cpus,
	double radius,
	int ignore,
	std::ostream& file
);
void Calc_Hirshfeld_atom(
	cube& CubeHirsh,
	cube& CubeRho,
	cube& CubeSpherical,
	WFN& wavy,
	int cpus,
	double radius,
	int ignore,
	std::ostream& file
);

void properties_calculation(options& opt);

void do_combine_mo(options& opt);

#include "wfn_class.h"
#include "cell.h"