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
	std::ofstream& file
);
void Calc_Static_Def(
	cube& CubeDEF,
	cube& CubeRho,
	cube& CubeSpher,
	WFN& wavy,
	int cpus,
	double radius,
	std::ofstream& file
);
void Calc_Spherical_Dens(
	cube& CubeSpher,
	WFN& wavy,
	int cpus,
	double radius,
	std::ofstream& file
);
void Calc_Rho(
	cube& CubeRho,
	WFN& wavy,
	int cpus,
	double radius,
	std::ofstream& file
);
void Calc_Rho_spherical_harmonics(
	cube& CubeRho,
	WFN& wavy,
	int cpus,
	std::ofstream& file
);
void Calc_MO_spherical_harmonics(
	cube& CubeRho,
	WFN& wavy,
	int cpus,
	int MO,
	std::ofstream& file
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
	std::ofstream& file,
	bool test
);
void Calc_ESP(
	cube& CubeESP,
	WFN& wavy,
	int cpus,
	double radius,
	bool no_date,
	std::ofstream& file
);
void Calc_MO(
	cube& CubeMO,
	int mo,
	WFN& wavy,
	int cpus,
	double radius,
	std::ofstream& file
);
void Calc_Hirshfeld(
	cube& CubeHDEF,
	cube& CubeRho,
	WFN& wavy,
	int cpus,
	double radius,
	int ignore,
	std::ofstream& file
);
void Calc_Hirshfeld(
	cube& CubeHDEF,
	cube& CubeRho,
	cube& CubeSpherical,
	WFN& wavy,
	int cpus,
	double radius,
	int ignore,
	std::ofstream& file
);
void Calc_Hirshfeld_atom(
	cube& CubeHirsh,
	cube& CubeRho,
	cube& CubeSpherical,
	WFN& wavy,
	int cpus,
	double radius,
	int ignore,
	std::ofstream& file
);

#include "wfn_class.h"
#include "cell.h"