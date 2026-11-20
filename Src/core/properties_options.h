#pragma once
#include "nos_types.h"

struct properties_options
{
	bool rho = false;
	bool eli = false;
	bool esp = false;
	bool elf = false;
	bool lap = false;
	bool rdg = false;
	bool hdef = false;
	bool def = false;
	bool hirsh = false;
	bool s_rho = false;
	//ELI-D(aa), ELI-D(bb) and the triplet ELI-D as three cubes from one spin-resolved pass
	bool eli_spin = false;
	bool all_mos = false;
	//Fukui functions f+/f-/f0 and the dual descriptor, frozen-orbital approximation
	bool fukui = false;
	//rho isosurface value (au) coloured by the ESP and written as obj; 0 = off
	double esp_isosurface = 0.0;
	double resolution = 0.1;
	double radius = 2.0;
	double integral_accuracy = -1;
	double promol_nci_rcut1 = 0.95;
	double promol_nci_rcut2 = 0.75;
	double promol_nci_rho_abs_max = 0.5;
	double promol_nci_rdg_max = 1.0;
	double promol_nci_colour_max = 0.015; // VMD/Olex2 colour range on sign(l2)rho in a.u., symmetric
	double promol_nci_iso = 0.5; // RDG of the _nci.obj surface
	//The _values.dat writer schedules grid points dynamically, so row order (not the values) varies between runs; forces single-threaded for golden files
	bool promol_nci_single_threaded = false;
	std::array<int, 3> NbSteps = { 0, 0, 0 };
	std::array<double, 6> MinMax = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	ivec MO_numbers;
	//-ibo_cube <list>: IBOs to write as cubes (ibo_selection syntax)
	std::string ibo_cube;
	//-nbo_cube / -nlmo_cube <list>: in-house NBOs / NLMOs to write as cubes (nbo_cube_wfn syntax)
	std::string nbo_cube, nlmo_cube;
	//-nbo_lewis <json>: the Lewis structure ($CHOOSE) those NBOs come from, e.g. an NRT run's leading structure
	std::filesystem::path nbo_lewis;
	int hirsh_number = 0;
	bool calc() const {
		return rho || eli || eli_spin || esp || elf || lap || rdg || hdef || def || hirsh || s_rho || all_mos || fukui || esp_isosurface > 0 || MO_numbers.size() > 0 || !ibo_cube.empty() || !nbo_cube.empty() || !nlmo_cube.empty();
	}
	size_t n_grid_points() const {
		size_t result = static_cast<size_t>(NbSteps[0]) * NbSteps[1] * NbSteps[2];
		return result;
	}
};
