#pragma once

#include <string.h>

class WFN;
struct bond {
	char label_1[10];
	char label_2[10];
	char label_3[10];
	bool success;
	bool dens;
	bool esp;
	std::string filename;
};

bond do_bonds(WFN &wavy, int mode_general, int mode_sel, bool mode_leng, bool mode_res, double res[], bool cub, double boxsize[], int atom1, int atom2, int atom3, const bool& debug, const bool& bohr, int runnumber, bool rho, bool rdg, bool eli, bool lap);
int autobonds(bool debug, WFN& wavy, const std::filesystem::path& inputfile, const bool& bohr);

#include "wfn_class.h"
