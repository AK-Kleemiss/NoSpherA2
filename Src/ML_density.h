#pragma once
#include <vector>
#include <string>
#include "convenience.h"
#include "npy.h"
#include "properties.h"
#include "JKFit.h"

class ML_density {
public:
	static void calc_diff(const std::string xyz_File, options& opt);
	void gbw2DM(std::string& fn, std::ostream& file, bool& debug);
	void cubeDiffDaniel();
private:
	static cube data_2_Cube(const std::string& xyz_File, const std::string& npy_Coeffs, int cpus);
	static void calc_cube(cube& CubeRho, WFN& wavy, int cpus, double radius, int& exp_coef, vec data, std::ostream& file);
};