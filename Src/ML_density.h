#pragma once
#include <vector>
#include <string>
#include "convenience.h"
#include "npy.h"


class ML {
public:
	static void calc_diff(const std::string xyz_File);
	void gbw2DM(std::string& fn, std::ostream& file, bool& debug);
	void cubeDiffDaniel();
private:
	static cube data_2_Cube(const std::string& xyz_File, const std::string& npy_Coeffs);
	static cube calc_cube(vec data, WFN& dummy, int& exp_coef);
};