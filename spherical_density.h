#pragma once
#ifndef __c_spherical_density__
#define __c_spherical_density__

#include <vector>

const double PI = 3.14159265358979323846;  /* pi */

class Thakkar {
private:
	int atomic_number;
	int first_ex();
	int previous_element_coef();
public:
	Thakkar(int g_atom_number) { atomic_number = g_atom_number; };
	double get_radial_density(double dist);
	int get_max_l();
	double get_max_alpha();
	std::vector<double> get_min_alpha_vector();
	double get_min_alpha();
};



#endif