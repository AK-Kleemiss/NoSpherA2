#pragma once
#ifndef __c_spherical_density__
#define __c_spherical_density__

#include <vector>
#include <iostream>

class Thakkar {
private:
	int atomic_number;
	const int first_ex();
	const int previous_element_coef();
public:
	Thakkar(const int g_atom_number) { atomic_number = g_atom_number; };
	const double get_radial_density(double dist);
	const int get_max_l();
	const double get_max_alpha();
	const std::vector<double> get_min_alpha_vector();
	const double get_min_alpha();
	const double get_form_factor(const double k_vector, std::ofstream &file, bool debug = false);
};



#endif