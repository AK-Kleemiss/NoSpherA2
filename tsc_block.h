#pragma once

#ifndef __TSC_BLOCK_H__
#define __TSC_BLOCK_H__

#include <iostream>
#include <vector>
#include <string>
#include <complex>

class tsc_block {
private:
	std::vector<std::vector<std::complex<double> > > sf; //[#scatterer] [#reflection]
	std::vector<std::string> scatterer; //Labels of reflections the correpsonding entry of sf belonds to
	std::vector<std::vector<int> > index; //[3] [index]

public:
	std::vector<std::complex<double> > get_sf_for_scatterer(const int nr) { if (nr >= 0 && nr < scatterer.size()) return sf[nr]; else return {}; };
	std::string get_scatterer(const int nr) { if (nr >= 0 && nr < scatterer.size()) return scatterer[nr]; else return "ERROR"; };
	int get_number_of_scatterers() { return scatterer.size(); };
	tsc_block(
		std::vector<std::vector<std::complex<double> > > given_sf,
		std::vector<std::string> given_scatterer, 
		std::vector<std::vector<int> > given_index) {
		sf = given_sf;
		scatterer = given_scatterer;
		index = given_index;
	};
	tsc_block() {};
	~tsc_block() {
		sf.clear();
		scatterer.clear();
		index.clear();
	};
};


#endif