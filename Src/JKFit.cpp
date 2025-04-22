#include "pch.h"
#include "JKFit.h"

#include "auxiliary_basis.hpp"

// Every BasisSet ist defined for all 118 Elements, basis sets, that do not cover all elements, are padded with 0s
BasisSet::BasisSet(const std::array<std::vector<primitive>, 6>& data) {
	for (int i = 0; i < 6; i++) {
		_data[i] = data[i];
	}
	for (int i = 6; i < 118; i++) {
		_data[i] = {};
	}
}


//// Every BasisSet ist defined for all 118 Elements, basis sets, that do not cover all elements, are padded with 0s
const std::array<std::vector<primitive>, 118>& BasisSet::get_data(){
	return _data;
}

const std::vector<primitive>& BasisSet::operator[](int element) const {
	//Test if the element is empty
	err_checkf(_data[element].size() != 0, "Basis set for element " + std::to_string(element) + " is not defined!", std::cout);
	return _data[element];
}

BasisSet& BasisSetLibrary::get_basis_set(std::string basis_name) {
	//Check if the supplied basis name is contained in part of a given basis set name
	std::string found_basis = "";
	if (basis_name == "") {
		std::cout << "No Basis Name Supplied! Aborting!!!" << std::endl;
		exit(1);
	}
	//Cast basis_name to lowercase
	std::transform(basis_name.begin(), basis_name.end(), basis_name.begin(), ::tolower);

	for (auto const& [key, val] : basisSets) {
		if (key.find(basis_name) != std::string::npos) {
			found_basis = key;
			break;
		}
	}
	err_checkf(found_basis != "", "Basis set " + basis_name + " not defined in BasisSetLibrary!", std::cout);
	std::cout << "Using basis set: " << found_basis << std::endl;
	return basisSets[found_basis];
}

BasisSetLibrary::BasisSetLibrary() {
	basisSets = built_in_basis_sets;
	basisSets["TESTING"] = BasisSet(TESTING);
}

bool BasisSetLibrary::check_basis_set_exists(std::string basis_name) {
	//Check if the supplied basis name is in the basisSets map
	return basisSets.find(basis_name) != basisSets.end();
}

const std::array<std::vector<primitive>, 6> TESTING =
{
	{
		{
			{0, 0, 13.087376, 1.0},
			{0, 0, 1.185515, 1.0},
			{0, 0, 0.368163, 1.0},
			{0, 1, 2.288385, 1.0},
			{0, 1, 1.311828, 1.0},
			{0, 2, 1.875349, 1.0},
		}, //H
		{
			{0, 0, 19.248269, 1.0},
			{0, 0, 3.746353 , 1.0},
			{0, 0, 0.729165, 1.0},
			{0, 1, 4.523992, 1.0},
			{0, 1, 1.648567, 1.0},
			{0, 2, 1.472906, 1.0},
		}, //He
		{
			{0, 0, 1.1, 1.0},
			{0, 1, 1.2, 1.0},
			{0, 2, 1.3, 1.0},
			{0, 3, 1.4, 1.0},
			{0, 4, 1.5, 1.0},
			{0, 5, 1.6, 1.0},
			{0, 6, 1.7, 1.0},
		}, //Li
		{}, //Be
		{}, //B
		{
			{0, 0, 120.498651, 1.0},
			{0, 0, 45.116782, 1.0},
			{0, 0, 16.892505, 1.0},
			{0, 0, 6.324846, 1.0},
			{0, 0, 2.368132, 1.0},
			{0, 0, 0.886670, 1.0},
			{0, 1, 13.216186, 1.0},
			{0, 1, 3.884909, 1.0},
			{0, 1, 1.141972, 1.0},
			{0, 1, 0.335684, 1.0},
			{0, 2, 3.750089, 1.0},
			{0, 2, 1.207332, 1.0},
			{0, 2, 0.388698, 1.0},
			{0, 3, 1.344106, 1.0},
			{0, 4, 0.769479, 1.0},
		},


	}
};
