#pragma once
#include <array>
#include <vector>
#include "pch.h"
#include "convenience.h"


extern const std::array<std::vector<primitive>, 6> TESTING;

// Templated BasisSet class definition

class BasisSet {
public:
    BasisSet(const std::array<std::vector<primitive>, 118>& data) : _data(data) {};
    BasisSet(const std::array<std::vector<primitive>, 6>& data);
    //Empty constructor
    BasisSet() = default;

    const std::array<std::vector<primitive>, 118>& get_data();
    const std::vector<primitive>& operator[](int element) const;

private:
    std::array<std::vector<primitive>, 118> _data;
};


// BasisSetLibrary class definition
class BasisSetLibrary {
public:
    // Constructor
    BasisSetLibrary();

    //Access basis set
    BasisSet& get_basis_set(std::string basis_name);
	bool check_basis_set_exists(std::string basis_name);

private:
    std::unordered_map<std::string, BasisSet> basisSets;
};
