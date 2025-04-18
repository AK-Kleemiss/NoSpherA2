#pragma once
#include <array>
#include <vector>
#include "pch.h"
#include "convenience.h"


extern const std::array<std::vector<primitive>, 6> TESTING;



extern const std::array<std::vector<primitive>, 118> def2_SVP_JKfit;
extern const std::array<std::vector<primitive>, 118> def2_SVP_RIfit;

extern const std::array<std::vector<primitive>, 118> HGBSP3_7;
extern const std::array<std::vector<primitive>, 118> def2_universal_JKfit;
extern const std::array<std::vector<primitive>, 118> def2_universal_Jfit;
extern const std::array<std::vector<primitive>, 118> def2_qzvppd_rifit;
extern const std::array<std::vector<primitive>, 118> combo_basis_fit; //Combination of cc-pvqz-jkfit, def2-qzvppd-rifit and HGBSP3-5 for the few remaining elements

extern const std::array<std::vector<primitive>, 35> CC_TZVP_JKfit;
extern const std::array<std::vector<primitive>, 35> CC_QZVP_JKfit;
extern const std::array<std::vector<primitive>, 118> CC_PVQZ_F12_OPTRI;
extern const std::array<std::vector<primitive>, 35> CC_PV5Z_JKfit;


// Templated BasisSet class definition

class BasisSet {
public:
    BasisSet(const std::array<std::vector<primitive>, 1>& data);
    BasisSet(const std::array<std::vector<primitive>, 6>& data);
    BasisSet(const std::array<std::vector<primitive>, 35>& data);
    BasisSet(const std::array<std::vector<primitive>, 86>& data);
    BasisSet(const std::array<std::vector<primitive>, 118>& data);
    BasisSet(const std::array<std::vector<primitive>, 118>& data) : _data(data) {};
    //Empty constructor
    BasisSet() {};

    const std::array<std::vector<primitive>, 118>& get_data();
    const std::vector<primitive>& operator[](int element) const;

private:
    std::array<std::vector<primitive>, 118> _data;
};

#include "auxiliary_basis.hpp"


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
