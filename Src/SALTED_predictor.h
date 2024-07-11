#pragma once

#include "convenience.h"
#include "JKFit.h"

#include "SALTED_math.h"
#include "SALTED_io.h"
#include "SALTED_utilities.h"
#include "SALTED_equicomb.h"


class SALTEDPredictor {
public:
    SALTEDPredictor(const WFN& wavy, const options &opt);

    vec gen_SALTED_densities(time_point& start, time_point& end_SALTED);
private:
    Config config;
    const WFN& wavy;
    const options& opt;

    int natoms;

    std::vector<std::string> atomic_symbols{};
    std::unordered_map<std::string, std::vector<int>> atom_idx{};
    std::unordered_map<std::string, int> natom_dict{}, lmax{}, nmax{};
    cvec4 v1, v2;
    void setup_atomic_environment();
    
    vec weights{};
    std::unordered_map<std::string, vec2> Vmat{}, power_env_sparse{};
    std::unordered_map<std::string, int> Mspe{};
    std::unordered_map<int, std::vector<int64_t>> vfps{};
    void read_model_data();

    vec predict();

};
