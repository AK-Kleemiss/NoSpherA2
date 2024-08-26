#pragma once
#include <filesystem>

#include "convenience.h"
#include "JKFit.h"

#include "SALTED_math.h"
#include "SALTED_io.h"
#include "SALTED_utilities.h"
#include "SALTED_equicomb.h"
#include "SALTED_math.h"



class SALTEDPredictor {
public:
    SALTEDPredictor(const WFN& wavy, const options &opt);

    const std::string get_dfbasis_name();
    vec gen_SALTED_densities();
    const std::string get_h5_filename() const {
        return config.h5_filename;
    };
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
    std::unordered_map<int, vec> wigner3j{};
    std:: unordered_map<std::string, vec> av_coefs{};
    void read_model_data();
#if has_RAS
    void read_model_data_h5();
#endif

    vec predict();

};
