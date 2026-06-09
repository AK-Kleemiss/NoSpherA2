#pragma once
#include <filesystem>

#include "convenience.h"
#include "SALTED_io.h"
#include "wfn_class.h"

class SALTEDPredictor
{
public:
    SALTEDPredictor(const WFN &wavy, options &opt_in);
    SALTEDPredictor() = default;

    const std::string get_dfbasis_name() const;
    vec gen_SALTED_densities();
    const std::filesystem::path get_salted_filename() const
    {
        return config.salted_filename;
    };
    WFN wavy;
    void shrink_intermediate_vectors();
    const bool basis_set_loaded() const { return bbasis_set_loaded; };

private:
    bool bbasis_set_loaded = false;
    Config config;
    int natoms;
    std::filesystem::path SALTED_DIR;
    std::filesystem::path coef_file;
    bool debug;
    std::vector<std::string> atomic_symbols{};
    std::unordered_map<std::string, std::vector<int>> atom_idx{};
    
    std::unordered_map<std::string, int> natom_dict{}, lmax{}, nmax{};
    cvec4 v1, v2;
    void setup_atomic_environment();

    vec weights{};
    std::unordered_map<std::string, dMatrix2> Vmat{};
    std::unordered_map<std::string, int> Mspe{};
    std::unordered_map<int, std::vector<int64_t>> vfps{};
    std::unordered_map<int, vec> wigner3j{};
    std::unordered_map<std::string, dMatrix2> power_env_sparse{};
    std::unordered_map<std::string, vec> av_coefs{};
    std::unordered_map<int, int> featsize{};
    void read_model_data();

    vec predict();
};

