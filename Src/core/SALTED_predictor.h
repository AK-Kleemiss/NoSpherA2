#pragma once
#include <filesystem>

#include "convenience.h"
#include "SALTED_io.h"
#include "SALTED_utilities.h"
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
    // Set when atoms were moved to the spherical Thakkar fill. The charge
    // constraint needs it: with a mixed ML/Thakkar system the split of a net
    // charge between the two regions is undefined.
    bool spherical_fill_used = false;
    // EEQ estimate of the charge sitting on the spherically filled atoms, and
    // how many there are. Used to say how wrong the neutral-fill assumption is.
    double filled_eeq_charge = 0.0;
    // The part of that charge the spherical fill can actually carry (an ion
    // must be tabulated for the element). The ML target is shifted by exactly
    // this, so predicted + filled still sums to the right number of electrons.
    double applied_fill_charge = 0.0;
    int n_filled = 0;
    Config config;
    int natoms;
    std::filesystem::path SALTED_DIR;
    std::filesystem::path coef_file;
    bool debug;
    std::vector<std::string> atomic_symbols{};
    std::unordered_map<std::string, std::vector<int>> atom_idx{};
    
    std::unordered_map<std::string, int> natom_dict{}, lmax{}, nmax{};
    SALTEDDescriptors v1, v2;
    // Both hyperparameter sets identical: v2 is conj(v1), never filled, equicomb conjugates on read
    bool v2_is_conj_of_v1 = false;
    // Lambda blocks held at once; 0 means all
    int lam_group_limit = 0;
    void setup_atomic_environment();

    vec weights{};
    std::unordered_map<std::string, dMatrix2> Vmat{};
    // Projector shapes for every species the model knows, present or not: the flat
    // weight vector is laid out over all of them, so absent widths shift the offsets
    std::unordered_map<std::string, std::array<size_t, 2>> proj_dims{};
    // Held open for the whole prediction, indexed by (species, lambda)
    std::unique_ptr<SALTED_BINARY_FILE> model_file{};
    std::unordered_map<std::string, SALTED_BINARY_FILE::block_ref> feat_index{};
    std::unordered_map<std::string, SALTED_BINARY_FILE::block_ref> proj_index{};
    std::unordered_set<std::string> model_species{};
    std::unordered_map<std::string, int> Mspe{};
    std::unordered_map<int, std::vector<int64_t>> vfps{};
    std::unordered_map<int, vec> wigner3j{};
    std::unordered_map<std::string, dMatrix2> power_env_sparse{};
    std::unordered_map<std::string, vec> av_coefs{};
    std::unordered_map<int, int> featsize{};
    void read_model_data();
    // Fetch / drop the model matrices of one lambda; read_model_data() only indexes
    void load_model_lambda(const int lam);
    void free_model_lambda(const int lam);

    vec predict();
};

