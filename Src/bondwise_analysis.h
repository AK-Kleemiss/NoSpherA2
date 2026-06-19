#pragma once

#include <string.h>
#include <vector>

class WFN;

// Average an atom-centred AO matrix over the 48 operations of O_h.
// Shells must be contiguous. Cartesian matrices use constants::type_vector;
// spherical matrices use the libcint real-spherical ordering and phases.
void symmetrize_atomic_matrix_oh(dMatrix2& matrix, const ivec& shell_angular_momenta,
    bool spherical = false);

struct bond {
    std::string label_1;
    std::string label_2;
    std::string label_3;
    std::string filename;
    bool success;
    bool dens;
    bool esp;
};

//structrue to contain multiple NAO atom sets and the overlap and density matrices
class Roby_information {
private:
    // Structure to hold the NAOs
    struct NAOResult {
        vec eigenvalues;  // Occupancies
        vec eigenvectors; // Coefficients
        int atom_index = -1;
        ivec matrix_elements;
        vec sub_DM;
        vec sub_OM;
    };

    struct bond_index_result {
        double covalent;
        double ionic;
        double total;
        double percent_covalent_Pyth;
        double percent_covalent_Arakai;
        double pair_population;
        double population_first;
        double population_second;
        std::pair<int, int> atom_indices;
        std::pair<int, int> atom_element_nr;
    };

    struct group_bond_index_result {
        double covalent;
        double ionic;
        double total;
        double percent_covalent_Pyth;
        double percent_covalent_Arakai;
        double pair_population;
        double population_first;
        double population_second;
        int group_index_first;
        int group_index_second;
        ivec atoms_first;
        ivec atoms_second;
    };

    std::vector<NAOResult> NAOs;
    dMatrix2 overlap_matrix;
    dMatrix2 density_matrix;
    dMatrix2 total_NAOs;
    std::vector<dMatrix2> projection_matrices;
    std::vector<dMatrix2> overlap_matrices;
    std::vector<bond_index_result> RGBI;
    std::vector<group_bond_index_result> RGBI_groups;
    NAOResult calculateAtomicNAO(const dMatrix2& D_full, const dMatrix2& S_full,
        const std::vector<int>& atom_indices, const ivec& shell_angular_momenta = {},
        bool spherical = false);
    double projection_matrix_and_expectation(const ivec& indices, const ivec& eigvals = {}, const ivec& eigvecs = {}, dMatrix2* given_NAO = nullptr, dMatrix2* proj_out = nullptr);
    void computeAllAtomicNAOs(WFN& wavy, bool symmetrize);
    ivec find_eigenvalue_pairs(const vec& eigvals, const double tolerance = 1E-4);
    void transform_Ionic_eigenvectors_to_Ionic_orbitals(dMatrix2& EVC,
        const vec& eigvals,
        const ivec& pairs,
        const int index_a,
        const int index_b,
        const ivec& pair_matrix_indices);
    dMatrix2 build_group_PAS(const ivec& group_atom_indices, const ivec& bond_bf_indices, const dMatrix2& P_G);
    void transform_group_Ionic_orbitals(dMatrix2& EVC, const vec& eigvals,
        const ivec& pairs, const ivec& group_a_atoms, const ivec& group_b_atoms,
        const ivec& bond_bf_indices, const dMatrix2& P_GA, const dMatrix2& P_GB);
    void computeGroupAnalysis(const ivec2& group_defs, const vec& atom_pops, const ivec& atom_charges);
    std::map<char, dMatrix2> make_covalent_from_ionic(
        const dMatrix2& theta_I,
        const vec& eigvals,
        const ivec& pairs);
    double Roby_population_analysis(ivec atoms);
public:
    Roby_information() = default;
    ~Roby_information() = default;
    Roby_information(const Roby_information&) = default;
    Roby_information(WFN& wavy, const ivec3& group_sets = {}, bool symmetrize = true);




};

bond do_bonds(WFN& wavy, int mode_general, int mode_sel, bool mode_leng, bool mode_res, double res[], bool cub, double boxsize[], int atom1, int atom2, int atom3, const bool& debug, const bool& bohr, int runnumber, bool rho, bool rdg, bool eli, bool lap);
int autobonds(bool debug, WFN& wavy, const std::filesystem::path& inputfile, const bool& bohr);

void bondwise_laplacian_plots(std::filesystem::path &wfn_name);
void ELI_analysis(const WFN &wavy, const options &opt);

#include "wfn_class.h"
