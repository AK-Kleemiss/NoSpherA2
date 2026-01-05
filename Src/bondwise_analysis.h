#pragma once

#include <string.h>
#include <vector>

class WFN;
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
        std::pair<int, int> atom_indices;
    };

    std::vector<NAOResult> NAOs;
    dMatrix2 overlap_matrix;
    dMatrix2 density_matrix;
    dMatrix2 total_NAOs;
    std::vector<dMatrix2> projection_matrices;
    std::vector<dMatrix2> overlap_matrices;
    std::vector<bond_index_result> RGBI;
    NAOResult calculateAtomicNAO(const dMatrix2& D_full, const dMatrix2& S_full, const std::vector<int>& atom_indices);
    double projection_matrix_and_expectation(const ivec& indices, const ivec& eigvals = {}, const ivec& eigvecs = {}, dMatrix2* given_NAO = nullptr);
    void computeAllAtomicNAOs(WFN& wavy);
    ivec find_eigenvalue_pairs(const vec& eigvals, const double tolerance = 1E-4);
    void transform_Ionic_eigenvectors_to_Ionic_orbitals(dMatrix2& EVC,
        const vec& eigvals,
        const ivec& pairs,
        const int index_a,
        const int index_b,
        const ivec& pair_matrix_indices);
    std::map<char, dMatrix2> make_covalent_from_ionic(
        const dMatrix2& theta_I,
        const vec& eigvals,
        const ivec& pairs);
    double Roby_population_analysis(ivec atoms);
public:
    Roby_information() = default;
    ~Roby_information() = default;
    Roby_information(const Roby_information&) = default;
    Roby_information(WFN& wavy);




};

bond do_bonds(WFN& wavy, int mode_general, int mode_sel, bool mode_leng, bool mode_res, double res[], bool cub, double boxsize[], int atom1, int atom2, int atom3, const bool& debug, const bool& bohr, int runnumber, bool rho, bool rdg, bool eli, bool lap);
int autobonds(bool debug, WFN& wavy, const std::filesystem::path& inputfile, const bool& bohr);

#include "wfn_class.h"
