#pragma once
#include "convenience.h"
#include "wfn_class.h"
#include "cell.h"
#include "integrator.h"

//Interaction energies between every pair of molecules in contact in a crystal. The density of each unique molecule is
//predicted or fitted once; a symmetry image is made by moving the atoms and rotating the aux coefficients with them.
namespace crystal_energies {
    struct symop { ivec2 rot; vec trans; };
    struct pair { int A, B, op; ivec n; double distance; vec centroid_A, centroid_B, t; vec2 R; std::string symop; DensityFitting::INTERACTION E; };
    struct job { std::filesystem::path cif, output; double cutoff = 3.8; pathvec structures, coefficients; };
    //c' = D c for a density moved by the orthogonal matrix R, D in the m = -l..l order of the aux basis
    vec2 real_sh_rotation(const int l, const vec2& R);
    //Moves the atoms of aux to R x + t (bohr) and rotates the coefficients with them
    void transform(WFN& aux, vec& coef, const vec2& R, const vec& t);
    std::vector<symop> symops(const cell& c);
    std::string symop_string(const symop& op, const ivec& n);
    //Cartesian cell matrix in bohr, columns a b c
    vec2 cell_matrix(const cell& c);
    vec2 inverse3(const vec2& M);
    //Coefficients of one molecule: read from coef_file, predicted with -SALTED, or an RI fit of the wavefunction with -ri_fit
    vec fitted_coefficients(const WFN& wavy, const std::filesystem::path& coef_file, WFN& aux, options& opt);
    //Every pair of molecules with two atoms closer than cutoff (Angstrom), each pair once with A as given and B = op(B) + n
    std::vector<pair> contacts(const cell& c, const std::vector<WFN>& molecules, const double cutoff);
    job read_job(const std::filesystem::path& file);
    void write_table(const std::vector<pair>& pairs, const pathvec& structures, std::ostream& out);
    int run(options& opt);
}
