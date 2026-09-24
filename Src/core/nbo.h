#pragma once
#include "convenience.h"
#include "nao.h"
#include "nbo_run.h"
#include <filesystem>

class WFN;

//Natural bond orbitals, their second-order perturbative interaction energies and natural resonance
//theory, computed in house.  The reference this is measured against is NBO 7 through nbo_run.h, and
//the results are filled into that header's structures so compare_nbo_results() and
//tests/nbo_reference/compare_nbo.py apply unchanged.
//
//The input is a FILE47 - the very archive NoSpherA2 already writes for gennbo - which is what makes
//this a comparison of two analyses rather than of two input preparations: overlap, density and Fock
//matrix are read from the same file the reference run consumed, in the same AO order and the same
//phase convention.  WFN::write_nbo() is the only place that knows how to build them, and it is
//already checked by Tr(P S) against the electron count.
//
//The search follows Foster and Weinhold, JACS 102 (1980) 7211 and the NBO 7 manual: core blocks
//first, then a threshold ladder alternating one-centre (lone pair) and two-centre (bond) blocks of
//the successively depleted density, until the accepted orbitals hold the whole electron pair count.
//NRT follows Glendening and Weinhold, J. Comput. Chem. 19 (1998) 593, as a convex quadratic
//program on the simplex; it lives in nrt.cpp.

//One FILE47, unpacked.  density and fock carry one block for a closed shell and two - alpha then
//beta - for an open one.
struct NboInput {
    int n = 0;
    bool open_shell = false;
    std::vector<NAOBasisFunction> ao;  //atom, l, shell inside (atom, l), component 0..2l
    dMatrix2 overlap;
    std::vector<dMatrix2> density;
    std::vector<dMatrix2> fock;
};

//One accepted Lewis or non-Lewis orbital in the NAO basis.
struct NboFunction {
    ivec centers;                      //0-based atom indices, one or two
    std::string type;                  //CR, LP, BD, BD*, LV, RY
    int multiplicity = 1;              //sigma = 1, the second and third bond of a multiple bond
    double occupancy = 0.0;
    double energy = 0.0;
    vec coefficients;                  //length n_nao, the orbital in the NAO basis
    //per centre polarisation, parallel to centers: |c_A|^2 and the l-resolved share inside it
    vec center_weight;
    vec2 center_lchar;                 //[centre][l], summing to 1 per centre
};

struct NboLewis {
    dMatrix2 gamma;                    //density in the NAO basis, trace = electron count
    std::vector<NboFunction> orbitals; //Lewis set first, then the non-Lewis set
    int n_lewis = 0;                   //how many of them are Lewis (CR, LP, BD)
    double rho_nl = 0.0;               //electrons outside the Lewis set
    double threshold = 0.0;            //occupancy threshold the ladder stopped at
    //bond multiplicities of the leading Lewis structure: topo(a,b) for a != b, lone pairs
    //(cores included) on the diagonal.  This is the starting point of the resonance search.
    ivec2 topo;
};

struct NboOptions {
    double occupancy_threshold = 1.90; //start of the ladder, NBO's default
    double e2_threshold_kcal = 0.5;    //what enters the printed E2 table
    //What enters the resonance search, NBO's NRTE2.  Calibrated, not inherited: ethane's strongest
    //interaction is 3.52 kcal (sigma(C-H) -> RY C) and its hyperconjugative sigma -> sigma* is
    //2.80 kcal, yet NBO finds ethane's six resonance structures, while water's strongest is 1.01 kcal
    //and NBO finds water one structure.  The gate has to sit between the two.
    double nrt_e2_kcal = 2.0;
    int nrt_max_arrows = 2;            //depth of the arrow-driven candidate generation
    //Bond-length screen for resonance candidates, looser than the search's 1.3: a resonance structure
    //may put a bond where the parent has none at all.  Ozone's third-largest reference structure at
    //23.63 % is the ring, O 1- O 3 at 2.24 A against a covalent-radius sum of 1.32 A.
    double nrt_bond_scale = 1.75;
    int nrt_max_candidates = 4000;
    bool nrt_max_set = false;          //-nrt_max given: obey the number, skip the size guard
    double nrt_weight_floor = 5.0e-5;  //weights below this are dropped from the reported set
    bool nrt = false;
    bool nrt_exhaustive = false;       //enumerate every feasible topology instead of arrows
    bool nrt_ion = true;               //allow candidates that move a bond to a lone pair (NRTION)
    bool nrt_symmetry = true;          //collapse candidates related by an automorphism
    bool nrt_components = true;        //split into connected components of the delocalisation graph
    ivec nrt_subspace;                 //1-based atoms the search may alter, empty = all of them
    int threads = 0;                   //0: OpenMP default
    std::filesystem::path file47;      //read this instead of writing a fresh one
    bool keep_file47 = false;
    bool debug = false;
};

/** Read overlap, density, Fock matrix and the AO map out of a FILE47. */
NboInput read_file47(const std::filesystem::path& file);

/** Density in the NAO basis: C^T S P S C, whose trace is the electron count. */
dMatrix2 nao_density(const dMatrix2& P, const dMatrix2& S, const dMatrix2& C);

/** A one-electron operator already given over the AOs, transformed to the NAO basis: C^T F C. */
dMatrix2 nao_operator(const dMatrix2& F, const dMatrix2& C);

/**
 * Atom pairs a two-centre orbital may be searched over: |R_AB| <= scale * (r_cov,A + r_cov,B).
 * Without this screen the search buys electron pairs with bonds between atoms that are not
 * neighbours - two lone pairs on different chlorines of TiCl4 span a two-centre block whose leading
 * eigenvalue is 1.96, which beats the third lone pair's 1.80 and takes its place.  NRT uses the same
 * screen to decide which pairs a candidate resonance structure may put a bond on.
 */
bvec2 bondable_pairs(const std::vector<atom>& atoms, double scale = 1.3);

/**
 * The NBO search on one spin's NAO basis.  n_pairs is the number of orbitals the Lewis set has to
 * fill (electrons / 2 for a closed shell, the spin's electron count for one spin of an open one),
 * and scale is the occupancy a full orbital carries, 2 or 1.
 */
NboLewis nbo_search(const NAOResult& nao, const dMatrix2& gamma, const bvec2& bondable, int n_pairs,
                    double scale, const NboOptions& options);

/**
 * Second-order perturbative donor-acceptor energies over the accepted orbitals.  lewis is not
 * const because the diagonal of the Fock matrix in the NBO basis is where the orbital energies
 * come from - this is the only place they are known - and they are filled in here.
 */
std::vector<NboE2Entry> nbo_e2(NboLewis& lewis, const dMatrix2& fock_nao, double threshold_kcal);

/**
 * Natural resonance theory on one spin's NAO basis, appended to nrt: the candidate resonance
 * structures, the convex quadratic program that weights them, and the bond orders, valencies and
 * topology matrices the converged weights imply.  lewis supplies the parent structure and the NAO
 * density, e2 the delocalisation graph the a-priori screens work on (its donor/acceptor indices must
 * still be 1-based into lewis.orbitals, i.e. it has to be passed before analyse_spin renumbers it),
 * bondable which pairs a candidate may put a bond on, and scale the occupancy of a full orbital.
 *
 * D(w) here is the plain Frobenius norm of Gamma - sum_a w_a Gamma_a; NBO's printed D carries a
 * normalisation this code does not reproduce, so the quantities to compare a reference against are
 * the weights by rank, the bond orders, the valencies and the retained-structure count.
 */
void native_nrt(NboNrt& nrt, const NAOResult& nao, const NboLewis& lewis,
                const std::vector<NboE2Entry>& e2, const bvec2& bondable,
                const NboOptions& options, const std::string& spin, double scale,
                std::ostream& log);

//The whole in-house analysis, in the shape the NBO 7 reference is stored in.  wavy is not const
//because writing the FILE47 is not: WFN::write_nbo() renormalises the stored basis on the way.
NboResults native_nbo(WFN& wavy, const NboOptions& options, std::ostream& log);

/** Print the tables of a native result in NBO's layout. */
void print_nbo(const NboResults& results, std::ostream& out);

/** The resonance tables of a native result, in NBO's layout. Called by print_nbo; separate so a
 *  test can print a hand-built result without an SCF. Prints nothing unless results.nrt.present. */
void print_nrt(const NboResults& results, std::ostream& out);
