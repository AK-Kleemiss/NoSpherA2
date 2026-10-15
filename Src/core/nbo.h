#pragma once
#include "convenience.h"
#include "nao.h"
#include "nbo_run.h"
#include <filesystem>

class WFN;

//Natural bond orbitals, second-order E2 and natural resonance theory, filled into nbo_run.h's structures
//so compare_nbo_results() applies unchanged.  The input is the FILE47 written for gennbo, so overlap,
//density and Fock share the reference run's AO order and phase convention.
//Search: Foster & Weinhold, JACS 102 (1980) 7211.  NRT: Glendening & Weinhold, J. Comput. Chem. 19
//(1998) 593, a convex quadratic program on the simplex, in nrt.cpp.

//density and fock: one block for a closed shell, alpha then beta for an open one
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
	//bond multiplicities of the leading Lewis structure, lone pairs (cores included) on the diagonal
	ivec2 topo;
};

struct NboOptions {
	double occupancy_threshold = 1.90; //start of the ladder, NBO's default
	double e2_threshold_kcal = 0.5;    //what enters the printed E2 table
	//NBO's NRTE2: between ethane's sigma -> sigma* (2.8 kcal, six NBO structures) and water's strongest
	//E2 (1.0 kcal, one structure)
	double nrt_e2_kcal = 2.0;
	int nrt_max_arrows = 2;            //depth of the arrow-driven candidate generation
	//looser than the search's 1.3: a resonance structure may bond atoms the parent does not, e.g.
	//ozone's O1-O3 ring bond at 2.24 A against a radius sum of 1.32 A
	double nrt_bond_scale = 1.75;
	int nrt_max_candidates = 4000;
	bool nrt_max_set = false;          //-nrt_max given: obey the number, skip the size guard
	double nrt_target = -1.0;          //-nrt_target: grow from weighted structures until D(w) <= it; < 0 = one search
	double nrt_screen = -1.0;          //-nrt_screen: a walk of nrt_screen_min+ structures builds those whose approximate reduced cost is below it; < 0 = all
	int nrt_screen_min = 200;
	int nrt_stall = 0;                 //-nrt_stall: a build stops once mixing has made no new 0.1 % low of the largest change for this many sweeps; 0 = off
	bool nrt_mirror = true;            //off with -nrt_no_mirror: one-move arrows that separate no charge, orbits completed against the candidates
	bool nrt_rotate = false;           //-nrt_rotate: an image of a converged structure starts from its orbitals, rotated; moves bond orders by up to 2e-4
	double nrt_weight_floor = 5.0e-5;  //weights below this are dropped from the reported set
	bool nrt = false;
	bool nrt_exhaustive = false;       //enumerate every feasible topology instead of arrows
	bool nrt_ion = true;               //allow candidates that move a bond to a lone pair (NRTION)
	bool nrt_symmetry = true;          //collapse candidates related by an automorphism
	bool nrt_components = true;        //split into connected components of the delocalisation graph
	ivec nrt_subspace;                 //1-based atoms the search may alter, empty = all of them
	//NBO's $CHOOSE: the search places exactly this Lewis structure, valence lone pairs on the diagonal and
	//bond orders off it; empty = the search picks.  native_nbo fills it from nrt_parent.
	ivec2 parent_lewis;
	std::filesystem::path nrt_parent;  //-nrt_parent <file.json>, read by read_lewis_structure
	int threads = 0;                   //0: OpenMP default
	int search_threads = 0;            //Lewis search only, 0: threads
									   //-fba runs it serial: its regions are too fine to share a loaded machine
	std::filesystem::path file47;      //read this instead of writing a fresh one
	bool keep_file47 = false;
	bool debug = false;
};

/** Read overlap, density, Fock matrix and the AO map out of a FILE47. */
NboInput read_file47(const std::filesystem::path& file);

/** Density in the NAO basis: C^T S P S C, whose trace is the electron count. */
dMatrix2 nao_density(const dMatrix2& P, const dMatrix2& S, const dMatrix2& C);

/** A one-electron AO operator in the NAO basis: C^T F C. */
dMatrix2 nao_operator(const dMatrix2& F, const dMatrix2& C);

/**
 * Atom pairs a two-centre orbital may be searched over: |R_AB| <= scale * (r_cov,A + r_cov,B).
 * Without it two-centre blocks between non-neighbours outbid genuine lone pairs (TiCl4).
 * NRT uses it for the bonds a candidate may place.  A pair with a third atom on the line between them
 * (angle above 150 degrees) is never bondable.
 */
bvec2 bondable_pairs(const std::vector<atom>& atoms, double scale = 1.3);

/**
 * NBO search on one spin's NAO basis.  n_pairs: orbitals the Lewis set fills (electrons / 2 closed
 * shell, the spin's electron count open shell); scale: occupancy of a full orbital, 2 or 1.
 */
NboLewis nbo_search(const NAOResult& nao, const dMatrix2& gamma, const bvec2& bondable, int n_pairs,
					double scale, const NboOptions& options, std::ostream& log = std::cout);

/**
 * -nrt_parent: read a "NoSpherA2 Lewis structure" JSON, or the "parent_lewis" object of a .native.nbo.json,
 * into topo (valence lone pairs on the diagonal, bond orders off it).  n_pairs counts the core pairs too;
 * wavefunction is the file name the structure is checked against.  Returns "" or the one-line reason it is
 * refused; warnings go to log.
 */
std::string read_lewis_structure(const std::filesystem::path& file, const NAOResult& nao, int n_pairs,
								 const std::string& wavefunction, ivec2& topo, std::ostream& log);

/** Second-order donor-acceptor energies.  Fills the orbital energies, the NBO-basis Fock diagonal, into lewis. */
std::vector<NboE2Entry> nbo_e2(NboLewis& lewis, const dMatrix2& fock_nao, double threshold_kcal);

/**
 * Natural resonance theory on one spin's NAO basis, appended to nrt.  e2's donor/acceptor indices must
 * still be 1-based into lewis.orbitals, i.e. passed before analyse_spin renumbers them.
 * D(w) is the plain Frobenius norm of Gamma - sum_a w_a Gamma_a, summed in squares over the spins of an open
 * shell; NBO 7's D, sqrt((D^2 / scale^2 + 0.001) / n_NAO) per spin, goes to nrt.spin_d.  With xyz (bohr, one row
 * per atom) the symmetry images are those of the geometry and an image's orbitals are its source's rotated,
 * without them the bondable graph's images are built one by one.
 */
void native_nrt(NboNrt& nrt, const NAOResult& nao, const NboLewis& lewis,
				const std::vector<NboE2Entry>& e2, const bvec2& bondable,
				const NboOptions& options, const std::string& spin, double scale,
				std::ostream& log, const vec2* xyz = nullptr);

//Exits on a wavefunction WFN::write_nbo() cannot turn into a FILE47: no contracted basis (.wfn/.wfx) or
//Cartesian d/f shells. -fba calls it before its NBO thread starts, where the exit would tear down RGBI.
void refuse_unsupported_nbo_source(const WFN& wavy, std::ostream& log);

//wavy is not const: WFN::write_nbo() renormalises the stored basis.
NboResults native_nbo(WFN& wavy, const NboOptions& options, std::ostream& log);

/** Print the tables of a native result in NoSpherA2's table layout; the JSON keeps gennbo's labels. */
void print_nbo(const NboResults& results, std::ostream& out);

/** The resonance tables alone, for a hand-built result; prints nothing unless results.nrt.present. */
void print_nrt(const NboResults& results, std::ostream& out);
