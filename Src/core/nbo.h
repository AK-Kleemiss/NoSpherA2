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
	std::vector<dMatrix2> lcaomo;      //$LCAOMO per spin block, row k = MO k in the AO basis, zero-padded to n rows
};

//One accepted Lewis or non-Lewis orbital in the NAO basis.
struct NboFunction {
	ivec centers;                      //0-based atom indices, one to max_centres
	std::string type;                  //CR, LP, BD, BD*, LV, RY; kC, kCn, kC* for k = 3, 4 (-nbo_centres)
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
	std::vector<ivec> multi;           //one sorted atom list per k-centre (k >= 3) Lewis orbital, outside topo
	ivec printed;                      //printed NBO number k + 1 is orbitals[printed[k]]
};

//-nbo_cube / -nlmo_cube: what native_nbo leaves behind for orbital pictures, closed shell only
struct NboOrbitalSpace {
	NAOResult nao;                     //nao.C: AO x NAO in the FILE47 AO order
	dMatrix2 overlap;
	dMatrix2 lcaomo;
	NboLewis lewis;
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
	//-nrt_tie: the fit pins the residual, not the weights (Gamma_parent + Gamma_ab = Gamma_a + Gamma_b for fixed
	//orbitals). Among weights whose D(w)^2 lies within this fraction of the optimum take the most parent weight,
	//then the most on single changes, then doubles, ... 0 = off (the first optimum the minimiser reaches)
	double nrt_tie = 1.0e-6;
	//-nrt_ref: multi-reference NRT. After a fit every structure of at least this fraction of the leading weight
	//becomes a reference, its arrow walk joins the pool and the fit repeats, until no new reference appears (at most
	//8 cycles); -nrt_tie then ranks by pair moves to the nearest reference. A parent the fit leaves light (benzene's
	//search puts a lone pair on C1, the Kekule structures carry the weight) no longer decides the pool: the leading
	//structure nearest the hybrid re-anchors the walk's charge screen each cycle; the leading one (ties to the nearest
	//the hybrid) is reported as the parent. A walk of nrt_screen_min+
	//structures builds those of approximate reduced cost below max(nrt_screen, 0) first, the rest stay deferred and
	//are re-priced after each refit (-nrt_screen off builds all). 0 = off
	double nrt_ref = 0.5;
	bool nrt = false;
	bool nrt_exhaustive = false;       //enumerate every feasible topology instead of arrows
	bool nrt_ion = true;               //allow candidates that move a bond to a lone pair (NRTION)
	bool nrt_symmetry = true;          //collapse candidates related by an automorphism
	bool nrt_components = true;        //split into connected components of the delocalisation graph
	ivec nrt_subspace;                 //1-based atoms the search may alter, empty = all of them
	//NBO's $CHOOSE: the search places exactly this Lewis structure, valence lone pairs on the diagonal and
	//bond orders off it; empty = the search picks.  native_nbo fills it from nrt_parent.
	ivec2 parent_lewis;
	std::vector<ivec> parent_multi;    //its k-centre orbitals, sorted atom lists, one per pair
	//-nbo_centres k: the ladder also places k'-centre two-electron bonds, 3 <= k' <= k (at most 4), over connected
	//sets of bondable atoms, once one- and two-centre blocks fail at a threshold; kept only if every centre carries
	//10 % of the orbital and the search without that unit loses 0.05 e of rho(NL); with -nrt_atoms every centre
	//lies in that subspace. 2 = off
	int max_centres = 2;
	//-nbo_tm_pairs n: a d-block metal holds at most n Lewis pairs (lone pairs + bonds); 0 = its valence NAO count, 6 for
	//ns + (n-1)d, NBO 7's footing. 9, the default, is the 18-electron space: under 6 Ni(CO)4 Ni-C is 0.64 and Fe(CO)5
	//0.81/0.68, back-donation invisible; under 9 1.37 and 1.29/1.15 with a lower D(w) in 3 of 3 (Pt(PH3)2ZnCl2 Pt-P
	//0.63 -> 1.03, D(w) 1.69 -> 1.50). The relax pass lifts the other atoms' caps first; a metal's only when nothing
	//else fills the pairs, and the log says so. Before this the relax pass lifted every cap at once, so Cr(CO)6 got 9
	//pairs, Pt(PH3)2AlCl3 7 and Ni(CO)4 6, M-L moving 2x (nrt_tm_reference validation, 24 Oct)
	int tm_pairs = 9;
	std::vector<ivec> skip_multi;      //k-centre sets the ladder may not place (the per-unit gate's leave-one-out)
	//-nrt_e2_avg <percent>: E2 belongs to one Lewis structure; rerun search + E2 under $CHOOSE for every NRT
	//structure degenerate with the leading one (weight within this many percent of it) and average the channels
	//with the renormalised weights. 0 = off. A plain weight cut lets in high-rho(NL) structures whose E2 is
	//meaningless (formate: 650 kcal/mol from a 1.9 e rho(NL) structure).
	double nrt_e2_avg_percent = 0.0;
	//-nrt_e2_avg_min <percent>: absolute weight floor instead of the window, to reach structures that change the
	//metal core when the near-degenerate ones differ only in remote ligand resonance. 0 = off.
	double nrt_e2_avg_min_percent = 0.0;
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
 * refused; warnings go to log.  multi: where "multi_centre" / "three_centre" units go; null refuses them.
 */
std::string read_lewis_structure(const std::filesystem::path& file, const NAOResult& nao, int n_pairs,
								 const std::string& wavefunction, ivec2& topo, std::ostream& log,
								 std::vector<ivec>* multi = nullptr);

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
NboResults native_nbo(WFN& wavy, const NboOptions& options, std::ostream& log, NboOrbitalSpace* space = nullptr);

/**
 * -nbo_cube / -nlmo_cube: a copy of wavy whose MO slot k holds selected orbital k as primitive coefficients, for
 * Calc_MO.  spec: NBO numbers as print_nbo prints them, "lewis", or labels with ':' for the blank ("LP:Pd1",
 * "BD*:Pd1-P6"), comma-separated.  An orbital reaches the primitives through the MOs, t = C_MO S x, so wavy must be
 * the closed-shell wavefunction the FILE47 was written from, with every MO; the lost norm 1 - |t|^2 is logged.
 * NLMOs are the occupied orbitals of maximum overlap with their Lewis NBOs: the NBOs projected on the occupied
 * space (Gamma / 2) and Loewdin-orthogonalised, the closed-shell limit of Reed and Weinhold's NLMOs
 * (J. Chem. Phys. 83 (1985) 1736).  picked[k] is the printed number and label of slot k.
 */
WFN nbo_cube_wfn(const WFN& wavy, const NboResults& r, const NboOrbitalSpace& space, const std::string& spec,
				 bool nlmo, std::vector<std::pair<int, std::string>>& picked, std::ostream& log);

/** Print the tables of a native result in NoSpherA2's table layout; the JSON keeps gennbo's labels. */
void print_nbo(const NboResults& results, std::ostream& out);

/** The resonance tables alone, for a hand-built result; prints nothing unless results.nrt.present. */
void print_nrt(const NboResults& results, std::ostream& out);
