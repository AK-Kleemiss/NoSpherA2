#pragma once

#include <filesystem>
#include <string>
#include <vector>

//Drives the external NBO 7 (FILE47 -> gennbo) and parses its output: the reference an in-house
//NBO implementation has to reproduce, checked with compare_nbo_results().

struct NboAtomPopulation {
	std::string element;
	int index = 0;         //1-based, as NBO prints it
	double charge = 0.0;
	double core = 0.0;
	double valence = 0.0;
	double rydberg = 0.0;
	double total = 0.0;
	double spin_density = 0.0;  //only filled for an open-shell run
	bool has_spin_density = false;
};

struct NboNao {
	int index = 0;
	std::string element;
	int center = 0;
	std::string lang;       //s, px, dxy, ...
	std::string type;       //Cor, Val, Ryd
	std::string shell;      //"1s", "2p", ...
	double occupancy = 0.0;
	//An open-shell run prints Occupancy/Spin in the spin-summed table, so the energy comes from the alpha table.
	double energy = 0.0;
	double spin_density = 0.0;
	bool has_spin_density = false;
};

//One hybrid of an NBO: its polarisation weight, coefficient and %s/%p/%d/%f content.
struct NboHybrid {
	std::string element;
	int center = 0;
	double weight_percent = 0.0;   //polarisation, empty (0) for one-centre NBOs
	double coefficient = 0.0;
	double s = 0.0;
	double p = 0.0;
	double d = 0.0;
	double f = 0.0;
	double sp_exponent() const { return s >= 1.0 && p >= 1.0 && s + p >= 95.0 ? p / s : 0.0; }
};

struct NboOrbital {
	int index = 0;
	std::string type;           //CR, LP, BD, BD*, LV, RY, RY*, ...
	std::string description;    //"BD ( 1) C  2- C  5", verbatim
	std::string spin;           //"", "alpha", "beta"
	std::vector<int> centers;   //1-based atom numbers
	double occupancy = 0.0;
	double energy = 0.0;
	std::vector<NboHybrid> hybrids;
};

struct NboE2Entry {
	int donor_index = 0;
	int acceptor_index = 0;
	std::string donor;
	std::string acceptor;
	std::string spin;
	double energy_kcal = 0.0;   //E(2)
	double e_diff = 0.0;        //E(NL)-E(L), a.u.
	double fij = 0.0;           //F(L,NL), a.u.
};

struct NboResonanceWeight {
	int structure = 0;              //the structure the weight vector names, not the table's RS column
	int rank = 0;                   //position in the printed table, which is sorted by weight
	double weight_percent = 0.0;    //the 2-decimal table value
	double weight_fraction = 0.0;   //the 5-decimal value of the NRTDTL weight vector, 0 without NRTDTL
	int idxres = 0;                 //NBO's internal structure id, only under NRTDTL
	std::string changes;            //"Added(Removed)" column, verbatim
	std::string spin;
};

struct NboBondOrder {
	int atom1 = 0;
	int atom2 = 0;
	double total = 0.0;
	double covalent = 0.0;
	double ionic = 0.0;
	//Upper triangle and diagonal of the symmetric table only; on the diagonal covalent and ionic print "---" and stay 0.
	bool diagonal = false;
	std::string spin;           //"", "alpha", "beta", "composite"
};

//What NRT derives from the converged weights per atom; the valency is the atom's bond-order sum.
struct NboValency {
	int atom = 0;
	std::string element;
	double valency = 0.0;
	double covalency = 0.0;
	double electrovalency = 0.0;
	double electron_count = 0.0;
	std::string spin;
};

//One line of the NRT search table: how the candidate set and the objective function moved.
struct NboNrtCycle {
	int cycle = 0;
	int structures_used = 0;
	int structures_found = 0;
	double d_w = 0.0;
	int kmax = 0, choose = 0, ion = 0, e2 = 0, sym = 0;
	double dbmax = 0.0, dbrms = 0.0;
	std::string spin;
};

//Integer bond-topology matrix, as printed.
struct NboTopo {
	std::string spin;
	std::vector<std::vector<int>> matrix;
};

//A candidate resonance structure, printed only under NRTDTL: everything examined, not only what survived.
struct NboNrtCandidate {
	int structure = 0;
	int idxres = 0;
	double rho_nl = 0.0;
	std::string spin;
	std::vector<std::vector<int>> topo;
};

//One step of the QP minimisation, printed only under NRTDTL.
struct NboQpIteration {
	int iteration = 0;
	int structures = 0;
	double d_w = 0.0;
	std::string kkt;
	double rho_nl = 0.0;
	int added = 0;              //structure that entered at this step, 0 if none
	std::string spin;
};

struct NboNrt {
	bool present = false;
	int structures_used = 0;
	int structures_found = 0;
	double d_w = 0.0;
	double d_0 = 0.0;
	double search_seconds = 0.0;
	double gram_seconds = 0.0;
	double minimize_seconds = 0.0;
	double other_seconds = 0.0;
	//Thresholds NBO echoes back. The delocalisation-list one (NRTE2) decides which E2 interactions enter
	//the resonance search, so a reference is only reproducible with it recorded.
	double parent_threshold_percent = 0.0;
	double deloc_threshold_kcal = 0.0;
	int max_search_cycles = 0;
	int initial_topo = 0, initial_nls = 0, initial_nbi = 0, initial_sym = 0;
	std::string symmetry;       //"Dih symmetry, 32 symmetry operator(s), ...", verbatim
	std::vector<NboResonanceWeight> weights;
	std::vector<NboBondOrder> bond_orders;
	std::vector<NboValency> valencies;
	std::vector<NboNrtCycle> cycles;
	std::vector<NboTopo> leading_topo;
	//NRTDTL only.
	std::vector<NboNrtCandidate> candidates;
	std::vector<NboQpIteration> qp_iterations;
	std::vector<std::string> arrows;          //"ARROWS generates N new structures from ..." lines
	//Other search notes (candidate budget, exhaustive count, dropped half-arrow intermediates), kept out of arrows
	std::vector<std::string> notes;
	std::vector<std::string> symmetry_forms;  //"Symmetry equivalent resonance forms" block
	std::string nrtstr_keylist;               //the $NRTSTR block NBO writes back, verbatim
};

struct NboResults {
	std::string name;
	std::string source;         //the wavefunction the archive came from
	std::string version;        //NBO banner
	std::string keywords;       //what went into the $NBO keylist
	//keywords NBO echoed back, in printed order: what it actually recognised, and all there is from output alone
	std::string keywords_reported;
	bool open_shell = false;
	double file47_seconds = 0.0;
	double nbo_seconds = 0.0;      //wall clock around the gennbo process, wrapper-measured
	//NBO's closing "completed in X CPU seconds (Y wall seconds)"; NBO 7.0.9 is serial, so CPU is one thread's time.
	double nbo_cpu_seconds = 0.0;
	double nbo_reported_wall_seconds = 0.0;
	//Native-route stages, summed over both spins; NRT reports its own phases in the nrt block.
	double nao_seconds = 0.0;
	double search_seconds = 0.0;
	double e2_seconds = 0.0;
	//E2 printing thresholds NBO echoes: below them nothing is printed, so a missing interaction looks like a small one.
	double e2_threshold_kcal = 0.0;
	double e2_intermolecular_threshold_kcal = 0.0;
	std::vector<NboAtomPopulation> npa;   //spin-summed table
	std::vector<NboNao> nao;              //spin-summed table
	//Per-spin NAO sets of an open shell; the spin-summed table hides per-spin errors that cancel.
	//Empty for a closed shell.
	std::vector<NboNao> nao_alpha;
	std::vector<NboNao> nao_beta;
	std::vector<NboOrbital> orbitals;
	std::vector<NboE2Entry> e2;
	NboNrt nrt;
};

/** Bumped when a parse or JSON change alters what a stored reference means, so a reader can refuse
	an older file. 2: NRT "RS" column read as a rank, composite alpha+beta valency table kept under its
	own label, open-shell NAO Spin column not read as Energy. */
constexpr int NBO_JSON_PARSER_VERSION = 2;

/** Parse a .nbo output file produced by NBO 6/7 into structured results. */
NboResults parse_nbo_output(const std::filesystem::path& nbo_file);

/** Write the results as JSON, the format the reference dataset is stored in. */
void write_nbo_json(const NboResults& results, const std::filesystem::path& json_file);

struct NboTolerances {
	double charge = 2.0e-3;
	double occupancy = 2.0e-3;
	double energy = 2.0e-3;         //a.u.
	double hybrid_percent = 0.5;    //percentage points of s/p/d character
	double e2_kcal = 0.1;           //kcal/mol
	double bond_order = 5.0e-3;
	double weight_percent = 0.5;    //percentage points of NRT weight
};

struct NboQuantityAgreement {
	std::string quantity;
	int compared = 0;
	int missing = 0;       //present in the reference, absent from the candidate
	int failed = 0;
	double max_deviation = 0.0;
	std::string worst;     //what the largest deviation was on
};

struct NboComparison {
	bool ok = false;
	std::vector<NboQuantityAgreement> quantities;
	std::vector<std::string> messages;
	std::string report() const;
};

/**
 * Per-quantity agreement between a reference result and a candidate one.
 * The reference drives: every reference entry must have a counterpart.
 */
NboComparison compare_nbo_results(const NboResults& reference, const NboResults& candidate,
	const NboTolerances& tolerances = NboTolerances());

struct NboRunOptions {
	std::filesystem::path wavefunction;
	std::filesystem::path work_dir;      //empty: next to the wavefunction
	std::string keywords;                //$NBO keylist additions, e.g. "NRT NRTE2=5"
	std::string executable;              //empty: auto (WSL ~/nbo7/gennbo on Windows, gennbo on PATH)
	std::filesystem::path json_out;      //empty: <stem>.nbo.json
	bool debug = false;
};

/** Write the .47, run the external NBO, parse the output, write the JSON. Returns NBO's exit status. */
int run_nbo(const NboRunOptions& options, std::ostream& log);
