#pragma once

#include <filesystem>
#include <string>
#include <vector>

/*
 * Driving the external NBO 7 program and reading its output back into structures.
 *
 * This is the reference side of the in-house NBO work: NoSpherA2 writes a FILE47,
 * hands it to the licensed gennbo, and parses everything back. The structures below
 * define what an in-house implementation has to reproduce, and compare_nbo_results()
 * is the check it gets measured with.
 */

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
	double energy = 0.0;    //"Spin" column in the spin-summed block of an open-shell run
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
	int structure = 0;
	double weight_percent = 0.0;
	std::string spin;
};

struct NboBondOrder {
	int atom1 = 0;
	int atom2 = 0;
	double total = 0.0;
	double covalent = 0.0;
	double ionic = 0.0;
	std::string spin;           //"", "alpha", "beta", "composite"
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
	//Thresholds NBO echoes back for this run. The delocalisation-list one is what NRTE2 sets and
	//it decides which E2 interactions enter the resonance search, so a reference is only
	//reproducible with it recorded.
	double parent_threshold_percent = 0.0;
	double deloc_threshold_kcal = 0.0;
	std::vector<NboResonanceWeight> weights;
	std::vector<NboBondOrder> bond_orders;
};

struct NboResults {
	std::string name;
	std::string source;         //the wavefunction the archive came from
	std::string version;        //NBO banner
	std::string keywords;       //what went into the $NBO keylist
	bool open_shell = false;
	double file47_seconds = 0.0;
	double nbo_seconds = 0.0;
	//The two E2 printing thresholds NBO echoes; below them the table is simply not printed, so a
	//comparison that does not know them cannot tell a missing interaction from a small one.
	double e2_threshold_kcal = 0.0;
	double e2_intermolecular_threshold_kcal = 0.0;
	std::vector<NboAtomPopulation> npa;   //spin-summed table
	std::vector<NboNao> nao;              //spin-summed table
	std::vector<NboOrbital> orbitals;
	std::vector<NboE2Entry> e2;
	NboNrt nrt;
};

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

/** Write the .47, run the external NBO, parse the output, write the JSON. Returns 0 on success. */
int run_nbo(const NboRunOptions& options, std::ostream& log);
