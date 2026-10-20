#pragma once
#include <occ/io/occ_input.h>
#include <occ/qm/wavefunction.h>
#include <iosfwd>
#include <string>

//The unrelaxed CCSD (or CCSD(T)) one-particle density of a closed-shell SCF reference from the Lambda
//equations (Stanton & Bartlett, JCP 98, 7029 (1993)), returned as natural orbitals
struct CcsdDensityOptions
{
	int n_frozen = -1;              //-1: occ's chemical core, 0: all electrons correlated
	std::string backend = "exact";  //exact | df
	std::string aux_basis;          //df and the local MP2 screening; empty resolves occ's correlation fitting basis
	double max_memory_gb = 1.0;
	double tol = 1e-5;              //max |D (t_new - t_old)| of the CCSD and Lambda amplitudes (ORCA's STol, TightSCF)
	int max_cycle = 200;
	bool triples = false;           //CCSD(T): the (T) terms in the Lambda equations and in the density
	bool local = false;             //cluster-in-molecule: one cluster per Pipek-Mezey LMO, densities assembled
	double cluster_radius = 8.0;    //Angstrom between LMO centres for a cluster candidate
	double domain_threshold = 0.01; //Mulliken population that puts an atom's PAOs into the candidates' domain
	double pair_threshold = 1e-5;   //|semicanonical MP2 pair energy| (Eh) with the central LMO to join its cluster
	double lno_threshold = 1e-6;    //MP2 virtual natural-orbital occupation kept in the cluster (local natural orbitals)
};

struct CcsdDensityResult
{
	occ::qm::Wavefunction natural;  //natural orbitals in mo.C, occupation/2 in mo.occupation, alpha density in mo.D
	occ::Vec occupations;           //natural occupation numbers, descending; they sum to the electron count
	double e_corr = 0.0;
	double e_triples = 0.0;         //(T) energy, from the same arrays that enter the density
	int n_frozen = 0;
	int lambda_iterations = 0;
	double t1_norm = 0.0, t2_norm = 0.0, l1_norm = 0.0, l2_norm = 0.0;
};

CcsdDensityResult ccsd_natural_orbitals(const occ::qm::Wavefunction &scf, const CcsdDensityOptions &opts, std::ostream &log);

//the same density from one local cluster per Pipek-Mezey LMO (cluster-in-molecule, Li & Piecuch), each solved
//by the code above and summed in the AO basis; clusters spanning the molecule give ccsd_natural_orbitals back
CcsdDensityResult local_ccsd_natural_orbitals(const occ::qm::Wavefunction &scf, const CcsdDensityOptions &opts, std::ostream &log);

//the relaxed MP2 density (MP2 Lagrangian with the Z-vector orbital response, PySCF grad/mp2.py) as natural orbitals;
//e_corr is E(MP2 corr), lambda_iterations the Z-vector iterations
CcsdDensityResult mp2_natural_orbitals(const occ::qm::Wavefunction &scf, const CcsdDensityOptions &opts, std::ostream &log);

//-occ with method = "mp2", "ccsd" or "ccsd(t)": occ's RHF, then the density; ccsd_* keys of the toml's [scf] table
CcsdDensityResult ccsd_density_from_occ_input(const std::string &toml_path, occ::io::OccInput config, std::ostream &log);

//every -occ toml goes through here: MP2, CCSD and CCSD(T) give their natural orbitals, anything else is occ's SCF
occ::qm::Wavefunction wavefunction_from_occ_input(const std::string &toml_path, const occ::io::OccInput &config, std::ostream &log);
