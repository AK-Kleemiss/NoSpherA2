#pragma once
#include <occ/io/occ_input.h>
#include <occ/qm/wavefunction.h>
#include <iosfwd>
#include <string>

//The unrelaxed CCSD one-particle density of a closed-shell SCF reference: occ's CCSD for the T amplitudes,
//the Lambda equations for the left-hand state, and the density they define (Stanton & Bartlett, JCP 98, 7029
//(1993)), returned as natural orbitals so every grid, tsc and property path reads it like any wavefunction.
struct CcsdDensityOptions
{
	int n_frozen = -1;              //-1: occ's chemical core, 0: all electrons correlated
	std::string backend = "exact";  //exact | df | thc, as occ's CCSD
	std::string aux_basis;          //df/thc; empty resolves occ's correlation fitting basis
	double max_memory_gb = 1.0;
	double tol = 1e-10;             //CCSD energy change and Lambda amplitude change
	int max_cycle = 200;
};

struct CcsdDensityResult
{
	occ::qm::Wavefunction natural;  //natural orbitals in mo.C, occupation/2 in mo.occupation, alpha density in mo.D
	occ::Vec occupations;           //natural occupation numbers, descending; they sum to the electron count
	double e_corr = 0.0;
	int n_frozen = 0;
	int lambda_iterations = 0;
	double t1_norm = 0.0, t2_norm = 0.0, l1_norm = 0.0, l2_norm = 0.0;
};

CcsdDensityResult ccsd_natural_orbitals(const occ::qm::Wavefunction &scf, const CcsdDensityOptions &opts, std::ostream &log);

//-occ with method = "ccsd": occ's restricted HF, then ccsd_natural_orbitals. The CCSD settings come from the toml's
//[scf] table (ccsd_frozen_core, ccsd_backend, ri-basis, ccsd_max_memory_gb), which occ's reader does not take.
//CCSD(T) is refused: the triples correction has no density here, and occ alone would hand back the HF one.
CcsdDensityResult ccsd_density_from_occ_input(const std::string &toml_path, occ::io::OccInput config, std::ostream &log);

//every -occ toml goes through here: CCSD gives its natural orbitals, CCSD(T) is refused, anything else is occ's SCF
occ::qm::Wavefunction wavefunction_from_occ_input(const std::string &toml_path, const occ::io::OccInput &config, std::ostream &log);
