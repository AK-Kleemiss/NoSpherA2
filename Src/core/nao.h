#pragma once
#include "convenience.h"
#include "atoms.h"

class WFN;

//NAO/NPA after Reed, Weinstock and Weinhold, J. Chem. Phys. 83 (1985) 735: symmetry-averaged (atom, l)
//block diagonalisation gives pre-NAOs, split into the minimal basis (NMB, core + valence) and Rydberg set (NRB);
//occupancy-weighted symmetric orthogonalisation (OWSO) of the NMB, Schmidt projection of the NRB out of it,
//OWSO of the NRB, then re-diagonalisation per (atom, l). Occupancies sum to Tr(PS) by construction.

enum class NAOClass { Core, Valence, Rydberg };

//One AO of the internal spherical basis, in the order Int_Params/libcint produce (atom-major,
//shells sorted by l inside an atom, 2l+1 components in libcint's m = -l..+l order).
struct NAOBasisFunction {
	int atom = 0;
	int l = 0;
	int shell = 0;  //index of the shell inside its (atom, l) group, in basis order
	int m = 0;      //0 .. 2l
};

struct NAO {
	int atom = 0;
	int l = 0;
	int m = 0;      //0 .. 2l, libcint component order
	int shell = 0;  //index inside (atom, l), sorted by decreasing occupancy
	int n = 0;      //principal quantum number label, n = l + 1 + shell
	NAOClass type = NAOClass::Rydberg;
	double occupation = 0.0;
};

struct NAOAtom {
	int index = 0;
	int Z = 0;
	std::string label;
	double Z_eff = 0.0;  //Z minus the core electrons an ECP replaces
	double core = 0.0;
	double valence = 0.0;
	double rydberg = 0.0;
	double population = 0.0;
	double charge = 0.0;
};

struct NAOResult {
	//nao x nao, column i holds the AO coefficients of NAO i; C^T S C = 1
	dMatrix2 C;
	std::vector<NAO> orbitals;
	std::vector<NAOAtom> atoms;
	double population = 0.0;
	double core = 0.0;
	double valence = 0.0;
	double rydberg = 0.0;
};

struct NPAResult {
	//The total and spin populations share one NAO basis for an unrestricted case.
	NAOResult total;
	bool spin_resolved = false;
	NAOResult alpha;
	NAOResult beta;
	//per atom alpha - beta population, empty unless spin_resolved
	vec spin_population;
};

//AO map of the internal spherical basis, from Int_Params so it matches the overlap and the stored density matrix
std::vector<NAOBasisFunction> spherical_ao_map(const WFN &wavy);

//in spherical_ao_map order
dMatrix2 ao_overlap(const WFN &wavy);

//Shells per l = 0..3 in the natural minimal basis of element Z, and how many of them are core
void natural_minimal_shells(const int Z, int (&n_shell)[4], int (&n_core)[4]);

//P and S both in spherical_ao_map order
NAOResult build_naos(const dMatrix2 &P, const dMatrix2 &S, const std::vector<NAOBasisFunction> &ao,
					 const std::vector<atom> &atoms, const ivec &ecp_electrons);

//spin-resolved when the reader kept a beta density
NPAResult natural_population_analysis(const WFN &wavy);

void print_npa(const NPAResult &result, std::ostream &out);
