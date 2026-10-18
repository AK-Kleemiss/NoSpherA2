#pragma once
#include "convenience.h"

class WFN;

//Intrinsic atomic orbitals and intrinsic bond orbitals after Knizia, J. Chem. Theory Comput. 9 (2013) 4834.
//The IAOs span the occupied space exactly and are built against the MINAO reference basis; the IBOs are the
//occupied orbitals rotated (2x2 Jacobi sweeps) to maximise sum_k sum_A (q_A^k)^4, q_A^k the IAO population of
//orbital k on atom A. Closed-shell only.

enum class IBOKind { Core, LonePair, Bond, Delocalised };

struct IBOResult {
	ivec mos;           //the occupied MOs that were localised, in WFN MO order
	dMatrix2 U;         //IBO k = sum_i U(i, k) MO mos[i]
	dMatrix2 iao_pop;   //(atom, IBO), each column sums to 1
	dMatrix2 mulliken;  //(atom, IBO), Mulliken population of the IBO as ORCA prints it
	vec charge;         //IAO partial charges
	std::vector<IBOKind> kind;
	ivec2 centres;      //per IBO the atoms by decreasing IAO population (top two, one for one-centre IBOs)
	vec energy;         //diagonal Fock element, sum_i U(i, k)^2 e_i
	double functional = 0.0, functional_start = 0.0;
	int sweeps = 0;
};

IBOResult intrinsic_bond_orbitals(const WFN &wavy);
void print_ibo(const IBOResult &r, const WFN &wavy, std::ostream &file);
//"6,7" (1-based IBO numbers as printed), "all", "C2:C3" or "2:3" (every bond-like IBO between two atoms), comma-separated
ivec ibo_selection(const IBOResult &r, const std::string &spec);
//Copy of wavy whose occupied MO slots r.mos[k] hold IBO k as primitive coefficients, for Calc_MO
WFN ibo_wfn(const WFN &wavy, const IBOResult &r);
