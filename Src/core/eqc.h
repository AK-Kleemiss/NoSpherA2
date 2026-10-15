#pragma once
//Rahm and Hoffmann's EQC decomposition (github.com/martinrahm/X-analysis) of A + B + ... -> C:
//E = nX + Vnn - Eee per wavefunction, Vnn over the full Z (ECP or not, as cclib), Eee the remainder.
//Mode 1 (-eqc gbw -eqc_frag): parent re-converged in occ from the gbw, fragments from its density blocks.
//Mode 2 (-eqc_wfn): terms from files; wfn, gbw and molden take E from the ORCA .out.
#include <array>
#include <string>
#include <vector>

struct options;

namespace eqc
{
	//cclib's factor, so the eV numbers compare with X-analysis digit for digit
	constexpr double hartree2eV = 27.21138505;

	//The terms of one wavefunction, in hartree
	struct terms
	{
		double E = 0.0;    //total SCF energy
		double nX = 0.0;   //sum over occupied spin orbitals of occupation * orbital energy
		double Vnn = 0.0;  //nuclear repulsion with the full atomic numbers
		int n = 0;         //electrons, ECP cores included
		double Eee() const { return Vnn - (E - nX); }
	};

	//Products minus reactants, in hartree; Q, the covalency index and dEeeE as X-analysis prints them
	struct reaction
	{
		double dE = 0.0, dnX = 0.0, dVnn = 0.0, dEee = 0.0;
		double dV() const { return dVnn - dEee; }  //Delta(Vnn - Eee)
		double Q = 0.0;          //2 Delta(nX)/Delta E - 1
		double covalency = 0.0;  //per cent; the ionicity index is 100 minus this
		double dEeeE = 0.0;      //per cent, |sum Eee / sum E| of the product minus that of the reactants
		int n = 0;               //electrons of the product, what the per-electron lines divide by
	};

	//X-analysis' covalency index for x = Delta(nX) and v = Delta(Vnn - Eee): 100|x|/(|v| + |x|)
	//when x < 0, 100 v/(v - x) when x > 0, NaN at x = 0 where its formula divides by zero
	double covalency(double x, double v);

	reaction react(const std::vector<terms> &reactants, const terms &product);

	//Sum over pairs of Z_i Z_j / r_ij, positions in bohr
	double nuclear_repulsion(const std::vector<int> &Z, const std::vector<std::array<double, 3>> &pos);

	//The last "FINAL SINGLE POINT ENERGY" of an ORCA output, NaN when there is none
	double orca_energy(const std::string &out_file);

	//Runs the analysis -eqc asked for; the tables go to std::cout and <stem>.eqc_log
	int run(const options &opt);
}
