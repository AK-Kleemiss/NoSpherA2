#pragma once
#include "convenience.h"
#include "wfn_class.h"
#include "cell.h"
#include "integrator.h"
#include "gaussian_atom.h"

//Interaction energies between every pair of molecules in contact in a crystal. The density of each unique molecule is
//predicted or fitted once; a symmetry image is made by moving the atoms and rotating the aux coefficients with them.
namespace crystal_energies {
	struct symop { ivec2 rot; vec trans; };
	//w neighbours of A have this energy (A,B and its inverse image), contact: two atoms closer than the cutoff
	struct pair { int A, B, op, w = 1; ivec n; double distance; bool contact = true; vec centroid_A, centroid_B, t; vec2 R; std::string symop; DensityFitting::Interaction_Energy E; };
	struct job { std::filesystem::path cif, output; double cutoff = 3.8, lattice = 0.0; pathvec structures, coefficients; };
	//M in e bohr, volume in bohr^3, Z molecules per cell, E = -2 pi |M|^2 / (3 V Z) Hartree per molecule
	struct dipole { vec M; double volume = 0.0, E = 0.0; int Z = 0; };
	//c' = D c for a density moved by the orthogonal matrix R, D in the m = -l..l order of the aux basis
	vec2 real_sh_rotation(const int l, const vec2& R);
	//Moves the atoms of aux to R x + t (bohr) and rotates the coefficients with them
	void transform(WFN& aux, vec& coef, const vec2& R, const vec& t);
	std::vector<symop> symops(const cell& c);
	std::string symop_string(const symop& op, const ivec& n);
	//Cartesian cell matrix in bohr, columns a b c
	vec2 cell_matrix(const cell& c);
	vec2 inverse3(const vec2& M);
	//Every pair of molecules with two atoms closer than cutoff or centroids closer than radius (Angstrom), each pair once
	//with A as given and B = op(B) + n
	std::vector<pair> contacts(const cell& c, const std::vector<WFN>& molecules, const double cutoff, const double radius = 0.0);
	/*unit-cell dipole from the molecular dipoles about their centroids (de Leeuw, Perram, Smith 1980): E takes a
	lattice sum over a sphere in vacuum to the tin-foil (Ewald) one, zero outside polar space groups
	*/
	dipole cell_dipole(const cell& c, const std::vector<Gaussian_Molecule>& molecules);
	job read_job(const std::filesystem::path& file);
	//the contact pairs, the cell dipole and for radius > 0 the lattice energy 1/2 sum_B E_AB over centroids within radius
	void write_table(const std::vector<pair>& pairs, const pathvec& structures, const dipole& D, std::ostream& out, const double radius = 0.0);
	int run(options& opt);
}
