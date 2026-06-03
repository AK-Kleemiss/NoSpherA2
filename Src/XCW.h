#pragma once
#include "convenience.h"
#include "GridManager.h"
#include "integration_params.h"
#include "scattering_factors.h"
#include "cell.h"
#include "basis_set.h"
#include <occ/qm/hf.h>

class XCW {
public:
	// Constructor
	XCW(const options& opt_in)
	{
		construct(opt_in);
	};


	// Calculates F_calc without DW factors (=1)
	void calc_F_calc_fast();
	// Does the XCW fitting routine
	void run_XCW_fitting();

private:

	// Structures
	struct ao_data {
		std::vector<primitive> prims;
		d3 pos;
		int m;
	};
	struct anom_atom {
		std::string identifier;
		cdouble dispersion;
	};


	// Constructor of the XCW class
	void construct(const options& opt_in);


	// Methods for evaluating essential crystallographic measures
	// Evaluates Debye-Waller factors (currently working up to U_anis = 2)
	void eval_DW();
	// Evaluates the rotational contribution to the phase factors
	void eval_phase();
	// Evaluates the translational contribution to the phase factors
	void eval_translation_phase();
	// Parses the anomalous dispersion information from a CIF style .txt file
	void parse_anom_atoms(std::vector<anom_atom>& anom_atoms);
	// Calculates direct corrections of the anomalous dispersion onto F_calc
	void eval_anom_disp();
	// Evaluates the scaling factor for |F_calc| by least squares fitting
	void eval_scale();
	// Calculates the chi^2 quality criterion
	double calc_chi2();


	// Methods for evaluating elements used for building the XCW perturbation potential
	// Creates primitive vectors from the basis set for calculating the XCW integrals
	void create_prims(std::vector<ao_data>& ao_data_shells);
	// Calculates the Hirshfeld weighted overlap integrals with phase factor
	cvec2 calculateXCWintegral(GridManager& grid_manager, const ao_data& mu_data, const ao_data& nu_data, const ivec& asym_atom_list, vec2& k_pt, bool& equal, vec2& mu_vals);
	// Helper function for flattening the I tensor
	size_t tri_index(int mu, int nu) const noexcept;
	// Helper function for flattening the I tensor
	size_t flattened_idx(int r, int mu, int nu) const noexcept;
	// Evaluates the I tensor (using DW factors & phase factors)
	void eval_I(std::vector<ao_data>& ao_data_shells);
	// Calculates F_calc using the I tensor and adds anomalous dispersion corrections
	void calc_F_calc(const dMatrixRef2& D);
	// Calculates the perturbation matrix elements
	void calc_perturb(occ::Mat& perturb);


	// Methods for performing the SCF procedure
	// Sets up a molecule object from the asym_atoms
	void setup_SCF_mol(occ::core::Molecule& mol);
	// Sets up the basis set with a previously generated molecule and basis set from JKFit, where the Olex2 basis sets are now located
	void setup_basis(occ::core::Molecule& mol, std::string& basis_set_name);
	// Executes a single SCF solver
	void do_SCF(const double& lambda, double& alpha, occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::qm::Wavefunction& last_wfn);
	// Executes a single SCF iteration
	bool SCF_iteration(occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& dm_old, occ::Mat& dm, const int& iter, double& ehf_last, const double& lambda, double& ehf, double& e_diff, double& alpha, const double& next_alpha, bool alpha_set);
	// Checks convergence for SCF cylce
	bool SCF_convergence_check(double& e_diff, occ::qm::SCF<occ::qm::HartreeFock>& scf);
	// Takes the SCF object from occ and creates the tscb file
	void create_tscb(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& lambda);


	// Methods for computing temporary variables, matrices and useful lists
	// Sets the rotational contribution to the phase factors
	void set_phase(cvec2& phase_in) { phase_fact = phase_in; }
	// Sets the Debye-Waller factors
	void set_DW(vec2& DW_in) { DW_fact = DW_in; }
	// Sets the translational contribution to the phase factors
	void set_translation_phase(cvec2& phase_in) { translation_phase = phase_in; }
	// Converts a matrix from fractional to reciprocal coordinates
	void U_frac2U_rec();
	// Converts a matrix from reciprocal to cartesian coordinates
	void U_rec2U_cart();
	// Generates a list that links the symmetry operations to symmetry-generated reflexes for given reflex r
	ivec generate_asym_lookup(const int r);
	
	
	// Available lists and variables that are often used
	int nmo;
	int ncen;
	int nr;
	int nr_small;
	int n_params;
	double scale;
	ivec asym_atom_list;
	vec U_iso;
	svec labels;
	bvec needs_grid;
	vec2 k_pt;
	vec2 DW_fact;
	cvec anom_corr;
	cvec F_calc;
	cvec2 phase_fact;
	cvec2 translation_phase;
	cvec I;
	std::vector<asym_atom> asym_atoms;
	std::vector<scattering_data> obs;
	hkl_list hkl;
	hkl_list hkl_enlarged;
	const options* opt;
	WFN wave;
	Int_Params params;
	cell unit_cell;
	occ::core::Molecule mol;
	occ::qm::AOBasis occ_basis_set;
	std::ofstream XCW_log;
};