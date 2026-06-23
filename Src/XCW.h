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
		: settings(loadSettings())
	{
		construct(opt_in);
	};


	// Calculates F_calc without DW factors (=1), needs rework
	//void calc_F_calc_fast();
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
	struct cryst_info {
		int nr_small;
		int nr;
		int n_params;
		double F_scale;
		double inv_scale;
		int ncen;
		int nmo;
		double GooF;
		double chi2;
		vec U_iso;
	};
	struct SCF_settings {
		double quant_diff;
		bool conv_quant_diff = false;
		double max_diis_error;
		bool conv_max_diis_error = false;
		double gradient;
		bool conv_gradient = false;
		double RMSP_diff;
		bool conv_RMSP_diff = false;
		double MaxP_diff;
		bool conv_MaxP_diff = false;
		double diis_stop_damping;
		bool apply_shift = true;
		double diis_stop_shift;
		bool apply_damping = true;
		std::string basis_set_name;
		bool grown;
		int n_params;
		int refine_against;
		occ::qm::SpinorbitalKind hf_type;
		double alpha;
		double level_shift;
		int num_xcw_steps;
		double xcw_step_size;
		int max_scf_iterations;
		int charge;
		int multiplicity;

		void clear() {
			conv_quant_diff = false;
			conv_max_diis_error = false;
			conv_gradient = false;
			conv_RMSP_diff = false;
			conv_MaxP_diff = false;
			apply_shift = true;
			apply_damping = true;
		}
		bool convergence_check() const {
			if (conv_quant_diff == true && conv_max_diis_error == true && conv_gradient == true && conv_RMSP_diff == true && conv_MaxP_diff == true) {
				return true;
			}
			else {
				return false;
			}
		}
		void update(const double& diis_error, std::ostream& file, double& alpha) {
			if (diis_error < diis_stop_damping && apply_damping == true) {
				apply_damping = false;
				print_centered_message("***Turned off damping***", 76, file);
				alpha = 0;
			}
			if (diis_error < diis_stop_shift && apply_shift == true) {
				apply_shift = false;
				print_centered_message("***Turned off level shift***", 76, file);
			}
		}
	};


	// Constructor of the XCW class
	void construct(const options& opt_in);
	// Loads the convergence settings
	SCF_settings loadSettings();


	// Mathematical helper functions
	// Helper function for coordinate transformation without building full tensors
	void get_voigt_index(const ivec& indices, int& ADP_idx);
	// Helper function for flattening the I tensor
	size_t tri_index(int mu, int nu) const noexcept;
	// Helper function for flattening the I tensor
	size_t flattened_idx(int r, int mu, int nu) const noexcept;
	// Converts the ADP matrix (just U) from cif format into reciprocal space
	void U_cif2U_star();
	// Converts all ADP tensors from reciprocal space into real space
	void U_star2U_cart();


	// Various helper functions
	// Generates a list that links the symmetry operations to symmetry-generated reflexes for given reflex r
	ivec generate_asym_lookup(const int r);
	// Sets up a molecule object from the asym_atoms
	void setup_SCF_mol(occ::core::Molecule& mol);
	// Sets up the basis set with a previously generated molecule and basis set from JKFit, where the Olex2 basis sets are now located
	void setup_basis(occ::core::Molecule& mol, std::string& basis_set_name, occ::qm::AOBasis& occ_basis_set);


	// Methods used for calculation of the I tensor
	// Combined method used to save memory, calculates both the I tensor and the correction for F_calc from anomalous dispersion
	void eval_I_anom_disp(std::vector<ao_data>& ao_data_shells, bool read);
	// Evaluates the I tensor
	// Needs: DW factors, phase factors, translational phase factors & XCW integrals
	void eval_I(std::vector<ao_data>& ao_data_shells, cvec2& DW_fact, cvec2& phase_fact, cvec2& translation_phase);
	// Evaluates Debye-Waller factors
	void eval_DW(cvec2& DW_fact);
	// Evaluates the rotational contribution to the phase factors
	void eval_phase(cvec2& phase_fact);
	// Evaluates the translational contribution to the phase factors
	void eval_translation_phase(cvec2& translation_phase);
	// Calculates the XCW_integrals
	//Needs: Primitve vectors
	cvec2 calculateXCWintegral(GridManager& grid_manager, const ao_data& mu_data, const ao_data& nu_data, const ivec& asym_atom_list, vec2& k_pt, bool& equal, vec2& mu_vals);
	// Creates primitive vectors from the basis set for calculating the XCW integrals
	void create_prims(std::vector<ao_data>& ao_data_shells, occ::qm::AOBasis& occ_basis_set);
	// Helper function that does eveything needed for XCW procedure
	occ::qm::HartreeFock setup_XCW_procedure(bool read, bool safe);


	// Methods used for calculation of F_calc
	// Calculates F_calc
	// Needs: I tensor, density matrix, anomalous dispersion correction
	void calc_F_calc(const dMatrix2& D);
	// Calculates direct corrections of the anomalous dispersion onto F_calc
	// Needs: Anomalous dispersion information from CIF
	void eval_anom_disp(cvec2& DW_fact, cvec2& phase_fact, cvec2& translation_phase);
	// Parses the anomalous dispersion information from a CIF style .txt file
	void parse_anom_atoms(std::vector<anom_atom>& anom_atoms);


	// Methods used to determine various crystallographic parameters
	// Evaluates the scaling factor for |F_calc| by least squares fitting
	void eval_scale();
	// Calculates quality criteria like GooF and chi^2
	void calc_criteria();


	// Calculates the perturbation matrix elements
	// Needs: I tensor, density matrix, F_calc, F_obs, sigma_obs, n_params & n_reflections
	void calc_perturb(occ::Mat& perturb, const occ::qm::SCF<occ::qm::HartreeFock>& scf);


	// Methods used for XCW fitting by run_XCW_fitting()
	// Executes a single SCF solver
	void do_SCF(const double& lambda, double& alpha, occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::qm::Wavefunction& last_wfn, bool& has_guess);
	// Executes a single SCF iteration
	bool SCF_iteration(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& lambda, double& alpha, double& e_diff_mem, double& quant, double& last_quant, occ::Mat& dm_last);
	// Checks convergence for SCF cylce
	bool SCF_convergence_check(const double& e_diff, const double& gradient, occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& dm_last);
	// Computes the orbital gradient for usage as a convergence criterion
	double compute_orbital_gradient(const occ::qm::SCF<occ::qm::HartreeFock>& scf);
	// Computes convergence criteria related to the density matrix (RMSP and MaxP)
	void get_density_criteria(double& RMSP_diff, double& maxP_diff, const occ::Mat& dm, const occ::Mat& dm_last);
	// Takes care of dynamic damping (very primitive)
	double dynamic_damping(const occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& current_alpha, const double& e_diff, double& e_diff_mem);
	// Applies level shift to fock matrix
	void apply_level_shift(const occ::Mat& C_old, const occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& F_diis);
	// Takes the SCF object from occ and creates the tscb file
	void create_tscb(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& lambda);
	// Builds the density matrix to use for structure factor calculations
	void build_effective_dm(const occ::qm::SCF<occ::qm::HartreeFock>& scf, dMatrix2& dm_ref, const occ::Mat& dm_old);


	// Available lists and variables that are often used
	ivec asym_atom_list;
	vec2 k_pt;
	// First F_calc (with anomalous dispersion), then anomalous correction
	cvec2 F_calc;
	cvec I;
	std::vector<asym_atom> asym_atoms;
	std::vector<scattering_data> obs;
	hkl_list hkl;
	hkl_list hkl_enlarged;
	const options* opt;
	WFN dummy_wave;
	cell unit_cell;
	std::ofstream XCW_log;
	SCF_settings settings;
	cryst_info cryst;
};