#pragma once
#include "convenience.h"
#include "GridManager.h"
#include "integration_params.h"
#include "scattering_factors.h"
#include "cell.h"
#include "basis_set.h"
#include "xcw_halting.h"
#include <occ/qm/hf.h>
#include "i_tensor_stream.h"

class XCW {
public:

	// Constructor
	XCW(const options& opt_in)
		: settings(loadSettings(opt_in.xcw_settings_path))
	{
		construct(opt_in);
	};


	// Calculates F_calc without DW factors (=1), needs rework
	//void calc_F_calc_fast();
	
	// Does the XCW fitting routine
	void run_XCW_fitting();

private:


	// STRUCTURES

	// Data for contracted basis function
	struct ao_data {
		std::vector<primitive> prims;
		d3 pos;
		int m;
	};

	// Data for anomalous dispersion correction
	struct anom_atom {
		std::string identifier;
		cdouble dispersion;
	};

	// Miscellaneous crystallographic data 
	struct cryst_info {
		int nr_small;
		int nr;
		int n_params;
		double F_scale;
		double inv_scale;
		int ncen;
		int nmo;
		double GooF1;
		double GooF2;
		double weighted_GooF1;
		double weighted_GooF2;
		vec U_iso;

		void grow_U_iso(const std::vector<asym_atom>& asym_atoms, const ivec3& symmetry_linking_list) {
			for (int i = 0; i < asym_atoms.size(); ++i) {
				if (asym_atoms[i].grown) {
					continue;
				}
				for (int j = 0; j < symmetry_linking_list[i].size(); j++) {
					U_iso[symmetry_linking_list[i][j][0]] = U_iso[i];
				}
			}
		};
	};

	/* Data concerning convergence criteria, molecule information, basis set, etc.
	   Also keeps track of the current state of convergence */
	struct SCF_settings {
		double quant_diff;
		double current_quant_diff;
		bool conv_quant_diff = false;
		double max_diis_error;
		double current_max_diis_error;
		bool conv_max_diis_error = false;
		double gradient;
		double current_gradient;
		bool conv_gradient = false;
		double RMSP_diff;
		double current_RMSP_diff;
		bool conv_RMSP_diff = false;
		double MaxP_diff;
		double current_MaxP_diff;
		bool conv_MaxP_diff = false;
		double diis_stop_damping;
		bool method_apply_shift = true;
		bool apply_shift = true;
		double diis_stop_shift;
		bool method_apply_damping = true;
		bool apply_damping = true;
		std::string basis_set_name;
		bool grown;
		int n_params;
		int refine_against;
		int XWR_type;
		occ::qm::SpinorbitalKind hf_type;
		double alpha;
		double level_shift;
		double xcw_start_value;
		int num_xcw_steps;
		double xcw_step_size;
		int max_scf_iterations;
		int charge;
		int multiplicity;
		bool safe_tensor;
		bool read_tensor;
		bool read_first_guess;
		// Largest I tensor held resident, in MB. Above it the tensor goes to disk
		// and is read back a window of reflections at a time; 0 means no limit,
		// which is the original behaviour. Set with `i_tensor_mb <n>` in the XCW
		// settings, or `stream` for the default budget.
		size_t i_tensor_max_mb;

		// Clears the convergence flags
		void clear() {
			conv_quant_diff = false;
			conv_max_diis_error = false;
			conv_gradient = false;
			conv_RMSP_diff = false;
			conv_MaxP_diff = false;
			apply_shift = true;
			apply_damping = true;
		}

		// Performs the convergence check and sets flags accordingly
		bool convergence_check() const {
			if (conv_quant_diff == true && conv_max_diis_error == true && conv_gradient == true && conv_RMSP_diff == true && conv_MaxP_diff == true) {
				return true;
			}
			else {
				return false;
			}
		}

		// Updates the SCF routine in regards to damping and level shift
		void update(std::ostream& file, double& alpha) {
			if (current_max_diis_error < diis_stop_damping && apply_damping == true) {
				apply_damping = false;
				print_centered_message("***Turned off damping***", 84, file);
			}
			if (current_max_diis_error < diis_stop_shift && apply_shift == true) {
				apply_shift = false;
				print_centered_message("***Turned off level shift***", 84, file);
			}
		}

	};


	// FUNCTIONS

	// Constructor of the XCW class
	void construct(const options& opt_in);

	// Loads the convergence settings
	SCF_settings loadSettings(const std::filesystem::path& settings_path);

	// Helper function that transforms indices into voigt notation
	void get_voigt_index(const ivec& indices, int& ADP_idx);

	// Helper function for flattening the I tensor
	size_t tri_index(int mu, int nu) const noexcept;
	// Helper function for flattening the I tensor
	size_t flattened_idx(int r, int mu, int nu) const noexcept;

	// Converts the ADP matrix (just U) from cif format into reciprocal space
	void U_cif2U_star();
	// Converts all ADP tensors from reciprocal space into real space
	void U_star2U_cart();

	// Generates a list that links the symmetry operations to symmetry-generated reflexes for given reflex r
	ivec generate_asym_lookup(const int r);

	// Sets up a molecule object from the asym_atoms
	void setup_SCF_mol(occ::core::Molecule& mol);

	// Sets up the basis set with a previously generated molecule and basis set from JKFit, where the Olex2 basis sets are now located
	void setup_basis(occ::core::Molecule& mol, std::string& basis_set_name, occ::qm::AOBasis& occ_basis_set);

	// Combined method used to save memory, calculates both the I tensor and the correction for F_calc from anomalous dispersion
	void eval_I_anom_disp(std::vector<ao_data>& ao_data_shells, bool read);

	// Evaluates the I tensor
	void eval_I(std::vector<ao_data>& ao_data_shells, cvec2& DW_fact, cvec2& phase_fact, cvec2& translation_phase, double& time_taken, long long& screen_counter, long long& skipped_grids);

	// Evaluates Debye-Waller factors
	void eval_DW(cvec2& DW_fact);

	// Evaluates the rotational contribution to the phase factors
	void eval_phase(cvec2& phase_fact);

	// Evaluates the translational contribution to the phase factors
	void eval_translation_phase(cvec2& translation_phase);


	// Creates primitive vectors from the basis set for calculating the XCW integrals
	void create_prims(std::vector<ao_data>& ao_data_shells, occ::qm::AOBasis& occ_basis_set);

	// Combined function that sets up the XCW procedure, evaluates I tensor (or loads it from file), sets up the Hartree-Fock object and evaluates anomalous dispersion correction
	occ::qm::HartreeFock setup_XCW_procedure(bool read_tensor, bool save_tensor);

	// I tensor storage: held resident, or written to disk and read back a window
	// of reflections at a time. See decide_i_storage.
	void decide_i_storage();
	void open_i_stream_for_reading();
	std::filesystem::path i_tensor_path() const;

	// Calculates F_calc
	void calc_F_calc(const dMatrix2& D);

	// Calculates direct corrections of the anomalous dispersion onto F_calc
	void eval_anom_disp(cvec2& DW_fact, cvec2& phase_fact, cvec2& translation_phase);

	// Parses the anomalous dispersion information from a CIF style .txt file
	void parse_anom_atoms(std::vector<anom_atom>& anom_atoms);

	// Evaluates the scaling factor for |F_calc| by least squares fitting
	void eval_scale();

	// Calculates quality criteria like GooF and chi^2. When
	// h2 weighting is set, both are computed with an additional
	// 1/|H|^2 weighting (XCW_plan.md sec. 6.2, residual self-energy
	// criterion) instead of the traditional unweighted sums.
	void calc_criteria();

	// Builds (once) the per-reflection 1/|H|^2 cache used by calc_criteria/
	// calc_perturb when h2 weighting is set. No-op otherwise.
	void ensure_inv_H2_weights();

	// Distributional (Gaussian) halting criterion (see xcw_halting.h and
	// tests/P1_test/XCW_plan.md). Computes standardized residuals z_h from
	// the current F_calc/obs/F_scale, evaluates the Anderson-Darling
	// statistic and supporting diagnostics, logs them, and stores the
	// result for the final lambda* recommendation. Only called when
	// opt->xcw_gaussian_halt is set.
	void evaluate_gaussian_halting(const double lambda);

	// Prints the full per-lambda Gaussian-halting table (XCW_log only), then
	// calls report_halting_progress_estimate(true). Called once at the end
	// of run_XCW_fitting().
	void report_gaussian_halting_summary();

	// Prints the recommended lambda* = argmin A^2 so far (subject to the
	// binned-trend test), a scan-boundary warning if that argmin sits at
	// the last evaluated lambda, and -- fitting the A^2(lambda) trend so
	// far with a small family of polynomial models and picking the best by
	// AIC (see xcw_halting.h) -- an extrapolated estimate of where the
	// true minimum likely lies, with the fit's residual/quality diagnostics
	// for every candidate model tried. Called periodically during the scan
	// (is_final=false, every 5 lambda steps) and once more at the end
	// (is_final=true, from report_gaussian_halting_summary). Uses whatever
	// is in gaussian_halt_history_ at call time, so periodic calls are
	// naturally based on partial data.
	void report_halting_progress_estimate(bool is_final);

	// Builds the (once-cached) ordered list of Miller indices matching the
	// index r used for obs[r]/F_calc[0][r] (see generate_asym_lookup),
	// needed to look up per-reflection resolution for the binned trend
	// test.
	void ensure_hkl_ordered();

	// Calculates the perturbation matrix elements
	void calc_perturb(occ::Mat& perturb, const occ::qm::SCF<occ::qm::HartreeFock>& scf);

	// Executes a single SCF solver (for specific lambda step)
	void do_SCF(const double& lambda, double& alpha, occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::qm::Wavefunction& last_wfn, bool& has_guess);

	// Executes a single SCF iteration
	bool SCF_iteration(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& lambda, double& alpha, double& e_diff_mem, double& quant, double& last_quant, occ::Mat& dm_last);

	// Checks convergence for SCF cycle
	bool SCF_convergence_check(occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& dm_last);

	// Computes the orbital gradient for usage as a convergence criterion
	double compute_orbital_gradient(const occ::qm::SCF<occ::qm::HartreeFock>& scf);

	// Computes convergence criteria related to the density matrix (RMSP and MaxP)
	void get_density_criteria(double& RMSP_diff, double& maxP_diff, const occ::Mat& dm, const occ::Mat& dm_last);

	// Takes care of dynamic damping
	double dynamic_damping(const occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& current_alpha, const double& e_diff, double& e_diff_mem);

	// Applies level shift to fock matrix
	void apply_level_shift(const occ::Mat& C_old, const occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& F_diis);

	// Takes the SCF object from occ and creates the tscb file
	void create_tscb(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& lambda);

	// Builds the density matrix to use for structure factor calculations
	void build_effective_dm(const occ::qm::SCF<occ::qm::HartreeFock>& scf, dMatrix2& dm_ref, const occ::Mat& dm_old);


	// OBJECTS
	ivec asym_atom_list;
	vec2 k_pt;
	// First F_calc (with anomalous dispersion), then anomalous correction
	cvec2 F_calc;
	// Held resident only while it fits settings.i_tensor_max_mb; otherwise empty
	// and i_file_ carries the tensor. Read through i_block(r) either way.
	cvec I;
	i_tensor_file i_file_;
	bool i_streamed_ = false;
	int i_window_ = 0;
	// The packed (mu, nu) run of reflection r, from the loaded window or from the resident tensor
	const cdouble* i_block(const int r) const
	{
		return i_streamed_ ? i_file_.block(r)
			: I.data() + static_cast<size_t>(r) * (static_cast<size_t>(cryst.nmo) * (cryst.nmo + 1) / 2);
	}
	std::vector<asym_atom> asym_atoms;
	std::vector<scattering_data> obs;
	hkl_list hkl;
	hkl_list hkl_enlarged;
	// Ordered snapshot of `hkl` (see ensure_hkl_ordered), i.e. hkl_ordered_[r]
	// is the Miller index of reflection r as used for obs[r]/F_calc[0][r].
	std::vector<i3> hkl_ordered_;
	// 1/|H_r|^2 per reflection, see ensure_inv_H2_weights.
	vec inv_H2_;
	// Per-lambda Gaussian halting diagnostics, see evaluate_gaussian_halting.
	std::vector<GaussianHaltEntry> gaussian_halt_history_;
	const options* opt;
	WFN dummy_wave;
	cell unit_cell;
	std::ofstream XCW_log;
	SCF_settings settings;
	cryst_info cryst;
};