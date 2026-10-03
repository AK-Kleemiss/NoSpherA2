#pragma once
#include "convenience.h"
#include "GridManager.h"
#include "integration_params.h"
#include "scattering_factors.h"
#include "cell.h"
#include "basis_set.h"
#include "xcw_halting.h"
#include "extinction.h"
#include <occ/qm/hf.h>
#include <occ/qm/second_order_scf.h>
#include "i_tensor_stream.h"
#include "stored_eri.h"
#include <thread>
#include <deque>

// Applies M to the U/C/D tensors of one atom in their Voigt storage: T'_{ij..} = sum_pq.. M_pi M_qj .. T_pq..
void transform_ADPs(vec2& ADPs, const vec2& M);

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
		// R1 = sum ||F_obs| - scale |F_calc|| / sum |F_obs|, over the fit set (Tonto's R1(gt))
		double R1;
		// The same numbers over every reflection, next to the fit-set ones above
		double R1_all;
		double GooF1_all;
		double GooF2_all;
		double weighted_GooF1_all;
		double weighted_GooF2_all;
		int n_fit;
		vec U_iso;

		void grow_U_iso(const std::vector<asym_atom>& asym_atoms, const ivec3& symmetry_linking_list) {
			for (int i = 0; i < symmetry_linking_list.size(); i++) {
				for (int j = 0; j < symmetry_linking_list[i].size(); j++) {
					if (symmetry_linking_list[i][j].size() != 0) {
						U_iso[j] = U_iso[i];
					}
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
		//The four damping/shift values default to off: `fast_conv` only clears the apply
		//flags, yet run_lambda copies alpha and level_shift and SCF_iteration mixes the
		//density with alpha regardless of the flag
		double diis_stop_damping = 0;
		//`slow_conv` chosen: the unperturbed first step still runs the normal schedule, see run_XCW_fitting
		bool slow_conv = false;
		//`soscf`: second-order steps as soon as the DIIS error is below trah_.start_threshold, see soscf_step
		bool soscf = false;
		//`check_hessian`: finite-difference check of the Hessian-vector product on the first
		//second-order step, reported in XCW.log
		bool check_hessian = false;
		//`i_sigma <x>`: only reflections with I/sigma(I) >= x enter chi^2 and the scale, as in
		//Tonto. Under the reader's sigma(F) = sigma(I)/2F that is F/sigma(F) >= 2x, so the
		//default 2 is SHELX's I > 2 sigma(I) and F > 4 sigma(F) at once
		double i_sigma_cutoff = 2;
		bool apply_shift = true;
		bool method_apply_shift = true;
		double diis_stop_shift = 0;
		bool apply_damping = true;
		bool method_apply_damping = true;
		std::string basis_set_name;
		//`df_basis <name>`: density fitting of the Fock build with this auxiliary basis
		std::string df_basis_name;
		//`guess_basis <name>`: the first lambda starts from a Hartree-Fock converged in this
		//(smaller) basis by OCC's own driver, its density projected into the orbital basis
		std::string guess_basis_name;
		bool grown = false;
		//`extinction <shelx|bc_gaussian|bc_lorentzian> [iso|aniso] [fixed] [start value]`
		extinction::model extinction_model = extinction::model::none;
		bool extinction_aniso = false;
		bool extinction_refine = true;
		double extinction_start = 1e-4;
		//`wavelength <lambda>`, in Angstrom; overrides _diffrn_radiation_wavelength from the CIF
		double wavelength = 0.0;
		int n_params;
		int refine_against;
		int XWR_type;
		occ::qm::SpinorbitalKind hf_type;
		double alpha = 0;
		double level_shift = 0;
		double xcw_start_value;
		int num_xcw_steps;
		double xcw_step_size;
		int max_scf_iterations;
		int charge;
		int multiplicity;
		bool read_tensor;
		bool read_first_guess;
		bool nbo_output = false;
		// Largest I tensor held resident, in MB. Above it the tensor goes to disk
		// and is read back a window of reflections at a time; 0 means no limit,
		// which is the original behaviour. Set with `i_tensor_mb <n>` in the XCW
		// settings, or `stream` for the default budget.
		size_t i_tensor_max_mb;
		// `i_float` in the settings file: hold the I tensor in single precision. The device
		// computes it in single anyway, so this stores what was computed rather than a
		// widened copy of it.
		bool i_tensor_single = false;
		bool i_tensor_double = false;
		// `I_tensor <path>` in the settings file: where the streamed tensor lives. Written
		// there, and reused from there when it is already the right size for this problem,
		// so that trying another refinement setting does not rebuild it.
		std::filesystem::path i_tensor_file_path;
		// `save <path>`: write the tensor there for a later `read <path>`, on a thread, so
		// the refinement starts at once instead of waiting for 100 GB to reach the disk.
		std::filesystem::path i_tensor_save_path;

		// Clears the convergence flags
		void clear() {
			conv_quant_diff = false;
			conv_max_diis_error = false;
			conv_gradient = false;
			conv_RMSP_diff = false;
			conv_MaxP_diff = false;
			apply_shift = method_apply_shift;
			apply_damping = method_apply_damping;
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
				alpha = 0;
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

	// Helper function for flattening the I tensor
	size_t tri_index(int mu, int nu) const noexcept;
	// Helper function for flattening the I tensor

	// Converts the ADP matrix (just U) from cif format into reciprocal space
	void U_cif2U_star();
	// Converts all ADP tensors from reciprocal space into real space
	void U_star2U_cart();
	// Rotates the ADP tensors copied onto grown atoms by their linking symmetry operation
	void rotate_grown_ADPs();

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
	occ::qm::HartreeFock setup_XCW_procedure(bool read_tensor);

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

	// Evaluates the scaling factor for |F_calc| by least squares fitting, and with it the
	// extinction coefficients when a model is being refined
	void eval_scale();
	// The closed-form weighted least-squares scale alone, with the extinction shape as it is
	void solve_scale();

	// Calculates quality criteria like GooF and chi^2. When
	// h2 weighting is set, both are computed with an additional
	// 1/|H|^2 weighting (XCW_plan.md sec. 6.2, residual self-energy
	// criterion) instead of the traditional unweighted sums.
	void calc_criteria();

	// Builds (once) the per-reflection 1/|H|^2 cache used by calc_criteria/
	// calc_perturb when h2 weighting is set. No-op otherwise.
	void ensure_inv_H2_weights();

	// Reads the wavelength (settings file, else _diffrn_radiation_wavelength in the CIF),
	// sizes the coefficient vector and builds the per-reflection extinction geometry. No-op
	// when no extinction model was asked for. Called once from construct.
	void setup_extinction(const std::filesystem::path& cif);
	// Recomputes y_r, sqrt(y_r), dI/d|Fc|^2 and d|Fc_ext|/d|Fc| from the current F_calc and
	// the current coefficients. No-op when no model is active.
	void update_extinction();
	// One Levenberg-damped Gauss-Newton step on the coefficients against the same weighted
	// residual eval_scale minimises for the scale, with the scale held at its current value.
	// False when the step was negligible or had to be rejected.
	bool refine_extinction_step();
	// "ext 0.000123" or the six tensor components, for the per-lambda log
	std::string extinction_report() const;
	// The parameters the criteria divide by: the settings file's `params` plus the extinction
	// coefficients, but only while those are actually being refined
	int n_params() const {
		return settings.n_params + static_cast<int>(settings.extinction_refine ? ext_p_.size() : 0);
	}
	// The extinction shape of reflection r, 1 where no model is active
	double ext_y(const int r) const { return ext_y_.empty() ? 1.0 : ext_y_[r]; }
	double ext_sqrt_y(const int r) const { return ext_sqrt_y_.empty() ? 1.0 : ext_sqrt_y_[r]; }
	// dI/d|Fc|^2 and d(sqrt(y)|Fc|)/d|Fc|, the chain factors the gradient needs
	double ext_g(const int r) const { return ext_g_.empty() ? 1.0 : ext_g_[r]; }
	double ext_m(const int r) const { return ext_m_.empty() ? 1.0 : ext_m_[r]; }
	// a_{r,p} of x_r = sum_p a_{r,p} P_p; 1 for the isotropic models, which have one P
	double ext_a(const int r, const size_t p) const {
		return ext_a_.empty() ? 1.0 : ext_a_[static_cast<size_t>(r) * ext_p_.size() + p];
	}
	double ext_x(const int r) const {
		double x = 0.0;
		for (size_t p = 0; p < ext_p_.size(); p++) x += ext_a(r, p) * ext_p_[p];
		return x;
	}

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
	// Sum_r Re(pre_r I_r(mu,nu)) over the fit set, the walk calc_perturb and the Hessian-vector
	// product share; out is nmo x nmo, symmetric, without any prefactor
	void contract_I(occ::Mat& out, const cvec& pre);

	// Executes a single SCF solver (for specific lambda step)
	void small_basis_guess(occ::qm::SCF<occ::qm::HartreeFock>& scf);
	bool do_SCF(const double& lambda, double& alpha, occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::qm::Wavefunction& last_wfn, bool& has_guess, bool write_result = true);

	// Executes a single SCF iteration
	bool SCF_iteration(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& lambda, double& alpha, double& e_diff_mem, double& quant, double& last_quant, occ::Mat& dm_last);

	// Checks convergence for SCF cycle
	bool SCF_convergence_check(occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& dm_last);

	// The Roothaan step and the DIIS of occ's SCF with their matrix products on MKL
	void solve_orbitals(occ::qm::SCF<occ::qm::HartreeFock>& scf, const occ::Mat& F) const;
	occ::Mat diis_update(occ::qm::SCF<occ::qm::HartreeFock>& scf);
	// CDIIS over the last diis_subspace_ Fock matrices and their commutators, see cdiis_extrapolate
	occ::Mat cdiis_extrapolate(const occ::Mat& F, const occ::Mat& E);
	std::deque<occ::Mat> diis_F_, diis_E_;
	static constexpr size_t diis_subspace_ = 8;
	occ::qm::ADIIS adiis_;
	occ::qm::EDIIS ediis_;
	// The best orbitals of the running lambda step and their E + lambda chi^2, the functional
	// the perturbed Fock matrix descends. A step climbing more than rescue_rise_ above that
	// best has lost the SCF (seen on the
	// Fe(phen)2(SCN)2 UHF singlet at lambda 0.04, right after the level shift went off):
	// rescue_scf restarts it from these orbitals with DIIS cleared and the level shift and
	// damping back on for ten times longer, up to three times per step
	occ::qm::MolecularOrbitals best_mo_;
	double best_quant_ = 0;
	int rescues_ = 0;
	static constexpr double rescue_rise_ = 1.0;
	bool rescue_scf(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double quant, double& alpha);

	// The way out of a plateau the Roothaan/DIIS map does not leave (Fe(phen)2(SCN)2 UHF
	// singlet at lambda 0.04: E and chi^2 flat for 130 iterations, orbital gradient stuck at
	// 1.5e-2, damping and rescue cycling): trust-region augmented-Hessian steps (TRAH) on the
	// occupied-virtual rotations of E + lambda chi^2, the functional the perturbed Fock matrix
	// is the gradient of. Each macro step solves the augmented-Hessian eigenproblem by Davidson
	// micro-iterations with the exact Hessian-vector product (Fock response plus the response
	// of the perturbation, scale included), the level shift set so the step fits the trust
	// radius, the orbitals moved by the Cayley transform, and the trust radius updated by
	// occ::qm::trust_radius_update from the model error the step revealed; a step that raises
	// the functional is re-solved at the smaller radius from the retained subspace, and two
	// rejections at trust_min hand the orbitals back to DIIS. It keeps the occupation, so it
	// cannot swap orbitals. Entered when the orbital gradient has not halved in
	// trah_.patience iterations, or with `soscf` once the DIIS error is below its
	// start_threshold; it stays on for the rest of the lambda step, and the radius it earned
	// carries into the next one. Micro-iterations do not count towards
	// max_iter; L-BFGS with a diagonal Hessian wandered for 100+ iterations on the same case
	// (E +-5e-6 Eh, a halved step every 2-3 iterations) where the curvature of chi^2 is stiff.
	bool soscf_ = false;
	int soscf_patience_iter_ = 0;
	double soscf_patience_grad_ = 0;
	// The patience, radius and noise knobs, and the policy that moves the radius, are
	// occ's: the same algorithm runs for plain -occ jobs out of second_order_scf.h, and
	// two copies of a convergence heuristic drift.
	occ::qm::SecondOrderSettings trah_;
	double soscf_trust_ = trah_.trust_first;
	// Consecutive rejected steps taken at the smallest radius.
	int soscf_floored_ = 0;
	std::vector<occ::Vec> trah_B_, trah_HB_;
	occ::Vec soscf_kappa_, soscf_grad_, soscf_hdiag_;
	occ::Mat soscf_C_;
	double soscf_phi_ = 0, soscf_pred_ = 0;
	bool soscf_boundary_ = false;
	int trah_micro_total_ = 0;
	void soscf_reset();
	void soscf_step(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const double phi);
	void trah_solve(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const bool extend);
	occ::Vec hessian_vector(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const occ::Vec& v);
	double rebuild_at(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const occ::Mat& C_from, const occ::Vec& kappa);
	occ::Vec gradient_at(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const occ::Mat& C_from, const occ::Vec& kappa);
	void check_hessian(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda);
	occ::Vec orbital_rotation_gradient(const occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Vec& diagonal_hessian) const;
	void rotate_orbitals(occ::qm::SCF<occ::qm::HartreeFock>& scf, const occ::Mat& C_from, const occ::Vec& kappa) const;

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
	static void flip_high_m_phases(occ::qm::Wavefunction& w);

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
	//A tensor built in single precision is held, streamed and saved in single: half the
	//memory and half the traffic of the two walks per iteration, and the values are the
	//ones the GEMM produced either way. i_float / i_double in the settings override the
	//choice the build precision makes.
	std::vector<std::complex<float>> I32;
	bool i_float_ = false;
	//A copy of the resident tensor on the device does both SCF walks there
	bool i_on_device_ = false;
	// The background writer for `save <path>`. Joined, never detached: a thread still
	// running at exit is how the GPU warm-up bug of 939268f happened, and this one holds a
	// FILE* and reads the resident tensor.
	std::thread i_writer_;
	//Incremental Fock build: the two-electron part and the density it was built from
	occ::Mat G_last_, D_last_build_;
	int last_full_build_ = 0;
	double next_full_build_error_ = 0.0;
	//Two-electron integrals over the screened-in pairs, built once per run when they fit in
	//memory: the Fock build then contracts them instead of recomputing every quartet per iteration
	stored_eri eri_;
	bool eri_on_device_ = false;
	occ::Mat eri_fock(const occ::qm::MolecularOrbitals& mo, bool screen) const;
	std::string i_writer_error_;
	void start_i_save();
	void finish_i_save();
	size_t i_budget(const char*& source, bool& automatic) const;
	static constexpr const char* i_tensor_default = "I_tensor_stream.bin";
	i_tensor_file i_file_;
	bool i_streamed_ = false;
	int i_window_ = 0;
	// The packed (mu, nu) run of reflection r, from the loaded window or from the resident tensor
	const cdouble* i_block(const int r) const
	{
		return i_streamed_ ? i_file_.block(r) : I.data() + static_cast<size_t>(r) * i_compact_;
	}
	const std::complex<float>* i_block32(const int r) const
	{
		return i_streamed_ ? i_file_.block32(r) : I32.data() + static_cast<size_t>(r) * i_compact_;
	}
	//Only the (mu, nu) pairs the overlap screening kept are stored, i_compact_ of the
	//nmo (nmo + 1) / 2, in the order these two lists give
	ivec i_pair_mu_, i_pair_nu_;
	size_t i_compact_ = 0;
	std::vector<asym_atom> asym_atoms;
	std::vector<scattering_data> obs;
	hkl_list hkl;
	hkl_list hkl_enlarged;
	ivec3 original_rotations;
	// Symmetry operations the structure factors are summed over: all of them, or one per coset
	// of the subgroup a grown cluster is closed under (cell::grown_subgroup)
	ivec sym_ops_;
	GridManager tsc_grids;
	// Ordered snapshot of `hkl` (see ensure_hkl_ordered), i.e. hkl_ordered_[r]
	// is the Miller index of reflection r as used for obs[r]/F_calc[0][r].
	std::vector<i3> hkl_ordered_;
	// 1/|H_r|^2 per reflection, see ensure_inv_H2_weights.
	vec inv_H2_;
	// The refined extinction coefficient, or the six Voigt components X11 X22 X33 X12 X13 X23
	// of the anisotropic tensor. Empty when no model is active, which is what every extinction
	// branch tests on.
	vec ext_p_;
	// Per-reflection geometry, built once: 0.001 lambda^3/sin(2 theta), cos(2 theta), and (for
	// the anisotropic models only) the nr_small x 6 coefficients a_{r,p}
	vec ext_c_, ext_cos2t_, ext_a_;
	// Per-iteration shape, see update_extinction. ext_dyc_ is dy/dt * c_r * |Fc_r|^2, so that
	// dy/dP_p = ext_dyc_[r] * a_{r,p}.
	vec ext_y_, ext_sqrt_y_, ext_g_, ext_m_, ext_dyc_;
	// Reflection r is in the fit set (I/sigma(I) >= i_sigma_cutoff), see construct
	bvec fit_mask_;
	// The criterion the SCF descends (XWR_type x refine_against), over the fit set or over all
	double criterion(bool all) const;
	// Per-lambda Gaussian halting diagnostics, see evaluate_gaussian_halting.
	std::vector<GaussianHaltEntry> gaussian_halt_history_;
	const options* opt;
	WFN dummy_wave;
	cell unit_cell;
	std::ofstream XCW_log;
	SCF_settings settings;
	cryst_info cryst;
};
