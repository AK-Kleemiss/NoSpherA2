#pragma once
#include "convenience.h"
#include "cell.h"
#include "SCF_log_writer.h"
#include "stored_eri.h"
#include <occ/qm/hf.h>
#include <occ/qm/second_order_scf.h>
#include <deque>

// A Hartree-Fock SCF on OCC's integrals with DIIS, damping, level shift, a rescue for steps that run
// away and second-order (TRAH) steps for what DIIS cannot finish. On its own it solves the plain
// system; a child adds a model term to the energy through the hooks below (XCW_solver: lambda chi^2),
// and the rescue, TRAH and the convergence check then work on E plus that term.
class SCF_solver {

public:

	SCF_solver() = default;
	SCF_solver(const occ::core::Molecule& mol, const occ::qm::AOBasis& basis, const options::SCF_settings& settings, SCF_log_writer& writer);
	virtual ~SCF_solver() = default;
	SCF_solver(SCF_solver&&) = default;
	SCF_solver& operator=(SCF_solver&&) = default;

	// Creates an OCC molecule from an asym_atom vector with charge and multiplicity
	static occ::core::Molecule setup_SCF_mol(const std::vector<asym_atom>& atoms, const int charge, const int multiplicity);

	// The basis set from the library (JKFit, where the Olex2 basis sets are now located) on the atoms of mol
	static occ::qm::AOBasis setup_basis(const occ::core::Molecule& mol, const std::string& basis_set_name);

	// Builds the stored two-electron integrals when they fit in avail less a reserve (0: unknown), and puts them on the device if asked to
	void setup_eri(const occ::qm::HartreeFock& hf, const size_t avail);
	// Gives the device back
	void release_device();
	// A density-fitted Fock build for the TRAH Hessian products alone
	void set_hessian_df(occ::qm::HartreeFock&& hf) { hessian_df_ = std::make_unique<occ::qm::HartreeFock>(std::move(hf)); }
	// The two-electron integrals are held: setup_eri found them room
	bool has_eri() const { return static_cast<bool>(eri_); }

	// One complete SCF, from the orbitals in wfn when use_guess, else from a fresh guess; the result
	// is left in wfn. False when it did not converge
	bool run(occ::qm::HartreeFock& hf, occ::qm::Wavefunction& wfn, const bool use_guess);

	// Mutable, a caller may change the schedule between runs
	options::SCF_settings settings;

protected:

	using scf_t = occ::qm::SCF<occ::qm::HartreeFock>;

	// The model term, nothing in the plain SCF. evaluate_model sees the current density before every
	// Fock build (timed from an SCF iteration, not from a TRAH rebuild); model_energy is the term the
	// functional adds to E and add_model_fock its derivative, added to the Fock matrix
	virtual void evaluate_model(const scf_t& scf, const bool timed) {}
	virtual double model_energy() const { return 0.0; }
	virtual void add_model_fock(occ::Mat& F) const {}
	// The model's part of the Hessian-vector product for the density change dD, added to dF
	virtual void add_model_response(const scf_t& scf, const occ::Mat& dD, occ::Mat& dF) {}
	// Without a model term the TRAH preconditioner is the diagonal alone
	virtual bool model_active() const { return false; }
	// Around the finite-difference Hessian check, which evaluates the model at other orbitals
	virtual void save_model_state() {}
	virtual void restore_model_state() {}
	// The SCF log: its head, the line of every iteration and what a converged run reports
	virtual void begin_log();
	virtual void log_iteration(const scf_t& scf);
	virtual void end_log(const scf_t& scf) {}
	// The functional the SCF descends, as the log names it
	virtual std::string functional_name() const { return "E"; }

	// Prefix of the timers and of the log lines that name the procedure, and what the log calls one run
	std::string label_ = "SCF", step_name_ = "run";
	int nmo_ = 0;
	SCF_log_writer* writer_ = nullptr;

private:

	// Two-electron integrals over the screened-in pairs, built once when they fit in memory: the Fock
	// build then contracts them instead of recomputing every quartet per iteration
	stored_eri eri_;
	bool eri_on_device_ = false;
	// Two-electron part of the Hessian products only: the gradient and energy keep the exact Fock
	// build, so an approximate Hessian can cost TRAH steps but never accuracy
	std::unique_ptr<occ::qm::HartreeFock> hessian_df_;

	occ::core::Molecule mol_;
	occ::qm::AOBasis basis_;

	// The SCF loop of one run
	bool do_SCF(double& alpha, scf_t& scf, occ::qm::Wavefunction& last_wfn, bool& has_guess);

	// Initial guess from a Hartree-Fock converged in the smaller guess basis
	void small_basis_guess(scf_t& scf);

	// Executes a single SCF iteration, true once converged
	bool SCF_iteration(scf_t& scf, double& alpha, double& quant_diff_mem, double& last_phi, occ::Mat& dm_last);

	// Checks convergence for SCF cycle
	bool SCF_convergence_check(scf_t& scf, occ::Mat& dm_last);

	// The Roothaan step and the DIIS of OCC's SCF with their matrix products on MKL
	void solve_orbitals(scf_t& scf, const occ::Mat& F) const;
	occ::Mat diis_update(scf_t& scf);
	// CDIIS over the last diis_subspace_ Fock matrices and their commutators, see cdiis_extrapolate
	occ::Mat cdiis_extrapolate(const occ::Mat& F, const occ::Mat& E);
	std::deque<occ::Mat> diis_F_, diis_E_;
	static constexpr size_t diis_subspace_ = 8;
	occ::qm::ADIIS adiis_;

	// The best orbitals of this run and their functional value. A step that climbs more than
	// rescue_rise_ above that best has lost the SCF: rescue_scf restarts it from these orbitals with
	// DIIS cleared, a stronger level shift and damping held longer, up to three times per run and
	// only from a cold start; past that, or warm-started, TRAH takes over from them (see soscf_)
	occ::qm::MolecularOrbitals best_mo_;
	// G_last_ at best_mo_.D, handed back by a rescue instead of building it again; reuse_best_G_ tells the next iteration
	occ::Mat best_G_;
	bool reuse_best_G_ = false;
	double best_phi_ = 0;
	int rescues_ = 0;
	// this run started from orbitals of an earlier one: a lost step goes straight to TRAH
	bool warm_start_ = false;
	// the stabilisers as the settings gave them, restored at the start of every run
	struct { double stop_shift, stop_damping, level_shift; } base_{ -1, -1, -1 };
	static constexpr double rescue_rise_ = 1.0;
	bool rescue_scf(scf_t& scf, const double phi, double& alpha);

	// Second-order steps (TRAH, trust-region augmented Hessian) for SCFs that DIIS leaves on a
	// plateau, where the energy stays flat and the orbital gradient stops falling. They rotate the
	// occupied into the virtual orbitals so as to lower the functional, E plus the model term, whose
	// gradient is the Fock matrix the SCF builds. Each macro step solves the augmented-Hessian
	// eigenproblem by Davidson micro-iterations with the exact Hessian-vector product (the Fock
	// response plus the response of the model term), with a level shift that makes the step fit the
	// trust radius. The orbitals are moved by the Cayley transform, and occ::qm::trust_radius_update
	// adapts the radius to how well the quadratic model predicted the change. A step that raises the
	// functional is solved again at the smaller radius from the subspace already built, and two
	// rejections at the smallest radius give the orbitals back to DIIS. The occupation is kept, so
	// TRAH cannot swap orbitals.
	// It starts at the first iteration of a warm start, when the orbital gradient has not halved in
	// trah_.patience iterations (patience_requested with `check_hessian`), with `check_hessian` once
	// the DIIS error is below trah_.start_threshold, or from best_mo_ when a step is lost and no
	// rescue is left. It stays on for the rest of the run, and the trust radius it reached carries
	// into the next run. The micro-iterations are preconditioned by up to trah_pcg_ conjugate-gradient
	// steps on the Hessian without its two-electron part and do not count towards max_iter.
	bool soscf_ = false;
	int soscf_patience_iter_ = 0;
	double soscf_patience_grad_ = 0;
	// The patience, radius and noise settings are OCC's, shared with the second-order SCF of plain -occ jobs
	occ::qm::SecondOrderSettings trah_;
	double soscf_trust_ = trah_.trust_first;
	// Consecutive rejected steps taken at the smallest radius
	int soscf_floored_ = 0;
	std::vector<occ::Vec> trah_B_, trah_HB_;
	occ::Vec soscf_kappa_, soscf_grad_, soscf_hdiag_;
	occ::Mat soscf_C_;
	double soscf_phi_ = 0, soscf_pred_ = 0;
	bool soscf_boundary_ = false;
	int trah_micro_total_ = 0, trah_pcg_total_ = 0;
	// Fock builds of this run: SCF iterations, Hessian-vector products, rebuilds at the orbitals
	struct { int iter, hessian, rebuild; } builds_{};
	static constexpr double unjudged_rise_ = 1e-6;
	// conjugate-gradient steps of the preconditioner at most, and the residual fraction they stop at
	static constexpr int trah_pcg_ = 10;
	static constexpr double trah_pcg_tol_ = 1e-2;
	// the Fock matrix of this iteration was built in full from the density of the current orbitals
	bool soscf_fresh_ = false;
	void soscf_reset();
	void soscf_step(scf_t& scf, const double phi);
	void trah_solve(scf_t& scf, const bool extend);
	occ::Vec trah_precondition(scf_t& scf, const occ::Vec& rhs, const double theta);
	occ::Vec hessian_vector(scf_t& scf, const occ::Vec& v, const bool two_electron = true);
	double rebuild_at(scf_t& scf, const occ::Mat& C_from, const occ::Vec& kappa);
	occ::Vec gradient_at(scf_t& scf, const occ::Mat& C_from, const occ::Vec& kappa);
	void check_hessian(scf_t& scf);
	occ::Vec orbital_rotation_gradient(const scf_t& scf, occ::Vec& diagonal_hessian) const;
	void rotate_orbitals(scf_t& scf, const occ::Mat& C_from, const occ::Vec& kappa) const;

	// Computes the orbital gradient for usage as a convergence criterion
	double compute_orbital_gradient(const scf_t& scf);

	// Computes convergence criteria related to the density matrix (RMSP and MaxP)
	void get_density_criteria(double& RMSP_diff, double& maxP_diff, const occ::Mat& dm, const occ::Mat& dm_last);

	// Takes care of dynamic damping
	double dynamic_damping(const scf_t& scf, const double& current_alpha, const double& e_diff, double& e_diff_mem);

	// Applies level shift to fock matrix
	void apply_level_shift(const occ::Mat& C_old, const scf_t& scf, occ::Mat& F_diis);

	// The two-electron Fock part G(D), from the stored integrals when there are any; screen skips what a small density cannot reach
	occ::Mat two_electron_fock(scf_t& scf, const bool screen) const;

	//Incremental Fock build: the two-electron part and the density it was built from
	occ::Mat G_last_, D_last_build_;
	int last_full_build_ = 0;
	double next_full_build_error_ = 0.0;
	occ::Mat eri_fock(const occ::qm::MolecularOrbitals& mo, bool screen) const;

};
