#pragma once
#include "convenience.h"
#include "structure_factors.h"
#include "XCW_solver.h"
#include "integrator.h"
#include "xcw_halting.h"

// Drives the X-ray constrained wavefunction fit: owns the structure factors and the SCF solver and
// is the only place that talks to both. It runs the lambda scan, extrapolates the orbitals from one
// lambda to the next, evaluates the Gaussian halting criterion and writes the results of every step.
class XCW_wrapper : private SCF_coupling {

public:

	// Builds the structure factors from the options, which XCW_wrapper keeps to hand out their parts
	XCW_wrapper(options& opt_in);

	// Does the XCW fitting routine
	void run();

	const structure_factors& structure() const { return sf; }

	// Distributional (Gaussian) halting criterion (see xcw_halting.h and tests/P1_test/XCW_plan.md):
	// standardized residuals z_h from the current F_calc, F_obs and scale, their Anderson-Darling
	// statistic and supporting diagnostics, logged and kept for the final lambda* recommendation.
	// Only called when `gaussian_halt` is set.
	void evaluate_gaussian_halting(const double lambda);

	// Writes the files of a converged lambda step: tscb (or tshb), wfn, fchk, F_calc table, NBO and RGBI as asked
	void write_step_results(const occ::qm::Wavefunction& wfn, const double lambda);

	// Per-lambda Gaussian halting diagnostics, see evaluate_gaussian_halting.
	std::vector<GaussianHaltEntry> gaussian_halt_history_;

private:

	// Sets up the molecule, the basis set (and optionally the density fitting basis), creates the SCF solver and returns the Hartree-Fock object
	occ::qm::HartreeFock setup_system();
	// Fits target's Fock build with the auxiliary basis `name`, its three-index integrals held when they fit in half of avail
	void density_fitting(occ::qm::HartreeFock& target, const std::string& name, const char* what, const size_t avail);

	// Debye-Waller and phase factors, anomalous dispersion and the I tensor (or its grid)
	void setup_structure_factors(const occ::qm::HartreeFock& hf);
	// The stored two-electron integrals and the density-fitted Hessian products of the SCF
	void setup_scf_integrals(const occ::qm::HartreeFock& hf);
	// The orbitals of an earlier run at the start lambda, when the settings ask for them
	bool read_start_guess(occ::qm::Wavefunction& wfn);

	// One SCF at lambda from the orbitals in wfn, reported when it converged
	bool run_lambda(occ::qm::HartreeFock& hf, const double lambda, occ::qm::Wavefunction& wfn, const bool use_guess, const bool write_result);
	// Reaches a lambda that failed in smaller steps from the last converged one
	bool approach_lambda(occ::qm::HartreeFock& hf, const double lambda, const double min_step, double& last_lambda, occ::qm::Wavefunction& last_wfn, occ::qm::Wavefunction& wfn);
	// Adds a step to the scan while the Gaussian halting minimum lies beyond it
	void extend_scan(const double lambda, int& planned_steps, int& extra_steps, const int max_extra_steps);
	// What a converged lambda step reports: the halting statistic, the table row and the tscb
	void report_lambda(const double lambda, const occ::qm::Wavefunction& wfn, const bool write_result);

	// SCF_coupling: what the SCF asks of the structure factors, see XCW_solver.h
	void evaluate(const dMatrix2& dm_eff) override;
	void perturbation(occ::Mat& perturb, const bool unrestricted) override;
	fit_quality quality() const override;
	occ::Mat perturbation_response(const dMatrix2& dD_eff, const double lambda) override;
	void save_state() override;
	void restore_state() override;
	bool on_device() const override;

	// The full per-lambda Gaussian halting table (SCF log only), then the final estimate; once at the end of run()
	void report_gaussian_halting_summary();

	// The recommended lambda* = argmin A^2 so far (subject to the binned-trend test), a warning when
	// that minimum sits at the last evaluated lambda, and an estimate of where the minimum lies from
	// polynomial fits of the A^2(lambda) trend, chosen by AIC (see xcw_halting.h). Every 5 lambda
	// steps during the scan (is_final false) and once at the end, from whatever history there is
	void report_halting_progress_estimate(bool is_final);

	static void flip_high_m_phases(occ::qm::Wavefunction& w);

	// -tshb: the converged density of a step as an RI fit, its scatterers keyed and ordered as the table's
	void write_tshb(WFN& wave, const tsc_block<int, cdouble>& table, const std::string& name);
	// The aux wavefunction of that fit, the restraints' partition grid and the table's atoms, made at the
	// first step and kept, since every step has the same atoms and basis
	WFN tshb_aux_{ e_origin::NOT_YET_DEFINED };
	DensityFitting::PartitionCache tshb_cache_;
	ivec tshb_atoms_;

	// TFVC, MBIS and EMBIS shares (sf.XG.ws) follow the density: refit between the solves of a step until they hold still (-1),
	// at every weight_every_-th evaluation of a solve, or (0) once at its first; at lambda 0 once, to the converged density.
	// weight_count_ counts the evaluations of the step, weight_rounds_ the solves the refits asked for
	int weight_every_ = -1, weight_count_ = 0, weight_rounds_ = 0;
	// The density and shares of the step before, for NOSPHERA2_XCW_WEIGHT_CHECK
	dMatrix2 weight_D_;
	vec2 weight_w_;
	// Refits the shares to the converged density of wfn; true while that moved F by more than NOSPHERA2_XCW_WEIGHT_FTOL of max |F|
	bool refit_moves(const occ::qm::Wavefunction& wfn, const double lambda);
	// How far the shares move F and R1(gt) over a step against what the density alone moves (NOSPHERA2_XCW_WEIGHT_CHECK)
	void weight_check(const occ::qm::Wavefunction& wfn, const double lambda);

	options* opt;
	structure_factors sf;
	SCF_log_writer writer;
	XCW_solver scf_solver;

};
