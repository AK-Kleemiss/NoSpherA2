#pragma once
#include "convenience.h"
#include "structure_factors.h"
#include "XCW_solver.h"
#include "GridManager.h"
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

	// Takes the converged wavefunction and creates the tscb file
	void create_tscb(const occ::qm::Wavefunction& wfn, const double& lambda);

	// Per-lambda Gaussian halting diagnostics, see evaluate_gaussian_halting.
	std::vector<GaussianHaltEntry> gaussian_halt_history_;

private:

	// Sets up the molecule, the basis set (and optionally the density fitting basis), creates the SCF solver and returns the Hartree-Fock object
	occ::qm::HartreeFock setup_system();
	// Fits target's Fock build with the auxiliary basis `name`, its three-index integrals held when they fit in half of avail
	void density_fitting(occ::qm::HartreeFock& target, const std::string& name, const char* what, const size_t avail);

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

	options* opt;
	structure_factors sf;
	SCF_log_writer writer;
	XCW_solver scf_solver;
	GridManager tsc_grids;
	// the grid route's tscb (XG.mode 1): the first step's tsc-route table gives the scatterers, the hkl rows,
	// the asym atom of each scatterer and the hkl_enlarged entry of each row; later steps fill it from grid_sfc
	ScattererLabels tscb_scat_;
	ivec2 tscb_hkl_;
	ivec tscb_atom_, tscb_entry_;
	// F_calc and scale as save_state left them
	cvec saved_F_calc_;
	double saved_scale_ = 0.0;

};
