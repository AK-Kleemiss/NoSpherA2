#pragma once
#include "convenience.h"
#include "structure_factors.h"
#include "SCF_wrapper.h"
#include "GridManager.h"
#include "xcw_halting.h"

// Drives the wavefunction fitting: owns the structure factors and the SCF and is the only place that talks to both
class XCW_solver : private SCF_coupling {

public:

	// Builds the structure factors from the options, which XCW_solver keeps to hand out their parts
	XCW_solver(options& opt_in);

	// Does the XCW fitting routine
	void run();

	const structure_factors& structure() const { return sf; }

	// Sum_r Re(pre_r I_r(mu,nu)) over the fit set, the walk perturbation and the Hessian-vector
	// product share; out is nmo x nmo, symmetric, without any prefactor
	void contract_I(occ::Mat& out, const cvec& pre);

	// Distributional (Gaussian) halting criterion (see xcw_halting.h and
	// tests/P1_test/XCW_plan.md). Computes standardized residuals z_h from
	// the current F_calc/obs/scale, evaluates the Anderson-Darling
	// statistic and supporting diagnostics, logs them, and stores the
	// result for the final lambda* recommendation. Only called when
	// `gaussian_halt` is set.
	void evaluate_gaussian_halting(const double lambda);

	// Takes the converged wavefunction and creates the tscb file
	void create_tscb(const occ::qm::Wavefunction& wfn, const double& lambda);

	// Properties related to the I tensor of sf, see structure_factors::eval_I (NOT THE I TENSOR ITSELF)
	structure_factors::I_tensor* I_tens;
	// Per-lambda Gaussian halting diagnostics, see evaluate_gaussian_halting.
	std::vector<GaussianHaltEntry> gaussian_halt_history_;

private:

	// Sets up the molecule, the basis set (and optionally the density fitting basis), creates the SCF_wrapper and returns the Hartree-Fock object
	occ::qm::HartreeFock setup_system();

	// What a converged lambda step reports: the halting statistic, the table row and the tscb
	void report_lambda(const double lambda, const occ::qm::Wavefunction& wfn, const bool write_result);

	// SCF_coupling: what the SCF asks of the structure factors, see SCF_wrapper.h
	void evaluate(const dMatrix2& dm_eff) override;
	void perturbation(occ::Mat& perturb, const bool unrestricted) override;
	fit_quality quality() const override;
	occ::Mat perturbation_response(const dMatrix2& dD_eff, const double lambda) override;
	void save_state() override;
	void restore_state() override;
	bool on_device() const override;

	// Prints the full per-lambda Gaussian-halting table (XCW_log only), then
	// calls report_halting_progress_estimate(true). Called once at the end
	// of run().
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

	static void flip_high_m_phases(occ::qm::Wavefunction& w);

	options* opt;
	structure_factors sf;
	SCF_log_writer writer;
	SCF_wrapper scf_solver;
	GridManager tsc_grids;
	// F_calc and scale as save_state left them
	cvec saved_F_calc_;
	double saved_scale_ = 0.0;

};
