#pragma once
#include "convenience.h"
#include "structure_factors.h"
#include "SCF_wrapper.h"
#include "GridManager.h"
#include "xcw_halting.h"

class XCW_solver {

public:

	XCW_solver() = default;
	XCW_solver(structure_factors& sf);

	// Does the XCW fitting routine
	void run();

	// Calculates the perturbation matrix elements
	void calc_perturb(occ::Mat& perturb, const occ::qm::SCF<occ::qm::HartreeFock>& scf);
	// Sum_r Re(pre_r I_r(mu,nu)) over the fit set, the walk calc_perturb and the Hessian-vector
	// product share; out is nmo x nmo, symmetric, without any prefactor
	void contract_I(occ::Mat& out, const cvec& pre);

	// Distributional (Gaussian) halting criterion (see xcw_halting.h and
	// tests/P1_test/XCW_plan.md). Computes standardized residuals z_h from
	// the current F_calc/obs/scale, evaluates the Anderson-Darling
	// statistic and supporting diagnostics, logs them, and stores the
	// result for the final lambda* recommendation. Only called when
	// `gaussian_halt` is set.
	void evaluate_gaussian_halting(const double lambda);

	// Takes the SCF object from occ and creates the tscb file
	void create_tscb(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& lambda);

	// The I tensor of sf, see structure_factors::eval_I_anom_disp
	structure_factors::I_tensor* I_tens;
	// Per-lambda Gaussian halting diagnostics, see evaluate_gaussian_halting.
	std::vector<GaussianHaltEntry> gaussian_halt_history_;

private:

	// Combined function that sets up the XCW procedure, evaluates I tensor (or loads it from file), sets up the Hartree-Fock object and evaluates anomalous dispersion correction
	occ::qm::HartreeFock setup_XCW_procedure();

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
	structure_factors* sf;
	SCF_wrapper scf_solver;
	GridManager tsc_grids;

};
