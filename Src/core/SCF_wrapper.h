#pragma once
#include "convenience.h"
#include "structure_factors.h"
#include "stored_eri.h"
#include <occ/qm/hf.h>
#include <occ/qm/second_order_scf.h>
#include <deque>

class XCW_solver;

class SCF_wrapper {

public:

	SCF_wrapper() = default;
	SCF_wrapper(XCW_solver& xcw_in, structure_factors& sf_in);

	// Sets up a molecule object from the asym_atoms
	void setup_SCF_mol(occ::core::Molecule& mol);

	// Sets up the basis set with a previously generated molecule and basis set from JKFit, where the Olex2 basis sets are now located
	void setup_basis(occ::core::Molecule& mol, std::string& basis_set_name, occ::qm::AOBasis& occ_basis_set);

	// Executes a single SCF solver (for specific lambda step)
	bool do_SCF(const double& lambda, double& alpha, occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::qm::Wavefunction& last_wfn, bool& has_guess, bool write_result = true);

	//Two-electron integrals over the screened-in pairs, built once per run when they fit in
	//memory: the Fock build then contracts them instead of recomputing every quartet per iteration
	stored_eri eri_;
	bool eri_on_device_ = false;

	std::ofstream SCF_log;

private:

	void small_basis_guess(occ::qm::SCF<occ::qm::HartreeFock>& scf);

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

	// Builds the density matrix to use for structure factor calculations
	void build_effective_dm(const occ::qm::SCF<occ::qm::HartreeFock>& scf, dMatrix2& dm_ref, const occ::Mat& dm_old);

	//Incremental Fock build: the two-electron part and the density it was built from
	occ::Mat G_last_, D_last_build_;
	int last_full_build_ = 0;
	double next_full_build_error_ = 0.0;
	occ::Mat eri_fock(const occ::qm::MolecularOrbitals& mo, bool screen) const;

	options* opt;
	structure_factors* sf;
	// The solver the SCF takes its perturbation from
	XCW_solver* xcw;
};
