#pragma once
#include "XCW_solver.h"

class BasisSet;

#ifdef NOSPHERA2_USE_DUSK
#include <dusk/scf.hpp>

// `scf_backend dusk`: the SCF of the XCW fit on the dusk library. occ and dusk order spherical functions alike
// (m = -l..l, p as y z x, libcint's phases; DuskBackendTests checks the overlap to 1e-12), so between the two only the
// storage changes, occ column-major and dusk row-major, and occ's density matrices are half of dusk's total D
namespace dusk_occ {
	dusk::Molecule molecule(const occ::core::Molecule& mol);
	// The shells of set on atoms in BasisSet::to_AOBasis's order, unit-normalised spherical functions
	bite::Basis basis(const BasisSet& set, const std::vector<occ::core::Atom>& atoms);
	dusk::Matrix to_dusk(const occ::Mat& A);
	occ::Mat to_occ(const dusk::Matrix& A);
	dusk::Wavefunction to_dusk(const occ::qm::Wavefunction& wfn);
	// The orbitals of r into wfn, whose basis and atoms are set; its energy is E without the perturbation
	void to_occ(const dusk::Result& r, occ::qm::Wavefunction& wfn);
	// What dusk::Options has of the SCF settings, see dusk_solver for what it has not
	dusk::Options scf_options(const options::SCF_settings& s);
}

// The XCW criterion as a perturbation of dusk's SCF: E_p = criterion, dE_p/dD from the coupling, strength lambda
class xcw_perturbation : public dusk::Perturbation {

public:

	explicit xcw_perturbation(SCF_coupling& coupling) : coupling_(&coupling) {}
	// The strength the SCF holds it at (it tells us): at 0 the perturbation matrix is not built, as in XCW_solver
	double lambda = 0;
	void strength_changed(double strength) override { lambda = strength; }
	// The fit quality at the last evaluate, and how many there were
	fit_quality fit{};
	int evaluations = 0;
	double evaluate(const dusk::Matrix& D, dusk::Matrix& F) override;
	// At the density of the last evaluate
	bool response(const dusk::Matrix& D, const dusk::Matrix& dD, dusk::Matrix& F1) override;
	void save() override { coupling_->save_F_calc(); }
	void restore() override { coupling_->restore_F_calc(); }

private:

	SCF_coupling* coupling_;
};

// XCW_solver's place in XCW_wrapper for `scf_backend dusk`: one dusk SCF for the whole lambda scan, warm-started from
// the orbitals the wrapper hands it; its J/K builder lives as long as the scan. The per-iteration table of SCF.log comes
// from dusk's iteration callback. occ-only: damp/shift and their stop values, check_hessian, hessian_df_basis and the
// GPU Fock build; save_eri/read_eri use dusk's own file (see XCW_wrapper::setup_system)
class dusk_solver {

public:

	dusk_solver(const occ::core::Molecule& mol, const occ::qm::AOBasis& occ_basis, const BasisSet& orbital, const dusk::Options& o,
		XCW_log_writer& writer, SCF_coupling& coupling);
	// One SCF at lambda from the orbitals in wfn when use_guess; the result is left in wfn. False when it did not converge
	bool run(const double lambda, occ::qm::Wavefunction& wfn, const bool use_guess);
	// NOSPHERA2_DUSK_FD_CHECK: at the density of wfn, dE_p/dD and the response of the perturbation against central
	// differences along random symmetric directions, as relative errors to out
	void fd_check(const occ::qm::Wavefunction& wfn, std::ostream& out);
	dusk::SCF scf;
	std::shared_ptr<xcw_perturbation> xcw;

private:

	occ::qm::AOBasis basis_;
	std::vector<occ::core::Atom> atoms_;
	XCW_log_writer* writer_;
	size_t index_ = 0;
};
#else
class dusk_solver {
public:
	bool run(const double, occ::qm::Wavefunction&, const bool) { return false; }
	void fd_check(const occ::qm::Wavefunction&, std::ostream&) {}
};
#endif
