//scf_backend dusk at the boundary to occ: the basis and its function order (the overlap of each library in the other's
//order), and molscf's plain Hartree-Fock against occ's, as ScfSolverTests does for SCF_solver
#include "pch.h"
#include <gtest/gtest.h>

#include "core/xcw/dusk_scf.h"
#include "core/xcw/SCF_solver.h"
#include "core/wfn/basis_set.h"
#include <occ/io/xyz.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>

#ifdef NOSPHERA2_USE_DUSK
namespace
{
	const char* water = "3\n\nO 0.000 0.000 0.117\nH 0.000 0.757 -0.469\nH 0.000 -0.757 -0.469\n";
	const char* methyl = "4\n\nC 0.000 0.000 0.000\nH 1.079 0.000 0.000\nH -0.540 0.935 0.000\nH -0.540 -0.935 0.000\n";

	struct pair_run {
		double occ_energy = 0, energy = 0, warm_energy = 0;
		bool converged = false, warm_converged = false;
		int warm_iterations = 0;
		double max_dD = 0;
	};

	//occ's SCF at an energy threshold of 1e-10 and molscf's at its tight preset, then molscf again from its own orbitals
	pair_run solve(const char* xyz, const std::string& basis_name, const occ::qm::SpinorbitalKind kind, const int multiplicity) {
		spdlog::set_level(spdlog::level::err);
		occ::core::Molecule mol = occ::io::molecule_from_xyz_string(xyz);
		mol.set_multiplicity(multiplicity);
		const occ::qm::AOBasis basis = SCF_solver::setup_basis(mol, basis_name);
		occ::qm::HartreeFock hf(basis);
		occ::qm::SCF<occ::qm::HartreeFock> reference(hf, kind);
		reference.set_charge_multiplicity(0, multiplicity);
		reference.convergence_settings.energy_threshold = 1e-10;
		pair_run r;
		r.occ_energy = reference.compute_scf_energy();
		molscf::Options o;
		o.reference = kind == occ::qm::SpinorbitalKind::Unrestricted ? molscf::Reference::Unrestricted : molscf::Reference::Restricted;
		o.convergence = molscf::Convergence::tight();
		o.integral_threshold = 1e-14;
		o.verbosity = 0;
		molscf::SCF scf(dusk::molecule(mol), dusk::basis(*BasisSetLibrary::get_basis_set(basis_name), mol.atoms()), o);
		const molscf::Result first = scf.run();
		r.converged = first.converged;
		r.energy = first.energy;
		occ::qm::Wavefunction wfn;
		wfn.basis = basis;
		wfn.atoms = mol.atoms();
		dusk::to_occ(first, wfn);
		occ::qm::MolecularOrbitals ref = reference.ctx.mo;
		ref.update_occupied_orbitals();
		ref.update_density_matrix();
		//the total density: occ puts a doublet's unpaired electron in beta (n_alpha 4, n_beta 5 for CH3), molscf in alpha
		r.max_dD = (wfn.mo.total_density_matrix() - ref.total_density_matrix()).cwiseAbs().maxCoeff();
		scf.set_guess(dusk::to_dusk(wfn));
		const molscf::Result warm = scf.run();
		r.warm_converged = warm.converged;
		r.warm_energy = warm.energy;
		r.warm_iterations = warm.iterations;
		return r;
	}
}

//The overlap molscf computes over the converted basis equals occ's, element by element: same functions, same order and phases
TEST(DuskBackendTests, OverlapMatchesOccThroughF)
{
	spdlog::set_level(spdlog::level::err);
	const occ::core::Molecule mol = occ::io::molecule_from_xyz_string(water);
	const occ::qm::AOBasis basis = SCF_solver::setup_basis(mol, "def2-tzvp");
	ASSERT_EQ(basis.l_max(), 3);
	occ::qm::HartreeFock hf(basis);
	const occ::Mat S_occ = hf.compute_overlap_matrix();
	const molscf::SCF scf(dusk::molecule(mol), dusk::basis(*BasisSetLibrary::get_basis_set("def2-tzvp"), mol.atoms()));
	ASSERT_EQ(scf.basis().nbf(), basis.nbf());
	const occ::Mat S = dusk::to_occ(scf.overlap());
	EXPECT_LT((S - S_occ).cwiseAbs().maxCoeff(), 1e-12);
	//and back: the storage change is its own inverse
	EXPECT_EQ(dusk::to_occ(dusk::to_dusk(S_occ)), S_occ);
}

TEST(DuskBackendTests, RestrictedHartreeFockMatchesOcc)
{
	const pair_run r = solve(water, "def2-tzvp", occ::qm::SpinorbitalKind::Restricted, 1);
	ASSERT_TRUE(r.converged);
	EXPECT_NEAR(r.energy, r.occ_energy, 1e-7);
	EXPECT_LT(r.max_dD, 1e-5);
}

TEST(DuskBackendTests, UnrestrictedHartreeFockMatchesOcc)
{
	const pair_run r = solve(methyl, "def2-svp", occ::qm::SpinorbitalKind::Unrestricted, 2);
	ASSERT_TRUE(r.converged);
	EXPECT_NEAR(r.energy, r.occ_energy, 1e-7);
	EXPECT_LT(r.max_dD, 1e-5);
}

//From converged orbitals, handed through occ's wavefunction and back, there is nothing left to do
TEST(DuskBackendTests, WarmRestartStaysAtTheSolution)
{
	const pair_run r = solve(water, "def2-svp", occ::qm::SpinorbitalKind::Restricted, 1);
	ASSERT_TRUE(r.warm_converged);
	EXPECT_NEAR(r.warm_energy, r.energy, 1e-8);
	EXPECT_LE(r.warm_iterations, 3);
}
#endif
