//The trust-region augmented-Hessian second-order step OCC's SCF hands over to once DIIS has
//reached the quadratic region (see <occ/qm/second_order_scf.h>), which is what -occ jobs from
//Olex2 get. A run that takes the second-order route has to end on the same solution as the plain
//DIIS run: for Hartree-Fock, where the Hessian-vector product is an exact Fock build of the
//density difference, for the unrestricted case, where the rotation is addressed block by block,
//and for a GGA functional, where the XC quadrature is not linear in the density and the response
//is differenced instead.
#include "pch.h"
#include <gtest/gtest.h>

#include "core/basis_set.h"
#include <occ/dft/dft.h>
#include <occ/io/xyz.h>
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <occ/qm/spinorbital.h>
#include <spdlog/spdlog.h>

namespace
{
	const char* water = "3\n\nO 0.000 0.000 0.117\nH 0.000 0.757 -0.469\nH 0.000 -0.757 -0.469\n";
	//Methyl radical: a doublet, so the unrestricted blocks have different occupations
	const char* methyl = "4\n\nC 0.000 0.000 0.000\nH 1.079 0.000 0.000\nH -0.540 0.935 0.000\nH -0.540 -0.935 0.000\n";

	occ::qm::AOBasis basis_for(const char* xyz, const char* name, const int mult)
	{
		occ::core::Molecule mol = occ::io::molecule_from_xyz_string(xyz);
		mol.set_multiplicity(mult);
		return BasisSetLibrary::get_basis_set(name)->to_AOBasis(mol.atoms());
	}

	//`engaged` reports whether the second-order step actually took over - without it a run that
	//quietly stayed with DIIS would pass the comparison against DIIS
	template <typename Proc>
	double run(Proc& proc, const occ::qm::SpinorbitalKind kind, const unsigned int mult,
	           const bool second_order, bool* engaged = nullptr)
	{
		occ::qm::SCF<Proc> scf(proc, kind);
		scf.second_order.settings.enabled = second_order;
		scf.second_order.settings.request_early = second_order;
		scf.set_charge_multiplicity(0, mult);
		scf.compute_initial_guess();
		const double energy = scf.compute_scf_energy();
		if (engaged)
			*engaged = scf.second_order.active();
		return energy;
	}
}

TEST(OccSecondOrderScf, RestrictedHartreeFockMatchesDiis)
{
	spdlog::set_level(spdlog::level::err);
	occ::qm::HartreeFock hf(basis_for(water, "def2-SVP", 1));
	const double diis = run(hf, occ::qm::SpinorbitalKind::Restricted, 1, false);
	bool engaged = false;
	const double trah = run(hf, occ::qm::SpinorbitalKind::Restricted, 1, true, &engaged);
	EXPECT_TRUE(engaged);
	EXPECT_NEAR(trah, diis, 1e-8);
}

TEST(OccSecondOrderScf, UnrestrictedHartreeFockMatchesDiis)
{
	spdlog::set_level(spdlog::level::err);
	occ::qm::HartreeFock hf(basis_for(methyl, "sto-3g", 2));
	const double diis = run(hf, occ::qm::SpinorbitalKind::Unrestricted, 2, false);
	bool engaged = false;
	const double trah = run(hf, occ::qm::SpinorbitalKind::Unrestricted, 2, true, &engaged);
	EXPECT_TRUE(engaged);
	EXPECT_NEAR(trah, diis, 1e-8);
}

//The differenced Fock response: PBE has no exact exchange, so every bit of the Hessian-vector
//product goes through the finite difference of the whole build
TEST(OccSecondOrderScf, GgaDftMatchesDiis)
{
	spdlog::set_level(spdlog::level::err);
	occ::numint::GridSettings grid;
	grid.max_angular_points = 170;
	grid.min_angular_points = 74;
	grid.radial_precision = 1e-7;
	occ::dft::DFT dft("pbe", basis_for(water, "def2-SVP", 1), grid);
	const double diis = run(dft, occ::qm::SpinorbitalKind::Restricted, 1, false);
	bool engaged = false;
	const double trah = run(dft, occ::qm::SpinorbitalKind::Restricted, 1, true, &engaged);
	EXPECT_TRUE(engaged);
	EXPECT_NEAR(trah, diis, 1e-7);
}
