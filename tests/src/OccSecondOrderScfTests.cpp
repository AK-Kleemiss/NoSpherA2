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

//The trust radius itself: the tests above only compare final energies, which a radius policy
//can get right while wasting Fock builds on the way. These drive occ::qm::trust_radius_update
//with the three numbers a macro step reports and check the properties the SCF relies on.
TEST(OccSecondOrderScf, TrustRadiusFollowsTheModelError)
{
	const occ::qm::SecondOrderSettings s;
	//a good step that the radius stopped may go twice as far next time
	EXPECT_NEAR(occ::qm::trust_radius_update(0.1, 0.1, -1e-3, -1e-3 + 1e-9, true, s), 0.2, 1e-12);
	//a good step that stopped short of the boundary leaves the radius alone: what stopped it was
	//the model, not the region, so its error says nothing about how far the region should reach.
	//Sizing the radius from that step instead collapsed it as the steps shrank towards
	//convergence, and P1 then spent seven macro steps per lambda climbing back out.
	EXPECT_EQ(occ::qm::trust_radius_update(0.5, 0.01, -1e-3, -1e-3 + 1e-9, false, s), 0.5);
	//a step whose model was badly wrong shrinks wherever it stopped: 80 % of the predicted
	//decrease missing asks for cbrt(0.1/0.8) = 0.5 of it
	EXPECT_NEAR(occ::qm::trust_radius_update(1.0, 0.5, -1.0, -0.2, false, s), 0.25, 1e-2);
	//a merely mediocre step keeps the radius it ran in. Sizing that band from the model error
	//instead held the iron case at a radius of 0.37 where the plain rule had reached 1.0, and it
	//cost ten macro steps
	EXPECT_EQ(occ::qm::trust_radius_update(0.4, 0.4, -1.0, -0.5, true, s), 0.4);
	//a shrink is bounded against the radius, so one bad step cannot collapse the region
	EXPECT_GE(occ::qm::trust_radius_update(0.4, 0.4, -1.0, -1e-9, true, s), 0.04);
	//every outcome stays inside the region's bounds
	EXPECT_LE(occ::qm::trust_radius_update(0.9, 0.9, -1.0, -1.0, true, s), s.trust_max);
	EXPECT_GE(occ::qm::trust_radius_update(1e-9, 1e-9, -1e-9, 1e-3, true, s), s.trust_min);
}

TEST(OccSecondOrderScf, RejectedStepAlwaysShrinksBelowWhatItTried)
{
	const occ::qm::SecondOrderSettings s;
	//the re-solve after a rejection runs in the subspace that produced the rejected step, so a
	//radius that did not fall below |kappa| would hand back the same step and the SCF would sit
	for (const double taken : {1e-3, 0.01, 0.1, 0.5, 1.0})
		for (const double predicted : {-1e-8, -1e-4, -1.0})
			for (const double actual : {1e-7, 1e-3, 1.0, 1e3})
				for (const bool boundary : {false, true})
					EXPECT_LT(occ::qm::trust_radius_update(1.0, taken, predicted, actual, boundary, s),
					          std::max(taken, s.trust_min * 1.000001))
						<< "taken " << taken << " predicted " << predicted << " actual " << actual;
	//nonsense from a failed solve must not produce a nonsense radius
	const double nan = std::numeric_limits<double>::quiet_NaN();
	EXPECT_GE(occ::qm::trust_radius_update(0.5, 0.5, -1e-3, nan, true, s), s.trust_min);
	EXPECT_LE(occ::qm::trust_radius_update(0.5, 0.5, -1e-3, nan, true, s), s.trust_max);
	EXPECT_GE(occ::qm::trust_radius_update(0.5, 0.0, 0.0, 0.0, false, s), s.trust_min);
}
