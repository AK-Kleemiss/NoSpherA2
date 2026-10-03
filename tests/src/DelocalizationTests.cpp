#include "pch.h"
#include "core/b2c.h"
#include "core/wfn_class.h"
#include "core/atoms.h"
#include "core/constants.h"
#include "core/cube.h"
#include "core/nos_math.h"
#include "core/properties.h"

#include <algorithm>
#include <filesystem>
#include <numeric>
#include <sstream>

//The delocalization indices are built from the basin overlap matrices alone, so they can be
//checked without an integration: a two-basin, two-electron bond has delta = 1 exactly, and the
//sum rule lambda(A) + sum_B delta(A,B)/2 = N(A) has to hold whatever the matrices say.
namespace
{
	//nmo orbitals, occ each, all alpha unless spins says otherwise; no primitives, the overlap
	//matrices are handed in ready
	WFN orbital_wfn(const vec &occ, const ivec &spin)
	{
		WFN wavy(e_origin::NOT_YET_DEFINED);
		wavy.push_back_atom("H", 0.0, 0.0, 0.0, 1);
		wavy.push_back_atom("H", 0.0, 0.0, 1.4, 1);
		for (size_t i = 0; i < occ.size(); i++)
			wavy.push_back_MO(static_cast<int>(i) + 1, occ[i], -0.5, spin[i]);
		return wavy;
	}

	//S[basin][mo] on the diagonal, off-diagonals zero: every orbital is split between the basins
	//by the fractions given and none of them overlaps another
	basin_overlaps diagonal_overlaps(const vec2 &fraction)
	{
		basin_overlaps ovl;
		ovl.nmo = static_cast<int>(fraction[0].size());
		ovl.mo_index.resize(ovl.nmo);
		std::iota(ovl.mo_index.begin(), ovl.mo_index.end(), 0);
		ovl.S.assign(fraction.size(), vec(ovl.triangle(), 0.0));
		for (size_t b = 0; b < fraction.size(); b++)
			for (int i = 0; i < ovl.nmo; i++)
				ovl.S[b][basin_overlaps::packed(i, i)] = fraction[b][i];
		return ovl;
	}
}

//The packing has to be symmetric and collision-free, or an off-diagonal overlap lands on the
//wrong pair and every index built from it is quietly wrong
TEST(Delocalization, PackedTriangleIsSymmetricAndDense)
{
	constexpr int n = 7;
	std::vector<bool> seen((size_t)n * (n + 1) / 2, false);
	for (int i = 0; i < n; i++)
		for (int j = 0; j <= i; j++) {
			const size_t p = basin_overlaps::packed(i, j);
			EXPECT_EQ(p, basin_overlaps::packed(j, i));
			ASSERT_LT(p, seen.size());
			EXPECT_FALSE(seen[p]);
			seen[p] = true;
		}
	for (bool s : seen) EXPECT_TRUE(s);
}

//H2 in a minimal basis: one doubly occupied orbital shared equally, delta(A,B) = 1
TEST(Delocalization, RestrictedTwoElectronBondIsOne)
{
	const WFN wavy = orbital_wfn({ 2.0 }, { 0 });
	const basin_overlaps ovl = diagonal_overlaps({ { 0.5 }, { 0.5 } });
	const delocalization_result r = delocalization_indices(wavy, ovl);
	ASSERT_EQ(r.di.size(), 1u);
	EXPECT_NEAR(r.di[0], 1.0, 1e-12);
	EXPECT_NEAR(r.lambda[0], 0.5, 1e-12);
	EXPECT_NEAR(r.population[0], 1.0, 1e-12);
	EXPECT_NEAR(r.identity_error, 0.0, 1e-12);
	//the sum rule, which is what makes the numbers mean anything
	EXPECT_NEAR(r.lambda[0] + 0.5 * r.di[0] - r.population[0], 0.0, 1e-12);
}

//The same bond written as separate alpha and beta orbitals has to give the same indices
TEST(Delocalization, UnrestrictedMatchesRestricted)
{
	const WFN wavy = orbital_wfn({ 1.0, 1.0 }, { 0, 1 });
	const basin_overlaps ovl = diagonal_overlaps({ { 0.5, 0.5 }, { 0.5, 0.5 } });
	const delocalization_result r = delocalization_indices(wavy, ovl);
	ASSERT_EQ(r.di.size(), 1u);
	EXPECT_NEAR(r.di[0], 1.0, 1e-12);
	EXPECT_NEAR(r.lambda[0], 0.5, 1e-12);
	EXPECT_NEAR(r.population[0], 1.0, 1e-12);
}

//An alpha electron on one atom and a beta electron on the other: opposite spins do not
//exchange, so nothing is shared and both electrons are localized
TEST(Delocalization, OppositeSpinsDoNotDelocalize)
{
	const WFN wavy = orbital_wfn({ 1.0, 1.0 }, { 0, 1 });
	const basin_overlaps ovl = diagonal_overlaps({ { 1.0, 0.0 }, { 0.0, 1.0 } });
	const delocalization_result r = delocalization_indices(wavy, ovl);
	ASSERT_EQ(r.di.size(), 1u);
	EXPECT_NEAR(r.di[0], 0.0, 1e-12);
	EXPECT_NEAR(r.lambda[0], 1.0, 1e-12);
	EXPECT_NEAR(r.population[0], 1.0, 1e-12);
	EXPECT_NEAR(r.lambda[0] + 0.5 * r.di[0] - r.population[0], 0.0, 1e-12);
}

//Basins that do not add up to the whole orbital are a quadrature that lost electrons, and the
//deviation of sum_A S^A from the identity is the only place that shows without a reference
TEST(Delocalization, IdentityErrorSeesLostDensity)
{
	const WFN wavy = orbital_wfn({ 2.0 }, { 0 });
	const basin_overlaps ovl = diagonal_overlaps({ { 0.4 }, { 0.5 } });
	const delocalization_result r = delocalization_indices(wavy, ovl);
	EXPECT_NEAR(r.identity_error, 0.1, 1e-12);
}

TEST(Delocalization, ReportNamesTheBasinsAndThePair)
{
	const WFN wavy = orbital_wfn({ 2.0 }, { 0 });
	const basin_overlaps ovl = diagonal_overlaps({ { 0.5 }, { 0.5 } });
	std::ostringstream out;
	report_delocalization(wavy, ovl, { "H0", "H1" }, out, 0.01);
	const std::string text = out.str();
	EXPECT_NE(text.find("Delocalization indices"), std::string::npos);
	EXPECT_NE(text.find("H0"), std::string::npos);
	EXPECT_NE(text.find("H1"), std::string::npos);
	EXPECT_NE(text.find("1.0000"), std::string::npos);
}

//The .wfn and .wfx readers have no record of where the beta orbitals begin and look for the
//place the orbital energies stop rising instead, so two exactly degenerate orbitals in a closed
//shell - the two pi lone pairs of OH- - make the reader label the second of them beta. An MO
//holding two electrons is a spatial orbital whatever that flag says; believing the flag gave
//m = 1 with a spin-orbital occupation of 2 and put every lambda and delta at exactly twice its
//value, while the population, which sees only the product m * occ, stayed right and hid it
TEST(Delocalization, DoublyOccupiedOrbitalsStayRestrictedWhateverTheSpinFlagSays)
{
	const basin_overlaps ovl = diagonal_overlaps({ { 0.5, 0.5 }, { 0.5, 0.5 } });
	const delocalization_result clean = delocalization_indices(orbital_wfn({ 2.0, 2.0 }, { 0, 0 }), ovl);
	const delocalization_result flagged = delocalization_indices(orbital_wfn({ 2.0, 2.0 }, { 0, 1 }), ovl);
	ASSERT_EQ(clean.di.size(), 1u);
	ASSERT_EQ(flagged.di.size(), 1u);
	//two pairs shared equally between two basins: N = 2, lambda = 1, delta = 2
	EXPECT_NEAR(clean.population[0], 2.0, 1e-12);
	EXPECT_NEAR(clean.lambda[0], 1.0, 1e-12);
	EXPECT_NEAR(clean.di[0], 2.0, 1e-12);
	EXPECT_NEAR(flagged.population[0], clean.population[0], 1e-12);
	EXPECT_NEAR(flagged.lambda[0], clean.lambda[0], 1e-12);
	EXPECT_NEAR(flagged.di[0], clean.di[0], 1e-12);
	EXPECT_NEAR(flagged.lambda[0] + 0.5 * flagged.di[0] - flagged.population[0], 0.0, 1e-12);
}

namespace
{
	//The QTAIM half of ELI_analysis: a density cube around the molecule, the basins the density's
	//own gradient carves out of it, and the populations and overlap matrices taken together on
	//the atomic quadrature grids. The synthetic matrices above never reach this path, which is
	//where the overlaps are actually made
	struct basin_integration
	{
		vec pop, volumes;
		double outside = 0.0;
		basin_overlaps ovl;
	};

	basin_integration integrate_qtaim(const WFN &wavy, const double spacing, const double radius, const int accuracy)
	{
		const std::vector<atom> atoms = wavy.get_atoms();
		const double pad = constants::ang2bohr(radius);
		d3 lo{ 1e30, 1e30, 1e30 }, hi{ -1e30, -1e30, -1e30 };
		for (const atom &a : atoms) {
			const d3 p = a.get_pos();
			for (int k = 0; k < 3; k++) {
				lo[k] = std::min(lo[k], p[k] - pad);
				hi[k] = std::max(hi[k], p[k] + pad);
			}
		}
		std::array<int, 3> n{};
		for (int k = 0; k < 3; k++) n[k] = static_cast<int>((hi[k] - lo[k]) / spacing) + 1;
		cube rho(n, 0, true);
		for (int k = 0; k < 3; k++) {
			rho.set_origin(k, lo[k]);
			rho.set_vector(k, k, spacing);
		}
		rho.calc_dv();
		std::ostringstream log;
		Calc_Rho(rho, wavy, radius, log, false);
		std::vector<d3> nuclei;
		for (const atom &a : atoms) nuclei.push_back(a.get_pos());
		const std::pair<cubei, std::vector<d4>> basins = topological_cube_analysis(&rho, atoms, false, true, 0.0, 1e-10, radius, 5e-3, &nuclei, &wavy);
		basin_integration out;
		out.pop = integrate_basins_on_atomic_grids(&rho, &basins.first, basins.second, wavy, accuracy, false, out.volumes, out.outside, nullptr, nullptr, 1, nullptr, &out.ovl);
		return out;
	}

	//S^A is an integral of phi_i phi_j over a piece of space, so it is positive semi-definite,
	//and the pieces add up to the whole, so S^A <= I: every eigenvalue lies in [0, 1] and
	//tr(S^2) <= tr(S). An off-diagonal element that carries too much weight breaks this while
	//leaving the diagonal - and with it the population - untouched, which is the one thing the
	//identity check cannot see
	void expect_overlaps_are_a_projection(const basin_overlaps &ovl, const double tol)
	{
		const int n = ovl.nmo;
		ASSERT_GT(n, 0);
		for (size_t b = 0; b < ovl.S.size(); b++) {
			vec A((size_t)n * n), W(n);
			double trace = 0.0, trace_sq = 0.0;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++) {
					const double s = ovl.at(static_cast<int>(b), i, j);
					A[(size_t)i * n + j] = s;
					trace_sq += s * s;
					if (i == j) trace += s;
				}
			EXPECT_LE(trace_sq, trace + tol) << "basin " << b + 1 << ": tr(S^2) above tr(S)";
			ASSERT_TRUE(try_make_Eigenvalues(A, W)) << "basin " << b + 1;
			EXPECT_GE(W.front(), -tol) << "basin " << b + 1 << " is not positive semi-definite";
			EXPECT_LE(W.back(), 1.0 + tol) << "basin " << b + 1 << " holds more than a whole orbital";
		}
	}

	//lambda(A) + sum_B delta(A,B)/2 = N(A) for every basin, and summed over the molecule
	//sum_A lambda(A) + sum_{A<B} delta(A,B) = N. Neither is a property of the code: they follow
	//from sum_A S^A = I, so a quadrature good enough to give the populations has to satisfy them
	void expect_sum_rule(const delocalization_result &r, const double total, const double tol)
	{
		double counted = 0.0;
		for (size_t b = 0; b < r.lambda.size(); b++) {
			double half = 0.0;
			for (size_t p = 0; p < r.pairs.size(); p++)
				if (r.pairs[p][0] == static_cast<int>(b) || r.pairs[p][1] == static_cast<int>(b)) half += 0.5 * r.di[p];
			EXPECT_NEAR(r.lambda[b] + half - r.population[b], 0.0, tol) << "basin " << b + 1 << " breaks the sum rule";
			counted += r.lambda[b] + half;
		}
		EXPECT_NEAR(counted, total, tol);
	}
}

//The whole path on a real wavefunction: OH- has five doubly occupied orbitals, two of them
//exactly degenerate, so it is both the smallest bonded case in tests/ and the one the spin-flag
//guess trips over. Everything here is physics the code does not get to choose
TEST(Delocalization, HydroxideIntegratesToPhysicalIndices)
{
	const std::filesystem::path wfn = nos_test_repo_root() / "tests" / "cytidine_tonto" / "OH.wfn";
	if (!std::filesystem::exists(wfn)) GTEST_SKIP() << "fixture missing: " << wfn.string();
	WFN wavy(wfn);
	ASSERT_EQ(wavy.get_ncen(), 2);
	const basin_integration in = integrate_qtaim(wavy, 0.25, 3.0, 3);
	ASSERT_EQ(in.pop.size(), 2u);
	ASSERT_EQ(in.ovl.nmo, 5);
	double total = 0.0;
	for (const double p : in.pop) total += p;
	EXPECT_NEAR(total, 10.0, 0.05) << "the quadrature lost electrons; nothing below is meaningful";

	expect_overlaps_are_a_projection(in.ovl, 5e-3);
	const delocalization_result r = delocalization_indices(wavy, in.ovl);
	EXPECT_LT(r.identity_error, 0.05);
	//the basins are labelled by their maxima, and the oxygen carries all but a fraction of an electron
	const int O = in.pop[0] > in.pop[1] ? 0 : 1, H = 1 - O;
	EXPECT_NEAR(r.population[O], in.pop[O], 1e-9) << "the AOM diagonal and the integrated population must be the same number";
	EXPECT_NEAR(r.population[H], in.pop[H], 1e-9);
	EXPECT_GT(r.population[O], 9.0);
	EXPECT_LT(r.population[H], 1.0);
	expect_sum_rule(r, total, 0.05);
	//a polar O-H bond shares about one pair; twice this was what the doubled m produced
	ASSERT_EQ(r.di.size(), 1u);
	EXPECT_GT(r.di[0], 0.5);
	EXPECT_LT(r.di[0], 1.2);
	//lambda cannot exceed the population it is a part of
	for (size_t b = 0; b < r.lambda.size(); b++) EXPECT_LE(r.lambda[b], r.population[b] + 1e-9);
}
