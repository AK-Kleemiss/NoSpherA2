#include "pch.h"

#include "core/eqc.h"

#include <cmath>
#include <vector>

// The reference numbers are X-analysis' own output (github.com/martinrahm/X-analysis, -m 2) for
// ORCA 6.1 HF/def2-SVP runs: H + H -> H2 and CH3 + F -> CH3F (homolytic), the two signs of
// Delta(nX) and so the two branches of the covalency index.
TEST(Eqc, CovalencyMatchesXAnalysisOnBothBranches)
{
	EXPECT_NEAR(eqc::covalency(-5.070243, 1.523598), 76.893621, 1e-5);   // H2: x < 0
	EXPECT_NEAR(eqc::covalency(0.145173, -3.454648), 95.967223, 1e-5);   // CH3F: x > 0
	EXPECT_TRUE(std::isnan(eqc::covalency(0.0, 1.0)));                   // X-analysis divides by zero
}

TEST(Eqc, ReactionTermsAreProductsMinusReactants)
{
	// Any split of the H2 Delta E and Delta(nX) over the reactants gives X-analysis' Q and covalency
	const double h = eqc::hartree2eV;
	eqc::terms a{ -0.5, -0.3, 0.0, 1 }, b{ -0.5, -0.3, 0.0, 1 }, p;
	p.E = a.E + b.E - 3.546645 / h;
	p.nX = a.nX + b.nX - 5.070243 / h;
	p.Vnn = 0.714285714;
	p.n = 2;
	const eqc::reaction r = eqc::react({ a, b }, p);
	EXPECT_NEAR(r.dE * h, -3.546645, 1e-9);
	EXPECT_NEAR(r.dnX * h, -5.070243, 1e-9);
	EXPECT_NEAR(r.dV() * h, 1.523598, 1e-9);  // Delta(Vnn - Eee) = Delta E - Delta(nX)
	EXPECT_NEAR(r.Q, 1.859177, 1e-6);
	EXPECT_NEAR(r.covalency, 76.893621, 1e-5);
	EXPECT_EQ(r.n, 2);
	EXPECT_NEAR(r.dEeeE, 100.0 * (std::abs(p.Eee() / p.E) - std::abs((a.Eee() + b.Eee()) / (a.E + b.E))), 1e-12);
}

TEST(Eqc, NuclearRepulsionSumsPairs)
{
	EXPECT_NEAR(eqc::nuclear_repulsion({ 1, 1 }, { {0, 0, 0}, {0, 0, 1.4} }), 1.0 / 1.4, 1e-14);
	// Z = 6, 1, 2 at the corners of a 3-4-5 triangle
	const double V = eqc::nuclear_repulsion({ 6, 1, 2 }, { {0, 0, 0}, {3, 0, 0}, {0, 4, 0} });
	EXPECT_NEAR(V, 6.0 / 3.0 + 12.0 / 4.0 + 2.0 / 5.0, 1e-14);
}
