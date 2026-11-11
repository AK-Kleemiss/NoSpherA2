#include "pch.h"
#include "core/b2c.h"
#include "core/ci_rdm2.h"
#include "core/wfn_class.h"

#include <cmath>
#include <numeric>

//The 2-RDM of a printed determinant CI: phases, normalisation and the exact delocalization indices built
//from it.  The strings are ORCA's ("2" both spins, "u" alpha, "d" beta, "0" empty), spin orbitals 1a 1b 2a 2b.
namespace
{
	CiRdm2 rdm_of(const std::vector<std::pair<std::string, double>> &dets)
	{
		CiRdm2 r;
		ci_rdm2_from_dets(dets, r);
		return r;
	}

	//sum_tv Gamma_ttvv = N(N-1)
	double trace(const CiRdm2 &r)
	{
		double s = 0.0;
		for (int t = 0; t < r.n; t++)
			for (int v = 0; v < r.n; v++) s += r.G[r.at(t, t, v, v)];
		return s;
	}

	//rho_2^ab(r,r) with orbital values phi at one point, active orbitals only
	double on_top(const CiRdm2 &r, const vec &phi)
	{
		double s = 0.0;
		for (int t = 0; t < r.n; t++)
			for (int u = 0; u < r.n; u++)
				for (int v = 0; v < r.n; v++)
					for (int w = 0; w < r.n; w++) s += r.Gab[r.at(t, u, v, w)] * phi[t] * phi[u] * phi[v] * phi[w];
		return s;
	}

	WFN orbital_wfn(const vec &occ)
	{
		WFN wavy(e_origin::NOT_YET_DEFINED);
		wavy.push_back_atom("H", 0.0, 0.0, 0.0, 1);
		wavy.push_back_atom("H", 0.0, 0.0, 1.4, 1);
		for (size_t i = 0; i < occ.size(); i++) wavy.push_back_MO(static_cast<int>(i) + 1, occ[i], -0.5, 0);
		return wavy;
	}

	basin_overlaps diagonal_overlaps(const vec2 &fraction)
	{
		basin_overlaps ovl;
		ovl.nmo = static_cast<int>(fraction[0].size());
		ovl.mo_index.resize(ovl.nmo);
		std::iota(ovl.mo_index.begin(), ovl.mo_index.end(), 0);
		ovl.S.assign(fraction.size(), vec(ovl.triangle(), 0.0));
		for (size_t b = 0; b < fraction.size(); b++)
			for (int i = 0; i < ovl.nmo; i++) ovl.S[b][basin_overlaps::packed(i, i)] = fraction[b][i];
		return ovl;
	}
}

//A closed-shell determinant: gamma is diagonal, Gamma traces to N(N-1), and the on-top density is rho_a rho_b
TEST(CiRdm2, SingleDeterminantIsMeanField)
{
	const CiRdm2 r = rdm_of({ { "20", 1.0 } });
	ASSERT_EQ(r.n, 2);
	EXPECT_NEAR(r.g1[0][0], 1.0, 1e-14);
	EXPECT_NEAR(r.g1[1][0], 1.0, 1e-14);
	EXPECT_NEAR(r.g1[0][3], 0.0, 1e-14);
	EXPECT_NEAR(trace(r), 2.0, 1e-14);
	const double a = 0.7, b = -0.3;
	EXPECT_NEAR(on_top(r, { a, b }), a * a * a * a, 1e-14);
}

//Two electrons in two orbitals, Ms = 0: the spatially symmetric singlet has on-top 2 a^2 b^2, the triplet
//none at all.  Alpha-first the singlet is |0a 1b| + |1a 0b|; ORCA orders "du" as 0b 1a, one
//swap from alpha-first, so the singlet reads [ud] +c, [du] -c.  A sign slip in the reordering swaps the two.
TEST(CiRdm2, OpenShellSingletAndTripletPhases)
{
	const double c = std::sqrt(0.5), a = 0.4, b = 0.9;
	const CiRdm2 singlet = rdm_of({ { "ud", c }, { "du", -c } });
	const CiRdm2 triplet = rdm_of({ { "ud", c }, { "du", c } });
	EXPECT_NEAR(on_top(singlet, { a, b }), 2.0 * a * a * b * b, 1e-14);
	EXPECT_NEAR(on_top(triplet, { a, b }), 0.0, 1e-14);
	for (const CiRdm2 *r : { &singlet, &triplet }) {
		EXPECT_NEAR(trace(*r), 2.0, 1e-14);
		EXPECT_NEAR(r->g1[0][0] + r->g1[1][0], 1.0, 1e-14);
		EXPECT_NEAR(r->g1[0][1], 0.0, 1e-14);
	}
}

//Same-spin pairs: a triplet Ms = 1 has no alpha-beta part, and Gamma still traces to N(N-1)
TEST(CiRdm2, HighSpinHasNoAlphaBetaPairs)
{
	const CiRdm2 r = rdm_of({ { "uu", 1.0 } });
	EXPECT_NEAR(trace(r), 2.0, 1e-14);
	for (const double x : r.Gab) EXPECT_EQ(x, 0.0);
}

//An idempotent 1-RDM handed through the exact path has a zero cumulant, so it has to repeat Mueller
TEST(CiRdm2, ExactDelocalizationOfADeterminantIsMueller)
{
	const CiRdm2 r = rdm_of({ { "20", 1.0 } });
	const WFN wavy = orbital_wfn({ 2.0, 0.0 });
	basin_overlaps ovl = diagonal_overlaps({ { 0.5 }, { 0.5 } });
	const delocalization_result m = delocalization_indices(wavy, ovl);
	//The AOM keeps only occupied orbitals; the empty active one has no slot and must not matter
	const delocalization_result e = delocalization_indices(wavy, ovl, &r);
	EXPECT_TRUE(e.exact);
	ASSERT_EQ(e.di.size(), 1u);
	EXPECT_NEAR(e.di[0], m.di[0], 1e-14);
	EXPECT_NEAR(e.lambda[0], m.lambda[0], 1e-14);
}

//Two-configuration H2: whatever the overlaps, lambda(A) + delta(A,B)/2 = N(A) has to hold exactly, and a
//basin split with S^A + S^B = 1 on both orbitals gives the closed form below
TEST(CiRdm2, TwoConfigurationSumRule)
{
	const double c1 = 0.95, c2 = -std::sqrt(1.0 - 0.95 * 0.95);
	const CiRdm2 r = rdm_of({ { "20", c1 }, { "02", c2 } });
	const WFN wavy = orbital_wfn({ 2.0 * c1 * c1, 2.0 * c2 * c2 });
	const double f0 = 0.5, f1 = 0.3;
	const basin_overlaps ovl = diagonal_overlaps({ { f0, f1 }, { 1.0 - f0, 1.0 - f1 } });
	const delocalization_result e = delocalization_indices(wavy, ovl, &r);
	ASSERT_EQ(e.di.size(), 1u);
	for (int b = 0; b < 2; b++) EXPECT_NEAR(e.lambda[b] + 0.5 * e.di[0] - e.population[b], 0.0, 1e-12);
	//rho_2(1,2) = 2 (c1 phi0 phi0 + c2 phi1 phi1)^2; with S^A_01 = 0 the cross term drops and
	//delta = 2 [N(A) N(B) - 2 c1^2 f0 (1 - f0) - 2 c2^2 f1 (1 - f1)]
	const double NA = 2.0 * (c1 * c1 * f0 + c2 * c2 * f1), NB = 2.0 - NA;
	const double expect = 2.0 * (NA * NB - 2.0 * c1 * c1 * f0 * (1.0 - f0) - 2.0 * c2 * c2 * f1 * (1.0 - f1));
	EXPECT_NEAR(e.di[0], expect, 1e-12);
}
