#include "pch.h"

#include "core/ibo.h"
#include "core/wfn_class.h"

//IAO/IBO checks on closed-shell ORCA wavefunctions. parent.molden is orca_2mkl's conversion of parent.gbw
//(ORCA 6.1.1), so both readers must give the same IBOs.

namespace {

	std::filesystem::path fixture(const std::string& dir, const std::string& file)
	{
		const auto p = nos_test_repo_root() / "tests" / dir / file;
		return std::filesystem::exists(p) ? p : std::filesystem::path{};
	}

	int count(const IBOResult& r, const IBOKind k)
	{
		return static_cast<int>(std::count(r.kind.begin(), r.kind.end(), k));
	}

	//charges sum to the molecular charge, populations to one per IBO, U orthogonal
	void check_invariants(const IBOResult& r, const double charge)
	{
		const int n = static_cast<int>(r.mos.size());
		EXPECT_NEAR(std::accumulate(r.charge.begin(), r.charge.end(), 0.0), charge, 1e-8);
		for (int k = 0; k < n; k++) {
			double iao = 0.0, mul = 0.0;
			for (size_t a = 0; a < r.iao_pop.extent(0); a++)
				iao += r.iao_pop(a, k), mul += r.mulliken(a, k);
			EXPECT_NEAR(iao, 1.0, 1e-8) << k;
			EXPECT_NEAR(mul, 1.0, 1e-8) << k;
			for (int l = 0; l < n; l++) {
				double d = 0.0;
				for (int i = 0; i < n; i++)
					d += r.U(i, k) * r.U(i, l);
				EXPECT_NEAR(d, k == l ? 1.0 : 0.0, 1e-10) << k << "," << l;
			}
		}
	}
}

//C2H4O: 3 cores, the two O lone pairs, C-C, two C-O and four C-H bonds
TEST(IboTests, EpoxideClassification)
{
	const auto p = fixture("epoxide_gbw", "epoxide.gbw");
	if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
	WFN wavy(p);
	const IBOResult r = intrinsic_bond_orbitals(wavy);
	ASSERT_EQ(r.mos.size(), 12u);
	check_invariants(r, 0.0);
	EXPECT_EQ(count(r, IBOKind::Core), 3);
	EXPECT_EQ(count(r, IBOKind::LonePair), 2);
	EXPECT_EQ(count(r, IBOKind::Bond), 7);
	for (int k = 0; k < 12; k++) {
		if (r.kind[k] == IBOKind::LonePair) EXPECT_EQ(r.centres[k][0], 0);
		if (r.kind[k] != IBOKind::Bond) continue;
		EXPECT_GT(r.iao_pop(r.centres[k][0], k) + r.iao_pop(r.centres[k][1], k), 0.85) << k;
	}
	EXPECT_GT(r.functional, r.functional_start);
	EXPECT_EQ(ibo_selection(r, "C2:C5").size(), 1u);
	EXPECT_EQ(ibo_selection(r, "1:2,O1:C5").size(), 2u);
	EXPECT_EQ(ibo_selection(r, "all,3").size(), 13u);
}

//The IBOs as primitive coefficients rotate the occupied space: the density at any point is unchanged and each
//IBO equals sum_i U(i, k) phi_i.
TEST(IboTests, IboWavefunctionKeepsTheDensity)
{
	const auto p = fixture("epoxide_gbw", "epoxide.gbw");
	if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
	WFN wavy(p);
	const IBOResult r = intrinsic_bond_orbitals(wavy);
	const WFN w = ibo_wfn(wavy, r);
	const int n = static_cast<int>(r.mos.size());
	for (const d3& x : {d3{0.5, 13.2, 1.6}, d3{0.0, 14.0, 3.0}, d3{-1.0, 13.5, 2.5}}) {
		vec phi(n);
		double rho = 0.0, rho_ibo = 0.0;
		for (int i = 0; i < n; i++) {
			phi[i] = wavy.computeMO(x, r.mos[i]);
			rho += 2.0 * phi[i] * phi[i];
		}
		for (int k = 0; k < n; k++) {
			double v = 0.0;
			for (int i = 0; i < n; i++)
				v += r.U(i, k) * phi[i];
			const double ibo = w.computeMO(x, r.mos[k]);
			EXPECT_NEAR(ibo, v, 1e-10) << k;
			rho_ibo += 2.0 * ibo * ibo;
		}
		EXPECT_NEAR(rho_ibo, rho, 1e-10);
	}
}

//CH3F, def2-SVP (d shells): gbw and the spherical molden of the same run give the same IBOs
TEST(IboTests, GbwAndMoldenAgree)
{
	const auto g = fixture("eqc_ch3f_pbe0", "parent.gbw"), m = fixture("eqc_ch3f_pbe0", "parent.molden");
	if (g.empty() || m.empty()) GTEST_SKIP() << "tests/eqc_ch3f_pbe0/parent.gbw or parent.molden not found";
	WFN wg(g), wm(m);
	const IBOResult a = intrinsic_bond_orbitals(wg), b = intrinsic_bond_orbitals(wm);
	check_invariants(a, 0.0);
	ASSERT_EQ(a.mos.size(), 9u);
	ASSERT_EQ(b.mos.size(), 9u);
	EXPECT_EQ(count(a, IBOKind::Core), 2);
	EXPECT_EQ(count(a, IBOKind::LonePair), 3);
	EXPECT_EQ(count(a, IBOKind::Bond), 4);
	EXPECT_EQ(ibo_selection(a, "C1:F2").size(), 1u);
	for (size_t at = 0; at < a.charge.size(); at++)
		EXPECT_NEAR(a.charge[at], b.charge[at], 1e-6) << at;
	EXPECT_NEAR(a.functional, b.functional, 1e-6);
	for (int k = 0; k < 9; k++) {
		EXPECT_EQ(a.kind[k], b.kind[k]) << k;
		EXPECT_EQ(a.centres[k], b.centres[k]) << k;
		EXPECT_NEAR(a.energy[k], b.energy[k], 1e-6) << k;
		for (size_t at = 0; at < a.charge.size(); at++)
			EXPECT_NEAR(a.mulliken(at, k), b.mulliken(at, k), 1e-5) << k << "," << at;
	}
}
