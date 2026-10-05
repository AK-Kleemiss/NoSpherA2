#include "pch.h"

#include "core/nao.h"
#include "core/wfn_class.h"

//NAO/NPA checks.  The numbers come from NBO 7.0.9 run on the same wavefunction.

namespace {

	std::filesystem::path fixture(const std::string& dir, const std::string& file)
	{
		const auto p = nos_test_repo_root() / "tests" / dir / file;
		return std::filesystem::exists(p) ? p : std::filesystem::path{};
	}

	//sum_i n_i = Tr(P S), the exact electron count, holds for any complete NAO set
	double trace_PS(const WFN& wavy)
	{
		const dMatrix2 P = wavy.get_dm();
		const dMatrix2 S = ao_overlap(wavy);
		double t = 0.0;
		for (size_t i = 0; i < P.extent(0); i++)
			for (size_t j = 0; j < P.extent(1); j++)
				t += P(i, j) * S(j, i);
		return t;
	}

	//Tr((P_alpha - P_beta) S), the number of unpaired electrons the density actually carries
	double trace_spin(const WFN& wavy)
	{
		const dMatrix2 P = wavy.get_dm(), Pb = wavy.get_dm_beta();
		const dMatrix2 S = ao_overlap(wavy);
		double t = 0.0;
		for (size_t i = 0; i < P.extent(0); i++)
			for (size_t j = 0; j < P.extent(1); j++)
				t += (P(i, j) - 2.0 * Pb(i, j)) * S(j, i);
		return t;
	}

	//largest |(C^T S C)_ij - delta_ij|
	double orthonormality_error(const dMatrix2& C, const dMatrix2& S)
	{
		const size_t n = C.extent(0);
		//SC first: two n^3 products instead of one n^4 loop
		std::vector<double> SC(n * n, 0.0);
		for (size_t i = 0; i < n; i++)
			for (size_t k = 0; k < n; k++) {
				const double s = S(i, k);
				if (s == 0.0) continue;
				for (size_t j = 0; j < n; j++)
					SC[i * n + j] += s * C(k, j);
			}
		double worst = 0.0;
		for (size_t a = 0; a < n; a++)
			for (size_t b = 0; b < n; b++) {
				double v = 0.0;
				for (size_t k = 0; k < n; k++)
					v += C(k, a) * SC[k * n + b];
				worst = std::max(worst, std::abs(v - (a == b ? 1.0 : 0.0)));
			}
		return worst;
	}
}

//The natural minimal basis depends on the element only: one s shell per period, p from period 2, d from 4,
//f from 6 (as in the epoxide NBO reference and PySCF's AOSHELL table).
TEST(NaoMinimalBasisTests, ShellCountsFollowThePeriodicTable)
{
	struct Case { int Z; int shell[4]; int core[4]; };
	const Case cases[] = {
		{  1, { 1, 0, 0, 0 }, { 0, 0, 0, 0 } },  //H   1s
		{  3, { 2, 0, 0, 0 }, { 1, 0, 0, 0 } },  //Li  [1s] 2s
		{  5, { 2, 1, 0, 0 }, { 1, 0, 0, 0 } },  //B   [1s] 2s 2p
		{  6, { 2, 1, 0, 0 }, { 1, 0, 0, 0 } },  //C
		{  8, { 2, 1, 0, 0 }, { 1, 0, 0, 0 } },  //O
		{ 11, { 3, 1, 0, 0 }, { 2, 1, 0, 0 } },  //Na  [1s2s2p] 3s
		{ 19, { 4, 2, 0, 0 }, { 3, 2, 0, 0 } },  //K   [..] 4s
		{ 21, { 4, 2, 1, 0 }, { 3, 2, 0, 0 } },  //Sc  [..] 4s 3d
		{ 26, { 4, 2, 1, 0 }, { 3, 2, 0, 0 } },  //Fe
		{ 31, { 4, 3, 1, 0 }, { 3, 2, 1, 0 } },  //Ga  [..3d] 4s 4p
	};
	for (const Case& c : cases) {
		int shell[4], core[4];
		natural_minimal_shells(c.Z, shell, core);
		for (int l = 0; l < 4; l++) {
			EXPECT_EQ(shell[l], c.shell[l]) << "Z = " << c.Z << ", l = " << l;
			EXPECT_EQ(core[l], c.core[l]) << "core, Z = " << c.Z << ", l = " << l;
		}
	}
}

//The AO map must follow the ordering the overlap is computed in: a shell's own diagonal block of S is the
//identity, and shells of different l on one atom do not overlap.
TEST(NaoBasisMapTests, EpoxideAoMapMatchesTheOverlap)
{
	const auto p = fixture("epoxide_gbw", "epoxide.gbw");
	if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
	WFN wavy(p);
	const std::vector<NAOBasisFunction> ao = spherical_ao_map(wavy);
	const dMatrix2 S = ao_overlap(wavy);
	ASSERT_EQ(S.extent(0), ao.size());
	for (size_t i = 0; i < ao.size(); i++)
		for (size_t j = 0; j < ao.size(); j++) {
			const bool same_shell = ao[i].atom == ao[j].atom && ao[i].l == ao[j].l && ao[i].shell == ao[j].shell;
			const bool same_atom_other_l = ao[i].atom == ao[j].atom && ao[i].l != ao[j].l;
			if (same_shell)
				EXPECT_NEAR(S(i, j), (ao[i].m == ao[j].m) ? 1.0 : 0.0, 1e-10) << i << "," << j;
			else if (same_atom_other_l)
				EXPECT_NEAR(S(i, j), 0.0, 1e-10) << i << "," << j;
		}
}

TEST(NaoBasisMapTests, OverlapCarriesTheDensitysPhaseConvention)
{
	const auto p = fixture("RGBI_groups", "nh3bh3.gbw");
	if (p.empty()) GTEST_SKIP() << "tests/RGBI_groups/nh3bh3.gbw not found";
	WFN wavy(p);
	const NPAResult npa = natural_population_analysis(wavy);
	EXPECT_NEAR(trace_PS(wavy), 18.0, 1e-9);
	EXPECT_NEAR(npa.total.population, 18.0, 1e-9);
	double charge_sum = 0.0;
	for (const NAOAtom& a : npa.total.atoms)
		charge_sum += a.charge;
	EXPECT_NEAR(charge_sum, 0.0, 1e-9);
}

TEST(NaoBasisMapTests, PrintedComponentsFollowLibcintAOOrder)
{
	const auto p = fixture("RGBI_groups", "nh3bh3.gbw");
	if (p.empty()) GTEST_SKIP() << "tests/RGBI_groups/nh3bh3.gbw not found";
	WFN wavy(p);
	const auto ao = spherical_ao_map(wavy);
	int shells = 0;
	for (int i = 0; i + 2 < ao.size(); i++) {
		if (ao[i].l != 1 || ao[i].m != 2) continue;
		EXPECT_EQ(ao[i + 1].m, 0);
		EXPECT_EQ(ao[i + 2].m, 1);
		shells++;
	}
	EXPECT_GT(shells, 0);
}

//A NAO's m is its AOs' m, which shell_label and lang_label print: the largest AO coefficient of a core or
//valence NAO sits on its own atom, l and m.
TEST(NaoBasisMapTests, NaoComponentIsTheComponentOfItsAOs)
{
	const auto p = fixture("RGBI_groups", "nh3bh3.gbw");
	if (p.empty()) GTEST_SKIP() << "tests/RGBI_groups/nh3bh3.gbw not found";
	WFN wavy(p);
	const auto ao = spherical_ao_map(wavy);
	const NPAResult npa = natural_population_analysis(wavy);
	int checked = 0;
	for (size_t k = 0; k < npa.total.orbitals.size(); k++) {
		const NAO& o = npa.total.orbitals[k];
		if (o.l == 0 || o.type == NAOClass::Rydberg) continue;
		size_t best = 0;
		for (size_t i = 1; i < ao.size(); i++)
			if (std::abs(npa.total.C(i, k)) > std::abs(npa.total.C(best, k))) best = i;
		EXPECT_EQ(ao[best].atom, o.atom) << "NAO " << k;
		EXPECT_EQ(ao[best].l, o.l) << "NAO " << k;
		EXPECT_EQ(ao[best].m, o.m) << "NAO " << k;
		checked++;
	}
	EXPECT_GE(checked, 6) << "the N and B 2p NAOs";
}

TEST(NaoEpoxideTests, OccupanciesSumToTheElectronCountAndTheTransformIsOrthogonal)
{
	const auto p = fixture("epoxide_gbw", "epoxide.gbw");
	if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
	WFN wavy(p);
	const NPAResult npa = natural_population_analysis(wavy);
	EXPECT_NEAR(npa.total.population, 24.0, 1e-9);
	EXPECT_NEAR(npa.total.population, trace_PS(wavy), 1e-9);
	EXPECT_NEAR(npa.total.core + npa.total.valence + npa.total.rydberg, npa.total.population, 1e-9);
	for (const NAO& o : npa.total.orbitals) {
		EXPECT_GT(o.occupation, -1e-6) << "negative occupancy";
		EXPECT_LT(o.occupation, 2.0 + 1e-6) << "occupancy above the Pauli limit";
	}
	EXPECT_LT(orthonormality_error(npa.total.C, ao_overlap(wavy)), 1e-9);
}

TEST(NaoEpoxideTests, NaturalChargesMatchNbo7)
{
	const auto p = fixture("epoxide_gbw", "epoxide.gbw");
	if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
	WFN wavy(p);
	const NPAResult npa = natural_population_analysis(wavy);
	const double nbo_charge[7] = { -0.56088, -0.04009, 0.16239, 0.15792, -0.01244, 0.14629, 0.14681 };
	const double nbo_core = 5.99986, nbo_valence = 17.94223, nbo_rydberg = 0.05791;
	ASSERT_EQ(npa.total.atoms.size(), 7u);
	for (size_t a = 0; a < 7; a++)
		EXPECT_NEAR(npa.total.atoms[a].charge, nbo_charge[a], 2.5e-2)
			<< npa.total.atoms[a].label << " " << a + 1;
	EXPECT_NEAR(npa.total.core, nbo_core, 1e-4);
	EXPECT_NEAR(npa.total.valence, nbo_valence, 1.5e-2);
	EXPECT_NEAR(npa.total.rydberg, nbo_rydberg, 1.5e-2);
}

TEST(NaoOpenShellTests, Nh3LiChargesAndSpinFollowNbo7)
{
	const auto p = fixture("RGBI_groups", "nh3li.gbw");
	if (p.empty()) GTEST_SKIP() << "tests/RGBI_groups/nh3li.gbw not found";
	WFN wavy(p);
	const NPAResult npa = natural_population_analysis(wavy);
	ASSERT_TRUE(npa.spin_resolved);
	ASSERT_EQ(npa.total.atoms.size(), 5u);
	const double nbo_charge[5] = { -1.20517, 0.39535, 0.39536, 0.39535, 0.01912 };
	const double nbo_spin[5] = { 0.03529, 0.00575, 0.00575, 0.00575, 0.94747 };
	double spin_sum = 0.0;
	for (size_t a = 0; a < 5; a++) {
		EXPECT_NEAR(npa.total.atoms[a].charge, nbo_charge[a], 3e-3) << "atom " << a + 1;
		EXPECT_NEAR(npa.spin_population[a], nbo_spin[a], 1e-3) << "spin, atom " << a + 1;
		spin_sum += npa.spin_population[a];
	}
	//the spin populations sum to Tr((P_alpha - P_beta) S) exactly, which is 1.00003 for this gbw as ORCA converged it
	EXPECT_NEAR(spin_sum, trace_spin(wavy), 1e-9);
	EXPECT_NEAR(npa.total.population, trace_PS(wavy), 1e-9);
	EXPECT_NEAR(npa.total.population, npa.alpha.population + npa.beta.population, 1e-9);
	for (size_t i = 0; i < npa.total.orbitals.size(); i++)
		EXPECT_NEAR(npa.total.orbitals[i].occupation,
		            npa.alpha.orbitals[i].occupation + npa.beta.orbitals[i].occupation, 1e-9);
	EXPECT_LT(orthonormality_error(npa.total.C, ao_overlap(wavy)), 1e-9);
	EXPECT_LT(orthonormality_error(npa.alpha.C, ao_overlap(wavy)), 1e-9);
	EXPECT_LT(orthonormality_error(npa.beta.C, ao_overlap(wavy)), 1e-9);
}

//A free H atom: one alpha electron, so spin population one and charge zero whatever the basis.
TEST(NaoOpenShellTests, HydrogenAtomCarriesOneUnpairedElectron)
{
	const auto p = fixture("ptb_H_file", "H.gbw");
	if (p.empty()) GTEST_SKIP() << "tests/ptb_H_file/H.gbw not found";
	WFN wavy(p);
	if (!wavy.get_is_unrestricted() || wavy.get_dm_beta().extent(0) == 0)
		GTEST_SKIP() << "H.gbw did not come back spin-resolved";
	const NPAResult npa = natural_population_analysis(wavy);
	ASSERT_TRUE(npa.spin_resolved);
	ASSERT_EQ(npa.total.atoms.size(), 1u);
	EXPECT_NEAR(npa.total.population, 1.0, 1e-9);
	EXPECT_NEAR(npa.total.atoms[0].charge, 0.0, 1e-9);
	EXPECT_NEAR(npa.spin_population[0], 1.0, 1e-9);
	EXPECT_NEAR(npa.alpha.population, 1.0, 1e-9);
	EXPECT_NEAR(npa.beta.population, 0.0, 1e-9);
}

//Tr(P S) checks that a reader and its overlap use the same AO basis.
TEST(NaoReaderConsistencyTests, EveryReaderConservesTheElectronCount)
{
	struct Case { const char* dir; const char* file; double tol; };
	const Case cases[] = {
		//A molden prints MO coefficients to about ten digits, so its density is only that precise.
		{ "molden_file", "F_open.molden",   1e-7 },
		{ "molden_file", "F_full.molden",   1e-7 },
		{ "molden_file", "Sc_full.molden",  1e-7 },
		{ "molden_file", "Ce_full.molden",  1e-7 },
		//High-l ORCA phase conventions are checked against the matching GBW.
		{ "CuF2_i_func/71", "calc_occupied.molden", 1e-3 },
		{ "CuF2_i_func/71", "calc.gbw",             1e-3 },
		//gbw controls
		{ "ECP_SF", "Au2Br2.gbw",       1e-9 },
		{ "RGBI_groups", "nh3li.gbw",   1e-9 },
		{ "ptb_H_file", "H.gbw",        1e-9 },
		{ "epoxide_gbw", "epoxide.gbw", 1e-9 },
	};
	for (const Case& c : cases) {
		const auto p = fixture(c.dir, c.file);
		if (p.empty()) { GTEST_LOG_(INFO) << "skipping absent " << c.dir << "/" << c.file; continue; }
		WFN wavy(p);
		//no contracted density or a cartesian basis is a refusal, which NaoRefusalTests and the CLI cover
		if (wavy.get_dm().extent(0) == 0 || wavy.get_d_f_switch()) {
			GTEST_LOG_(INFO) << "no spherical contracted density in " << c.file;
			continue;
		}
		double occ = 0.0;
		for (int i = 0; i < wavy.get_nmo(); i++)
			occ += wavy.get_MO_occ(i);
		EXPECT_NEAR(trace_PS(wavy), occ, c.tol) << c.dir << "/" << c.file;
	}
}

TEST(NaoReaderConsistencyTests, MoldenAndGbwOfTheSameCalculationAgree)
{
	const auto g = fixture("CuF2_i_func/71", "calc.gbw");
	const auto m = fixture("CuF2_i_func/71", "calc_occupied.molden");
	if (g.empty() || m.empty()) GTEST_SKIP() << "tests/CuF2_i_func/71 fixtures not found";
	WFN gbw(g), mol(m);
	ASSERT_EQ(gbw.get_origin(), e_origin::gbw);
	ASSERT_EQ(mol.get_origin(), e_origin::molden);
	EXPECT_NEAR(trace_PS(mol), trace_PS(gbw), 1e-4);
	const NPAResult a = natural_population_analysis(gbw), b = natural_population_analysis(mol);
	ASSERT_EQ(a.total.atoms.size(), b.total.atoms.size());
	for (size_t i = 0; i < a.total.atoms.size(); i++)
		EXPECT_NEAR(b.total.atoms[i].charge, a.total.atoms[i].charge, 5e-3) << "atom " << i + 1;
}

TEST(NaoReaderConsistencyTests, ASphericalFchkBasisIsNormalised)
{
	const auto p = fixture("alanine_occ", "alanine.owf.fchk");
	if (p.empty()) GTEST_SKIP() << "tests/alanine_occ/alanine.owf.fchk not found";
	WFN wavy(p);
	ASSERT_EQ(wavy.get_origin(), e_origin::fchk);
	//A cartesian fchk (NiP3_fchk/good.fchk) declares more functions than the spherical basis Int_Params rebuilds,
	//which no normalisation closes; this fixture is spherical.
	ASSERT_FALSE(wavy.get_d_f_switch()) << "fixture is no longer spherical";
	const dMatrix2 S = ao_overlap(wavy);
	ASSERT_GT(S.extent(0), size_t(0));
	int off = 0;
	double worst = 0.0;
	for (size_t i = 0; i < S.extent(0); i++) {
		const double d = std::abs(S(i, i) - 1.0);
		worst = std::max(worst, d);
		if (d > 1e-8) off++;
	}
	EXPECT_EQ(off, 0) << off << " of " << S.extent(0) << " AOs are not normalised, worst |S_ii - 1| = "
		<< worst << " - e_origin::fchk has fallen out of its normalisation branch again";
	//the basis Int_Params rebuilds from the shells must hold exactly the fchk's declared count
	EXPECT_EQ(S.extent(0), size_t(228)) << "the fchk declares 228 basis functions";
	//No Tr(P S) here: read_fchk leaves WFN::DM empty, the .47 writer builds its own density.
}

TEST(NaoPrintTests, NoOccupancyIsPrintedAsANegativeZero)
{
	const auto p = fixture("molden_file", "Sc_full.molden");
	if (p.empty()) GTEST_SKIP() << "tests/molden_file/Sc_full.molden not found";
	WFN wavy(p);
	ASSERT_FALSE(wavy.get_d_f_switch()) << "NPA needs a spherical basis; fixture is no longer spherical";
	const NPAResult r = natural_population_analysis(wavy);
	std::ostringstream os;
	print_npa(r, os);
	const std::string out = os.str();
	ASSERT_NE(out.find("Natural atomic orbital occupancies"), std::string::npos) << "no occupancy table was printed";
	EXPECT_EQ(out.find("-0.00000"), std::string::npos)
		<< "an occupancy printed as a negative zero: the sign is roundoff from the diagonalisation and it "
		   "flips when the molecule is translated, so it says nothing and reads as a negative population";
	EXPECT_NE(out.find("0.00000"), std::string::npos)
		<< "the nearly empty NAO rows have gone missing entirely, which is not what the clamp does";
}
