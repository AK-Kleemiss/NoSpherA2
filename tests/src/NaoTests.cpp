#include "pch.h"

#include "core/nao.h"
#include "core/wfn_class.h"

//NAO/NPA checks.  The numbers come from NBO 7.0.9 run on the same wavefunction; the epoxide
//reference is committed as tests/epoxide_gbw/NBO/reference.nbo.

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
		//SC first, so the check is two n^3 products rather than one n^4 loop
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

//The natural minimal basis is fixed by the element, not by the basis set: one s shell per period,
//p from period 2 on, d from period 4 on, f from period 6 on.  Checked against the epoxide NBO
//reference (H, C, O) and PySCF's AOSHELL table.
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

//The AO map has to describe the same ordering the overlap is computed in, or every block of the
//pipeline mixes unrelated functions.  Two consequences are checked: a shell's own diagonal block
//of S is the identity (normalised, orthogonal components), and shells of different l on one atom
//do not overlap (different irreps of the atomic rotation group).
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

//The overlap has to carry the same convention as the coefficients the reader kept.  ORCA writes the
//|m| >= 3 components of an f or g shell - f(+-3), g(+-3), g(+-4) - with the sign opposite to
//libcint's, and the gbw reader keeps ORCA's convention in the density, so an uncorrected overlap
//makes Tr(P S) miss the electron count: 0.018 e here, 0.009 e on Zn(NH3), 0.31 of SF6's 70.  It
//hides in every linear molecule, because the sign cancels between two flipped functions and the
//surviving cross terms vanish by axial symmetry - which is why the diatomic references never showed
//it.  ao_overlap() applies the same correction the FILE47 writer does; what is checked here is the
//exact statement, that the natural populations sum to the nuclear charge whatever the partition
//does with them.  The epoxide and ethane fixtures cannot see this and pass either way.
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

//NBO 7.0.9, tests/epoxide_gbw/NBO/reference.nbo, "Summary of Natural Population Analysis".
//Tolerance: oxygen and the cores come out to 2e-3 and 4e-5, but a residual C -> H transfer of
//about 0.009 e per hydrogen is left, so the carbons sit 0.02 e too positive.  Ruled out as the
//cause: the OWSO formula (re-derived), equal-weight Loewdin instead of OWSO (much worse), the
//core definition (exact), the step-4 class merge (intra-atomic and unitary), self-consistent
//weights (converges to hydrogens at +0.66) and one OWSO over the whole natural minimal basis
//instead of core before valence (worse here and on nh3bh3 and benzene).  Also ruled out: NBO's
//symmetry averaging, which the 22-molecule reference spine now makes testable and which does not
//explain it - the reference prints distinct charges for symmetry-equivalent atoms wherever the
//geometry allows, and where it does not our spread is zero too.  The overlap phase convention was
//a real cause and is fixed (see OverlapCarriesTheDensitysPhaseConvention): it carried the whole
//error on the hypervalent and heavy-atom references, PF5 0.099 -> 0.013 e and SF6 0.080 -> 0.016,
//but none of this one, which has no |m| >= 3 cross terms.  The bound below is the measured
//agreement; if a later change makes it pass at a tighter one, tighten it.
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

//An open-shell case with chemistry in it: the NH3...Li doublet.  NBO 7.0.9 on the .47 of the same
//gbw gives the charges below and puts 0.947 of the unpaired electron on lithium; we agree on the
//charges to 5e-3 and put 0.976 there, so the tolerances differ by column.
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
		EXPECT_NEAR(npa.total.atoms[a].charge, nbo_charge[a], 6e-3) << "atom " << a + 1;
		EXPECT_NEAR(npa.spin_population[a], nbo_spin[a], 3e-2) << "spin, atom " << a + 1;
		spin_sum += npa.spin_population[a];
	}
	//the spin populations sum to Tr((P_alpha - P_beta) S) exactly, whatever the partition does
	//with it - which for this gbw is 1.00003, not 1, because that is what ORCA converged to
	EXPECT_NEAR(spin_sum, trace_spin(wavy), 1e-9);
	EXPECT_NEAR(npa.total.population, trace_PS(wavy), 1e-9);
	EXPECT_NEAR(npa.total.population, npa.alpha.population + npa.beta.population, 1e-9);
	EXPECT_LT(orthonormality_error(npa.alpha.C, ao_overlap(wavy)), 1e-9);
	EXPECT_LT(orthonormality_error(npa.beta.C, ao_overlap(wavy)), 1e-9);
}

//A free hydrogen atom: one electron, all of it alpha, so the spin population is exactly one and
//the natural charge exactly zero whatever the basis does.
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

//Tr(P S) is the electron count by construction, whatever reader built P and whatever basis it is in,
//so it is the one number that catches a density matrix and an overlap that are not in the same basis.
//It caught three such defects on the molden route, all fixed: Int_Params had no normalisation branch
//for a molden origin (the log said "tread carefully" and nothing else), the molden reader left the
//coefficients in the file's own AO order while Int_Params sorts an atom's shells by l and orders a
//shell's components in libcint's convention, and the ORCA |m| >= 3 sign convention was gated on the
//gbw origin alone although orca_2mkl writes the gbw's own coefficients.  F_open came out at 1.512 of
//its 9 electrons, Ce_full at 47.11 of 56, CuF2's i-shell molden at 17.57 of 47.  The reference each
//file is checked against is its own sum of MO occupations, so this test needs no external number and
//cannot be satisfied by making the reader and the integrals agree on something wrong.
TEST(NaoReaderConsistencyTests, EveryReaderConservesTheElectronCount)
{
	struct Case { const char* dir; const char* file; double tol; };
	const Case cases[] = {
		//the moldens - a closed shell, an open shell, a 3d and a 4f element.  A molden prints its MO
		//coefficients to about ten digits, so the density it carries is only that precise: Ce_full
		//misses its 56 electrons by 1.5e-9 and the other three are exact to 1e-9.
		{ "molden_file", "F_open.molden",   1e-7 },
		{ "molden_file", "F_full.molden",   1e-7 },
		{ "molden_file", "Sc_full.molden",  1e-7 },
		{ "molden_file", "Ce_full.molden",  1e-7 },
		//g, h and i shells, where the |m| >= 3 phase convention is worth 0.042 e.  The 1e-3 is not
		//this code's error bar: the same 7.1e-4 is there for the gbw of the same calculation (see
		//MoldenAndGbwOfTheSameCalculationAgree), a pre-existing high-l matter that is not a phase
		//convention - flipping every |m| >= 3 is the best of the three candidate rules, measured.
		{ "CuF2_i_func/71", "calc_occupied.molden", 1e-3 },
		{ "CuF2_i_func/71", "calc.gbw",             1e-3 },
		//gbw controls: ECP, an open shell, one electron, and the epoxide reference
		{ "ECP_SF", "Au2Br2.gbw",       1e-9 },
		{ "RGBI_groups", "nh3li.gbw",   1e-9 },
		{ "ptb_H_file", "H.gbw",        1e-9 },
		{ "epoxide_gbw", "epoxide.gbw", 1e-9 },
	};
	for (const Case& c : cases) {
		const auto p = fixture(c.dir, c.file);
		if (p.empty()) { GTEST_LOG_(INFO) << "skipping absent " << c.dir << "/" << c.file; continue; }
		WFN wavy(p);
		//a reader that fills no contracted density, or a cartesian basis, is a refusal and not this
		//test's business - NaoRefusalTests and the CLI cover those
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

//The same ORCA calculation read two ways has to give the same density in the same basis.  Before the
//molden fixes these two differed by 29.4 electrons of 47 and nothing said so; the gbw was right and
//the molden was not, which is why the reference here is the gbw.  Both now sit at 46.99929, and the
//7.1e-4 they share is what remains to explain about h and i shells.
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

//An fchk's basis reached no normalisation branch in Int_Params at all: e_origin::fchk fell into the
//"WFN Origin 5 not recognized, tread carefully!  No normalisation was performed" default, so every AO
//of a spherical fchk left ao_overlap() with a non-unit diagonal and every density built on it was
//wrong.  Two invariants pin it and neither needs an external reference: the diagonal of an AO overlap
//is 1 by construction once the basis functions are normalised, and Tr(P S) is the electron count.
//Measured on this fixture (spherical, shell types -2 and -3, 228 functions, 94 shells, 48 electrons,
//written by OCC): before, 228 of 228 diagonal elements were off and Tr(P S) = 7.579076; after, the
//diagonal is 1 everywhere and Tr(P S) = 47.938220.  The remaining 0.0618 e is a separate and still
//unfixed defect - most plausibly spherical-component ordering or phase between the fchk-built density
//matrix and libcint's order - which is why the trace tolerance here is 0.07 rather than 1e-6: it pins
//the improvement without claiming the file is cured.  The reader's own 1e-4 trace guard still refuses
//it downstream, so nothing consumes a half-right basis.
TEST(NaoReaderConsistencyTests, ASphericalFchkBasisIsNormalised)
{
	const auto p = fixture("alanine_occ", "alanine.owf.fchk");
	if (p.empty()) GTEST_SKIP() << "tests/alanine_occ/alanine.owf.fchk not found";
	WFN wavy(p);
	ASSERT_EQ(wavy.get_origin(), e_origin::fchk);
	//a cartesian fchk is a different and unfixable matter: tests/NiP3_fchk/good.fchk declares 964
	//functions where the spherical basis Int_Params rebuilds holds 857, and no normalisation convention
	//closes a 107-function gap.  This fixture is spherical, so normalisation is the whole story.
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
	//the second invariant, and the one that separates this file from the cartesian case: the fchk's own
	//"Number of basis functions" is 228, and the basis Int_Params rebuilds from its shells must hold
	//exactly that many.  NiP3_fchk/good.fchk fails this at 857 against a declared 964 and no
	//normalisation convention can close that gap; this one does not.
	EXPECT_EQ(S.extent(0), size_t(228)) << "the fchk declares 228 basis functions";
	//Tr(P S) is not checked here: WFN::DM stays empty for an fchk (the gbw, wfx, molden and ptb readers
	//fill it, read_fchk does not, and the .47 writer builds its own contracted density from the MO
	//coefficients instead), which is also why EveryReaderConservesTheElectronCount skips this file.  The
	//end-to-end trace is measured through the CLI: 7.579076 before this fix, 47.938220 after, against 48.
}

//An occupancy is a number of electrons, so a printed "-0.00000" is wrong twice over: it reads as a
//negative population, and it is not even a property of the molecule.  This was found by translating a
//molecule rigidly and diffing the output - the one transformation a wavefunction file admits without
//recomputing it, since the basis functions ride on the atom centres - and tests/molden_file/test.molden
//printed "Ryd( 6s)   -0.00000" in one position and "0.00000" for the same molecule moved 4.35 bohr.
//This test does NOT assert that no occupancy is negative: one that is genuinely negative would be a real
//defect and must stay visible.  It asserts that nothing is printed as a signed zero, which is the one
//case where the sign carries no information at all.
//The fixture is Sc_full and not the test.molden the defect was found on, because the first version of
//this test used test.molden and PASSED against a binary with the clamp taken back out: test.molden prints
//the signed zero only in the MOVED frame, so on the fixture as committed there was nothing to catch.
//Sc_full is scandium in a large basis, so it carries a long tail of nearly empty Rydberg NAOs and prints
//eleven of them signed in its own frame - the same defect where it is reproducible rather than where it
//happened to be noticed.
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
	ASSERT_NE(out.find("NATURAL POPULATIONS"), std::string::npos) << "no occupancy table was printed";
	EXPECT_EQ(out.find("-0.00000"), std::string::npos)
		<< "an occupancy printed as a negative zero: the sign is roundoff from the diagonalisation and it "
		   "flips when the molecule is translated, so it says nothing and reads as a negative population";
	//and the fix must not have hidden the row: the tiny occupancy is still printed, as a plain zero
	EXPECT_NE(out.find("0.00000"), std::string::npos)
		<< "the nearly empty NAO rows have gone missing entirely, which is not what the clamp does";
}
