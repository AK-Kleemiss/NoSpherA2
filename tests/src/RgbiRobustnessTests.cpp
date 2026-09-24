#include "pch.h"
#include "core/wfn_class.h"
#include "core/bondwise_analysis.h"

#include <occ/core/parallel.h>

#include <cmath>
#include <sstream>
#include <string>

//RGBI robustness: the ANO route builds a free atom per centre and keeps a fixed number of its
//natural orbitals.  Both halves of that sentence had a defect that only showed as a number that
//moved between two identical runs, so these tests pin the reproducibility and the two invariants
//that make it hold - a free atom in its ground state, and a subspace of occupied orbitals only.

namespace {

	struct CoutCapture {
		std::ostringstream buffer;
		std::streambuf *old;
		CoutCapture() : old(std::cout.rdbuf(buffer.rdbuf())) {}
		~CoutCapture() { std::cout.rdbuf(old); }
		std::string str() const { return buffer.str(); }
	};

	//tests/TFVC/water.gbw is water plus a non-bonded helium: the one fixture in the tree whose
	//free atoms include a closed-shell period-1 element.
	std::filesystem::path water_he_fixture()
	{
		const auto p = nos_test_repo_root() / "tests" / "TFVC" / "water.gbw";
		return std::filesystem::exists(p) ? p : std::filesystem::path{};
	}

	std::string rgbi_ano_output(const bool EVs)
	{
		const auto p = water_he_fixture();
		if (p.empty())
			return {};
		CoutCapture cap;
		WFN wavy(p);
		Roby_information roby(wavy, {}, true, true, EVs, false);
		return cap.str();
	}

	//tests/RGBI holds a Tonto job: the water archive Tonto's own Roby-Gould analysis ran on, and
	//its output in tests/RGBI/stdout.  That output is the only external RGBI reference in the tree.
	std::filesystem::path tonto_water_archive()
	{
		const auto p = nos_test_repo_root() / "tests" / "RGBI" / "h2o.MOs,r";
		return std::filesystem::exists(p) ? p : std::filesystem::path{};
	}

	double value_after(const std::string &text, const std::string &key)
	{
		const size_t at = text.find(key);
		if (at == std::string::npos)
			return std::nan("");
		std::istringstream in(text.substr(at + key.size()));
		double value = std::nan("");
		in >> value;
		return value;
	}
}

//A non-bonded closed-shell atom has to get its own electrons back.  Helium's free atom was built
//as a triplet - 1s(1)2s(1), an excited state with two exactly degenerate natural orbitals - so the
//rank-1 atomic subspace cut through the degeneracy and kept whichever of the two the threaded
//eigensolver happened to return first: this population came out anywhere between 0.01 and 1.80.
TEST(RgbiRobustnessTests, NonBondedHeliumKeepsItsTwoElectrons)
{
	const std::string out = rgbi_ano_output(false);
	if (out.empty())
		GTEST_SKIP() << "tests/TFVC/water.gbw not found";
	EXPECT_NEAR(value_after(out, "Population of atom 3: "), 2.0, 5e-3);
	//and the whole-system projection then accounts for all but a thousandth of the 12 electrons
	EXPECT_NEAR(value_after(out, "Number of electrons in Roby Analysis:  "), 11.899, 5e-3);
}

//The defect's signature was a number that changed between two identical runs, so that is what is
//asserted: same process, same input, same result to every digit that is printed.
TEST(RgbiRobustnessTests, AnoPopulationsAreReproducible)
{
	const std::string first = rgbi_ano_output(false);
	if (first.empty())
		GTEST_SKIP() << "tests/TFVC/water.gbw not found";
	const std::string second = rgbi_ano_output(false);
	for (int atom = 0; atom < 4; atom++) {
		const std::string key = "Population of atom " + std::to_string(atom) + ": ";
		EXPECT_DOUBLE_EQ(value_after(first, key), value_after(second, key)) << atom;
	}
	EXPECT_DOUBLE_EQ(value_after(first, "Number of electrons in Roby Analysis:  "),
		value_after(second, "Number of electrons in Roby Analysis:  "));
}

//The rank of an atomic subspace is fixed by the element, and for hydrogen in a triple-zeta basis it
//is 1 of 14: the other 13 natural orbitals are empty.  Keeping an empty one adds an arbitrary
//direction out of a degenerate null space, which is the second way this analysis lost its
//reproducibility, so no kept orbital may have a vanishing occupation.
TEST(RgbiRobustnessTests, KeptAnoOrbitalsAreOccupied)
{
	const std::string out = rgbi_ano_output(true);
	if (out.empty())
		GTEST_SKIP() << "tests/TFVC/water.gbw not found";
	int kept = 0;
	std::istringstream in(out);
	std::string line;
	while (std::getline(in, line)) {
		if (line.find("  kept") == std::string::npos)
			continue;
		kept++;
		const double occupation = std::stod(line);
		EXPECT_GT(occupation, 1e-8) << line;
	}
	EXPECT_EQ(kept, 8); //O 1s 2s 2p(3), three one-orbital atoms: H, H, He
	//a subspace that still cut through a degeneracy would say so
	EXPECT_EQ(out.find("cuts through a degenerate"), std::string::npos);
}

//A centre the shipped minimal basis does not reach - cerium - takes occ's other guess route, and
//that route used to corrupt the heap for an unrestricted free atom and abort in a malloc inside
//libcint.  It is fixed by starting that one case from the core Hamiltonian instead, and the case
//has to stay that one case: starting *every* free atom there moved the light-atom populations
//above by up to 0.84 electrons, which is a different converged atom and not a better one.  So this
//pins both halves - the heavy atom completes and really computes its own free atom, and the light
//atoms are left to the automatic choice, which the goldens in BondwiseTests.cpp check.
TEST(RgbiRobustnessTests, CeriumFreeAtomRunsAndIsNotFallenBackOn)
{
	const auto p = nos_test_repo_root() / "tests" / "molden_file" / "Ce_full.molden";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/molden_file/Ce_full.molden not found";
	std::string out;
	{
		CoutCapture cap;
		WFN wavy(p);
		Roby_information roby(wavy, {}, true, true, false, false);
		out = cap.str();
	}
	//the free atom was computed, not replaced by the molecular local orbitals
	EXPECT_NE(out.find("ANO fallback summary: no atom-level fallbacks were needed."), std::string::npos);
	//and the projected population is a population: finite, positive, and short of Z only by what
	//the ANO cutoff leaves outside, which the run reports separately
	const double population = value_after(out, "Population of atom 0: ");
	EXPECT_TRUE(std::isfinite(population));
	EXPECT_GT(population, 0.5 * 58.0);
	EXPECT_LT(population, 58.0 + 1e-6);
}

//The external reference.  tests/RGBI/stdout is Tonto 26.01.05's own Roby-Gould output for the
//archive next to it, and it states the option set it used: NAOs, no spherical averaging, def2-SVP,
//Cartesian d functions.  Run with those options, every number NoSpherA2 prints has to be the number
//Tonto printed, to the two decimals Tonto prints.  Nothing else in the tree checks this analysis
//against anything but itself.
//
//  Tonto            n_A 9.45  n_B 1.49  n_AB 9.72  s_AB 1.22  Cov 0.94  Ion 0.21  Tot 0.96
//                   % Pythagorean 95.25  % Araki 86.01   populations O 9.45, H 1.49, total 9.99
TEST(RgbiRobustnessTests, TontoWaterRobyGouldNumbersAreReproduced)
{
	const auto p = tonto_water_archive();
	if (p.empty())
		GTEST_SKIP() << "tests/RGBI/h2o.MOs,r not found";
	std::string out;
	{
		CoutCapture cap;
		WFN wavy(p);
		//symmetrize = false (Tonto: "Use spherical averaging? F"), use_ano_basis = false ("Use NAOs? T")
		Roby_information roby(wavy, {}, false, false, false, false);
		out = cap.str();
	}
	EXPECT_NEAR(value_after(out, "Population of atom 0: "), 9.45, 5e-3);
	EXPECT_NEAR(value_after(out, "Population of atom 1: "), 1.49, 5e-3);
	EXPECT_NEAR(value_after(out, "Population of atom 2: "), 1.49, 5e-3);
	EXPECT_NEAR(value_after(out, "Number of electrons in Roby Analysis:  "), 9.99, 5e-3);

	//the two bond rows, which are one symmetry orbit and so have to agree with each other as well
	const double tonto[] = { 9.45, 1.49, 9.72, 1.22, 0.94, 0.21, 0.96, 95.25, 86.01 };
	const char *names[] = { "n_A", "n_B", "n_AB", "s_AB", "Cov", "Ion", "Tot", "%Pythagorean", "%Araki" };
	int rows = 0;
	std::istringstream in(out);
	std::string line;
	while (std::getline(in, line)) {
		if (line.find("O -  H") == std::string::npos)
			continue;
		vec numbers;
		std::istringstream cells(line.substr(line.find("O -  H") + 6));
		double v = 0.0;
		while (cells >> v)
			numbers.push_back(v);
		ASSERT_EQ(numbers.size(), 9u) << line;
		for (int i = 0; i < 9; i++)
			EXPECT_NEAR(numbers[i], tonto[i], 5e-3) << names[i] << " in: " << line;
		rows++;
	}
	EXPECT_EQ(rows, 2); //O-H1 and O-H2
}

//The free-atom SCFs run through occ, and occ parallelises with TBB. Pinning them to one thread for
//reproducibility is only half a fix: the guard that does it has to put the process back the way it
//found it, or the first RGBI analysis leaves every later occ user in the same binary serial. Before
//anything installs a tbb::global_control, occ::parallel::get_num_threads() already answers 1, so a
//guard that naively restores "what get_num_threads() said" installs a control at 1 and never lifts it.
TEST(RgbiRobustnessTests, TheAnalysisLeavesOccsThreadCountAsItFoundIt)
{
	const auto p = water_he_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/TFVC/water.gbw not found";

	//case 1: a caller who had asked for a thread count gets that count back
	occ::parallel::set_num_threads(3);
	(void)rgbi_ano_output(false);
	EXPECT_EQ(occ::parallel::get_num_threads(), 3);

	//case 2: a caller who never set one is left without a control at all, not pinned to 1
	occ::parallel::shutdown_tbb();
	occ::parallel::nthreads = 4; //bookkeeping says 4, no control installed - occ's own starting state
	(void)rgbi_ano_output(false);
	EXPECT_EQ(occ::parallel::get_tbb_control(), nullptr);
	EXPECT_EQ(occ::parallel::get_num_threads(), 4);

	occ::parallel::shutdown_tbb();
	occ::parallel::nthreads = 1;
}

//The reproducibility check that has teeth. tests/ECP_SF/Au2Br2.gbw is centrosymmetric: its two Au
//atoms are one orbit, as are the two Au-Br bonds and the two Au-P bonds, so the molecule's own
//symmetry is a reference that needs no second program. Those numbers used to disagree inside a single
//run - populations 16.040429 against 16.029264, bond totals apart by 0.019 and covalent percentages
//by 1.4 - because each free-atom Fock matrix reduced in whatever order TBB's work stealing produced
//and a free atom's open shell is degenerate enough for the last bit to pick a different member of the
//manifold. Symmetry-equivalent centres now agree to every digit, which is the only way a published
//bond index can be compared to anything.
TEST(RgbiRobustnessTests, SymmetryEquivalentGoldCentresAgreeToEveryDigit_full)
{
	if (const char *env = std::getenv("RUN_FULL_TEST"); !env || std::string(env) == "0" || std::string(env) == "false")
		GTEST_SKIP() << "Set RUN_FULL_TEST=1 for the 53-atom Au2Br2 Roby-Gould analysis";
	const auto p = nos_test_repo_root() / "tests" / "ECP_SF" / "Au2Br2.gbw";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/ECP_SF/Au2Br2.gbw not found";
	std::string out;
	{
		CoutCapture cap;
		WFN wavy(p);
		Roby_information roby(wavy, {}, true, true, false, false);
		out = cap.str();
	}

	//the two gold atoms
	const double au0 = value_after(out, "Population of atom 0: ");
	const double au1 = value_after(out, "Population of atom 1: ");
	ASSERT_TRUE(std::isfinite(au0) && std::isfinite(au1)) << "no populations in the output";
	EXPECT_DOUBLE_EQ(au0, au1);

	//and the two bond orbits. The inversion centre maps 0<->1 (Au), 2<->3 (Br), 4<->5 (P), so the
	//rows "0 - 2" and "1 - 3" are one bond, as are "0 - 4" and "1 - 5".
	auto bond_row = [&out](const std::string &pair) {
		vec numbers;
		std::istringstream in(out);
		std::string line;
		while (std::getline(in, line)) {
			//the pair has to be the first thing on the line, or "0 - 2" also matches "10 - 2"
			const size_t start = line.find_first_not_of(" \t");
			if (start == std::string::npos || line.compare(start, pair.size(), pair) != 0)
				continue;
			const size_t at = start;
			std::istringstream cells(line.substr(line.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ", at)));
			std::string element_a, dash, element_b;
			cells >> element_a >> dash >> element_b;
			double v = 0.0;
			while (cells >> v)
				numbers.push_back(v);
			break;
		}
		return numbers;
	};
	const char *orbits[][2] = { { "0 - 2", "1 - 3" }, { "0 - 4", "1 - 5" } };
	for (const auto &orbit : orbits) {
		const vec first = bond_row(orbit[0]), second = bond_row(orbit[1]);
		ASSERT_EQ(first.size(), 9u) << "no bond row " << orbit[0];
		ASSERT_EQ(second.size(), 9u) << "no bond row " << orbit[1];
		for (size_t i = 0; i < 9; i++)
			EXPECT_DOUBLE_EQ(first[i], second[i]) << "column " << i << " of " << orbit[0] << " vs " << orbit[1];
	}
}
