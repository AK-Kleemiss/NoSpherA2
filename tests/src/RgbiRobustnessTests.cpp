#include "pch.h"
#include "core/wfn_class.h"
#include "core/bondwise_analysis.h"

#include <occ/core/parallel.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
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

	//The nine numbers of one bond row of the RGBI table, or an empty vector if there is no such row.
	//The indices are right-aligned in their own fields ("   0 -   2   Au - Br  ..."), so they are read as
	//numbers: matching the literal "0 - 2" finds no row at all - a copy of this reader did exactly that
	//behind RUN_FULL_TEST and asserted on its own parser before comparing a single number - and "0 -"
	//would also match "10 -".
	vec bond_row(const std::string &text, const int a, const int b)
	{
		vec numbers;
		std::istringstream in(text);
		std::string line;
		while (std::getline(in, line)) {
			std::istringstream cells(line);
			int i = 0, j = 0;
			char dash = 0;
			if (!(cells >> i >> dash >> j) || dash != '-' || i != a || j != b)
				continue;
			std::string element_a, element_dash, element_b;
			if (!(cells >> element_a >> element_dash >> element_b) || element_dash != "-")
				continue;
			double v = 0.0;
			while (cells >> v)
				numbers.push_back(v);
			break;
		}
		return numbers;
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
	//Ce's free-atom SCF is the most expensive thing in this suite by four orders of magnitude: 2414355 ms
	//of a 96-core node in the 25 Sep 04:38 run, against 165 ms for the next slowest RGBI test. A gate
	//nobody waits for is a gate nobody runs, and this is the only test covering the f-element free-atom
	//path, so it stays - behind a switch, and named in the skip message rather than quietly dropped.
	if (!std::getenv("NOS_RGBI_SLOW_TESTS"))
		GTEST_SKIP() << "Ce's free-atom SCF took 2414 s on AKL007 (96 cores, 25 Sep 04:38); "
		                "set NOS_RGBI_SLOW_TESTS=1 to run it";
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
	//And the part this test used to pass over in silence: in the suite's own run, pinned to one thread,
	//Ce's free-atom SCF did not converge - occ's full 100 iterations with |dE|/E at 9.9e-10 and
	//max|FDS-SDF| stalled at 7.7e-5, logged at error level, the last density returned rather than
	//thrown, and used and cached as cerium's free-atom reference. The population above was still a
	//population and this test still passed, which is why the run now warns.
	//
	//Not asserted here, deliberately, and measured rather than assumed: whether it converges depends on
	//the threading. Unpinned with 8 threads (NOS_RGBI_NO_PIN=1, AKL007, 96 cores) the same binary and
	//fixture converged in 37.5 s and printed no warning; the pinned run took 2414 s and never did. An
	//expectation either way would be asserting the node's thread count. The warning itself is checked by
	//AFreeAtomThatRunsOutOfIterationsSaysSo below, on water, in about a second.
}

//occ does not throw when an SCF runs out of iterations - scf_impl.h logs one line at error level and
//returns the last energy - so RGBI used, and cached, free-atom densities that never converged, and the
//only fixture that reached that state costs 2414 s and answers differently depending on the thread pin.
//NOS_RGBI_FREE_ATOM_MAXITER caps the iterations so the reporting can be checked on water: both
//directions, because a warning that is always printed is not a warning.
TEST(RgbiRobustnessTests, AFreeAtomThatRunsOutOfIterationsSaysSo)
{
	if (water_he_fixture().empty())
		GTEST_SKIP() << "tests/TFVC/water.gbw not found";
	const auto count = [](const std::string &hay, const std::string &token) {
		size_t n = 0;
		for (size_t at = hay.find(token); at != std::string::npos; at = hay.find(token, at + 1))
			n++;
		return n;
		};
	const std::string token = "WARNING: the free-atom SCF of ";

	//the control arm, cold so it really runs the SCFs: O, H and He all converge and nothing is said
	clear_rgbi_free_atom_cache();
	const std::string quiet = rgbi_ano_output(false);
	ASSERT_FALSE(quiet.empty()) << "the analysis produced no output at all";
	EXPECT_EQ(count(quiet, token), 0u) << "water's free atoms converge, so there is nothing to warn about:\n"
		<< quiet;

	clear_rgbi_free_atom_cache();
#ifdef _WIN32
	_putenv_s("NOS_RGBI_FREE_ATOM_MAXITER", "1");
#else
	setenv("NOS_RGBI_FREE_ATOM_MAXITER", "1", 1);
#endif
	const std::string capped = rgbi_ano_output(false);
#ifdef _WIN32
	_putenv_s("NOS_RGBI_FREE_ATOM_MAXITER", "");
#else
	unsetenv("NOS_RGBI_FREE_ATOM_MAXITER");
#endif
	//and never leave a one-iteration density in the cache for whatever test runs next
	clear_rgbi_free_atom_cache();
	ASSERT_FALSE(capped.empty()) << "the capped arm produced no output at all";

	//three distinct free atoms, three warnings, each naming its element by Z and carrying both residuals
	EXPECT_EQ(count(capped, token), 3u) << "water plus helium has 3 distinct free atoms:\n" << capped;
	EXPECT_NE(capped.find("(Z=8,"), std::string::npos) << capped;
	EXPECT_NE(capped.find("(Z=1,"), std::string::npos) << capped;
	EXPECT_NE(capped.find("(Z=2,"), std::string::npos) << capped;
	EXPECT_EQ(count(capped, "did not converge in 1 iterations"), 3u) << capped;
	EXPECT_EQ(count(capped, "max|FDS-SDF|="), 3u) << capped;
	//and it reports rather than refuses: the analysis still finishes and still prints its populations
	EXPECT_TRUE(std::isfinite(value_after(capped, "Population of atom 0: "))) << capped;
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
//
//The test was gated behind RUN_FULL_TEST when it was written and in that state it never ran once - and
//it could not have passed: bond_row matched the literal string "0 - 2" against a line printed as
//"   0 -   2   Au - Br  ...", where each index is right-aligned in its own field, so it found no row and
//asserted on its own parser before it compared a single number. The indices are read as numbers now, and
//the gate is gone: fe628ab9's free-atom cache brought the whole analysis to 4 s at the command line and
//5.6 s in the test binary, which is not a gate's worth of time. Ungating it also puts the ECP heavy-atom
//RGBI path into the default suite, which is where 09a9c925 belongs - both deployed share binaries abort
//inside "Calculating ANOs for all atoms..." on this exact file, 4 runs out of 4, with malloc() and
//SIGSEGV heap diagnostics, while this branch finishes it twice with bit-identical output.
TEST(RgbiRobustnessTests, SymmetryEquivalentGoldCentresAgreeToEveryDigit)
{
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

	//the three orbits of the inversion centre: 0<->1 (Au), 2<->3 (Br), 4<->5 (P)
	for (const int a : { 0, 2, 4 }) {
		const double first = value_after(out, "Population of atom " + std::to_string(a) + ": ");
		const double second = value_after(out, "Population of atom " + std::to_string(a + 1) + ": ");
		ASSERT_TRUE(std::isfinite(first) && std::isfinite(second)) << "no population for atoms " << a << " and " << a + 1;
		EXPECT_DOUBLE_EQ(first, second) << "populations of the symmetry-equivalent atoms " << a << " and " << a + 1;
	}

	//The one quantity on this path that does NOT come out bit-identical between two centres of the same
	//orbit: the population left outside the ANO cutoff, 8.8788628 against 8.8788629 for the two gold
	//atoms and 3.16409 against 3.1640901 for the two phosphorus atoms, the same on both sides at 1 and
	//8 threads and in two repeats of each (AKL007, 25 Sep). That is a reduction-order residual in the
	//last printed digit of a number that is itself a difference of two large ones, so it is pinned at
	//1e-6 instead of being asserted equal - and pinned rather than ignored, because the populations and
	//all nine bond columns below ARE bit-identical, and a drift here would be the first sign of the
	//old defect coming back.
	for (const int a : { 0, 2, 4 }) {
		const std::string key = "Atomic projector population outside selected ANO cutoff of atom ";
		const double first = value_after(out, key + std::to_string(a) + ": ");
		const double second = value_after(out, key + std::to_string(a + 1) + ": ");
		ASSERT_TRUE(std::isfinite(first) && std::isfinite(second)) << "no outside-cutoff population for atoms " << a << " and " << a + 1;
		EXPECT_NEAR(first, second, 1e-6) << "outside-cutoff population of atoms " << a << " and " << a + 1;
	}

	//and the two bond orbits: the rows "0 - 2" and "1 - 3" are one bond, as are "0 - 4" and "1 - 5".
	const int orbits[][4] = { { 0, 2, 1, 3 }, { 0, 4, 1, 5 } };
	for (const auto &orbit : orbits) {
		const vec first = bond_row(out, orbit[0], orbit[1]), second = bond_row(out, orbit[2], orbit[3]);
		ASSERT_EQ(first.size(), 9u) << "no bond row " << orbit[0] << " - " << orbit[1];
		ASSERT_EQ(second.size(), 9u) << "no bond row " << orbit[2] << " - " << orbit[3];
		for (size_t i = 0; i < 9; i++)
			EXPECT_DOUBLE_EQ(first[i], second[i]) << "column " << i << " of " << orbit[0] << " - " << orbit[1]
				<< " vs " << orbit[2] << " - " << orbit[3];
	}
}

//A free atom's density depends on its basis and its electron count, not on which atom it is: two
//hydrogens of the same molecule have the same free-atom SCF to the last bit. The ANO route computed
//one per centre anyway, so tests/Fe_gbw/Fe.gbw ran 21 one-atom SCFs for 4 distinct answers - and with
//the reductions pinned to one thread for reproducibility, that repetition was most of the run: the
//analysis did not finish inside 900 s where the unpinned build took 405 s.
//
//The invariant is stated so that it holds whether or not an earlier test in the same binary has
//already warmed the cache: every centre must be accounted for, and the number of SCFs actually run
//must be smaller than the number of centres. tests/TFVC/water.gbw is O, two H and a non-bonded He -
//4 centres, 3 distinct free atoms - so a cold run is 3 SCFs and 1 hit, a warm one 0 and 4, and 4
//SCFs is the defect.
//
//Made red on purpose by returning before the cache lookup: the count goes to 4 SCFs and 0 hits.
TEST(RgbiRobustnessTests, IdenticalCentresShareOneFreeAtomScf)
{
	const auto p = water_he_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/TFVC/water.gbw not found";

	struct ScopedEnv {
		ScopedEnv() {
#ifdef _WIN32
			_putenv_s("NOS_RGBI_DEBUG", "1");
#else
			setenv("NOS_RGBI_DEBUG", "1", 1);
#endif
		}
		~ScopedEnv() {
#ifdef _WIN32
			_putenv_s("NOS_RGBI_DEBUG", "");
#else
			unsetenv("NOS_RGBI_DEBUG");
#endif
		}
	} debug_on;

	//Cold, or this counts the test order rather than the run: with the cache left warm by a neighbour this
	//arm reports 0 SCFs and 4 hits, which satisfies every bound below while proving nothing.
	clear_rgbi_free_atom_cache();
	const std::string out = rgbi_ano_output(false);
	ASSERT_FALSE(out.empty()) << "the analysis produced no output at all";

	const auto count = [&out](const std::string &token) {
		size_t n = 0;
		for (size_t at = out.find(token); at != std::string::npos; at = out.find(token, at + 1))
			n++;
		return n;
		};
	//FREEATOM-CACHED is a prefix of nothing, but FREEATOM-START and FREEATOM both match a
	//FREEATOM-START line, so count the completed-SCF lines by their own distinct opening.
	const size_t scfs = count("FREEATOM-START");
	const size_t hits = count("FREEATOM-CACHED");

	EXPECT_EQ(scfs + hits, 4u) << "water.gbw has 4 centres and each must be answered once:\n" << out;
	//Exact, now that the cache is cleared above: O, H and He are three distinct free atoms and the second
	//hydrogen is the one repeat. A bound like "fewer than 4" was satisfied by a warm cache doing none.
	EXPECT_EQ(scfs, 3u) << "water plus helium has 3 distinct free atoms, so 3 SCFs and no more:\n" << out;
	EXPECT_EQ(hits, 1u) << "the second hydrogen is the only repeat and it must come from the cache:\n" << out;

	//The half that makes this a check and not a stopwatch. NOS_RGBI_NO_FREEATOM_CACHE recomputes every
	//centre in this same process, so the cached bond table has something to be identical to. Without it
	//the cache offers no evidence at all: it is barely faster. On the one fixture it exists for,
	//tests/Fe_gbw/Fe.gbw, all eight arms of job 582380 - four configurations, two repeats - print md5
	//1672c4eae0c8 and 20 rows, and the 17 repeats the cache removes are worth 0.93 s unpinned and 3.37 s
	//pinned, the non-Fe SCFs summing to 0.356 s against 1.290 s. Do not quote the whole-run seconds for
	//this: over those two repeats the cached arm was 6.5 s faster and then 3.0 s slower. So the eighth
	//fixture is verified and the cache's justification is entirely this identity.
#ifdef _WIN32
	_putenv_s("NOS_RGBI_NO_FREEATOM_CACHE", "1");
#else
	setenv("NOS_RGBI_NO_FREEATOM_CACHE", "1", 1);
#endif
	const std::string uncached = rgbi_ano_output(false);
#ifdef _WIN32
	_putenv_s("NOS_RGBI_NO_FREEATOM_CACHE", "");
#else
	unsetenv("NOS_RGBI_NO_FREEATOM_CACHE");
#endif
	ASSERT_FALSE(uncached.empty()) << "the uncached arm produced no output at all";

	const auto count_in = [](const std::string &hay, const std::string &token) {
		size_t n = 0;
		for (size_t at = hay.find(token); at != std::string::npos; at = hay.find(token, at + 1))
			n++;
		return n;
		};
	EXPECT_EQ(count_in(uncached, "FREEATOM-START"), 4u)
		<< "the off switch is not off: fewer than 4 SCFs with the cache disabled";
	EXPECT_EQ(count_in(uncached, "FREEATOM-CACHED"), 0u)
		<< "a disabled cache still reported a hit";

	//Compare the tables, not the whole stream: the debug lines differ by construction.
	const auto table = [](const std::string &s) {
		const size_t from = s.find("Atom Nr");
		return from == std::string::npos ? std::string() : s.substr(from);
		};
	ASSERT_FALSE(table(out).empty()) << "no bond table in the cached run";
	EXPECT_EQ(table(out), table(uncached))
		<< "the cache changed the answer, which is the only way it can be wrong\ncached:\n"
		<< table(out) << "\nuncached:\n" << table(uncached);
}

//The cache removes the repeats; it does not make the remaining SCFs any faster, and on
//tests/Fe_gbw/Fe.gbw the repeats it removes are worth 0.93 s of a 404.0 s run, reproduced to 0.6 % over two
//repeats - the earlier claim that it cost about 1300 s was the pin measured across two jobs, not this cache.
//The distinct SCFs run one behind the other and do not depend on each other, so NOS_RGBI_PARALLEL_FREEATOM
//runs them up front and concurrently. On this fixture that can reach only the 0.93 s, because 402.653 s of
//the run is one Fe; the fixture where it could pay is one with many distinct heavy centres. The risk it buys is
//the only one worth testing for: occ's SCF is not documented re-entrant, and a free-atom density that
//comes out subtly different under concurrency would be invisible in a timing table.
//
//So the assertion is not that it is faster - on 3 light atoms it cannot be - but that the answer did
//not move. Two things about how it is written are the whole difficulty, and both were found by making
//the check red on purpose rather than by reading it:
//
//The serial arm runs with the cache OFF. Written the obvious way round - serial first with the cache on,
//then the parallel arm - the serial arm fills the cache, the warm pass finds every centre already
//answered, no SCF ever runs concurrently, and the check compares a concurrent run against itself. It
//passed a deliberately perturbed free-atom density that way, because the perturbation lived on the SCF
//path and no SCF ran.
//
//And the comparison is over the 14-digit debug lines, not only the printed bond table. A bond table
//prints three decimals, so perturbing one element of the atomic density by a part in 1e7 changes no
//printed digit and the table - and any md5 taken over it - is identical. That sensitivity floor applies
//to every md5-identical claim in this harness: they say "agrees to the digits it prints", which is not
//"agrees bit for bit". sumD, normD and the SCF energy are printed at setprecision(14) under
//NOS_RGBI_DEBUG, so the check compares those and can see what the table cannot.
TEST(RgbiRobustnessTests, WarmingTheFreeAtomCacheInParallelDoesNotChangeTheAnswer)
{
	const auto p = water_he_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/TFVC/water.gbw not found";

	const auto set_env = [](const char *name, const char *value) {
#ifdef _WIN32
		_putenv_s(name, value);
#else
		if (*value)
			setenv(name, value, 1);
		else
			unsetenv(name);
#endif
		};
	const auto table = [](const std::string &s) {
		const size_t from = s.find("Atom Nr");
		return from == std::string::npos ? std::string() : s.substr(from);
		};
	//Every completed free-atom SCF, by its numbers and not by which atom it sat on: the label is dropped
	//so that the two H of water collapse to one entry, which is the whole claim of the cache, and so that
	//the warm pass's representative compares equal to whichever H the serial loop reached first.
	const auto scf_numbers = [](const std::string &s) {
		std::set<std::string> found;
		for (size_t at = s.find("FREEATOM "); at != std::string::npos; at = s.find("FREEATOM ", at + 1)) {
			const size_t eol = s.find('\n', at);
			const std::string line = s.substr(at, eol == std::string::npos ? eol : eol - at);
			if (line.find(" E=") == std::string::npos)
				continue;                         //a -START or -WARM line, not a completed SCF
			const size_t z = line.find(" Z=");    //drop "FREEATOM <label>"
			if (z != std::string::npos)
				found.insert(line.substr(z));
		}
		return found;
		};

	//The cache off, so these are four SCFs actually run, in sequence, and nothing is left behind for the
	//parallel arm to find.
	set_env("NOS_RGBI_DEBUG", "1");
	set_env("NOS_RGBI_NO_FREEATOM_CACHE", "1");
	const std::string serial = rgbi_ano_output(false);
	set_env("NOS_RGBI_NO_FREEATOM_CACHE", "");
	ASSERT_FALSE(table(serial).empty()) << "no bond table from the serial run";

	//Cold on purpose. The cache lives for the process, so without this the warm pass finds every atom of
	//this fixture already answered by an earlier test in this same binary, runs no SCF at all, and the
	//fourteen-digit comparison below has nothing to compare. That is how this check failed once it
	//acquired neighbours, having passed while it ran alone: 0 SCFs in the parallel arm, and the product
	//behaving correctly the whole time.
	clear_rgbi_free_atom_cache();
	set_env("NOS_RGBI_PARALLEL_FREEATOM", "1");
	const std::string parallel = rgbi_ano_output(false);
	set_env("NOS_RGBI_PARALLEL_FREEATOM", "");
	set_env("NOS_RGBI_DEBUG", "");
	ASSERT_FALSE(table(parallel).empty()) << "no bond table from the parallel run";

	//Proof the flag did something, so that an equal answer is not equal because nothing ran differently.
	EXPECT_NE(parallel.find("FREEATOM-WARM 3 distinct of 4 centres"), std::string::npos)
		<< "the parallel warm pass did not run, so this check compared two identical code paths";
	EXPECT_EQ(serial.find("FREEATOM-WARM"), std::string::npos)
		<< "the serial run warmed the cache in parallel without being asked to";

	const auto serial_scfs = scf_numbers(serial);
	const auto parallel_scfs = scf_numbers(parallel);
	ASSERT_EQ(serial_scfs.size(), 3u)
		<< "4 centres of 3 distinct free atoms did not produce 3 distinct SCF results in the serial arm";
	ASSERT_EQ(parallel_scfs.size(), 3u)
		<< "the parallel arm did not run 3 SCFs, so concurrency was never exercised here";
	EXPECT_EQ(serial_scfs, parallel_scfs)
		<< "a free-atom SCF gave a different answer when run concurrently - the density, its norm or its "
		   "energy differs in the 14 printed digits";

	EXPECT_EQ(table(serial), table(parallel))
		<< "running the free-atom SCFs concurrently changed the answer\nserial:\n"
		<< table(serial) << "\nparallel:\n" << table(parallel);
}

//The determinism guard's own effect, asserted directly instead of through the digits it protects.
//
//TheAnalysisLeavesOccsThreadCountAsItFoundIt above checks that the guard cleans up. Nothing checked
//that it ever engaged, and it stopped: occ declares `inline int nthreads = 1` and installs no
//tbb::global_control until set_num_threads is called, so a guard that skipped the call when
//get_num_threads() already read 1 installed nothing, and every free-atom SCF ran its reductions across
//every core. The damage was 1.3e-13 in norm(D) - two chemically identical hydrogens of this same fixture
//disagreed in 7 of 50 processes - which no bond table printed at three decimals can see, and none of the
//eight md5 fixtures of this harness did.
//
//So the observable is the pin, not the digits: a check on the digits would have gone red about one time
//in seven, and a gate that fails one run in seven gets muted. This one is exact. NOS_RGBI_NO_PIN turns
//the pin off in this same binary and makes it red immediately, which is how it was confirmed to be able
//to fail at all.
TEST(RgbiRobustnessTests, EveryFreeAtomScfRunsWithOccPinnedForReal)
{
	const auto p = water_he_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/TFVC/water.gbw not found";

	struct ScopedEnv {
		ScopedEnv() {
#ifdef _WIN32
			_putenv_s("NOS_RGBI_DEBUG", "1");
			_putenv_s("NOS_RGBI_NO_FREEATOM_CACHE", "1");
#else
			setenv("NOS_RGBI_DEBUG", "1", 1);
			setenv("NOS_RGBI_NO_FREEATOM_CACHE", "1", 1);
#endif
		}
		~ScopedEnv() {
#ifdef _WIN32
			_putenv_s("NOS_RGBI_DEBUG", "");
			_putenv_s("NOS_RGBI_NO_FREEATOM_CACHE", "");
#else
			unsetenv("NOS_RGBI_DEBUG");
			unsetenv("NOS_RGBI_NO_FREEATOM_CACHE");
#endif
		}
	} debug_uncached;

	//Uncached on purpose: a cache hit runs no SCF, so the cached run would only report on the
	//centres that happened to miss. Every centre has to be seen pinned, including the repeats.
	const std::string out = rgbi_ano_output(false);
	ASSERT_FALSE(out.empty()) << "the analysis produced no output at all";

	size_t starts = 0, pinned = 0;
	for (size_t at = out.find("FREEATOM-START"); at != std::string::npos;
		at = out.find("FREEATOM-START", at + 1)) {
		const size_t eol = out.find('\n', at);
		const std::string line = out.substr(at, eol == std::string::npos ? eol : eol - at);
		starts++;
		if (line.find("pinned=1") != std::string::npos)
			pinned++;
	}
	ASSERT_EQ(starts, 4u)
		<< "4 centres with the cache disabled must run 4 free-atom SCFs, so this test is looking at "
		"nothing if it does not see 4 of them:\n" << out;
	EXPECT_EQ(pinned, starts)
		<< pinned << " of " << starts << " free-atom SCFs began with a tbb::global_control admitting "
		"one thread. The rest ran their reductions concurrently, and the free-atom density they return "
		"is then reproducible only to the digits a bond table prints:\n" << out;

	//And the guard still has to put occ back: an assertion that it engaged is worth nothing if the
	//way it engaged leaves the rest of the binary pinned to one thread.
	EXPECT_EQ(occ::parallel::get_tbb_control(), nullptr)
		<< "the analysis left a TBB control behind, so every later occ user in this process is serial";
}

//RGBI refuses a basis whose shells go beyond h, because its O_h symmetrization has no transform for
//them. It used to refuse from inside symmetrize_atomic_matrix_oh(), which is reached only after the
//overlap matrix is built and, on the ANO route, after a free-atom SCF per element: -rgbi on
//tests/CuF2_i_func/71/calc.gbw ran for 467.3 s before printing it. The decision needs the shell
//types and nothing else, so it is now taken from the basis before any of that work. err_checkf()
//exits the process, so what is tested here is the question the guard asks - highest_shell_angular_
//momentum() over the whole molecule - rather than the exit; a test that reproduces the refusal
//end to end would have to carry a 670-MO fixture and wait for it.
namespace {
	//type is the NoSpherA2 basis_set_entry convention, l + 1: an s shell is 1, an i shell is 7.
	atom shell_atom(const std::string &label, const int Z, const double z, const ivec &l_values)
	{
		atom a(label, {}, 1, 0.0, 0.0, z, Z);
		for (int s = 0; s < static_cast<int>(l_values.size()); s++)
			a.push_back_basis_set(1.0 + 0.1 * s, 1.0, l_values[s] + 1, s);
		return a;
	}

	WFN wfn_of(const std::vector<atom> &atoms)
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		for (const atom &a : atoms) w.push_back_atom(a);
		return w;
	}
}

TEST(RgbiRobustnessTests, TheUnsupportedShellQuestionIsAskedOfEveryAtomsBasis)
{
	//s through h is what the symmetrization supports, and the whole point of asking early is that the
	//answer must not depend on how far the analysis got.
	EXPECT_EQ(highest_shell_angular_momentum(wfn_of({ shell_atom("H", 1, 0.0, {0}) })), 0);
	EXPECT_EQ(highest_shell_angular_momentum(wfn_of({ shell_atom("Cu", 29, 0.0, {0, 1, 2, 3, 4, 5}) })), 5);

	//The defect shape this guards against: a check that inspects the first atom only. Copper carries
	//the i shells in CuF2_i_func and fluorine does not, so an i shell on any centre but the first has
	//to be seen.
	const WFN i_on_the_second = wfn_of({ shell_atom("F", 9, 0.0, {0, 1, 2}),
										 shell_atom("Cu", 29, 3.5, {0, 1, 2, 3, 4, 5, 6}) });
	EXPECT_EQ(highest_shell_angular_momentum(i_on_the_second), 6)
		<< "an i shell on the second atom was not seen, so RGBI would run the whole overlap and only "
		"then refuse - which is the defect this replaced";

	//A basis with no shells at all has no highest l. The Roby_information constructor refuses that
	//case separately (a plain .wfn), and this must not turn into a 0 that looks supported.
	EXPECT_EQ(highest_shell_angular_momentum(wfn_of({ atom("H", {}, 1, 0.0, 0.0, 0.0, 1) })), -1);
}

//A refusal that sends the reader back to the file it just refused is worse than no advice at all, and
//that is what this one did: -rgbi on tests/molden_file/f_ref.wfx printed "f_ref.wfx carries none ... run
//RGBI on a .wfx, .fchk, .molden, .gbw or a Tonto archive instead". Measured on this binary, 6 reader
//cells: .gbw and .molden complete, .wfx and .wfn are refused for an empty basis, .fchk is refused for a
//missing contracted density matrix, and the Tonto fixture in tests/cytidine_tonto cannot be read at all
//because its companion stdout file is not in the tree. So two of the four formats the message named
//cannot run this analysis and one of them was the input.
TEST(RgbiRobustnessTests, TheRefusalNeverSuggestsTheFormatItIsRefusing)
{
	//the whole list, for the case where the refused file is neither of them
	EXPECT_EQ(rgbi_supported_input_phrase(".wfn"), "a .gbw or a .molden");
	//and each of the two working formats drops itself, because a .gbw whose basis did not survive the
	//reader must not be answered with "use a .gbw"
	EXPECT_EQ(rgbi_supported_input_phrase(".gbw"), "a .molden");
	EXPECT_EQ(rgbi_supported_input_phrase(".molden"), "a .gbw");

	//The property, over every extension the matrix covers, written so that adding a format to the list
	//cannot reintroduce the defect: whatever the phrase says, it does not say the input.
	for (const char *ext : { ".wfn", ".wfx", ".ffn", ".fchk", ".molden", ".gbw", "GBW", "wfx" }) {
		const std::string phrase = rgbi_supported_input_phrase(ext);
		std::string lower = ext;
		std::transform(lower.begin(), lower.end(), lower.begin(),
			[](unsigned char c) { return (char)std::tolower(c); });
		if (!lower.empty() && lower.front() != '.')
			lower.insert(lower.begin(), '.');
		EXPECT_EQ(phrase.find(lower), std::string::npos)
			<< "the refusal of a " << ext << " offers a " << lower << ": " << phrase;
		EXPECT_FALSE(phrase.empty()) << "a refusal with no alternative at all, for " << ext;
	}
}

//An exactly octahedral molecule has ONE bond orbit: all six Te-F bonds are the same bond, so the six
//rows of the table must agree in every printed digit. On this branch before the fix they did not - they
//came out as three pairs, s_AB 0.192 / 0.194 / 0.192 and Cov. 0.459 / 0.449 / 0.455 - and the log said
//why: "the atomic subspace of rank 13 cuts through a degenerate occupation (0.45375369)". Te's ANO
//occupations group as 3, 3, 3, 2, 3, so the boundaries are at 11 and 14 and the fixed rank of 13 kept
//two members of a threefold set. A degenerate set spans one subspace and which vectors inside it the
//diagonalizer hands back is arbitrary, so a rank that cuts one makes the atomic projector itself
//arbitrary; it then no longer commutes with the molecule's symmetry and bonds the symmetry makes
//identical come out different. calculateAtomicNAO now extends the rank to the end of the set it would
//have cut.
//
//WHAT WOULD MAKE THIS TEST WRONG rather than red: TeF6/def2-TZVP is octahedral to the last digit of the
//input geometry (1.815 A along each axis, ORCA NoUseSym), so any spread at all is the code's. The three
//cheaper members of the same bisection - SF6/def2-SVP (l <= 2), SF6/def2-TZVP (l <= 3) and
//SF6/def2-QZVP (l = 4) - all came out exactly Oh on this binary without the fix, which is what makes
//the ECP/degenerate-rank path and not the spherical transforms the thing under test here.
//
//The binary without the extension IS the red run: on it this fixture prints s_AB 0.192 / 0.194 / 0.192
//and Cov. 0.459 / 0.449 / 0.455, so the five exact columns below fail 2 + 2 + 2 (measured on AKL007,
//binary b34902b06d12, 25 Sep). With the extension the same binary prints six identical rows in those
//five columns (8c137772530d, same node, same fixture).
TEST(RgbiRobustnessTests, OctahedralTeF6HasOneBondOrbitNotThree)
{
	const auto p = nos_test_repo_root() / "tests" / "RGBI_groups" / "tef6_tzvp.gbw";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/RGBI_groups/tef6_tzvp.gbw not found";
	std::string out;
	{
		CoutCapture cap;
		WFN wavy(p);
		Roby_information roby(wavy, {}, true, true, false, false);
		out = cap.str();
	}

	//The six fluorines are one orbit, so one population, and they now print one number: 9.8018244 six
	//times. Before the atomic reference was averaged over all rotations instead of over O_h they read
	//9.8018244 for the four equatorial ones and 9.8018242 for the two along z - the D4h pattern, at the
	//last printed digit. The pin stays a tolerance of 1e-7 rather than becoming EXPECT_DOUBLE_EQ only
	//because these are PARSED PRINTED digits: two numbers that differ below the seventh decimal can still
	//straddle a rounding boundary and print apart, which would be a flake and not a defect. 1e-7 is that
	//printing width, and the residual it allows is one order of magnitude below the one that was there.
	const double f0 = value_after(out, "Population of atom 1: ");
	ASSERT_TRUE(std::isfinite(f0)) << "no population for atom 1";
	for (int a = 2; a <= 6; a++) {
		const double f = value_after(out, "Population of atom " + std::to_string(a) + ": ");
		ASSERT_TRUE(std::isfinite(f)) << "no population for atom " << a;
		EXPECT_NEAR(f0, f, 1e-7) << "population of fluorine " << a << " against fluorine 1";
	}

	//ALL NINE columns now agree in every printed digit, and the assertion below is the equality it used to
	//be unable to make. Two earlier states of this same line are worth keeping, because each was a real
	//measurement: all nine split before the degenerate-rank extension, and five of nine agreed after it
	//while the four ionic columns kept a 1e-3 residual with the D4h pattern (Ion. -0.432 on the four
	//equatorial bonds against -0.431 on the two along z). That residual was the atomic reference being
	//averaged over the 48 operations of O_h instead of over all rotations: an O_h average of a shell with
	//l >= 2 still leaves more than one invariant - e_g and t_2g stay separate - so what survives is the
	//part of the atom's own anisotropy that happens to line up with the Cartesian axes of the input file,
	//and the two fluorines on z do not lie in the frame the way the four on x and y do. Averaging over
	//SO(3) instead leaves exactly one number per pair of shells of equal l, by Schur's lemma, and the
	//pattern is gone: measured max |O_h average - exact average| in a fluorine block is 6.519e-05 on the
	//four equatorial and 2.833e-04 on the two axial ones, which is where the 1e-3 in Ion. came from.
	//
	//WHERE IT IS NOT: two things this comment used to claim, both since measured and dropped.
	//  * Not the pair metric's rank. Every one of the six bonds keeps 81 of 81 eigenvalues on both
	//    metrics, the smallest kept is 4.29e-03 - two and a half decades above the floor - and the table
	//    is byte-identical with NOS_RGBI_PINV_CUTOFF swept from 1E-4 to 1E-8.
	//  * Not the theta split either. The two pi subspaces of ONE bond come out at 85.109 and 85.113 deg
	//    where C4v site symmetry makes them exactly degenerate, and this comment called that "the pair
	//    path's own asymmetry". SF6/def2-QZVP - the same geometry with no ECP - splits the same way, at
	//    80.138 against 80.144, while agreeing in all nine columns of all six bonds. A split that is
	//    present where the result is exact is not the cause of a result that is not, and 1/cos(85 deg) is
	//    11.5, so 6e-3 deg is ~1e-5 in the underlying ratio.
	//Where it was, and the state of the four-corner option matrix now. {sym, no_sym} x {ANO, NAO} on both
	//molecules, measured on AKL007 25 Sep, printed spread over the six bonds per column:
	//                   before (O_h average)                     after (exact rotational average)
	//  no_sym + NAO     0 in all nine, both molecules            0 in all nine, unchanged
	//  sym    + NAO     4 + 2, SF6 Pyth. 33.926 / 34.260         0 in all nine, both molecules
	//  sym    + ANO     SF6 0, TeF6 1e-3 in Ion. (this test)     0 in all nine, both molecules
	//  no_sym + ANO     2 + 2 + 2, SF6 s_AB 0.659/0.632/0.661    0 in all nine, both molecules
	//All four corners are exact now, and the last row took a second change and an argument that had been
	//got wrong here. This comment used to say that row was "expected to" split because -rgbi_no_sym asks
	//for no average and "nothing is entitled to fix it". That reads the flag as applying to something it
	//does not: on the ANO route the matrix being averaged is not the molecule's, it is a free atom's own
	//density from occ's atomic SCF, and a free atom is spherically symmetric. The average there is a
	//property of the reference, not an approximation imposed on the molecule, and skipping it let the
	//SCF's arbitrary choice of m components put the axes of the FILE into every bond. So the exact
	//rotational average now runs on the ANO route whatever the flag says, the flag still governs the
	//molecular route it was named for, and a run that asks for it on the ANO route is told in one printed
	//line that it changed nothing there. OctahedralTeF6IsExactlyOhOnTheAnoPathWithoutSymmetrization is
	//that row's own test. The assertion below is EXPECT_DOUBLE_EQ on all nine columns.
	//
	//WHAT WOULD MAKE THIS FAIL: any change that lets the atom's orientation back into its own reference -
	//an average applied per bond rather than per atom, the O_h route reinstated for a spherical basis, a
	//degenerate ANO set truncated mid-shell. It uses no external program and no reference numbers: it only
	//asserts that a symmetry the molecule has survives into the output.
	const vec first = bond_row(out, 0, 1);
	ASSERT_EQ(first.size(), 9u) << "no bond row 0 - 1";
	for (int b = 2; b <= 6; b++) {
		const vec row = bond_row(out, 0, b);
		ASSERT_EQ(row.size(), 9u) << "no bond row 0 - " << b;
		for (size_t i = 0; i < 9; i++)
			EXPECT_DOUBLE_EQ(first[i], row[i]) << "column " << i << " of Te-F bond 0 - " << b
				<< " against 0 - 1: an octahedral molecule has one bond orbit, and on the DEFAULT option "
				"combination all nine columns reproduce it since the atomic reference is averaged over all "
				"rotations instead of over O_h";
	}
}

//The invariant that LOCALISED the residual the test above used to pin, and the one that keeps it honest
//now that it is gone: with the NAO orbital basis and the symmetrization of the free-atom matrix BOTH OFF,
//an octahedral molecule's six bonds come out identical in every one of the nine columns, to the last bit -
//EXPECT_DOUBLE_EQ, not a tolerance. It is the arm that carried the diagnosis, because it said where the
//defect was NOT: the pair decomposition, the projector construction, the reductions and the printing are
//all exactly Oh-covariant on their own, so what was left could only be the atomic reference. It was: the
//reference was averaged over O_h instead of over all rotations. Keeping this test after the fix is not
//redundant with the one above - it is the only arm that exercises the whole pair machinery with NO atomic
//reference construction in the way at all, so a future break in one of the two is attributable.
//
//WHAT WOULD MAKE THIS FAIL, and it is worth saying because a passing symmetry test is easy to trust too
//much: any change that makes the NAO path's per-bond work depend on the bond's orientation - a reduction
//order tied to atom index, a cutoff applied per bond rather than per orbit, an eigenvector phase leaking
//into a population. It does NOT test the numbers themselves against anything external; it tests that the
//molecule's own symmetry survives, which needs no second program and no reference at all. TeF6 also
//carries an ECP and f functions, so the invariant covers the two input kinds that break things most.
//
//This is also not an obscure corner of the option space: it is the same setting
//TontoWaterRobyGouldNumbersAreReproduced runs, because it is what Tonto itself uses ("Use spherical
//averaging? F", "Use NAOs? T"). The one combination that reproduces the published reference numbers is
//the one that is exactly Oh here.
//
//Made red on purpose by flipping the symmetrization on: 16 assertions fail, the six rows splitting in
//s_AB (0.577 against 0.585) and in all four ionic columns.
TEST(RgbiRobustnessTests, OctahedralTeF6IsExactlyOhWithoutTheAtomicReference)
{
	const auto p = nos_test_repo_root() / "tests" / "RGBI_groups" / "tef6_tzvp.gbw";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/RGBI_groups/tef6_tzvp.gbw not found";
	std::string out;
	{
		CoutCapture cap;
		WFN wavy(p);
		//no group sets, NO symmetrization, NO ANO basis, no eigenvalues, no theta table
		Roby_information roby(wavy, {}, false, false, false, false);
		out = cap.str();
	}

	const double f0 = value_after(out, "Population of atom 1: ");
	ASSERT_TRUE(std::isfinite(f0)) << "no population for atom 1";
	for (int a = 2; a <= 6; a++)
		EXPECT_DOUBLE_EQ(f0, value_after(out, "Population of atom " + std::to_string(a) + ": "))
			<< "population of fluorine " << a << " against fluorine 1 on the NAO path without symmetrization";

	const vec first = bond_row(out, 0, 1);
	ASSERT_EQ(first.size(), 9u) << "no bond row 0 - 1";
	for (int b = 2; b <= 6; b++) {
		const vec row = bond_row(out, 0, b);
		ASSERT_EQ(row.size(), 9u) << "no bond row 0 - " << b;
		for (size_t i = 0; i < 9; i++)
			EXPECT_DOUBLE_EQ(first[i], row[i]) << "column " << i << " of Te-F bond 0 - " << b << " against "
				"0 - 1: this path reproduces the octahedron exactly, so any difference at all is new";
	}
}

//The fourth corner of {sym, no_sym} x {ANO, NAO}, and the last one that was not exact: -rgbi_no_sym with
//the ANO basis. It split TeF6's six bonds 2 + 2 + 2 - worst 2.893 in the Pythagorean index, 0.659 against
//0.632 in s_AB on SF6 - and the reason is not the molecule. On this route the matrix that gets averaged
//is a FREE ATOM's own density from occ's atomic SCF, not the molecular density, and a free atom is
//spherically symmetric: a single-determinant SCF on an open-shell atom only appears anisotropic because
//it has to put its electrons in particular m components. Leaving that unaveraged carried the axes of the
//INPUT FILE into every bond touching the atom, which is why the split follows x, y and z. So the exact
//rotational average now runs on the ANO route whatever symmetrize says, and -rgbi_no_sym keeps its
//meaning on the molecular route it was named for.
//
//WHAT WOULD MAKE THIS FAIL: putting the flag back in front of the free-atom average, or any change that
//lets an atomic reference depend on the orientation of the file - the O_h average reinstated for a
//spherical basis, an average applied per bond instead of per atom, a degenerate ANO set cut mid-shell.
//WHAT WOULD MAKE IT WRONG rather than red: nothing about the fixture. TeF6/def2-TZVP is octahedral to the
//last digit of the input geometry, so a spread in symmetry-equivalent bonds is the code's by construction.
//It asserts no reference numbers and needs no second program.
//
//Made red on purpose, this exact call on the binary before the change (AKL007, 0ab63b56dc00, 8 threads,
//4 s): the six rows come out as three pairs and the pairs are the Cartesian axes.
//     bonds 1,2   17.493  9.792  27.092  0.193  0.453  -0.430  0.624  52.525  51.608
//     bonds 3,4   17.493  9.798  27.100  0.191  0.448  -0.451  0.636  49.632  49.766
//     bonds 5,6   17.493  9.793  27.093  0.193  0.452  -0.431  0.624  52.443  51.556
//Eight of the nine columns split - worst 2.893 in Pyth. (52.525 against 49.632), 0.021 in Ion. - and the
//populations split 6.1e-03 (9.7922512 / 9.798316 / 9.7930333), which is what the EXPECT_NEAR above
//catches at 1e-7. Only column 0, the total Te population, survived. Running THIS test against a binary
//with the one condition put back fails 32 assertions: 28 of the 45 column equalities and 4 of the 5
//population checks, the Note check passing because the printed line is a separate change. SF6/def2-SVP
//through def2-QZVP split the same way on the same corner (s_AB 0.659 / 0.632 / 0.661), so it is not the
//ECP and not the f functions.
//And the fixed numbers are not merely self-consistent: this corner now prints the same nine numbers as
//the DEFAULT corner of the same file to every digit (17.476 9.802 27.083 0.194 0.449 -0.428 0.621 52.359
//51.502), which is the statement the printed Note makes - on the ANO route the flag decides nothing.
TEST(RgbiRobustnessTests, OctahedralTeF6IsExactlyOhOnTheAnoPathWithoutSymmetrization)
{
	const auto p = nos_test_repo_root() / "tests" / "RGBI_groups" / "tef6_tzvp.gbw";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/RGBI_groups/tef6_tzvp.gbw not found";
	std::string out;
	{
		CoutCapture cap;
		WFN wavy(p);
		//no group sets, NO symmetrization, ANO basis, no eigenvalues, no theta table
		Roby_information roby(wavy, {}, false, true, false, false);
		out = cap.str();
	}

	//A switch that decides nothing has to say so out loud, or it is the same defect as a switch nobody
	//reads: this run asked for -rgbi_no_sym and the ANO reference is averaged anyway.
	EXPECT_NE(out.find("-rgbi_no_sym does not change the atomic reference on the ANO route"), std::string::npos)
		<< "a run that asks for no symmetrization on the ANO route must be told that the free-atom "
		   "reference is spherical regardless";

	const double f0 = value_after(out, "Population of atom 1: ");
	ASSERT_TRUE(std::isfinite(f0)) << "no population for atom 1";
	for (int a = 2; a <= 6; a++)
		EXPECT_NEAR(f0, value_after(out, "Population of atom " + std::to_string(a) + ": "), 1e-7)
			<< "population of fluorine " << a << " against fluorine 1 on the ANO path without symmetrization";

	const vec first = bond_row(out, 0, 1);
	ASSERT_EQ(first.size(), 9u) << "no bond row 0 - 1";
	for (int b = 2; b <= 6; b++) {
		const vec row = bond_row(out, 0, b);
		ASSERT_EQ(row.size(), 9u) << "no bond row 0 - " << b;
		for (size_t i = 0; i < 9; i++)
			EXPECT_DOUBLE_EQ(first[i], row[i]) << "column " << i << " of Te-F bond 0 - " << b << " against "
				"0 - 1: the free-atom ANO reference is spherically averaged whether or not -rgbi_no_sym "
				"was given, so this corner reproduces the octahedron like the other three";
	}
}
