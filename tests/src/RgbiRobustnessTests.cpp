#include "pch.h"
#include "core/tuning.h"
#include "core/wfn_class.h"
#include "core/bondwise_analysis.h"

#include <occ/core/parallel.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
#include <sstream>
#include <string>

//RGBI on the ANO route builds a free atom per centre and keeps a fixed number of its natural orbitals. Its
//reproducibility rests on two invariants: the free atom is in its ground state, and the kept subspace holds
//occupied orbitals only.

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

	//The nine numbers of one bond row of the RGBI table, or an empty vector if there is no such row. The indices
	//are right-aligned in their own fields ("   0 -   2   Au - Br  ..."), so they are read as numbers: the
	//literal "0 - 2" matches no row, and "0 -" would also match "10 -".
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

//A non-bonded closed-shell atom gets its own electrons back. A triplet free He (1s(1)2s(1)) has two exactly
//degenerate natural orbitals, and a rank-1 subspace then keeps whichever the threaded eigensolver returns first.
TEST(RgbiRobustnessTests, NonBondedHeliumKeepsItsTwoElectrons)
{
	const std::string out = rgbi_ano_output(false);
	if (out.empty())
		GTEST_SKIP() << "tests/TFVC/water.gbw not found";
	EXPECT_NEAR(value_after(out, "Population of atom 3: "), 2.0, 5e-3);
	//and the whole-system projection then accounts for all but a thousandth of the 12 electrons
	EXPECT_NEAR(value_after(out, "Number of electrons in Roby Analysis:  "), 11.899, 5e-3);
}

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

//The rank of an atomic subspace is fixed by the element; for H in a triple-zeta basis it is 1 of 14 and the
//other 13 natural orbitals are empty. Keeping an empty one adds an arbitrary direction from a degenerate null
//space, so no kept orbital may have a vanishing occupation.
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

//Cerium is beyond the shipped minimal basis, and occ's other guess route corrupts the heap for an unrestricted
//free atom, so that one case starts from the core Hamiltonian. Only that case: starting every free atom there
//converges the light atoms to a different state, which the goldens in BondwiseTests.cpp check.
TEST(RgbiRobustnessTests, CeriumFreeAtomRunsAndIsNotFallenBackOn)
{
	const auto p = nos_test_repo_root() / "tests" / "molden_file" / "Ce_full.molden";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/molden_file/Ce_full.molden not found";
	//The most expensive test here by orders of magnitude and the only one on the f-element free-atom path, so it
	//runs behind a switch and is named in the skip message.
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
	//Not asserted: pinned to one thread the Ce SCF runs out of iterations, unpinned it converges, so an expectation
	//either way would assert the node's thread count. The run warns instead, which
	//AFreeAtomThatRunsOutOfIterationsSaysSo checks on water.
}

//occ does not throw when an SCF runs out of iterations (scf_impl.h logs at error level and returns the last
//energy), so RGBI warns itself. NOS_RGBI_FREE_ATOM_MAXITER caps the iterations so the warning is checked on
//water, in both directions: a warning that is always printed is not a warning.
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
	set_tuning("NOS_RGBI_FREE_ATOM_MAXITER", "1");
	const std::string capped = rgbi_ano_output(false);
	set_tuning("NOS_RGBI_FREE_ATOM_MAXITER", nullptr);
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

//tests/RGBI/stdout is Tonto 26.01.05's Roby-Gould output for the archive next to it, run with the options it
//states (NAOs, no spherical averaging, def2-SVP, Cartesian d). With those options every printed number must be
//Tonto's to the two decimals Tonto prints; nothing else checks this analysis against an external program.
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

//The free-atom SCFs are pinned to one TBB thread, and the guard must restore the process as it found it or
//every later occ user in the binary stays serial. Before any tbb::global_control exists get_num_threads()
//already answers 1, so restoring "what get_num_threads() said" would install a control at 1 for good.
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

//Inversion-related Au, Br and P centres must retain matching RGBI populations and bond indices.
TEST(RgbiRobustnessTests, SymmetryEquivalentGoldCentresAgreeWithinPrintedPrecision)
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
		//The reported populations can differ by one unit in their sixth decimal.
		EXPECT_NEAR(first, second, 1.1e-6) << "populations of the symmetry-equivalent atoms " << a << " and " << a + 1;
	}

	//The omitted population is a difference of two projector populations.
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

//A free atom's density depends only on its basis and electron count, so identical centres share one SCF.
//water.gbw is O, two H and a non-bonded He: 4 centres and 3 distinct free atoms, so 3 SCFs and 1 cache hit.
TEST(RgbiRobustnessTests, IdenticalCentresShareOneFreeAtomScf)
{
	const auto p = water_he_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/TFVC/water.gbw not found";

	struct ScopedEnv {
		ScopedEnv() {
			set_tuning("NOS_RGBI_DEBUG", "1");
		}
		~ScopedEnv() {
			set_tuning("NOS_RGBI_DEBUG", nullptr);
		}
	} debug_on;

	//Cold, or the counts depend on test order: a cache warmed by a neighbour gives 0 SCFs and 4 hits.
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
	EXPECT_EQ(scfs, 3u) << "water plus helium has 3 distinct free atoms, so 3 SCFs and no more:\n" << out;
	EXPECT_EQ(hits, 1u) << "the second hydrogen is the only repeat and it must come from the cache:\n" << out;

	//NOS_RGBI_NO_FREEATOM_CACHE recomputes every centre in this process, so the cached bond table has something
	//to be identical to; that identity, not speed, is what justifies the cache.
	set_tuning("NOS_RGBI_NO_FREEATOM_CACHE", "1");
	const std::string uncached = rgbi_ano_output(false);
	set_tuning("NOS_RGBI_NO_FREEATOM_CACHE", nullptr);
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

//The distinct free-atom SCFs are independent, so NOS_RGBI_PARALLEL_FREEATOM runs them up front and
//concurrently. occ's SCF is not documented re-entrant, so the assertion is that the answer did not move, not
//that it is faster. The serial arm runs with the cache off: with it on, the serial arm fills the cache, the
//warm pass runs no SCF, and the parallel arm is compared with itself. The comparison is over the 14-digit
//NOS_RGBI_DEBUG lines (sumD, normD, SCF energy), as a three-decimal bond table hides a 1e-7 change in an
//atomic density.
TEST(RgbiRobustnessTests, WarmingTheFreeAtomCacheInParallelDoesNotChangeTheAnswer)
{
	const auto p = water_he_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/TFVC/water.gbw not found";

	const auto set_env = [](const char *name, const char *value) { set_tuning(name, *value ? value : nullptr); };
	const auto table = [](const std::string &s) {
		const size_t from = s.find("Atom Nr");
		return from == std::string::npos ? std::string() : s.substr(from);
		};
	//Every completed free-atom SCF by its numbers, label dropped: the two H of water collapse to one entry, and
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

	//The cache off, so these are SCFs actually run in sequence, and nothing is left for the parallel arm to find.
	set_env("NOS_RGBI_DEBUG", "1");
	set_env("NOS_RGBI_NO_FREEATOM_CACHE", "1");
	const std::string serial = rgbi_ano_output(false);
	set_env("NOS_RGBI_NO_FREEATOM_CACHE", "");
	ASSERT_FALSE(table(serial).empty()) << "no bond table from the serial run";

	//Cold on purpose: the cache lives for the process, so an earlier test would leave every atom answered, the
	//warm pass would run no SCF, and the 14-digit comparison would have nothing to compare.
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

//The determinism guard's effect, asserted directly. occ declares `inline int nthreads = 1` and installs no
//tbb::global_control until set_num_threads is called, so a guard that skips the call when get_num_threads()
//already reads 1 installs nothing and the reductions run across every core. That moves norm(D) only at
//roundoff and only in some runs, below every printed digit, so the observable is the pin itself;
//NOS_RGBI_NO_PIN turns it off and makes this red.
TEST(RgbiRobustnessTests, EveryFreeAtomScfRunsWithOccPinnedForReal)
{
	const auto p = water_he_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/TFVC/water.gbw not found";

	struct ScopedEnv {
		ScopedEnv() {
			set_tuning("NOS_RGBI_DEBUG", "1");
			set_tuning("NOS_RGBI_NO_FREEATOM_CACHE", "1");
		}
		~ScopedEnv() {
			set_tuning("NOS_RGBI_DEBUG", nullptr);
			set_tuning("NOS_RGBI_NO_FREEATOM_CACHE", nullptr);
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

//RGBI refuses shells beyond h, as its O_h symmetrization has no transform for them. The decision needs only
//the shell types, so it is taken from the basis before the overlap and the free-atom SCFs. err_checkf() exits
//the process, so this tests the question the guard asks, highest_shell_angular_momentum() over the molecule.
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

	//The defect shape guarded against is a check of the first atom only: in CuF2_i_func copper carries the i
	//shells and fluorine does not.
	const WFN i_on_the_second = wfn_of({ shell_atom("F", 9, 0.0, {0, 1, 2}),
										 shell_atom("Cu", 29, 3.5, {0, 1, 2, 3, 4, 5, 6}) });
	EXPECT_EQ(highest_shell_angular_momentum(i_on_the_second), 6)
		<< "an i shell on the second atom was not seen, so RGBI would run the whole overlap and only "
		"then refuse - which is the defect this replaced";

	//A basis with no shells at all has no highest l. The Roby_information constructor refuses that
	//case separately (a plain .wfn), and this must not turn into a 0 that looks supported.
	EXPECT_EQ(highest_shell_angular_momentum(wfn_of({ atom("H", {}, 1, 0.0, 0.0, 0.0, 1) })), -1);
}

//Only .gbw, .molden and a pure-shell .fchk reach RGBI: .wfx and .wfn carry no basis, a cartesian .fchk no
//contracted density matrix, and the Tonto fixture in tests/cytidine_tonto lacks its companion stdout.
TEST(RgbiRobustnessTests, TheRefusalNeverSuggestsTheFormatItIsRefusing)
{
	//the whole list, for the case where the refused file is none of them
	EXPECT_EQ(rgbi_supported_input_phrase(".wfn"), "a .gbw, a .molden or a .fchk");
	//and each working format drops itself, because a .gbw whose basis did not survive the
	//reader must not be answered with "use a .gbw"
	EXPECT_EQ(rgbi_supported_input_phrase(".gbw"), "a .molden or a .fchk");
	EXPECT_EQ(rgbi_supported_input_phrase(".molden"), "a .gbw or a .fchk");
	EXPECT_EQ(rgbi_supported_input_phrase(".fchk"), "a .gbw or a .molden");

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

//An exactly octahedral molecule has one bond orbit, so the six Te-F rows agree in every printed digit. Te's
//ANO occupations group as 3, 3, 3, 2, 3, so a fixed rank of 13 would cut a threefold degenerate set; which
//vectors of a degenerate set the diagonalizer returns is arbitrary, so the projector would break the
//symmetry. calculateAtomicNAO extends the rank to the end of the set it would cut. TeF6/def2-TZVP is
//octahedral to the last digit of the input geometry, so any spread is the code's; SF6 from def2-SVP to
//def2-QZVP stays Oh without the extension, so this fixture is the one on the ECP/degenerate-rank path.
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

	//The six fluorines are one orbit. A tolerance of 1e-7 rather than equality because these are parsed printed
	//digits: two values that differ below the seventh decimal can straddle a rounding boundary and print apart.
	const double f0 = value_after(out, "Population of atom 1: ");
	ASSERT_TRUE(std::isfinite(f0)) << "no population for atom 1";
	for (int a = 2; a <= 6; a++) {
		const double f = value_after(out, "Population of atom " + std::to_string(a) + ": ");
		ASSERT_TRUE(std::isfinite(f)) << "no population for atom " << a;
		EXPECT_NEAR(f0, f, 1e-7) << "population of fluorine " << a << " against fluorine 1";
	}

	//All nine columns agree exactly because the atomic reference is averaged over SO(3), not over the 48
	//operations of O_h: an O_h average of a shell with l >= 2 keeps e_g and t_2g apart, so the part of the atom's
	//anisotropy that lines up with the input file's axes survives and the axial fluorines differ from the
	//equatorial ones (the D4h pattern). By Schur's lemma the SO(3) average leaves one number per pair of shells of
	//equal l. It fails if an atom's orientation leaks back into its reference: an average per bond instead of per
	//atom, the O_h route for a spherical basis, a degenerate ANO set truncated mid-shell.
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

//On the NAO route without symmetrization no atomic reference is constructed, so the six bonds agree to the
//last bit in all nine columns: the pair decomposition, projectors, reductions and printing are exactly
//Oh-covariant on their own. Kept beside the test above so a break is attributable to one of the two. It fails
//on per-bond work that depends on the bond's orientation: a reduction order tied to atom index, a cutoff per
//bond rather than per orbit, an eigenvector phase leaking into a population. TeF6 carries an ECP and f
//functions. This is Tonto's own setting, as in TontoWaterRobyGouldNumbersAreReproduced.
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

//The ANO route without symmetrization: the averaged matrix is a free atom's density, and a free atom is
//spherical, while an unaveraged open-shell SCF carries the input file's x, y and z into every bond. So the
//exact rotational average runs whatever -rgbi_no_sym says, and this corner prints the default corner's
//numbers. It fails if the flag gates the free-atom average again, or on any orientation dependence of an
//atomic reference (see OctahedralTeF6HasOneBondOrbitNotThree).
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

//"Number of cutoff electrons" is the difference of two sums of order the electron count, so with nothing cut
//off it is a cancellation residual of arbitrary sign; Co2 is a fixture where it is only that residual. The
//sign is asserted, not zero, as a real ANO cutoff gives a legitimate positive number.
TEST(RgbiRobustnessTests, NoRgbiRunReportsANegativeNumberOfCutoffElectrons)
{
	const auto p = nos_test_repo_root() / "tests" / "molden_file" / "Co2.molden";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/molden_file/Co2.molden not found";
	std::string out;
	{
		CoutCapture cap;
		WFN wavy(p);
		Roby_information roby(wavy, {}, true, true, false, false);
		out = cap.str();
	}
	const std::string key = "Number of cutoff electrons:";
	const double cut = value_after(out, key);
	ASSERT_TRUE(std::isfinite(cut)) << "no cutoff line was printed - an ANO run is the only one that "
		"prints it, so either the default basis changed or the run did not get that far";
	EXPECT_GE(cut, 0.0) << "a negative count of cutoff electrons: nothing can be omitted a negative number "
		"of times, so this is either the cancellation residual leaking its arbitrary sign or a real defect "
		"in the cutoff projection";
	//And the printed token, because IEEE has -0.0 >= 0.0, so the check above passes on a printed "-0" -
	//the other half of what the clamp exists to prevent. The token and not the whole line, because a
	//legitimate small value prints as 3.55271e-15 and that minus belongs to the exponent.
	std::istringstream rest(out.substr(out.find(key) + key.size()));
	std::string token;
	ASSERT_TRUE(static_cast<bool>(rest >> token)) << "the cutoff line carries no value";
	EXPECT_NE(token.front(), '-') << "the number of cutoff electrons is printed with a leading minus: "
		<< token;
}
