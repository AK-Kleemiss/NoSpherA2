#include "pch.h"
#include "core/wfn_class.h"
#include "core/bondwise_analysis.h"

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
