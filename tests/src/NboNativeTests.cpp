#include "pch.h"
#include "core/nbo.h"
#include "core/wfn_class.h"

#include <cmath>
#include <sstream>

//The in-house NBO route, checked by the sum rules it cannot break and still be right.
//
//Nbo47Tests drives parse_nbo_output(), so it checks the parser against stored gennbo output; NrtTests
//drives native_nrt() on a synthetic NAO basis, so it checks the algebra.  Nothing called
//native_nbo() at all - the one function that takes a real wavefunction and produces the tables a
//colleague reads.  That is the gap this file closes, and it closes it without a reference: every
//assertion here is an identity the analysis has to satisfy whatever the right answer turns out to be.
//
//Why these and not a golden: two of the defects found in this code were columns read or written one
//position off (the NRT RS column, an open-shell NAO table printing Spin where Energy was expected),
//and a golden file agrees with whatever the code did on the day it was written.  A sum that has to
//come out at the wavefunction's own electron count cannot go stale.

namespace {

	std::filesystem::path ethane_fixture()
	{
		const auto p = nos_test_repo_root() / "tests" / "TFVC" / "ethane.gbw";
		return std::filesystem::exists(p) ? p : std::filesystem::path{};
	}

	std::filesystem::path fchk_fixture()
	{
		const auto p = nos_test_repo_root() / "tests" / "alanine_occ" / "alanine.owf.fchk";
		return std::filesystem::exists(p) ? p : std::filesystem::path{};
	}

}  // namespace

//Ethane: 8 atoms, closed shell, no ECP, so the electron count the analysis must reproduce is the
//wavefunction's own MO occupation sum and the charges must add to the molecular charge.
//
//  1. every atom's core + valence + Rydberg is its total.  Not a tautology: the three come from the
//     NAO classifier and the total from the atom's population, so an NAO put in no class or in two
//     breaks this and nothing else in the output would show it.
//  2. the totals add to count_nr_electrons().  This is Tr(gamma) redistributed over atoms, so it
//     fails if the NAO transform is not norm-preserving or if an atom's block is missed.
//  3. the charges add to the molecular charge, which checks Z_eff - population in the other
//     direction and would catch ECP bookkeeping applied to a molecule that has none.
//  4. the NAO table and the population table are filled by two separate loops over two separate
//     structures; their electron counts have to agree.
//  5. no NAO may hold less than zero or more than two electrons - they are eigenvalues of a
//     one-particle density matrix in an orthonormal basis.
TEST(NboNativeTests, TheNativePopulationsObeyTheirSumRules)
{
	const auto p = ethane_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/TFVC/ethane.gbw not found";

	WFN wavy(p);
	ASSERT_FALSE(wavy.get_has_ECPs()) << "the charge sum rule below assumes Z_eff == Z";
	const double electrons = wavy.count_nr_electrons();
	ASSERT_GT(electrons, 0.0);

	NboOptions options;  //no NRT: this is about the populations, and NRT has its own tests
	std::ostringstream log;
	const NboResults r = native_nbo(wavy, options, log);

	ASSERT_EQ(static_cast<int>(r.npa.size()), wavy.get_ncen());
	double total_sum = 0.0, charge_sum = 0.0;
	for (const NboAtomPopulation &a : r.npa) {
		EXPECT_NEAR(a.core + a.valence + a.rydberg, a.total, 1.0e-8)
			<< "atom " << a.index << " " << a.element << ": the NAO classes do not partition it";
		total_sum += a.total;
		charge_sum += a.charge;
	}
	EXPECT_NEAR(total_sum, electrons, 1.0e-6);
	EXPECT_NEAR(charge_sum, static_cast<double>(wavy.get_charge()), 1.0e-6);

	ASSERT_FALSE(r.nao.empty());
	double nao_sum = 0.0;
	for (const NboNao &n : r.nao) {
		EXPECT_GE(n.occupancy, -1.0e-8) << "NAO " << n.index << " holds a negative population";
		EXPECT_LE(n.occupancy, 2.0 + 1.0e-8) << "NAO " << n.index << " holds more than two electrons";
		nao_sum += n.occupancy;
	}
	EXPECT_NEAR(nao_sum, total_sum, 1.0e-6)
		<< "the NAO table and the population table disagree about how many electrons there are";
}

//Every accepted orbital names centres that exist, and a two-centre orbital names two different atoms.
//
//What this test deliberately does NOT assert: that the orbital occupancies add to the electron count.
//That looks like the obvious invariant and it is not one here.  nbo.cpp OWSO-orthonormalises the Lewis
//set, but the BD* antibonds are built as c_B h_A - c_A h_B and are explicitly not orthogonal to it,
//and the RY complement is taken with Q = I - M (MtM)^-1 Mt against that non-orthonormal M.  Measured
//over 15 inputs the residual runs from -0.388 e (22 electrons, ECP) to +0.905 e (116 electrons, open
//shell), i.e. it is a property of the construction, not a failure.  Asserting it would give a check
//that is red on 9 of 15 inputs from the day it lands, and a check that is always red says nothing.
TEST(NboNativeTests, EveryAcceptedOrbitalNamesRealCentres)
{
	const auto p = ethane_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/TFVC/ethane.gbw not found";

	WFN wavy(p);
	NboOptions options;
	std::ostringstream log;
	const NboResults r = native_nbo(wavy, options, log);

	ASSERT_FALSE(r.orbitals.empty());
	//An occupancy is <phi|gamma|phi> for a normalised phi, so it is bounded by the largest eigenvalue
	//of the spin-summed density matrix whether or not phi is orthogonal to the rest of the set.  This
	//holds on all 15 inputs measured; the SUM does not, which is why it is not asserted above.
	for (const NboOrbital &o : r.orbitals) {
		EXPECT_GE(o.occupancy, -1.0e-8) << "orbital " << o.index << " " << o.type << " is negative";
		EXPECT_LE(o.occupancy, 2.0 + 1.0e-8) << "orbital " << o.index << " " << o.type << " is over two";
	}

	//every orbital names centres that exist, and a two-centre one names two different atoms
	for (const NboOrbital &o : r.orbitals) {
		ASSERT_FALSE(o.centers.empty()) << "orbital " << o.index << " " << o.type << " has no centre";
		for (const int c : o.centers) {
			EXPECT_GE(c, 1);
			EXPECT_LE(c, static_cast<int>(wavy.get_ncen()));
		}
		if (o.centers.size() == 2)
			EXPECT_NE(o.centers[0], o.centers[1]) << "orbital " << o.index << " bonds an atom to itself";
	}
}

//The .47 writer already computed two numbers that say its archive is wrong for a .fchk source - all
//228 AOs unnormalised in the overlap, Tr(P*S) = 7.579076 electrons against the 48 the occupied
//orbitals hold - and printed neither unless a progress log had been passed, which -nbo_native does
//only under -debug. Without this refusal the run exits 0 and prints an NPA table whose charges sum
//to +40.42 on a neutral molecule and whose largest NAO occupancy is 2.16 electrons.
//
//std::cout is redirected onto stderr inside the child because err_checkf writes the message to the
//stream it is handed and a death test matches stderr; without that the assertion could only check
//the exit code, which is the weaker half of what is being asserted.
TEST(NboNativeTests, AnArchiveThatDoesNotDescribeTheWavefunctionIsRefusedByName)
{
	const auto p = fchk_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/alanine_occ/alanine.owf.fchk not found";
	const auto out = std::filesystem::temp_directory_path() / "nos_nbo_native_refusal.47";

	EXPECT_EXIT(
		{
			std::cout.rdbuf(std::cerr.rdbuf());
			WFN wavy(p);
			wavy.write_nbo(out, false);
		},
		::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE),
		"does not describe this wavefunction");

	std::error_code ec;
	std::filesystem::remove(out, ec);
}

//The .47 writer's refusal on tests/alanine_occ/alanine.owf.fchk blames "the wrong normalisation
//convention", and for that file it is right: its shells are spherical, 228 declared and 228 built,
//and the gbw conversion factor brings Tr(P*S) from 7.579 to 47.938 of 48. For
//tests/NiP3_fchk/good.fchk the message is misleading. That file declares 964 basis functions because
//its 44 d and 21 f shells are cartesian; the libcint basis Int_Params builds counts 2l+1 per shell
//and holds 857. No constant can map one space onto the other - five candidate conventions were
//measured against Tr(P*S) before the counts were compared, and the best of them still left 367 AOs
//with a non-unit overlap diagonal, which is 220 spherical d plus 147 spherical f exactly.
//
//So the reader says it, and this pins both directions: the cartesian file must name its two counts,
//and the spherical file must stay silent, because a warning that fires on every fchk is a warning
//nobody reads. Made red on purpose by asserting 857 on the spherical file: it fails with the
//reader's silence printed beside it.
TEST(NboNativeTests, ACartesianFchkNamesTheTwoCountsAndASphericalOneStaysQuiet)
{
	const auto cartesian = nos_test_repo_root() / "tests" / "NiP3_fchk" / "good.fchk";
	const auto spherical = fchk_fixture();
	if (!std::filesystem::exists(cartesian) || spherical.empty())
		GTEST_SKIP() << "tests/NiP3_fchk/good.fchk or tests/alanine_occ/alanine.owf.fchk not found";

	std::ostringstream loud;
	WFN cart;
	ASSERT_TRUE(cart.read_fchk(cartesian, loud, false));
	const std::string said = loud.str();
	EXPECT_NE(said.find("964 basis functions"), std::string::npos) << said;
	EXPECT_NE(said.find("857 spherical functions"), std::string::npos) << said;
	EXPECT_NE(said.find("cartesian"), std::string::npos) << said;

	std::ostringstream quiet;
	WFN sph;
	ASSERT_TRUE(sph.read_fchk(spherical, quiet, false));
	EXPECT_EQ(quiet.str().find("WARNING"), std::string::npos)
		<< "a spherical fchk needs no warning, and got: " << quiet.str();
}
