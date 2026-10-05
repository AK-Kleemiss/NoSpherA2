#include "pch.h"
#include "core/nbo.h"
#include "core/wfn_class.h"

#include <cmath>
#include <sstream>

//native_nbo() on a real wavefunction, checked by sum rules rather than a golden: a column read one
//position off agrees with a golden written from the same code, but not with the electron count.

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

//Ethane, closed shell, no ECP: the totals must reproduce the MO occupation sum and the charges the
//molecular charge.  core + valence + Rydberg == total catches an NAO the classifier puts in no class or
//two; the NAO table and the population table come from separate loops; NAO occupancies are eigenvalues
//of a one-particle density in an orthonormal basis, so lie in [0, 2].
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

//Lewis, antibonds and Rydbergs together are an orthonormal basis of the NAO space, so the occupancies
//sum to Tr(gamma), the electron count.
TEST(NboNativeTests, EveryAcceptedOrbitalNamesRealCentres)
{
	const auto p = ethane_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/TFVC/ethane.gbw not found";

	WFN wavy(p);
	const double electrons = wavy.count_nr_electrons();
	NboOptions options;
	std::ostringstream log;
	const NboResults r = native_nbo(wavy, options, log);

	ASSERT_FALSE(r.orbitals.empty());
	double occupancy_sum = 0.0;
	for (const NboOrbital &o : r.orbitals) occupancy_sum += o.occupancy;
	EXPECT_NEAR(occupancy_sum, electrons, 1.0e-6);
	//<phi|gamma|phi> for normalised phi is bounded by gamma's largest eigenvalue, orthogonal or not
	for (const NboOrbital &o : r.orbitals) {
		EXPECT_GE(o.occupancy, -1.0e-8) << "orbital " << o.index << " " << o.type << " is negative";
		EXPECT_LE(o.occupancy, 2.0 + 1.0e-8) << "orbital " << o.index << " " << o.type << " is over two";
	}

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

//A .fchk archive whose AOs are unnormalised in the overlap (Tr(P*S) far from the electron count) must
//be refused by name instead of yielding an NPA table with wrong charges.  std::cout goes to stderr in
//the child: err_checkf writes to the stream it is handed and a death test matches stderr.
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

//A cartesian fchk declares more functions than the 2l+1 per shell Int_Params builds, and no
//normalisation constant maps one onto the other, so the reader names both counts; a spherical fchk
//stays silent, since a warning on every fchk is read by nobody.
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
