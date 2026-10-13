#include "pch.h"
#include "core/nbo.h"
#include "core/wfn_class.h"

#include <cmath>
#include <fstream>
#include <map>
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

//A density whose Tr(P*S) is far from the electron count must be refused by name instead of yielding an
//NPA table with wrong charges. Without read_fchk's DM the .47 writer rebuilds the density from the
//primitives, where the first MO is stretched by 10 %.
//std::cout goes to stderr in the child: err_checkf writes to the stream it is handed and a death test
//matches stderr.
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
			wavy.set_dm(dMatrix2());
			for (int i = 0; i < wavy.get_nex(); i++)
				wavy.set_MO_coef(0, i, 1.1 * wavy.get_MO_coef(0, i));
			wavy.write_nbo(out, false);
		},
		::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE),
		"does not describe this wavefunction");

	std::error_code ec;
	std::filesystem::remove(out, ec);
}

//The same closed-shell fchk with the DM read_fchk builds is a valid archive (a refusal exits the process).
//F = S C^T eps C S belongs to that density only if C is orthonormal in the archive's S: Tr(P F) = sum n_m eps_m.
TEST(NboNativeTests, ASphericalFchkWritesAnArchive)
{
	const auto p = fchk_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/alanine_occ/alanine.owf.fchk not found";
	const auto out = std::filesystem::temp_directory_path() / "nos_nbo_native_fchk.47";
	WFN wavy(p);
	wavy.write_nbo(out, false);
	ASSERT_TRUE(std::filesystem::exists(out));
	const NboInput in = read_file47(out);
	ASSERT_EQ(in.fock.size(), 1u);
	double tr = 0.0, sum = 0.0;
	for (int i = 0; i < in.n; i++)
		for (int j = 0; j < in.n; j++)
			tr += in.density[0](i, j) * in.fock[0](j, i);
	for (int m = 0; m < wavy.get_nmo(); m++)
		sum += wavy.get_MO_occ(m) * wavy.get_MO_energy(m);
	EXPECT_NEAR(tr, sum, 1e-6 * std::abs(sum));
	std::error_code ec;
	std::filesystem::remove(out, ec);
}

//-nrt_target 0 on ethane: the second cycle walks from the weighted structures and states the Frank-Wolfe floor of
//its fit at the warm start, which the fit cannot go below
TEST(NboNativeTests, AGrownNrtCycleStatesAFloorItsFitRespects)
{
	const auto p = ethane_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/TFVC/ethane.gbw not found";
	WFN wavy(p);
	NboOptions options;
	options.nrt = true;
	options.nrt_target = 0.0;
	std::ostringstream log;
	const NboResults r = native_nbo(wavy, options, log);
	int floors = 0;
	for (const std::string& n : r.nrt.notes) {
		const size_t b = n.find("bound before the fit ");
		if (b == std::string::npos) continue;
		floors++;
		EXPECT_LE(std::stod(n.substr(b + 21)), std::stod(n.substr(n.find("D(w) = ") + 7)) + 1e-6) << n;
	}
	EXPECT_GT(floors, 0);
}

//CH3 doublet, `occ scf` B3LYP(VWN5)/def2-TZVP at the geometry of tests/nbo_reference_v2/data/ch3, where NBO 7
//splits each C-H bond 0.404/0.096 (covalent/ionic) in the spin holding the radical and 0.432/0.068 in the other,
//with D 0.0103 and 0.0065. OCC puts the radical in beta, NBO 7's ORCA run in alpha.
TEST(NboNativeTests, MethylRadicalSplitsItsBondsPerSpinAsNbo7Does)
{
	const auto p = nos_test_repo_root() / "tests" / "ch3_occ" / "ch3.owf.fchk";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/ch3_occ/ch3.owf.fchk not found";
	WFN wavy(p);
	NboOptions options;
	options.nrt = true;
	std::ostringstream log;
	const NboResults r = native_nbo(wavy, options, log);
	ASSERT_TRUE(r.open_shell);
	EXPECT_NEAR(r.npa[0].spin_density, -1.08333, 2e-3);
	int bonds = 0;
	for (const NboBondOrder& b : r.nrt.bond_orders) {
		if (b.diagonal || b.total == 0.0) continue;
		const double cov = b.spin == "beta" ? 0.404 : 0.432;
		EXPECT_NEAR(b.covalent, cov, 3e-3) << b.spin << " C1-H" << b.atom2;
		EXPECT_NEAR(b.covalent + b.ionic, 0.5, 1e-9);
		bonds++;
	}
	EXPECT_EQ(bonds, 6);
	ASSERT_EQ(r.nrt.spin_d.size(), 2u);
	EXPECT_GT(r.nrt.spin_d[1].d_0, r.nrt.spin_d[0].d_0) << "the radical's spin carries the larger D";
}

//HI from `occ scf` B3LYP(VWN5)/def2-TZVP: the fchk's "Nuclear charges" give iodine 25, Z minus its 28-electron ECP
//core. The atom keeps Z = 53 with the core as its ECP count, so NBO sees iodine and its 18 remaining core electrons.
TEST(NboNativeTests, AnEcpFchkKeepsTheElementAndItsCore)
{
	const auto p = nos_test_repo_root() / "tests" / "hi_occ" / "hi.owf.fchk";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/hi_occ/hi.owf.fchk not found";
	WFN wavy(p);
	EXPECT_EQ(wavy.get_atom_charge(0), 53);
	EXPECT_TRUE(wavy.get_has_ECPs());
	EXPECT_EQ(wavy.get_atom_ECP_electrons(0), 28);
	EXPECT_NEAR(wavy.count_nr_electrons(), 26.0, 1e-6);
	NboOptions options;
	std::ostringstream log;
	const NboResults r = native_nbo(wavy, options, log);
	ASSERT_EQ(r.npa.size(), 2u);
	EXPECT_EQ(r.npa[0].element, "I");
	EXPECT_NEAR(r.npa[0].core, 18.0, 1e-3);
	EXPECT_NEAR(r.npa[0].charge + r.npa[1].charge, 0.0, 1e-4);
	EXPECT_LT(std::abs(r.npa[0].charge), 0.5);
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

namespace {

	std::filesystem::path write_text(const std::string& name, const std::string& text)
	{
		const auto p = std::filesystem::temp_directory_path() / name;
		std::ofstream(p) << text;
		return p;
	}

	//the Lewis set of a native result as -nrt_parent JSON, bare 1-based indices
	std::string lewis_json(const NboResults& r)
	{
		std::map<int, int> lp;
		std::map<std::pair<int, int>, int> bd;
		for (const NboOrbital& o : r.orbitals)
			if (o.type == "LP") lp[o.centers[0]]++;
			else if (o.type == "BD") bd[{ std::min(o.centers[0], o.centers[1]), std::max(o.centers[0], o.centers[1]) }]++;
		std::string s = "{\"format\": \"NoSpherA2 Lewis structure\", \"version\": 1, \"lone_pairs\": {";
		for (const auto& [a, n] : lp) s += (s.back() == '{' ? "\"" : ", \"") + std::to_string(a) + "\": " + std::to_string(n);
		s += "}, \"bonds\": [";
		for (const auto& [ab, n] : bd)
			s += (s.back() == '[' ? "[" : ", [") + std::to_string(ab.first) + ", " + std::to_string(ab.second) + ", " +
				 std::to_string(n) + "]";
		return s + "]}";
	}

	//O1 H2 H3 with one O core NAO: water's 5 pairs are 1 core + 2 lone pairs + 2 bonds
	NAOResult water_naos()
	{
		NAOResult nao;
		const int Z[3] = { 8, 1, 1 };
		for (int a = 0; a < 3; a++) {
			NAOAtom at;
			at.index = a;
			at.Z = Z[a];
			nao.atoms.push_back(at);
			NAO o;
			o.atom = a;
			o.type = a == 0 ? NAOClass::Core : NAOClass::Valence;
			nao.orbitals.push_back(o);
		}
		return nao;
	}

}  // namespace

//$CHOOSE with the search's own structure must change nothing: the constrained ladder looks only at the
//blocks the structure names, and the free ladder's choice was one of them at every step.
TEST(NboNativeTests, NrtParentFedTheSearchsOwnStructureReproducesTheWeights)
{
	const auto p = ethane_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/TFVC/ethane.gbw not found";
	WFN wavy(p);
	NboOptions free;
	free.nrt = true;
	std::ostringstream log0;
	const NboResults own = native_nbo(wavy, free, log0);
	ASSERT_TRUE(own.nrt.present);

	NboOptions chosen = free;
	chosen.nrt_parent = write_text("nos_nrt_parent_own.json", lewis_json(own));
	std::ostringstream log1;
	const NboResults again = native_nbo(wavy, chosen, log1);
	EXPECT_NE(log1.str().find("NRT parent from"), std::string::npos) << log1.str();
	EXPECT_NE(log1.str().find("the same structure"), std::string::npos) << log1.str();
	ASSERT_EQ(again.nrt.weights.size(), own.nrt.weights.size());
	for (size_t i = 0; i < own.nrt.weights.size(); i++)
		EXPECT_NEAR(again.nrt.weights[i].weight_fraction, own.nrt.weights[i].weight_fraction, 1.0e-10) << i;
	ASSERT_EQ(again.nrt.bond_orders.size(), own.nrt.bond_orders.size());
	for (size_t i = 0; i < own.nrt.bond_orders.size(); i++)
		EXPECT_NEAR(again.nrt.bond_orders[i].total, own.nrt.bond_orders[i].total, 1.0e-10) << i;
	std::filesystem::remove(chosen.nrt_parent);
}

//Ethane's hyperconjugative structure H(a)+ C=C H(b)-: the Lewis set must carry it and the log name the change.
TEST(NboNativeTests, NrtParentPlacesADifferentStructureAndTheLogSaysWhatChanged)
{
	const auto p = ethane_fixture();
	if (p.empty())
		GTEST_SKIP() << "tests/TFVC/ethane.gbw not found";
	WFN wavy(p);
	NboOptions free;
	free.nrt = true;
	std::ostringstream log0;
	const NboResults own = native_nbo(wavy, free, log0);
	int c1 = 0, c2 = 0, ha = 0, hb = 0;
	for (const NboOrbital& o : own.orbitals)
		if (o.type == "BD" && o.hybrids[0].element == "C" && o.hybrids[1].element == "C") {
			c1 = o.centers[0];
			c2 = o.centers[1];
		}
	ASSERT_GT(c1, 0);
	std::string bonds = "[\"C" + std::to_string(c1) + "\", \"C" + std::to_string(c2) + "\", 2]";
	for (const NboOrbital& o : own.orbitals) {
		if (o.type != "BD" || o.hybrids[0].element == o.hybrids[1].element) continue;
		const int c = o.hybrids[0].element == "C" ? o.centers[0] : o.centers[1];
		const int h = o.centers[0] + o.centers[1] - c;
		if (c == c1 && !ha) { ha = h; continue; }
		if (c == c2 && !hb) { hb = h; continue; }
		bonds += ", [\"C" + std::to_string(c) + "\", \"H" + std::to_string(h) + "\", 1]";
	}
	ASSERT_GT(ha * hb, 0);
	NboOptions chosen = free;
	chosen.nrt_parent = write_text("nos_nrt_parent_hyper.json",
		"{\"atoms\": 8, \"lone_pairs\": {\"H" + std::to_string(hb) + "\": 1}, \"bonds\": [" + bonds + "]}");
	std::ostringstream log1;
	const NboResults r = native_nbo(wavy, chosen, log1);
	int cc = 0, lp_hb = 0, bd_ha = 0;
	for (const NboOrbital& o : r.orbitals) {
		if (o.type == "BD" && o.centers == std::vector<int>{ std::min(c1, c2), std::max(c1, c2) }) cc++;
		if (o.type == "LP" && o.centers[0] == hb) lp_hb++;
		if (o.type == "BD" && (o.centers[0] == ha || o.centers[1] == ha)) bd_ha++;
	}
	EXPECT_EQ(cc, 2);
	EXPECT_EQ(lp_hb, 1);
	EXPECT_EQ(bd_ha, 0);
	EXPECT_TRUE(r.nrt.present);
	const std::string said = log1.str();
	EXPECT_NE(said.find("NRT parent from " + chosen.nrt_parent.string()), std::string::npos) << said;
	EXPECT_NE(said.find("4 entries differ"), std::string::npos) << said;
	std::filesystem::remove(chosen.nrt_parent);
}

TEST(NboNativeTests, NrtParentFilesAreReadOrRefusedWithTheOffendingEntry)
{
	const NAOResult nao = water_naos();
	const std::string head = "{\"format\": \"NoSpherA2 Lewis structure\", \"version\": 1, \"wavefunction\": \"water.gbw\", "
							 "\"atoms\": 3, \"lone_pairs\": {\"O1\": 2}, \"bonds\": ";
	const auto read = [&](const std::string& text, ivec2& topo, std::string* warned = nullptr) {
		const auto f = write_text("nos_nrt_parent_case.json", text);
		std::ostringstream log;
		const std::string why = read_lewis_structure(f, nao, 5, "water.gbw", topo, log);
		if (warned) *warned = log.str();
		std::filesystem::remove(f);
		return why;
	};
	ivec2 topo;
	ASSERT_EQ(read(head + "[[\"O1\", \"H2\", 1], [3, \"o1\", 1]]}", topo), "");
	EXPECT_EQ(topo, (ivec2{ { 2, 1, 1 }, { 1, 0, 0 }, { 1, 0, 0 } }));
	//a whole .native.nbo.json, and a structure written for another file is used with a warning
	std::string warned;
	ASSERT_EQ(read("{\"name\": \"w\", \"parent_lewis\": {\"wavefunction\": \"other.gbw\", \"lone_pairs\": {\"1\": 2}, "
				   "\"bonds\": [[\"O1\", \"H2\", 1], [\"O1\", \"H3\", 1]]}}", topo, &warned), "");
	EXPECT_NE(warned.find("other.gbw"), std::string::npos) << warned;
	const std::vector<std::pair<std::string, std::string>> bad = {
		{ head + "[[\"O1\", \"H2\", 1], [\"O1\", \"H4\", 1]]}", "atom \"H4\" is outside 1..3" },
		{ head + "[[\"O1\", \"H2\", 1], [\"O1\", \"C3\", 1]]}", "atom 3 is H, not C" },
		{ head + "[[\"O1\", \"H2\", 1], [\"H2\", \"O1\", 1]]}", "[\"H2\",\"O1\",1] is listed twice" },
		{ head + "[[\"O1\", \"H2\", 4], [\"O1\", \"H3\", 1]]}", "the order must be 1, 2 or 3" },
		{ head + "[[\"O1\", \"H2\", 1]]}", "1 bond pairs + 2 lone pairs + 1 core pairs hold 8 electrons, "
										   "the wavefunction has 10: 1 pair(s) too few" },
		{ head + "[[\"O1\", \"H2\", 1], [\"O1\", \"H3\", 1]], \"three_centre\": [[\"O1\", \"H2\", \"H3\"]]}", "three_centre" },
		{ "{\"alpha\": {}, \"beta\": {}}", "closed-shell" },
		{ "{\"atoms\": 4, \"bonds\": []}", "\"atoms\" is 4, the wavefunction has 3" },
		{ "{\"lonepairs\": {}}", "unknown key \"lonepairs\"" },
		{ "{\"bonds\": [[\"O1\", \"H2\", 1]", "not valid JSON" },
	};
	for (const auto& [text, expected] : bad) {
		const std::string why = read(text, topo);
		EXPECT_NE(why.find(expected), std::string::npos) << text << "\n  gave: " << why;
		EXPECT_EQ(why.find('\n'), std::string::npos) << "not one line: " << why;
	}
}

//Ni(CO)4, B3LYP/def2-TZVP (def2-SVP does not show the second trap): NBO 7 puts one dative Ni-C bond
//and d10 lone pairs on Ni. The search used to accept a Ni d pair plus a C lone pair as Ni-C bonds
//(four electrons, two centres), and NRT then built the
//parent with a bond slot that took a d pair before the lone pairs were placed, so the parent it fitted
//was not the structure NBO found and dropped out of the weights.
TEST(NboNativeTests, NickelTetracarbonylKeepsOneDativeBondAndItsParent)
{
	const auto p = nos_test_repo_root() / "tests" / "NiCO4_gbw" / "ni_co_4.gbw";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/NiCO4_gbw/ni_co_4.gbw not found";
	WFN wavy(p);
	const double electrons = wavy.count_nr_electrons();
	NboOptions options;
	options.nrt = true;
	std::ostringstream log;
	const NboResults r = native_nbo(wavy, options, log);
	int ni_c = 0;
	double lewis = 0.0;
	for (const NboOrbital &o : r.orbitals) {
		if (o.type == "CR" || o.type == "LP" || o.type == "BD") lewis += o.occupancy;
		if (o.type == "BD" && o.centers.size() == 2) {
			const std::string a = r.npa[o.centers[0] - 1].element, b = r.npa[o.centers[1] - 1].element;
			if ((a == "Ni" && b == "C") || (a == "C" && b == "Ni")) ni_c++;
		}
	}
	EXPECT_EQ(ni_c, 1) << log.str();
	ASSERT_TRUE(r.nrt.present);
	//D(0) is the parent's residual: n_NAO D(0)^2 = rho_NL(parent) + 0.001 in NBO's units for an idempotent density
	ASSERT_EQ(r.nrt.spin_d.size(), 1u);
	const double d0 = r.nrt.spin_d[0].d_0;
	EXPECT_NEAR(static_cast<double>(r.nao.size()) * d0 * d0 - 0.001, electrons - lewis, 1.0e-3) << "NRT built a different parent than NBO found";
	//Td: one Ni-C and one C-O bond order however the walk's budget cut the candidates
	std::map<std::string, std::pair<double, double>> range;
	for (const NboBondOrder &b : r.nrt.bond_orders) {
		if (b.diagonal || b.total < 0.1) continue;
		std::string a = r.npa[b.atom1 - 1].element, c = r.npa[b.atom2 - 1].element;
		if (a > c) std::swap(a, c);
		auto it = range.emplace(a + "-" + c, std::make_pair(b.total, b.total)).first;
		it->second = { std::min(it->second.first, b.total), std::max(it->second.second, b.total) };
	}
	for (const auto &x : range)
		EXPECT_NEAR(x.second.first, x.second.second, 1e-6) << x.first;
	//every isometry holds in the NAO basis with ORCA's f phases, and images start from their source's orbitals rotated
	const auto noted = [&](const std::string &s) {
		return std::any_of(r.nrt.notes.begin(), r.nrt.notes.end(), [&](const std::string &n) { return n.find(s) != std::string::npos; });
	};
	EXPECT_FALSE(noted("miss the NAO density"));
	EXPECT_TRUE(noted("symmetry image's orbitals, rotated"));
	//Ni-C is barely determined: weights within 0.1 % of the optimal D(w) move it by more than 0.1 either way
	for (const NboBondOrder &b : r.nrt.bond_orders) {
		if (b.diagonal || b.total < 0.1) continue;
		EXPECT_LE(b.range_lo, b.total + 1e-9);
		EXPECT_GE(b.range_hi, b.total - 1e-9);
		if (r.npa[b.atom1 - 1].element == "Ni") EXPECT_GT(b.range_hi - b.range_lo, 0.2);
	}
}
