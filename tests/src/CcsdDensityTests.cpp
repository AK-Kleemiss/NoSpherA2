//CCSD densities from both sources NoSpherA2 reads: occ's CCSD with the Lambda-equation density (ccsd_density.cpp)
//and ORCA's unrelaxed CCSD natural orbitals (.mdci.nat). Water cc-pVDZ, all electrons, the geometry of
//tests/ccsd_water/water_ccsd.inp. References: PySCF RHF 1e-12, CCSD 1e-11, ccsd_lambda 1e-10, make_rdm1; for
//CCSD(T) ccsd_t, ccsd_t_lambda_slow 1e-10 and ccsd_t_rdm_slow in its response form (for_grad=True).
#include "pch.h"
#include <gtest/gtest.h>

#include "core/ccsd_density.h"
#include "core/convenience.h"
#include "core/wfn_class.h"
#include <occ/io/xyz.h>
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>
#undef I

#include <algorithm>
#include <array>
#include <sstream>

namespace
{
	const char *water = "3\n\nO 0.000000 0.000000 0.117300\nH 0.000000 0.757200 -0.469200\nH 0.000000 -0.757200 -0.469200\n";
	constexpr double ref_e_corr = -0.2133274269;
	constexpr double ref_l1 = 0.0132555089, ref_l2 = 0.2125482514;
	constexpr std::array<double, 24> ref_occ{
		1.999916791935, 1.984591127484, 1.973055495266, 1.966977132722, 1.964441817473, 0.026576000090,
		0.024606534618, 0.016958989231, 0.011482604142, 0.005877291267, 0.005610198708, 0.005102431406,
		0.004501850550, 0.004391650529, 0.001153986147, 0.001135571676, 0.000733257552, 0.000666305288,
		0.000595994224, 0.000541860862, 0.000498173499, 0.000488628631, 0.000056697599, 0.000039609101};
	constexpr double ref_e_t = -0.0030587074;  //ORCA 6.1.0 on the same job: -0.003058698
	constexpr double ref_t_l1 = 0.0135673060, ref_t_l2 = 0.2202679718;
	constexpr std::array<double, 24> ref_t_occ{
		1.999915560472, 1.984126920375, 1.971987096386, 1.965109636192, 1.962228897254, 0.028145280442,
		0.025910267256, 0.017475185700, 0.011824173422, 0.006209180127, 0.005911724793, 0.005345731873,
		0.004639591047, 0.004522823917, 0.001290017319, 0.001240118379, 0.000812346188, 0.000788191094,
		0.000716511774, 0.000597192399, 0.000545450072, 0.000541000350, 0.000069385498, 0.000047717671};

	//the -occ path's SCF thresholds (1e-6) are looser than the references; the test pins the Lambda code, not the SCF
	occ::qm::Wavefunction tight_rhf()
	{
		spdlog::set_level(spdlog::level::err);
		const auto mol = occ::io::molecule_from_xyz_string(water);
		auto basis = occ::gto::AOBasis::load(mol.atoms(), "cc-pvdz");
		basis.set_pure(true);
		occ::qm::HartreeFock hf(basis);
		occ::qm::SCF<occ::qm::HartreeFock> scf(hf, occ::qm::SpinorbitalKind::Restricted);
		scf.convergence_settings.energy_threshold = 1e-12;
		scf.convergence_settings.commutator_threshold = 1e-9;
		scf.compute_scf_energy();
		return scf.wavefunction();
	}

	const CcsdDensityResult &occ_ccsd()
	{
		static const CcsdDensityResult result = [] {
			CcsdDensityOptions opts;
			opts.n_frozen = 0;
			opts.tol = 1e-9;  //PySCF comparisons at 1e-8, not ORCA's 1e-5
			std::ostringstream log;
			return ccsd_natural_orbitals(tight_rhf(), opts, log);
		}();
		return result;
	}

	std::filesystem::path orca_nat()
	{
		return nos_test_repo_root() / "tests" / "ccsd_water" / "water_ccsd.mdci.nat";
	}
}

TEST(CcsdDensity, OccLambdaDensityMatchesPyscf)
{
	if (!ensure_occ_data_path(nullptr)) GTEST_SKIP() << "OCC_DATA_PATH not set";
	const CcsdDensityResult &r = occ_ccsd();
	EXPECT_EQ(r.n_frozen, 0);
	EXPECT_NEAR(r.e_corr, ref_e_corr, 1e-8);
	EXPECT_NEAR(r.l1_norm, ref_l1, 1e-8);
	EXPECT_NEAR(r.l2_norm, ref_l2, 1e-8);
	ASSERT_EQ(r.occupations.size(), static_cast<long>(ref_occ.size()));
	for (size_t k = 0; k < ref_occ.size(); k++)
		EXPECT_NEAR(r.occupations(k), ref_occ[k], 1e-8) << "natural orbital " << k;
	EXPECT_NEAR(r.occupations.sum(), 10.0, 1e-10);
}

//.nat goes to the gbw reader, and the occupations it hands on are ORCA's CCSD ones, not the 2/0 of the reference
TEST(CcsdDensity, OrcaNaturalOrbitalsCarryTheCcsdOccupations)
{
	if (!std::filesystem::exists(orca_nat())) GTEST_SKIP() << "fixture missing: " << orca_nat().string();
	const WFN orca(orca_nat());
	ASSERT_EQ(orca.get_nmo(), static_cast<int>(ref_occ.size()));
	std::vector<double> occ;
	for (int i = 0; i < orca.get_nmo(); i++) occ.push_back(orca.get_MO_occ(i));
	std::sort(occ.rbegin(), occ.rend());
	double sum = 0.0;
	for (size_t k = 0; k < ref_occ.size(); k++)
	{
		EXPECT_NEAR(occ[k], ref_occ[k], 1e-6) << "natural orbital " << k;
		sum += occ[k];
	}
	EXPECT_NEAR(sum, 10.0, 1e-6);
}

//ORCA's property file says which density a .mdci.nat holds; a CCSD(T) run's is MDCI(SD), and the log must say so
TEST(CcsdDensity, OrcaTriplesRunIsReportedAsCcsdDensity)
{
	if (!std::filesystem::exists(orca_nat())) GTEST_SKIP() << "fixture missing: " << orca_nat().string();
	const std::filesystem::path dir = std::filesystem::temp_directory_path() / "nos_mdci_density_label";
	std::filesystem::create_directories(dir);
	std::filesystem::copy_file(orca_nat(), dir / "w.mdci.nat", std::filesystem::copy_options::overwrite_existing);
	const auto read_with = [&](const std::string &energy_method) {
		std::ofstream(dir / "w.property.txt") << "   &Method [&Type \"String\"] \"" << energy_method << "\"\n"
			<< "   &Method [&Type \"String\"] \"MDCI(SD)\"\n   &Level [&Type \"String\"] \"Unrelaxed density\"\n";
		WFN wfn;
		std::ostringstream log;
		wfn.read_known_wavefunction_format(dir / "w.mdci.nat", log);
		return log.str();
	};
	const std::string ccsd = read_with("MDCI(SD)"), ccsd_t = read_with("MDCI(SD(T))");
	EXPECT_NE(ccsd.find("ORCA density in w.mdci.nat: MDCI(SD), Unrelaxed density"), std::string::npos) << ccsd;
	EXPECT_EQ(ccsd.find("WARNING"), std::string::npos) << ccsd;
	EXPECT_NE(ccsd_t.find("without (T)"), std::string::npos) << ccsd_t;
	std::filesystem::remove_all(dir);
}

//the two independent CCSD implementations must give one density, and that density must not be the HF one
TEST(CcsdDensity, OccAndOrcaGiveTheSameDensity)
{
	if (!ensure_occ_data_path(nullptr)) GTEST_SKIP() << "OCC_DATA_PATH not set";
	if (!std::filesystem::exists(orca_nat())) GTEST_SKIP() << "fixture missing: " << orca_nat().string();
	const WFN occ_wfn(occ_ccsd().natural), orca(orca_nat()), hf(tight_rhf());
	ASSERT_EQ(occ_wfn.get_ncen(), 3);
	for (unsigned a = 0; a < 3; a++)
		for (unsigned x = 0; x < 3; x++)
			ASSERT_NEAR(occ_wfn.get_atom_coordinate(a, x), orca.get_atom_coordinate(a, x), 1e-5);

	auto at = [&](const unsigned a, const double dx, const double dy, const double dz) {
		return d3{ orca.get_atom_coordinate(a, 0) + dx, orca.get_atom_coordinate(a, 1) + dy, orca.get_atom_coordinate(a, 2) + dz };
	};
	std::vector<d3> points{ at(0, 0.02, 0.0, 0.03), at(0, 0.6, 0.0, 0.3), at(0, 0.0, 0.0, 1.5), at(1, 0.0, 0.05, 0.1), at(2, 0.0, -0.4, 0.4) };
	points.push_back(d3{ 0.5 * (points[0][0] + points[3][0]), 0.5 * (points[0][1] + points[3][1]), 0.5 * (points[0][2] + points[3][2]) });
	double hf_dev = 0.0;
	for (const d3 &p : points)
	{
		const double rho = orca.compute_dens(p);
		EXPECT_NEAR(occ_wfn.compute_dens(p) / rho, 1.0, 1e-5) << "at " << p[0] << " " << p[1] << " " << p[2];
		hf_dev = std::max(hf_dev, std::abs(hf.compute_dens(p) / rho - 1.0));
	}
	EXPECT_GT(hf_dev, 1e-3);
}

//clusters that span the molecule (every LMO in every cluster, every atom in every domain) must give the canonical
//density back: each cluster is then the full problem in rotated orbitals, and the partition sums to the whole
TEST(CcsdDensity, LocalClustersSpanningTheMoleculeReproduceCanonical)
{
	if (!ensure_occ_data_path(nullptr)) GTEST_SKIP() << "OCC_DATA_PATH not set";
	const occ::qm::Wavefunction scf = tight_rhf();
	for (const bool triples : { false, true })
	{
		CcsdDensityOptions opts;
		opts.n_frozen = 0;
		opts.tol = 1e-9;  //PySCF comparisons at 1e-8, not ORCA's 1e-5
		opts.triples = triples;
		opts.cluster_radius = 1e3;
		opts.domain_threshold = -1.0;
		opts.pair_threshold = -1.0;
		opts.lno_threshold = -1.0;
		opts.aux_basis = "def2-svp-rifit";  //screening only, all of it off here; occ's cc-pvdz-rifit lacks oxygen
		std::ostringstream log;
		const CcsdDensityResult canonical = ccsd_natural_orbitals(scf, opts, log);
		const CcsdDensityResult local = local_ccsd_natural_orbitals(scf, opts, log);
		ASSERT_EQ(local.occupations.size(), canonical.occupations.size());
		for (long k = 0; k < canonical.occupations.size(); k++)
			EXPECT_NEAR(local.occupations(k), canonical.occupations(k), 1e-8) << (triples ? "CCSD(T)" : "CCSD") << " natural orbital " << k;
		EXPECT_NEAR(local.e_corr, canonical.e_corr, 1e-9);
	}
}

//the (T) terms in Lambda and in the density; E(T) from the same arrays pins them against PySCF's and ORCA's (T),
//and the occupations must move off the CCSD ones (by up to 1.6e-3 here)
TEST(CcsdDensity, OccTriplesDensityMatchesPyscf)
{
	if (!ensure_occ_data_path(nullptr)) GTEST_SKIP() << "OCC_DATA_PATH not set";
	CcsdDensityOptions opts;
	opts.n_frozen = 0;
	opts.tol = 1e-9;  //PySCF comparisons at 1e-8, not ORCA's 1e-5
	opts.triples = true;
	std::ostringstream log;
	const CcsdDensityResult r = ccsd_natural_orbitals(tight_rhf(), opts, log);
	EXPECT_NEAR(r.e_corr, ref_e_corr, 1e-8);
	EXPECT_NEAR(r.e_triples, ref_e_t, 1e-9);
	EXPECT_NEAR(r.l1_norm, ref_t_l1, 1e-8);
	EXPECT_NEAR(r.l2_norm, ref_t_l2, 1e-8);
	ASSERT_EQ(r.occupations.size(), static_cast<long>(ref_t_occ.size()));
	double moved = 0.0;
	for (size_t k = 0; k < ref_t_occ.size(); k++)
	{
		EXPECT_NEAR(r.occupations(k), ref_t_occ[k], 1e-8) << "natural orbital " << k;
		moved = std::max(moved, std::abs(r.occupations(k) - ref_occ[k]));
	}
	EXPECT_GT(moved, 1e-3);
	EXPECT_NEAR(r.occupations.sum(), 10.0, 1e-10);
}

//ORCA's CCSD(T) density is AutoCI's unrelaxed one (MDCI keeps the CCSD density), exported by orca_2json and read as
//natural orbitals of the job's gbw: PySCF's CCSD(T) occupations to ORCA's STol 1e-5, and occ's CCSD(T) density in space
TEST(CcsdDensity, OrcaAutociTriplesDensityMatchesOcc)
{
	if (!ensure_occ_data_path(nullptr)) GTEST_SKIP() << "OCC_DATA_PATH not set";
	const std::filesystem::path json = nos_test_repo_root() / "tests" / "ccsd_water" / "water_ccsdt_autoci.json";
	if (!std::filesystem::exists(json)) GTEST_SKIP() << "fixture missing: " << json.string();
	WFN orca;
	std::ostringstream log;
	ASSERT_TRUE(orca.read_orca_json_density(json, log));
	EXPECT_NE(log.str().find("ORCA density autocipur"), std::string::npos) << log.str();
	ASSERT_EQ(orca.get_nmo(), static_cast<int>(ref_t_occ.size()));
	std::vector<double> occ;
	for (int i = 0; i < orca.get_nmo(); i++) occ.push_back(orca.get_MO_occ(i));
	std::sort(occ.rbegin(), occ.rend());
	for (size_t k = 0; k < ref_t_occ.size(); k++) EXPECT_NEAR(occ[k], ref_t_occ[k], 2e-6) << "natural orbital " << k;
	CcsdDensityOptions opts;
	opts.n_frozen = 0;
	opts.tol = 1e-9;
	opts.triples = true;
	const WFN occ_wfn(ccsd_natural_orbitals(tight_rhf(), opts, log).natural), occ_ccsd_wfn(occ_ccsd().natural);
	double ccsd_dev = 0.0;
	for (unsigned a = 0; a < 3; a++)
		for (const d3 &shift : { d3{ 0.02, 0.0, 0.03 }, d3{ 0.6, 0.0, 0.3 }, d3{ 0.0, 0.4, -0.4 } })
		{
			const d3 p{ orca.get_atom_coordinate(a, 0) + shift[0], orca.get_atom_coordinate(a, 1) + shift[1], orca.get_atom_coordinate(a, 2) + shift[2] };
			const double rho = orca.compute_dens(p);
			EXPECT_NEAR(occ_wfn.compute_dens(p) / rho, 1.0, 1e-5) << "atom " << a;
			ccsd_dev = std::max(ccsd_dev, std::abs(occ_ccsd_wfn.compute_dens(p) / rho - 1.0));
		}
	EXPECT_GT(ccsd_dev, 1e-4);
}
