//CCSD densities from both sources NoSpherA2 reads: occ's CCSD with the Lambda-equation density (ccsd_density.cpp)
//and ORCA's unrelaxed CCSD natural orbitals (.mdci.nat). Water cc-pVDZ, all electrons, the geometry of
//tests/ccsd_water/water_ccsd.inp. References: PySCF RHF 1e-12, CCSD 1e-11, ccsd_lambda 1e-10, make_rdm1.
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

TEST(CcsdDensity, TriplesHaveNoDensityAndAreRefused)
{
	occ::io::OccInput config;
	config.method.name = "ccsd(t)";
	std::ostringstream log;
	EXPECT_EXIT((void)wavefunction_from_occ_input("unused.toml", config, log),
		::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}
