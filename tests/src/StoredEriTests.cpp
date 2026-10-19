#include "pch.h"
#include <gtest/gtest.h>

#include "core/basis_set.h"
#include "core/stored_eri.h"
#if defined(NOSPHERA2_USE_GPU)
#include "core/itensor_gpu.h"
#endif
#include <occ/io/xyz.h>
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <occ/qm/spinorbital.h>
#include <spdlog/spdlog.h>

namespace
{
	const char* water = "3\n\nO 0.000 0.000 0.117\nH 0.000 0.757 -0.469\nH 0.000 -0.757 -0.469\n";
	//Two waters 14 A apart: the pairs across the gap fall below the Schwarz screen
	const char* two_waters =
		"6\n\nO 0.000 0.000 0.117\nH 0.000 0.757 -0.469\nH 0.000 -0.757 -0.469\n"
		"O 14.000 0.000 0.117\nH 14.000 0.757 -0.469\nH 14.000 -0.757 -0.469\n";
	const char* methyl = "4\n\nC 0.000 0.000 0.000\nH 1.079 0.000 0.000\nH -0.540 0.935 0.000\nH -0.540 -0.935 0.000\n";

	struct scf_case {
		occ::core::Molecule mol;
		occ::qm::AOBasis basis;
		std::unique_ptr<occ::qm::HartreeFock> hf;
		std::unique_ptr<occ::qm::SCF<occ::qm::HartreeFock>> scf;
		occ::qm::MolecularOrbitals guess, converged;
	};

	scf_case run(const char* xyz, const char* basis_name, const occ::qm::SpinorbitalKind kind, const int mult = 1) {
		spdlog::set_level(spdlog::level::err);
		scf_case c;
		c.mol = occ::io::molecule_from_xyz_string(xyz);
		c.mol.set_multiplicity(mult);
		c.basis = BasisSetLibrary::get_basis_set(basis_name)->to_AOBasis(c.mol.atoms());
		c.hf = std::make_unique<occ::qm::HartreeFock>(c.basis);
		c.scf = std::make_unique<occ::qm::SCF<occ::qm::HartreeFock>>(*c.hf, kind);
		c.scf->set_charge_multiplicity(0, mult);
		c.scf->compute_initial_guess();
		c.guess = c.scf->ctx.mo;
		c.scf->compute_scf_energy();
		c.converged = c.scf->wavefunction().mo;
		return c;
	}

	occ::Mat direct_fock(const scf_case& c, const occ::qm::MolecularOrbitals& mo) {
		return c.hf->compute_fock(mo, c.hf->compute_schwarz_ints());
	}

	double max_abs(const occ::Mat& m) { return m.cwiseAbs().maxCoeff(); }

	void expect_full_build(const scf_case& c, const stored_eri& eri, const double tol, const char* what) {
		const occ::Mat F_direct = direct_fock(c, c.converged), F_stored = eri.fock(c.converged);
		const double diff = max_abs(F_direct - F_stored);
		std::cout << what << ": nbf " << eri.nbf() << ", " << eri.npairs() << " of " << eri.nbf() * (eri.nbf() + 1) / 2
			<< " pairs kept, |F| " << max_abs(F_direct) << ", max |F_stored - F_direct| " << diff << std::endl;
		EXPECT_LT(diff, tol);
	}

	//G(D1) against G(D0) + G_screened(D1 - D0), the incremental step of the SCF loop
	void expect_incremental(const scf_case& c, const stored_eri& eri, const double tol, const char* what) {
		occ::qm::MolecularOrbitals diff = c.converged;
		diff.D = c.converged.D - c.guess.D;
		const occ::Mat F_full = eri.fock(c.converged), F_inc = eri.fock(c.guess) + eri.fock(diff, true);
		const double err = max_abs(F_full - F_inc);
		std::cout << what << " incremental: |D1 - D0| " << max_abs(diff.D) << ", max |F_inc - F_full| " << err << std::endl;
		EXPECT_LT(err, tol);
	}
}

TEST(StoredEriTests, DenseRestrictedMatchesTheDirectBuild)
{
	const scf_case c = run(water, "def2-SVP", occ::qm::SpinorbitalKind::Restricted);
	stored_eri eri;
	ASSERT_TRUE(eri.build(*c.hf, 0, std::cout));
	EXPECT_TRUE(eri.dense());
	expect_full_build(c, eri, 1e-11, "water def2-SVP RHF");
	expect_incremental(c, eri, 1e-9, "water def2-SVP RHF");
}

TEST(StoredEriTests, SparseRestrictedMatchesTheDirectBuild)
{
	const scf_case c = run(two_waters, "def2-SVP", occ::qm::SpinorbitalKind::Restricted);
	stored_eri eri;
	ASSERT_TRUE(eri.build(*c.hf, 0, std::cout));
	EXPECT_FALSE(eri.dense());
	expect_full_build(c, eri, 1e-10, "two waters def2-SVP RHF");
	expect_incremental(c, eri, 1e-9, "two waters def2-SVP RHF");
}

TEST(StoredEriTests, UnrestrictedMatchesTheDirectBuild)
{
	const scf_case c = run(methyl, "def2-SVP", occ::qm::SpinorbitalKind::Unrestricted, 2);
	stored_eri eri;
	ASSERT_TRUE(eri.build(*c.hf, 0, std::cout));
	expect_full_build(c, eri, 1e-11, "methyl def2-SVP UHF");
	expect_incremental(c, eri, 1e-9, "methyl def2-SVP UHF");
}

TEST(StoredEriTests, BudgetRefusesWhatDoesNotFit)
{
	const scf_case c = run(water, "sto-3g", occ::qm::SpinorbitalKind::Restricted);
	stored_eri eri;
	std::ostringstream log;
	EXPECT_FALSE(eri.build(*c.hf, 8, log));
	EXPECT_FALSE(static_cast<bool>(eri));
	EXPECT_NE(log.str().find("do not fit"), std::string::npos);
}

//save_eri/read_eri: the file gives back every array build() made, bit for bit, and a basis on another geometry is refused
TEST(StoredEriTests, SavedIntegralsReadBackExactlyAndOnlyForTheirBasis)
{
	const scf_case c = run(two_waters, "def2-SVP", occ::qm::SpinorbitalKind::Restricted);
	stored_eri built, read;
	ASSERT_TRUE(built.build(*c.hf, 0, std::cout));
	const auto path = std::filesystem::temp_directory_path() / "nos_stored_eri_test.bin";
	ASSERT_TRUE(built.save(path, *c.hf));
	std::ostringstream log;
	ASSERT_TRUE(read.load(path, *c.hf, 0, log)) << log.str();
	EXPECT_NE(log.str().find("read from"), std::string::npos) << log.str();
	ASSERT_EQ(read.npairs(), built.npairs());
	EXPECT_EQ(read.nbf(), built.nbf());
	EXPECT_EQ(read.pair_a(), built.pair_a());
	EXPECT_EQ(read.pair_b(), built.pair_b());
	EXPECT_EQ(read.first_pair(), built.first_pair());
	EXPECT_EQ(read.pair_index(), built.pair_index());
	EXPECT_EQ(read.schwarz(), built.schwarz());
	EXPECT_EQ(std::memcmp(read.data(), built.data(), built.bytes()), 0);
	const occ::Mat D = c.converged.D;
	occ::Mat Jb, Kb, Jr, Kr;
	built.JK(D, Jb, Kb);
	read.JK(D, Jr, Kr);
	EXPECT_EQ((Jb - Jr).cwiseAbs().maxCoeff(), 0.0);
	EXPECT_EQ((Kb - Kr).cwiseAbs().maxCoeff(), 0.0);
	//water alone: other atoms, other pairs
	const scf_case w = run(water, "def2-SVP", occ::qm::SpinorbitalKind::Restricted);
	std::ostringstream refused;
	EXPECT_FALSE(read.load(path, *w.hf, 0, refused));
	EXPECT_FALSE(static_cast<bool>(read));
	EXPECT_NE(refused.str().find("holds no two-electron integrals of this basis"), std::string::npos) << refused.str();
	//a file cut short is refused, not half read
	std::filesystem::resize_file(path, std::filesystem::file_size(path) - 8);
	EXPECT_FALSE(read.load(path, *c.hf, 0, refused));
	std::filesystem::remove(path);
}

#if defined(NOSPHERA2_USE_GPU)
TEST(StoredEriTests, DeviceContractionMatchesTheCpu)
{
	if (!itensor_gpu_available()) GTEST_SKIP() << "no device";
	for (const char* xyz : { water, two_waters }) {
		const scf_case c = run(xyz, "def2-SVP", occ::qm::SpinorbitalKind::Restricted);
		stored_eri eri;
		ASSERT_TRUE(eri.build(*c.hf, 0, std::cout));
		ASSERT_TRUE(eri_gpu_hold(eri.data(), eri.nbf(), eri.npairs(), eri.pair_a().data(), eri.pair_b().data(), eri.first_pair().data(), eri.pair_index().data()));
		occ::Mat J, K, Jd, Kd;
		eri.JK(c.converged.D, J, K);
		Jd.resize(eri.nbf(), eri.nbf());
		Kd.resize(eri.nbf(), eri.nbf());
		ASSERT_TRUE(eri_gpu_JK(c.converged.D.data(), Jd.data(), Kd.data()));
		eri_gpu_release();
		const double dj = max_abs(J - Jd), dk = max_abs(K - Kd);
		std::cout << (eri.dense() ? "dense" : "sparse") << " device: max |J - J_dev| " << dj << ", max |K - K_dev| " << dk << std::endl;
		EXPECT_LT(dj, 1e-11);
		EXPECT_LT(dk, 1e-11);
	}
}
#endif
