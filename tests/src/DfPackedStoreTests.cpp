#include "pch.h"
#include <gtest/gtest.h>

#include "core/basis_set.h"
#include <occ/io/xyz.h>
#include <occ/qm/hf.h>
#include <occ/qm/integral_engine_df.h>
#include <occ/qm/opmatrix.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>

//occ's packed three-centre DF store (occ/qm/df_packed_store.h) against the integral-direct builds on the same orbitals:
//stored whole, stored in part under a budget (J semi-direct, K direct) and not at all
namespace
{
	using Policy = occ::qm::IntegralEngineDF::Policy;
	//Two waters 14 A apart: the shell-pair screen drops the pairs across the gap, so the store is really packed
	const char* two_waters =
		"6\n\nO 0.000 0.000 0.117\nH 0.000 0.757 -0.469\nH 0.000 -0.757 -0.469\n"
		"O 14.000 0.000 0.117\nH 14.000 0.757 -0.469\nH 14.000 -0.757 -0.469\n";
	const char* methyl = "4\n\nC 0.000 0.000 0.000\nH 1.079 0.000 0.000\nH -0.540 0.935 0.000\nH -0.540 -0.935 0.000\n";

	struct df_case {
		occ::core::Molecule mol;
		occ::qm::AOBasis basis, aux;
		occ::qm::MolecularOrbitals mo;
	};

	df_case setup(const char* xyz, const occ::qm::SpinorbitalKind kind, const int mult = 1) {
		spdlog::set_level(spdlog::level::err);
		df_case c;
		c.mol = occ::io::molecule_from_xyz_string(xyz);
		c.mol.set_multiplicity(mult);
		c.basis = BasisSetLibrary::get_basis_set("def2-svp")->to_AOBasis(c.mol.atoms());
		c.aux = BasisSetLibrary::get_basis_set("def2-universal-jkfit")->to_AOBasis(c.mol.atoms());
		occ::qm::HartreeFock hf(c.basis);
		occ::qm::SCF<occ::qm::HartreeFock> scf(hf, kind);
		scf.set_charge_multiplicity(0, mult);
		scf.compute_initial_guess();
		c.mo = scf.ctx.mo;
		return c;
	}

	std::unique_ptr<occ::qm::IntegralEngineDF> engine(const df_case& c, const Policy p, const size_t limit = 0) {
		auto df = std::make_unique<occ::qm::IntegralEngineDF>(c.mol.atoms(), c.basis.shells(), c.aux.shells());
		df->set_integral_policy(p);
		if (limit) df->set_integral_store_memory_limit(limit);
		return df;
	}

	double max_abs(const occ::Mat& m) { return m.cwiseAbs().maxCoeff(); }

	void expect_same_jk(const df_case& c, occ::qm::IntegralEngineDF& df, const char* what) {
		auto direct = engine(c, Policy::Direct);
		const occ::Mat Jd = direct->coulomb(c.mo), Kd = direct->exchange(c.mo);
		const occ::Mat J = df.coulomb(c.mo), K = df.exchange(c.mo);
		ASSERT_EQ(J.rows(), Jd.rows());
		ASSERT_EQ(K.cols(), Kd.cols());
		std::cout << what << ": |J| " << max_abs(Jd) << " max |dJ| " << max_abs(J - Jd) << ", |K| " << max_abs(Kd) << " max |dK| " << max_abs(K - Kd) << std::endl;
		EXPECT_LT(max_abs(J - Jd), 1e-10);
		EXPECT_LT(max_abs(K - Kd), 1e-10);
		//a second build reuses the store (and, after the first exchange, its folded form)
		EXPECT_LT(max_abs(df.coulomb(c.mo) - Jd), 1e-10);
	}

	size_t full_store_bytes(occ::qm::IntegralEngineDF& df) {
		return occ::qm::PackedDFStore::bytes_for(df.ao_engine().aobasis(), df.ao_engine().shellpairs(), df.aux_engine().nbf());
	}
}

TEST(DfPackedStoreTests, RestrictedStoreMatchesDirectAndIsPacked)
{
	const df_case c = setup(two_waters, occ::qm::SpinorbitalKind::Restricted);
	auto df = engine(c, Policy::Stored);
	expect_same_jk(c, *df, "restricted");
	const auto& st = df->packed_store();
	ASSERT_FALSE(st.empty());
	EXPECT_TRUE(st.complete());
	EXPECT_TRUE(st.folded());
	const size_t nbf = c.basis.nbf(), dense = nbf * nbf * c.aux.nbf() * sizeof(double);
	EXPECT_LT(st.bytes(), dense / 2);
}

TEST(DfPackedStoreTests, UnrestrictedStoreMatchesDirect)
{
	const df_case c = setup(methyl, occ::qm::SpinorbitalKind::Unrestricted, 2);
	auto df = engine(c, Policy::Stored);
	expect_same_jk(c, *df, "unrestricted");
}

TEST(DfPackedStoreTests, GeneralSpinStoreMatchesDirect)
{
	df_case c = setup(two_waters, occ::qm::SpinorbitalKind::Restricted);
	//One spin-mixed spinor per spatial orbital, so the alpha-beta blocks of K are not zero
	const occ::Mat& C = c.mo.Cocc;
	const int nbf = C.rows(), nocc = C.cols();
	occ::Mat Cg = occ::Mat::Zero(2 * nbf, nocc);
	for (int i = 0; i < nocc; i++) {
		const double t = 0.3 + 0.1 * i;
		Cg.col(i).head(nbf) = std::cos(t) * C.col(i);
		Cg.col(i).tail(nbf) = std::sin(t) * C.col(i);
	}
	c.mo.kind = occ::qm::SpinorbitalKind::General;
	c.mo.Cocc = Cg;
	c.mo.D = Cg * Cg.transpose();
	auto df = engine(c, Policy::Stored);
	expect_same_jk(c, *df, "general");
	//K is Hermitian, and its alpha-beta block is there: the direct kernel once added it into beta-beta instead
	const occ::Mat K = df->exchange(c.mo);
	EXPECT_GT(max_abs(occ::qm::block::ab(K)), 1e-6);
	EXPECT_LT(max_abs(K - K.transpose()), 1e-12);
	EXPECT_LT(max_abs(engine(c, Policy::Direct)->exchange(c.mo) - K), 1e-10);
}

TEST(DfPackedStoreTests, ChooseUnderHalfTheBudgetStoresPartAndMatchesDirect)
{
	const df_case c = setup(two_waters, occ::qm::SpinorbitalKind::Restricted);
	auto probe = engine(c, Policy::Direct);
	const size_t full = full_store_bytes(*probe);
	auto df = engine(c, Policy::Choose, full / 2);
	expect_same_jk(c, *df, "half budget");
	const auto& st = df->packed_store();
	EXPECT_FALSE(st.empty());
	EXPECT_FALSE(st.complete());
	EXPECT_FALSE(st.folded());
	EXPECT_LE(st.bytes(), full / 2);
}

TEST(DfPackedStoreTests, ChooseUnderATinyBudgetStoresNothing)
{
	const df_case c = setup(methyl, occ::qm::SpinorbitalKind::Unrestricted, 2);
	auto df = engine(c, Policy::Choose, 1);
	expect_same_jk(c, *df, "tiny budget");
	EXPECT_TRUE(df->packed_store().empty());
}
