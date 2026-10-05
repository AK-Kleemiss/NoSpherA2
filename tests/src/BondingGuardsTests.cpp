#include "pch.h"
#include "core/wfn_class.h"
#include "core/constants.h"
#include "core/nos_math.h"

#include <cmath>
#include <filesystem>
#include <sstream>

//Inputs the bonding analyses must refuse rather than turn into a number: a shell type build_DM does not
//know, an index list longer than its matrix, and a wavefunction with fewer electrons than its nuclei carry.

namespace {

	double norm_s(const double a) { return std::pow(2.0 * a / constants::PI, 0.75); }

	//an s shell and a g shell on one centre, in the layout the wfn reader leaves behind
	WFN make_s_and_g()
	{
		WFN w(e_origin::wfn);
		w.push_back_atom("He", 0.0, 0.0, 0.0, 2);
		w.push_back_atom_basis_set(0, 1.0, 1.0, 1, 0);
		w.push_back_atom_basis_set(0, 0.7, 1.0, 5, 1);
		w.push_back_MO(1, 2.0, -0.9);
		double s[1] = { norm_s(1.0) };
		w.add_primitive(1, 1, 1.0, s);
		double zero[1] = { 0.0 };
		for (int type = 21; type <= 35; type++) //the fifteen cartesian g components
			w.add_primitive(1, type, 0.7, zero);
		w.set_exp_cutoff();
		return w;
	}
}

//build_DM's general-l primitive normalisation must reproduce the s/p/d/f constants exactly, or every
//number build_DM feeds moves.
TEST(BondingGuardsTests, AxialPrimNormReproducesTheHardcodedSpdfConstants)
{
	for (const double a : { 0.3, 1.0, 7.5, 120.0 }) {
		EXPECT_NEAR(constants::axial_prim_norm(0, a),
			std::pow(8 * std::pow(a, 3) / constants::PI3, 0.25), 1e-13 * std::abs(constants::axial_prim_norm(0, a)));
		EXPECT_NEAR(constants::axial_prim_norm(1, a),
			std::pow(128 * std::pow(a, 5) / constants::PI3, 0.25), 1e-13 * std::abs(constants::axial_prim_norm(1, a)));
		EXPECT_NEAR(constants::axial_prim_norm(2, a),
			std::pow(2048 * std::pow(a, 7) / (9 * constants::PI3), 0.25), 1e-13 * std::abs(constants::axial_prim_norm(2, a)));
		EXPECT_NEAR(constants::axial_prim_norm(3, a),
			std::pow(32768 * std::pow(a, 9) / (225 * constants::PI3), 0.25), 1e-13 * std::abs(constants::axial_prim_norm(3, a)));
	}
	//and it keeps going where the switch stopped: x^l exp(-a r^2) integrates to one
	for (int l = 0; l <= 6; l++) {
		const double a = 1.3, n = constants::axial_prim_norm(l, a);
		//<x^l e^{-ar^2} | x^l e^{-ar^2}> = (2l-1)!!/(4a)^l (PI/2a)^{3/2}
		const double self = constants::double_ft[std::max(2 * l - 1, 0)] / std::pow(4 * a, l) *
			std::pow(constants::PI / (2 * a), 1.5);
		EXPECT_NEAR(n * n * self, 1.0, 1e-12) << "l = " << l;
	}
}

//build_DM's shell walk knows s, p, d and f only; a g shell would leave `factor` at the previous shell's
//value and push no constants, shifting every later shell's norm_const by one. It refuses until the g
//component order in this file is settled.
TEST(BondingGuardsTests, BuildDmRefusesAGShellInsteadOfSlippingAShell)
{
	WFN w = make_s_and_g();
	EXPECT_FALSE(w.build_DM("unused-because-basis-is-loaded", false));
	EXPECT_EQ(w.get_DM_size(), 0);
}

//Every caller builds its own index list from a shell description, so an index past the matrix is named
//instead of read past its end.
TEST(BondingGuardsDeathTest, SubmatrixIndexPastTheMatrixIsNamed)
{
	dMatrix2 full(3, 3);
	std::fill(full.container().begin(), full.container().end(), 1.0);
	vec sub(4, 0.0);
	EXPECT_EXIT(get_submatrix(full, sub, ivec{ 0, 3 }),
		::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	vec one(1, 0.0);
	EXPECT_EXIT(get_submatrix(full, one, ivec{ -1 }),
		::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//No file format states an ECP unless -ECP is passed, so ECP atoms keep their full nuclear charge and
//anything filling orbitals from Z (RGBI's free-atom SCF) works with electrons the basis never described.
//Au2Br2: two def2 gold ECPs, 120 electrons; bromine is all-electron in def2, so only gold falls short.
TEST(BondingGuardsTests, AnExactDef2CoreShortfallIsDeclaredAsAnEcp)
{
	const auto p = nos_test_repo_root() / "tests" / "ECP_SF" / "Au2Br2.gbw";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/ECP_SF/Au2Br2.gbw not found";
	std::ostringstream log;
	WFN w(e_origin::NOT_YET_DEFINED);
	w.read_known_wavefunction_format(p, log, false);
	EXPECT_TRUE(w.get_has_ECPs());
	EXPECT_EQ(w.get_nr_ECP_electrons(), 120u);
}

//and a wavefunction that describes all of its electrons is left alone: the detection acts on an
//exact match against the def2 table, never on a shortfall it cannot name.
TEST(BondingGuardsTests, AnAllElectronWavefunctionIsNotGivenEcps)
{
	const auto p = nos_test_repo_root() / "tests" / "TFVC" / "water.gbw";
	if (!std::filesystem::exists(p))
		GTEST_SKIP() << "tests/TFVC/water.gbw not found";
	std::ostringstream log;
	WFN w(e_origin::NOT_YET_DEFINED);
	w.read_known_wavefunction_format(p, log, false);
	EXPECT_FALSE(w.get_has_ECPs());
	EXPECT_EQ(w.get_nr_ECP_electrons(), 0u);
}
