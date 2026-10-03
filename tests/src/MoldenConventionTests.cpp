#include "pch.h"

#include "core/wfn_class.h"

#include <filesystem>

// Which normalisation a molden file's [GTO] coefficients belong to is not written anywhere in
// the file, and getting it wrong changes the density by orders of magnitude rather than by a
// little. Two readings are in circulation:
//
//   bare      the tabulated coefficient multiplies exp(-a r^2) and the contracted function is
//             already normalised - what orca_2mkl writes and what this reader assumes
//   normprim  the coefficient multiplies an individually normalised primitive and the
//             contracted shell is then renormalised - the reading the format documentation
//             describes, and what DGrid applies
//
// The file itself settles it: only under the reading it was written in are its MO vectors
// normalised. For tests/molden_file/F2.molden (2 F atoms, one s and one p shell each, an
// ECP valence basis, 7 doubly occupied MOs) the bare reading gives every MO norm 1.00000 and
// 14.000000 electrons, the normprim reading 0.889..1.344 and 13.727. The reference values
// below come from an independent evaluator built on that arithmetic, not from this program:
// scratchpad molden_rho.py / f2_mo_norm.py of session 61956d9d.
//
// The reader passes this as written, so the 695x disagreement with DGrid 5.2 at the F nucleus is
// DGrid applying the normprim reading to a bare file, not a bug here. Neither side's reported
// number belongs to either reading: bare gives 0.021245 at the nucleus and normprim 102.045967,
// while the discrepancy was reported as 0.0897 against DGrid's 62.33.
TEST(MoldenConvention, F2ValenceDensityMatchesAnIndependentEvaluation)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "F2.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN w(f);
	ASSERT_EQ(w.get_ncen(), 2);
	// At a nucleus only the s shells contribute; the near cancellation of the four s
	// coefficients is what makes this value small, and it is where a wrong normalisation
	// shows up as a factor of thousands rather than as a small shift.
	EXPECT_NEAR(w.compute_dens({0.0, 0.0, 0.0}), 0.021245, 1e-5) << "density at the F nucleus";
	EXPECT_NEAR(w.compute_dens({1.41729459346535, 0.0, 0.0}), 0.238070, 1e-5) << "bond midpoint";
	EXPECT_NEAR(w.compute_dens({0.5, 0.3, -0.2}), 0.810097, 1e-5) << "off-axis point";
}
