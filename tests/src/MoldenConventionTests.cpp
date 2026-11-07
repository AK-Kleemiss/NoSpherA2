#include "pch.h"

#include "core/wfn_class.h"

#include <filesystem>

// Molden [GTO] coefficients are read as bare (multiplying exp(-a r^2), contraction already normalised,
// as orca_2mkl writes), not normprim (normalised primitives, as DGrid reads). Only the bare reading
// normalises F2.molden's MOs; references come from an independent evaluator using it.
TEST(MoldenConvention, F2ValenceDensityMatchesAnIndependentEvaluation)
{
	const std::filesystem::path f = nos_test_repo_root() / "tests" / "molden_file" / "F2.molden";
	if (!std::filesystem::exists(f)) GTEST_SKIP() << "missing fixture " << f.string();
	WFN w(f);
	ASSERT_EQ(w.get_ncen(), 2);
	// Only s shells contribute at the nucleus; their near cancellation makes a wrong normalisation a factor of thousands.
	EXPECT_NEAR(w.compute_dens({0.0, 0.0, 0.0}), 0.021245, 1e-5) << "density at the F nucleus";
	EXPECT_NEAR(w.compute_dens({1.41729459346535, 0.0, 0.0}), 0.238070, 1e-5) << "bond midpoint";
	EXPECT_NEAR(w.compute_dens({0.5, 0.3, -0.2}), 0.810097, 1e-5) << "off-axis point";
}
