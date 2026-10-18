#include "pch.h"

#include "core/spherical_density.h"
#include "core/constants.h"

// The delta_k series replaces the extrapolation rho_q = (1-|q|) rho_0 + |q| rho_ion for |q| > 1 by
// rho_q = rho_n - f * delta_{n+1}, rho_n = rho_1 - sum_{k=2..n} delta_k, between adjacent bound states.

namespace
{
	// The k the rest of NoSpherA2 hands to get_form_factor is 4 pi s a0 -- see
	// k_of_reflection in scattering_factors.cpp.
	double k_of_s(const double s)
	{
		return constants::bohr2ang(constants::FOUR_PI * s);
	}

	// Elements carrying delta_2 AND delta_3, so charge states up to 3 exist.
	const int kDeltaElements[] = {20 /*Ca*/, 24 /*Cr*/, 26 /*Fe*/,
								  29 /*Cu*/, 30 /*Zn*/, 35 /*Br*/};

	// f(k) -> electron count as k -> 0, but NOT at k = 0 itself: the Slater
	// form factor is singular exactly there (Fe(2.5) gives 113.8 instead of
	// 23.5). Harmless in practice, since k = 4 pi s a0 vanishes only for the
	// unmeasured 000 reflection, but it means the count has to be probed just
	// off zero. By k = 1e-8 the value has converged.
	const double kNearZero = 1e-8;

	// Tolerance on that probe. The Slater tables reproduce their own electron
	// count to ~1.5e-6 here, and a blend multiplies that by its weight (a
	// charge of 4.5 was measured 7e-6 out on the FALLBACK path, which never
	// touches a delta table). So this is baseline table precision, not the
	// delta series. It still sits four orders below the ~1.0 error that a
	// miscounted delta or an off-by-one in the chain would produce.
	const double kCountTol = 1e-4;

	// Numerically integrate 4 pi r^2 rho(r) on a log grid. Used to confirm the
	// r-space tables and the s-space tables describe the same atom.
	double integrate_density(const HE_Spherical_Atom &a)
	{
		const double lo = std::log(1e-8), hi = std::log(40.0);
		const int n = 400000;
		const double dl = (hi - lo) / (n - 1);
		double sum = 0.0;
		for (int i = 0; i < n; i++)
		{
			const double r = std::exp(lo + i * dl);
			const double w = (i == 0 || i == n - 1) ? 0.5 : 1.0;
			// dr = r dlnr
			sum += w * 4.0 * constants::PI * r * r * a.get_radial_density(r) * r * dl;
		}
		return sum;
	}
}

TEST(DeltaSeriesTests, PastPlusOneInterpolatesInsteadOfExtrapolating)
{
	for (const int Z : kDeltaElements)
	{
		for (const double q : {1.5, 2.0, 2.5, 3.0})
		{
			const HE_Spherical_Atom a(Z, q);
			EXPECT_TRUE(a.uses_delta_series())
				<< "Z=" << Z << " q=" << q;
			EXPECT_FALSE(a.is_extrapolating())
				<< "Z=" << Z << " q=" << q;
		}
	}
}

// Each delta_k holds exactly one electron, so f(k->0) = Z - q
TEST(DeltaSeriesTests, ElectronCountStaysExactPastPlusOne)
{
	for (const int Z : kDeltaElements)
	{
		for (const double q : {1.5, 2.0, 2.5, 3.0})
		{
			const HE_Spherical_Atom a(Z, q);
			EXPECT_NEAR(a.get_form_factor(kNearZero), Z - q, kCountTol)
				<< "Z=" << Z << " q=" << q;
		}
	}
}

// Relative to the local density: 1e-5 e/bohr^3 is nothing at r ~ 1 bohr but a sign error in the tail
TEST(DeltaSeriesTests, DISABLED_ProfileTheNegativeExcursion)
{
	for (const int Z : {20, 24, 26, 29, 30, 35})
	{
		const HE_Spherical_Atom a(Z, 2.5);
		const HE_Spherical_Atom neutral(Z, 0.0);
		double worst = 0.0, worst_r = 0.0, worst_rel = 0.0;
		for (double r = 1e-4; r < 30.0; r *= 1.02)
		{
			const double d = a.get_radial_density(r);
			if (d < worst)
			{
				worst = d;
				worst_r = r;
				const double scale = neutral.get_radial_density(r);
				worst_rel = (scale > 0.0) ? d / scale : 0.0;
			}
		}
		// What the OLD path would have done at the same charge:
		// rho = (1-q) rho_neutral + q rho_cation, with q = 2.5. HE(Z,1.0)
		// returns the tabulated cation and HE(Z,0.0) the neutral, so the old
		// blend can be reconstructed exactly for comparison.
		const HE_Spherical_Atom cation(Z, 1.0);
		double old_worst = 0.0, old_worst_r = 0.0;
		for (double r = 1e-4; r < 30.0; r *= 1.02)
		{
			const double d = -1.5 * neutral.get_radial_density(r)
							 + 2.5 * cation.get_radial_density(r);
			if (d < old_worst) { old_worst = d; old_worst_r = r; }
		}
		std::cout << "Z=" << Z << "  delta: min rho = " << worst
				  << " at r = " << worst_r
				  << " (ratio to neutral " << worst_rel << ")"
				  << "   |   extrapolated: min rho = " << old_worst
				  << " at r = " << old_worst_r << std::endl;
	}
}

TEST(DeltaSeriesTests, DISABLED_ProbeFormFactorNearZero)
{
	const HE_Spherical_Atom a(26, 2.5);
	for (const double k : {0.0, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1})
		std::cout << "  k=" << k << "  f=" << a.get_form_factor(k) << std::endl;
	std::cout << "  expected electron count = " << 26.0 - 2.5 << std::endl;
}

// The Slater cation base and the GTO differences do not cancel exactly, so a small dip survives
TEST(DeltaSeriesTests, ResidualNegativeExcursionStaysNegligible)
{
	for (const int Z : kDeltaElements)
	{
		const HE_Spherical_Atom a(Z, 2.5);
		EXPECT_GE(a.most_negative_density(), -1e-4) << "Z=" << Z;
	}
}

// 10x leaves margin, so the test is not brittle about the exact factor
TEST(DeltaSeriesTests, DeltaRouteBeatsTheExtrapolationItReplaces)
{
	for (const int Z : kDeltaElements)
	{
		const HE_Spherical_Atom a(Z, 2.5);
		const HE_Spherical_Atom neutral(Z, 0.0);
		const HE_Spherical_Atom cation(Z, 1.0);
		double old_worst = 0.0;
		for (double r = 1e-4; r < 30.0; r *= 1.02)
		{
			const double d = -1.5 * neutral.get_radial_density(r)
							 + 2.5 * cation.get_radial_density(r);
			old_worst = std::min(old_worst, d);
		}
		ASSERT_LT(old_worst, 0.0) << "Z=" << Z
			<< ": the old path did not dip, so there is nothing to improve on"
			   " and this test is no longer measuring what it claims";
		EXPECT_LT(std::abs(a.most_negative_density()), std::abs(old_worst) / 10.0)
			<< "Z=" << Z << " delta=" << a.most_negative_density()
			<< " extrapolated=" << old_worst;
	}
}

// An off-by-one log-grid index or r tables misaligned with the f tables break the count while each looks plausible
TEST(DeltaSeriesTests, IntegratedDensityMatchesTheElectronCount)
{
	for (const int Z : {20, 26, 29})
	{
		for (const double q : {1.5, 2.5})
		{
			const HE_Spherical_Atom a(Z, q);
			EXPECT_NEAR(integrate_density(a), Z - q, 5e-3)
				<< "Z=" << Z << " q=" << q;
		}
	}
}

// At exactly +1 the tabulated ion is returned; the delta route must start from it
TEST(DeltaSeriesTests, ContinuousAcrossPlusOne)
{
	for (const int Z : kDeltaElements)
	{
		const HE_Spherical_Atom at_one(Z, 1.0);
		const HE_Spherical_Atom just_past(Z, 1.0 + 1e-7);
		ASSERT_FALSE(at_one.uses_delta_series()) << "Z=" << Z;
		ASSERT_TRUE(just_past.uses_delta_series()) << "Z=" << Z;
		for (const double r : {0.05, 0.2, 0.5, 1.0, 2.0, 5.0})
		{
			EXPECT_NEAR(just_past.get_radial_density(r),
						at_one.get_radial_density(r), 1e-6)
				<< "Z=" << Z << " r=" << r;
		}
		for (const double s : {0.02, 0.25, 0.6, 1.2, 3.0})
		{
			EXPECT_NEAR(just_past.get_form_factor(k_of_s(s)),
						at_one.get_form_factor(k_of_s(s)), 1e-6)
				<< "Z=" << Z << " s=" << s;
		}
	}
}

TEST(DeltaSeriesTests, UpToPlusOneIsUnaffected)
{
	for (const int Z : kDeltaElements)
	{
		for (const double q : {0.25, 0.5, 1.0})
		{
			const HE_Spherical_Atom a(Z, q);
			EXPECT_FALSE(a.uses_delta_series()) << "Z=" << Z << " q=" << q;
			EXPECT_FALSE(a.is_extrapolating()) << "Z=" << Z << " q=" << q;
			EXPECT_NEAR(a.get_form_factor(kNearZero), Z - q, kCountTol)
				<< "Z=" << Z << " q=" << q;
		}
	}
}

// Carbon's delta_2 failed validation and rho_m sums from k=2, so delta_3 alone is unusable
TEST(DeltaSeriesTests, ElementWithoutAContiguousChainFallsBack)
{
	const HE_Spherical_Atom c(6, 2.5);
	EXPECT_FALSE(c.uses_delta_series());
	EXPECT_TRUE(c.is_extrapolating());
	// The fallback still conserves electrons exactly; only the shape suffers.
	EXPECT_NEAR(c.get_form_factor(kNearZero), 6.0 - 2.5, kCountTol);
}

// Anions past -1 are unbound, so there is no reference state to walk towards
TEST(DeltaSeriesTests, AnionsPastMinusOneDoNotUseTheSeries)
{
	for (const int Z : {8 /*O*/, 17 /*Cl*/})
	{
		const HE_Spherical_Atom a(Z, -2.0);
		EXPECT_FALSE(a.uses_delta_series()) << "Z=" << Z;
	}
}

TEST(DeltaSeriesTests, ChargeBeyondTheTabulatedRangeFallsBack)
{
	const HE_Spherical_Atom a(26, 4.5);   // Fe, tables stop at delta_3
	EXPECT_FALSE(a.uses_delta_series());
	EXPECT_TRUE(a.is_extrapolating());
	EXPECT_NEAR(a.get_form_factor(kNearZero), 26.0 - 4.5, kCountTol);
}
