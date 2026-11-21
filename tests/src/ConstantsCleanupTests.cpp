#include "pch.h"

#include "core/base/convenience.h"
#include "core/base/constants.h"

#include <array>
#include <bit>
#include <cstring>
#include <random>

namespace
{
	//master's factorial generators (constants.h at d624d8c4), verbatim
	constexpr std::array<size_t, constants::MAX_FACTORIAL> factorial_const_master_ref() {
		std::array<size_t, constants::MAX_FACTORIAL> ft{};

		ft[0] = 1;
		for (int i = 1; i < constants::MAX_FACTORIAL; i++) {
			ft[i] = ft[i - 1] * i;
		}
		return ft;
	}
	constexpr std::array<double, constants::MAX_FACTORIAL> factorial_double_master_ref() {
		std::array<double, constants::MAX_FACTORIAL> f{};
		f[0] = 1.0;
		for (int i = 1; i < constants::MAX_FACTORIAL; i++) f[i] = f[i - 1] * i;
		return f;
	}
}

TEST(ConstantsCleanupTests, FactorialTablesAreMastersBitForBit)
{
	const auto ft = factorial_const_master_ref();
	const auto ftd = factorial_double_master_ref();
	ASSERT_EQ(sizeof(ft), sizeof(constants::ft));
	ASSERT_EQ(sizeof(ftd), sizeof(constants::ftd));
	EXPECT_EQ(std::memcmp(ft.data(), constants::ft.data(), sizeof(ft)), 0);
	EXPECT_EQ(std::memcmp(ftd.data(), constants::ftd.data(), sizeof(ftd)), 0);
}

TEST(ConstantsCleanupTests, DensityAccuracyIsMastersValue)
{
	const double acc = constants::density_accuracy;
	EXPECT_EQ(std::bit_cast<uint64_t>(acc), std::bit_cast<uint64_t>(5.0e-5));
}

//sum_m Y_lm^2 = (2l+1)/(4 pi) bounds every |Y_lm| on the unit sphere, so the stored maxima may not undercut a sampled value
TEST(ConstantsCleanupTests, SphericalHarmonicMaxBoundsEveryValue)
{
	std::mt19937_64 rng(20261120);
	std::normal_distribution<double> g(0.0, 1.0);
	for (int l = 0; l <= 9; l++)
		for (int m = -l; m <= l; m++) {
			const double mx = constants::spherical_harmonic_max(l, m);
			EXPECT_GT(mx, 0.0) << l << " " << m;
			EXPECT_LE(mx, sqrt((2 * l + 1) * constants::INV_FOUR_PI) * (1 + 1e-12)) << l << " " << m;
			double seen = 0.0;
			for (int s = 0; s < 4000; s++) {
				double d[5] = { g(rng), g(rng), g(rng), 0.0, 0.0 };
				const double r = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
				for (int i = 0; i < 3; i++) d[i] /= r;
				seen = std::max(seen, std::abs(constants::spherical_harmonic(l, m, d)));
			}
			EXPECT_LE(seen, mx * (1 + 1e-12)) << l << " " << m;
		}
}
