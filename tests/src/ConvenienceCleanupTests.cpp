//Differential tests for the convenience/NoSpherA2 cleanup: master's code, copied verbatim, against the production one,
//compared bit for bit
#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"

#include <bit>
#include <random>

namespace
{
	//master d624d8c4, primitive::primitive(const SimplePrimitive &) in convenience.cpp
	struct primitive_master_ref
	{
		double exp_l_plus_3_2, norm_const, normalized_coefficient;
		primitive_master_ref(const SimplePrimitive &other)
		{
			const int type = other.type;
			const double exp = other.exp, coefficient = other.coefficient;
			exp_l_plus_3_2 = pow(other.exp, type + 1.5);
			norm_const = pow(
				pow(2, 7 + 4 * type) * pow(exp, 3 + 2 * type) / constants::PI / pow(doublefactorial(2 * type + 1), 2),
				0.25);
			normalized_coefficient = coefficient * norm_const;
		}
	};

	uint64_t bits(const double x) { return std::bit_cast<uint64_t>(x); }

	void expect_same(const SimplePrimitive &s)
	{
		const primitive_master_ref ref(s);
		const primitive p(s), q(s.center, s.type, s.exp, s.coefficient);
		for (const primitive *x : { &p, &q })
		{
			EXPECT_EQ(x->get_center(), s.center);
			EXPECT_EQ(x->get_type(), s.type);
			EXPECT_EQ(bits(x->get_exp()), bits(s.exp)) << s.type << " " << s.exp;
			EXPECT_EQ(bits(x->get_coef()), bits(s.coefficient)) << s.type << " " << s.exp;
			EXPECT_EQ(bits(x->get_exp_l_plus_3_2()), bits(ref.exp_l_plus_3_2)) << s.type << " " << s.exp;
			EXPECT_EQ(bits(x->normalization_constant()), bits(ref.norm_const)) << s.type << " " << s.exp;
			EXPECT_EQ(bits(x->get_normalized_coefficient()), bits(ref.normalized_coefficient)) << s.type << " " << s.exp;
		}
	}
}

TEST(ConvenienceCleanup, PrimitiveConstructorsMatchMaster)
{
	const double exps[] = { 0.0, 1e-300, 1e-12, 0.035, 0.5, 1.0, 2.0, 13.7, 5.0e3, 1.2e7, 1e150, -0.0, -0.7 };
	const double coefs[] = { 0.0, -0.0, 1.0, -0.4, 3.3e-9, 7.1e8 };
	for (int t = 0; t <= 7; t++)
		for (const double e : exps)
			for (const double c : coefs)
				expect_same({ t % 3, t, e, c, 0 });
	std::mt19937_64 rng(20261120);
	std::uniform_real_distribution<double> le(-4.0, 6.0), co(-2.0, 2.0);
	std::uniform_int_distribution<int> ty(0, 6);
	for (int i = 0; i < 20000; i++)
		expect_same({ i % 5, ty(rng), std::pow(10.0, le(rng)), co(rng), 0 });
}
