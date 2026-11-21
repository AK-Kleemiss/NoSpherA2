//Tests for the convenience/NoSpherA2 cleanup. PrimitiveConstructorsMatchMaster is differential: master's code, copied
//verbatim, against the production one, compared bit for bit. CellBoxFromCifHasTheCellsAxes is a geometry check of the
//moved, unchanged readxyzMinMax_fromCIF with tolerances, not a differential test.
#include "pch.h"

#include "core/base/convenience.h"
#include "core/base/constants.h"

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

//The step matrix of a cell cube: column j times NbSteps[j] is cell axis j in bohr, so the columns have the cell's
//lengths and the angles between them are the cell's angles; NbSteps covers each axis at the resolution
TEST(ConvenienceCleanup, CellBoxFromCifHasTheCellsAxes)
{
	const struct { const char *cif; double a, b, c, al, be, ga; } cases[] = {
		{ "tests/P1_test/P1_test.cif", 4.8819, 6.4182, 8.9900, 74.040, 83.655, 80.497 },
		{ "tests/epoxide_gbw/epoxide.cif", 4.633, 8.400, 6.577, 90.0, 100.37, 90.0 } };
	for (const auto &t : cases)
	{
		properties_options opts;
		opts.resolution = 0.13;
		vec2 cm(3, vec(3, 0.0));
		readxyzMinMax_fromCIF(nos_test_repo_root() / t.cif, opts, cm);
		const double len[3] = { t.a, t.b, t.c }, ang[3] = { t.al, t.be, t.ga };
		d3 col[3];
		for (int j = 0; j < 3; j++)
		{
			EXPECT_EQ(opts.NbSteps[j], (int)std::ceil(len[j] / 0.13)) << t.cif;
			for (int i = 0; i < 3; i++)
				col[j][i] = cm[i][j] * opts.NbSteps[j];
			EXPECT_NEAR(array_length(col[j]), constants::ang2bohr(len[j]), 1e-10) << t.cif << " axis " << j;
		}
		EXPECT_EQ(cm[1][0], 0.0);
		EXPECT_EQ(cm[2][0], 0.0);
		EXPECT_EQ(cm[2][1], 0.0);
		for (int k = 0; k < 3; k++)
		{
			const d3 &u = col[(k + 1) % 3], &v = col[(k + 2) % 3];
			EXPECT_NEAR(vec_dot(u, v) / (array_length(u) * array_length(v)), std::cos(ang[k] * constants::PI / 180.0), 1e-12) << t.cif << " angle " << k;
		}
		for (int i = 0; i < 3; i++)
		{
			EXPECT_EQ(opts.MinMax[i], 0.0);
			//x and y are converted with 0.529177249, z with constants::a0: equal to 1e-7 relative
			EXPECT_NEAR(opts.MinMax[3 + i], col[0][i] + col[1][i] + col[2][i], 1e-5) << t.cif << " corner " << i;
		}
	}
}
