#pragma once
#include "../pch.h"

struct SimplePrimitive {
	int center;     // Center index
	int type;       // Primitive type (0=s, 1=p, 2=d, etc.)
	double exp;     // Exponent
	double coefficient; // Coefficient
	int shell;
};

class primitive
{
private:
	int center, type;
	double exp, coefficient;
	double norm_const = -10;
	double exp_l_plus_3_2 = -10;
	double normalized_coefficient = -10;

public:
	void normalize()
	{
		coefficient *= normalization_constant();
	};
	void unnormalize()
	{
		coefficient /= normalization_constant();
	};
	double normalization_constant() const
	{
		// assuming type is equal to angular momentum
		return norm_const;
	}
	primitive() : center(0), type(0), exp(0.0), coefficient(0.0) {};
	primitive(int c, int t, double e, double coef);
	primitive(const SimplePrimitive& p);
	bool operator==(const primitive& other) const
	{
		return center == other.center &&
			type == other.type &&
			exp == other.exp &&
			coefficient == other.coefficient &&
			exp_l_plus_3_2 == other.exp_l_plus_3_2;
	};
	int get_center() const
	{
		return center;
	};
	int get_type() const
	{
		return type;
	};
	double get_exp_l_plus_3_2() const
	{
		return exp_l_plus_3_2;
	};
	double get_normalized_coefficient() const
	{
		return normalized_coefficient;
	};
	double get_exp() const
	{
		return exp;
	};
	double get_coef() const
	{
		return coefficient;
	};
	void set_center(const int& c)
	{
		center = c;
	};
	void set_type(const int& t)
	{
		type = t;
	};
	void set_exp(const double& e)
	{
		exp = e;
	};
	void set_coef(const double& c)
	{
		coefficient = c;
	};
	void set_norm_const(const double& nc)
	{
		norm_const = nc;
		normalized_coefficient = nc * coefficient;
	};
	double eval_gaussian(const double& r) const
	{
		return pow(r, type) * std::exp(-exp * r * r) * normalized_coefficient;
	};
	inline double eval_gaussian_unnormalized(const double& r2) const
	{
		return std::exp(-exp * r2) * coefficient;
	};
};

struct ECP_primitive : primitive
{
	int n;
	ECP_primitive() : primitive(), n(0) {}
	ECP_primitive(int c, int t, double e, double coef, int n) : primitive(c, t, e, coef), n(n) {}
};
