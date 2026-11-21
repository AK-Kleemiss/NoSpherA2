#include "pch.h"

#include "core/nbo.h"
#include "core/nao_eigen.h"
#include "core/constants.h"
#include "core/wfn_class.h"

#include <bit>
#include <cstring>
#include <random>

//to_eigen / to_dmatrix of nao_eigen.h against frozen copies of the anonymous-namespace helpers they replace
//(master d624d8c4: nbo.cpp, nrt.cpp's to_eigen is the same text, and nao.cpp), compared bit for bit.

namespace {
	using Eigen::MatrixXd;

	MatrixXd to_eigen_nbo_master_ref(const dMatrix2& m)
	{
		const int r = static_cast<int>(m.extent(0)), c = static_cast<int>(m.extent(1));
		MatrixXd out(r, c);
		for (int i = 0; i < r; i++)
			for (int j = 0; j < c; j++)
				out(i, j) = m(i, j);
		return out;
	}

	dMatrix2 to_dmatrix_nbo_master_ref(const MatrixXd& m)
	{
		dMatrix2 out(m.rows(), m.cols());
		for (int i = 0; i < m.rows(); i++)
			for (int j = 0; j < m.cols(); j++)
				out(i, j) = m(i, j);
		return out;
	}

	MatrixXd to_eigen_nao_master_ref(const dMatrix2 &m)
	{
		const int n = static_cast<int>(m.extent(0)), c = static_cast<int>(m.extent(1));
		MatrixXd out(n, c);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < c; j++)
				out(i, j) = m(i, j);
		return out;
	}

	dMatrix2 to_dmatrix_nao_master_ref(const MatrixXd &m)
	{
		dMatrix2 out(m.rows(), m.cols());
		for (int i = 0; i < m.rows(); i++)
			for (int j = 0; j < m.cols(); j++)
				out(i, j) = m(i, j);
		return out;
	}

	bool same_bits(const MatrixXd& a, const MatrixXd& b)
	{
		return a.rows() == b.rows() && a.cols() == b.cols() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
	}

	bool same_bits(const dMatrix2& a, const dMatrix2& b)
	{
		return a.extent(0) == b.extent(0) && a.extent(1) == b.extent(1) && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
	}

	const double special[] = { 0.0, -0.0, 1.0, -1.0, 5e-324, -5e-324, 2.2250738585072014e-308, 1e-300, 1.7976931348623157e308,
							   -1.7976931348623157e308, std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
							   std::numeric_limits<double>::quiet_NaN(), 0.1, -3.5e17 };
}

TEST(NboCleanupTests, DmatrixEigenCopiesMatchTheMasterHelpers)
{
	std::mt19937_64 rng(20261120);
	const int ns = static_cast<int>(sizeof(special) / sizeof(special[0]));
	for (int pass = 0; pass < 400; pass++) {
		const int r = pass < 64 ? pass / 8 : static_cast<int>(rng() % 40), c = pass < 64 ? pass % 8 : static_cast<int>(rng() % 40);
		dMatrix2 d(r, c);
		MatrixXd e(r, c);
		for (int i = 0; i < r; i++)
			for (int j = 0; j < c; j++) {
				const uint64_t u = rng();
				d(i, j) = (u & 3) == 0 ? special[(u >> 2) % ns] : (u & 3) == 1 ? std::bit_cast<double>(rng()) : std::ldexp(static_cast<double>(u >> 11), -40) - 4096.0;
				e(i, j) = std::bit_cast<double>(rng());
			}
		const MatrixXd te = to_eigen(d);
		EXPECT_TRUE(same_bits(te, to_eigen_nbo_master_ref(d))) << r << "x" << c;
		EXPECT_TRUE(same_bits(te, to_eigen_nao_master_ref(d))) << r << "x" << c;
		const dMatrix2 td = to_dmatrix(e);
		EXPECT_TRUE(same_bits(td, to_dmatrix_nbo_master_ref(e))) << r << "x" << c;
		EXPECT_TRUE(same_bits(td, to_dmatrix_nao_master_ref(e))) << r << "x" << c;
		EXPECT_TRUE(same_bits(to_dmatrix(te), d)) << r << "x" << c;
	}
}

namespace {
	std::vector<atom> atoms_at(const std::vector<std::pair<int, std::array<double, 3>>>& at)
	{
		WFN w;
		for (const auto& a : at)
			w.push_back_atom(constants::atnr2letter(a.first), constants::ang2bohr(a.second[0]), constants::ang2bohr(a.second[1]),
							 constants::ang2bohr(a.second[2]), a.first);
		return w.get_atoms();
	}
}

TEST(NboCleanupTests, BondablePairsNeverBondThroughAnAtomOnTheLine)
{
	//Ni-C 1.80 A, C-O 1.15 A: at the NRT scale 1.75 the Ni-O distance is inside the cut either way
	const bvec2 lin = bondable_pairs(atoms_at({ { 28, { 0, 0, 0 } }, { 6, { 0, 0, 1.8 } }, { 8, { 0, 0, 2.95 } } }), 1.75);
	EXPECT_TRUE(lin[0][1] && lin[1][0] && lin[1][2] && lin[2][1]);
	EXPECT_FALSE(lin[0][2] || lin[2][0]);
	const double a = 60.0 * constants::PI / 180.0;
	const bvec2 bent = bondable_pairs(atoms_at({ { 28, { 0, 0, 0 } }, { 6, { 0, 0, 1.8 } },
												 { 8, { 1.15 * std::sin(a), 0, 1.8 + 1.15 * std::cos(a) } } }), 1.75);
	EXPECT_TRUE(bent[0][2] && bent[2][0]);
	EXPECT_FALSE(bondable_pairs(atoms_at({ { 28, { 0, 0, 0 } }, { 6, { 0, 0, 1.8 } }, { 8, { 0, 0, 2.95 } } }), 1.3)[0][2]);
	for (int i = 0; i < 3; i++) EXPECT_FALSE(lin[i][i]);
}

TEST(NboCleanupTests, BondablePairsGiveALonelyAtomItsNearestNeighbour)
{
	//1.5 A is far beyond 1.3 (r_H + r_H); H2 is 10 A from H1 and 11.5 A from H0
	const bvec2 b = bondable_pairs(atoms_at({ { 1, { 0, 0, 0 } }, { 1, { 1.5, 0, 0 } }, { 1, { 11.5, 0, 0 } } }), 1.3);
	EXPECT_TRUE(b[0][1] && b[1][0]);
	EXPECT_TRUE(b[1][2] && b[2][1]);
	EXPECT_FALSE(b[0][2] || b[2][0]);
	//O-H 0.96 A is bondable, H...H 1.52 A in water is not, and nobody needs a fallback
	const double h = 104.5 * constants::PI / 360.0;
	const bvec2 w = bondable_pairs(atoms_at({ { 8, { 0, 0, 0 } }, { 1, { 0.96 * std::sin(h), 0, 0.96 * std::cos(h) } },
											  { 1, { -0.96 * std::sin(h), 0, 0.96 * std::cos(h) } } }));
	EXPECT_TRUE(w[0][1] && w[0][2]);
	EXPECT_FALSE(w[1][2]);
}
