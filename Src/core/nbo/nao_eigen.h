#pragma once
#include "convenience.h"

#include <Eigen/Dense>

//dMatrix2 <-> Eigen copies shared by nao.cpp, nbo.cpp and nrt.cpp
inline Eigen::MatrixXd to_eigen(const dMatrix2& m)
{
	const int r = static_cast<int>(m.extent(0)), c = static_cast<int>(m.extent(1));
	Eigen::MatrixXd out(r, c);
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			out(i, j) = m(i, j);
	return out;
}

inline dMatrix2 to_dmatrix(const Eigen::MatrixXd& m)
{
	dMatrix2 out(m.rows(), m.cols());
	for (int i = 0; i < m.rows(); i++)
		for (int j = 0; j < m.cols(); j++)
			out(i, j) = m(i, j);
	return out;
}
