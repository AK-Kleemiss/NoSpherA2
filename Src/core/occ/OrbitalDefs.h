//
// Created by lucas on 23/11/25.
//
#pragma once
#include <cmath>
#include <Eigen/Dense>

inline std::span<const double> occ_vec_span(occ::Vec &v)
{
    return std::span<const double>(v.data(), v.size());
}

inline std::span<const double> eigen_vec_span(const Eigen::VectorXd &v)
{
    return std::span<const double>(v.data(), v.size());
}
