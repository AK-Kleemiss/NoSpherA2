//
// Created by lucas on 11/30/25.
//
#ifndef NOSPHERA2_COMPAREWFNS_H
#define NOSPHERA2_COMPAREWFNS_H

#include <string>
#include <fmt/ranges.h> // Essential for printing tuples/vectors
#include <format>
#include "../Src/convenience.h"
#include <occ/qm/wavefunction.h>
#include "../Src/wfn_class.h"
struct VecSize
{
    int x{1};
    int y{1};
    [[nodiscard]] bool is_empty() const { return (x==1) && (y==1); }
};

std::string inline get_shape(const occ::Mat& mat)
{
    return std::format("({}, {})", mat.rows(), mat.cols());
}
bool compare_MOs(const std::vector<MO>& moNOS, const std::vector<MO>& moOCC, double tol = 1e-6);

template <typename T>
bool compare_vectors(std::vector<T> vecNOS, std::vector<T> vecOCC, std::string label,
    VecSize size = VecSize(), double tol = 1e-6, bool compare_diff_size = false, bool print_on_success = true)
{
    bool differ =false;
    if (vecNOS.size() != vecOCC.size())
    {
        fmt::println("[{}] sizes differ, [NOS]: {}, [OCC]: {}.\n", label, vecNOS.size(), vecOCC.size());
        if (compare_diff_size)
        {
            if (vecNOS.size() < vecOCC.size())
            {
                vecNOS.resize(vecOCC.size(), 0);
            } else
                vecOCC.resize(vecNOS.size(), 0);
        } else
            return false;
    }

    bool comp;
    for (int i=0; i<vecNOS.size(); i++)
    {
        if constexpr (std::is_integral_v<T>)  comp =vecNOS[i] != vecOCC[i];
        else comp =  abs(vecNOS[i] - vecOCC[i]) > tol;
        if (comp)
        {
            std::cout << "idx: " << i << " - "<< label<<"NOS: " << vecNOS[i] << " - " << label << "OCC: " <<vecOCC[i] << "\n";
            differ = true;
        }
    }
    if (differ)
        std::cout << "[" << label << "] from occ constructor and from file differ.\n";
    else
    {
        if (print_on_success)
            std::cout << "[" << label << "] from occ constructor and from file are equal.\n";
        return true;
    }
    if (size.is_empty())
    {
        size.x = vecNOS.size();
    }

    return false;
}

bool compare_Atoms(const std::vector<atom>& atomsNOS, const std::vector<atom>& atomsOCC, double tol = 1e-6);

struct ConversionTimes
{
    long file2occ;
    long file2nos;
    long occ2nos;
};
struct BenchmarkResults
{
    ConversionTimes times;
    std::optional<WFN> file2NOS;
    std::optional<WFN> occ2NOS;
};
BenchmarkResults conversion_benchmark(std::string filepath, bool get_wfns = false);

bool wfnequal(WFN& wfn1, WFN& wfn2);
#endif //NOSPHERA2_COMPAREWFNS_H
