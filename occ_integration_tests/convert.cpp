//
// Created by lucas on 11/21/25.
//
#include <algorithm>
#include <vector>
#include <ranges>
#include <fmt/core.h>
#include <fmt/ranges.h> // Essential for printing tuples/vectors
#include <cstdio>
#include <occ/qm/wavefunction.h>
#include <Eigen/Eigen>
#include <iostream>
#include <ostream>
#include "../Src/convenience.h"
#include "../Src/wfn_class.h"
#include <tuple>
#include <type_traits>
#include <fmt/base.h>
using namespace occ::qm;
using namespace std;
struct VecSize
{
    int x{1};
    int y{1};
    bool is_empty() const { return (x==1) && (y==1); }
};
std::string get_shape(const occ::Mat& mat)
{
    return std::format("({}, {})", mat.rows(), mat.cols());
}
template <typename T>
void compare_vectors(std::vector<T> vecNOS, std::vector<T> vecOCC, std::string label,
    VecSize size = VecSize(), double tol = 1e-1, bool compare_diff_size = false)
{
    bool differ =false;
    if (vecNOS.size() != vecOCC.size())
    {
        fmt::print("[{}] sizes differ, [NOS]: {}, [OCC]: {}.\n", label, vecNOS.size(), vecOCC.size());
        if (compare_diff_size)
        {
            if (vecNOS.size() < vecOCC.size())
            {
                vecNOS.resize(vecOCC.size(), 0);
            } else
                vecOCC.resize(vecNOS.size(), 0);
        } else
            return;
    }

    bool comp;
    for (int i=0; i<vecNOS.size(); i++)
    {
        if constexpr (std::is_integral_v<T>)  comp =vecNOS[i] != vecOCC[i];
        else comp =  abs(vecNOS[i]) - abs(vecOCC[i]) > tol;
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
        std::cout << "[" << label << "] from occ constructor and from file are equal.\n";
        return;
    }
    if (size.is_empty())
    {
        size.x = vecNOS.size();
    }

    if (differ) {
        ofstream logFile;
        logFile.open (fmt::format("{}_matrix", label));
        for (int x=0; x < size.x; x++)
        {
            for (int y=0; y< size.y; y++)
            {
                int idx = y + x*size.y;
                if (idx >= vecOCC.size()) break;
                logFile << fmt::format("  {}", vecOCC[idx]);
            }
            logFile << std::endl;
        }
        logFile.close();
    }
}

void compare_MOs(std::vector<MO>& moNOS, std::vector<MO>& moOCC, double tol = 1e-1)
{
    if (moNOS.size() != moOCC.size())
    {
        fmt::print("[MO] sizes differ, [NOS]: {}, [OCC]: {}.\n", moNOS.size(), moOCC.size());

    }
    for (int i=0; i< moOCC.size(); i++)
    {
        if (abs(moNOS[i].get_energy()) - abs(moOCC[i].get_energy()) > tol)
        {
            fmt::print("[MO] energies differ, [NOS]: {}, [OCC]: {}. \n", moNOS[i].get_energy(), moOCC[i].get_energy());
        }
        if (abs(moNOS[i].get_occ()) - abs(moOCC[i].get_occ()) > tol)
        {
            fmt::print("[MO] occupancies (occ) differ, [NOS]: {}, [OCC]: {}. \n", moNOS[i].get_occ(), moOCC[i].get_occ());
            return;
        }
        if (moNOS[i].get_op() != moOCC[i].get_op())
        {
            fmt::print("[MO] OPs differ, [NOS]: {}, [OCC]: {}. \n", moNOS[i].get_op(), moOCC[i].get_op());
            return;
        }
        const auto& moNOSC = moNOS[i].get_coefficients();
        const auto& moOCCC = moOCC[i].get_coefficients();
        if (moNOSC.size() != moOCCC.size())
        {
            fmt::print("[MO] coefficients differ in size, [NOS]: {}, [OCC]: {}. \n",
                moNOS[i].get_coefficients().size(), moOCC[i].get_coefficients().size());
            compare_vectors(moNOSC, moOCCC, fmt::format("moCOEFS[{}]", i), VecSize(), 1e-6,true);
            return;
        }
    }
    fmt::print("[MOs] from occ constructor and from file are equal.");
}


WFN wfn_from_nos(const std::string &filepath)
{
    return WFN(filepath);
}

int main()
{
    std::string filepathnos("/home/lucas/CLionProjects/NoSpherA2/tests/alanine_occ/alanine.owf.fchk");
    std::string filepathocc("/home/lucas/CLionProjects/NoSpherA2/tests/alanine_occ/alanine.owf.fchk");
    auto wfn_nos = WFN(filepathnos);
    Wavefunction wfn = Wavefunction::load(filepathocc);
    auto wfn_from_occ = WFN(wfn);
    auto centersvecocc = wfn_from_occ.get_centers();
    auto centersvecnos = wfn_nos.get_centers();
    compare_vectors(centersvecnos, centersvecocc, "centers");

    auto typesvecocc = wfn_from_occ.get_types();
    auto typesvecnos = wfn_nos.get_types();
    compare_vectors(typesvecnos, typesvecocc, "types",
        VecSize{.x = 21, .y = 20});

    auto expsvecocc = wfn_from_occ.get_exponents();
    auto expsvecnos = wfn_nos.get_exponents();
    compare_vectors(expsvecnos, expsvecocc, "exponents");

    auto mosOCC = wfn_from_occ.get_MOs_vec();
    auto mosNOS = wfn_nos.get_MOs_vec();
    compare_MOs(mosNOS, mosOCC);
    return 0;
}
