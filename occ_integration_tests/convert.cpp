//
// Created by lucas on 11/21/25.
//
#include <algorithm>
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
void compare_vectors(const std::vector<T>& vecNOS, const std::vector<T>& vecOCC, std::string label,
    VecSize size = VecSize(), double tol = 1e-1)
{
    bool differ =false;
    if (vecNOS.size() != vecOCC.size())
    {
        fmt::print("[{}] sizes differ, [NOS]: {}, [OCC]: {}.\n", label, vecNOS.size(), vecOCC.size());
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


WFN wfn_from_nos(const std::string &filepath)
{
    return WFN(filepath);
}

int main()
{
    std::string filepathnos("/home/lucas/CLionProjects/NoSpherA2/occ_integration_tests/alanine.wfn");
    std::string filepathocc("/home/lucas/CLionProjects/NoSpherA2/occ_integration_tests/alanine.owf.fchk");
    auto wfn_nos = WFN(filepathnos);
    Wavefunction wfn = Wavefunction::load(filepathocc);
    auto wfn_from_occ = WFN(wfn);
    auto centersvecocc = wfn_from_occ.get_centers();
    auto centersvecnos = wfn_nos.get_centers();

    auto typesvecocc = wfn_from_occ.get_types();
    auto typesvecnos = wfn_nos.get_types();

    auto expsvecocc = wfn_from_occ.get_exponents();
    auto expsvecnos = wfn_nos.get_exponents();
    compare_vectors(typesvecnos, typesvecocc, "types",
        VecSize{.x = 21, .y = 20});
    compare_vectors(expsvecnos, expsvecocc, "exponents");
    compare_vectors(centersvecnos, centersvecocc, "centers");
    return 0;
}
