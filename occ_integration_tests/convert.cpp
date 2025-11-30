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
#include <chrono>
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
bool compare_vectors(std::vector<T> vecNOS, std::vector<T> vecOCC, std::string label,
    VecSize size = VecSize(), double tol = 1e-12, bool compare_diff_size = false, bool print_on_success = true)
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
    return false;
}

bool compare_MOs(std::vector<MO>& moNOS, std::vector<MO>& moOCC, double tol = 1e-6)
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
            return false;
        }
        if (moNOS[i].get_op() != moOCC[i].get_op())
        {
            fmt::print("[MO] OPs differ, [NOS]: {}, [OCC]: {}. \n", moNOS[i].get_op(), moOCC[i].get_op());
            return false;
        }
        const auto& moNOSC = moNOS[i].get_coefficients();
        const auto& moOCCC = moOCC[i].get_coefficients();
        if (moNOSC.size() != moOCCC.size())
        {
            fmt::println("[MO] coefficients differ in size, [NOS]: {}, [OCC]: {}. \n",
                moNOS[i].get_coefficients().size(), moOCC[i].get_coefficients().size());
            compare_vectors(moNOSC, moOCCC, fmt::format("moCOEFS[{}]", i), VecSize(), 1e-6,true);
            return false;
        }
        if (!compare_vectors(moNOSC, moOCCC, fmt::format("moCOEFS[{}]", i), VecSize(), 1e-6,true, false)) return false;
    }
    fmt::println("[MOs] from occ constructor and from file are equal.");
    return true;
}
void compare_Atoms(std::vector<atom>& atomsNOS, std::vector<atom>& atomsOCC, double tol = 1e-6)
{
    if (atomsNOS.size() != atomsOCC.size())
    {
        fmt::print("[atoms] sizes differ, [NOS]: {}, [OCC]: {}.\n", atomsNOS.size(), atomsOCC.size());

    }
    for (int i=0; i< atomsOCC.size(); i++)
    {
        if (abs(atomsNOS[i].get_charge()) - abs(atomsOCC[i].get_charge()) > tol)
        {
            fmt::print("[Atoms] charges differ, [NOS]: {}, [OCC]: {}. \n", atomsNOS[i].get_charge(), atomsOCC[i].get_charge());
        }
        if (abs(atomsNOS[i].get_ECP_electrons()) - abs(atomsOCC[i].get_ECP_electrons()) > tol)
        {
            fmt::print("[Atoms] ECP electrons (occ) differ, [NOS]: {}, [OCC]: {}. \n", atomsNOS[i].get_ECP_electrons(), atomsOCC[i].get_ECP_electrons());
            return;
        }
        auto coordsNOS = atomsNOS[i].get_coords();
        auto coordsOCC = atomsOCC[i].get_coords();
        compare_vectors(std::vector(coordsNOS.begin(), coordsNOS.end()),
            std::vector(coordsOCC.begin(), coordsOCC.end()),
            fmt::format("Atoms", i), VecSize(), 1e-11,true, false);
    }
    fmt::println("[Atoms] from occ constructor and from file are equal.");
}

WFN wfn_from_nos(const std::string &filepath)
{
    return WFN(filepath);
}

int main()
{
    std::string filepath("/home/lucas/CLionProjects/NoSpherA2/occ_integration_tests/alanine.owf.fchk");
    auto start_nos = std::chrono::high_resolution_clock::now();
    auto wfn_nos = WFN(filepath);
    auto end_nos = std::chrono::high_resolution_clock::now();
    auto duration_nos = std::chrono::duration_cast<std::chrono::microseconds>(end_nos - start_nos);
    fmt::println("Conversion from file to NospherA2 took {}µs.", duration_nos.count());
    auto start_occ = std::chrono::high_resolution_clock::now();
    Wavefunction wfn = Wavefunction::load(filepath);
    auto end_occ = std::chrono::high_resolution_clock::now();
    auto duration_occ = std::chrono::duration_cast<std::chrono::microseconds>(end_occ - start_occ);

    auto start_occ2nos = std::chrono::high_resolution_clock::now();
    auto wfn_from_occ = WFN(wfn);
    auto end_occ2nos = std::chrono::high_resolution_clock::now();
    auto duration_occ2nos = std::chrono::duration_cast<std::chrono::microseconds>(end_occ2nos - start_occ2nos);

    fmt::println("Conversion from file to OCC took {}µs.", duration_occ.count());
    fmt::println("Conversion from OCC to NospherA2 took {}µs.", duration_occ2nos.count());
    fmt::println("Total, OCC and OCC2NOS took {}µs.", duration_occ2nos.count() + duration_occ.count());
    fmt::println("---------------------------------------------");

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

    auto atomsOCC = wfn_from_occ.get_atoms();
    auto atomsNOS = wfn_nos.get_atoms();
    compare_Atoms(atomsNOS, atomsOCC);

    auto mosOCC = wfn_from_occ.get_MOs_vec();
    auto mosNOS = wfn_nos.get_MOs_vec();
    if (!compare_MOs(mosNOS, mosOCC))
        return 1;
    return 0;
}
