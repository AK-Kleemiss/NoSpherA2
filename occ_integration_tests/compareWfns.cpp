//
// Created by lucas on 11/30/25.
//
#include "compareWfns.h"
#include <ranges>
#include <iostream>
#include "../Src/convenience.h"
#include <type_traits>
#include <fmt/base.h>
#include <fmt/core.h>

#include "../Src/constants.h"
class MO;
using namespace occ::qm;
using namespace std;

bool compare_MOs(const std::vector<MO>& mo1, const std::vector<MO>& mo2, double tol)
{
    if (mo1.size() != mo2.size())
    {
        fmt::print("[MO] sizes differ, [1]: {}, [1]: {}.\n", mo1.size(), mo2.size());
    }
    for (int i=0; i< mo1.size(); i++)
    {
        if (abs(mo1[i].get_energy()) - abs(mo2[i].get_energy()) > tol)
        {
            fmt::print("[MO] energies differ, [1]: {}, [1]: {}. \n", mo1[i].get_energy(), mo2[i].get_energy());
        }
        if (abs(mo1[i].get_occ()) - abs(mo2[i].get_occ()) > tol)
        {
            fmt::print("[MO] occupancies (occ) differ, [1]: {}, [1]: {}. \n", mo1[i].get_occ(), mo2[i].get_occ());
            return false;
        }
        if (mo1[i].get_op() != mo2[i].get_op())
        {
            fmt::print("[MO] OPs differ, [1]: {}, [1]: {}. \n", mo1[i].get_op(), mo2[i].get_op());
            return false;
        }
        const auto& mo1C = mo1[i].get_coefficients();
        const auto& mo2C = mo2[i].get_coefficients();
        if (mo1C.size() != mo2C.size())
        {
            fmt::println("[MO] coefficients differ in size, [1]: {}, [1]: {}. \n",
               mo1C.size(), mo2C.size());
            compare_vectors(mo1C, mo2C, fmt::format("moCOEFS[{}]", i), VecSize(), 1e-6,true);
            return false;
        }
        if (!compare_vectors(mo1C, mo2C, fmt::format("moCOEFS[{}]", i), VecSize(), 1e-6,true, false)) return false;
    }
    fmt::println("[MOs] from occ constructor and from file are equal.");
    return true;
}
bool compare_Atoms(const std::vector<atom>& atoms1, const std::vector<atom>& atoms2, double tol)
{
    if (atoms1.size() != atoms2.size())
    {
        fmt::print("[atoms] sizes differ, [1]: {}, [1]: {}.\n", atoms1.size(), atoms2.size());
        return false;
    }
    for (int i=0; i< atoms1.size(); i++)
    {
        if (abs(atoms1[i].get_charge()) - abs(atoms2[i].get_charge()) > tol)
        {
            fmt::print("[Atoms] charges differ, [1]: {}, [1]: {}. \n", atoms1[i].get_charge(), atoms2[i].get_charge());
            return false;
        }
        if (abs(atoms1[i].get_ECP_electrons()) - abs(atoms2[i].get_ECP_electrons()) > tol)
        {
            fmt::print("[Atoms] ECP electrons (occ) differ, [1]: {}, [1]: {}. \n", atoms1[i].get_ECP_electrons(), atoms2[i].get_ECP_electrons());
            return false;
        }
        auto coords1 = atoms1[i].get_coords();
        auto coords2 = atoms2[i].get_coords();
        if (!compare_vectors(std::vector(coords1.begin(), coords1.end()),
            std::vector(coords2.begin(), coords2.end()),
            fmt::format("Atoms", i), VecSize(), 1e-11,true, false))
            return false;
    }
    fmt::println("[Atoms] from occ constructor and from file are equal.");
    return true;
}

BenchmarkResults conversion_benchmark(std::string filepath, bool get_wfns)
{
    auto start_nos = std::chrono::high_resolution_clock::now();
    auto wfn_nos = WFN(filepath);
    auto end_nos = std::chrono::high_resolution_clock::now();
    auto duration_nos = std::chrono::duration_cast<std::chrono::microseconds>(end_nos - start_nos);
    auto start_occ = std::chrono::high_resolution_clock::now();
    Wavefunction wfn = Wavefunction::load(filepath);
    auto end_occ = std::chrono::high_resolution_clock::now();
    auto duration_occ = std::chrono::duration_cast<std::chrono::microseconds>(end_occ - start_occ);

    auto start_occ2nos = std::chrono::high_resolution_clock::now();
    auto wfn_from_occ = WFN(wfn);
    auto end_occ2nos = std::chrono::high_resolution_clock::now();
    auto duration_occ2nos = std::chrono::duration_cast<std::chrono::microseconds>(end_occ2nos - start_occ2nos);
    ConversionTimes times = {static_cast<long>(duration_occ.count()), static_cast<long>(duration_nos.count()), static_cast<long>(duration_occ2nos.count())};
    return get_wfns ? BenchmarkResults{times, wfn_nos, wfn_from_occ} : BenchmarkResults{times, std::nullopt, std::nullopt};
}

bool wfnequal(WFN& wfn1, WFN& wfn2)
{
    const auto centers1 = wfn1.get_centers();
    const auto centers2 = wfn2.get_centers();
    const bool centers = compare_vectors(centers1, centers2, "centers");

    const auto types1 = wfn1.get_types();
    const auto types2 = wfn2.get_types();
    bool types = compare_vectors(types1, types2, "types",
        VecSize{.x = 21, .y = 20});

    const auto exps1 = wfn1.get_exponents();
    const auto exps2 = wfn2.get_exponents();
    bool exps = compare_vectors(exps1, exps2, "exponents");

    const auto atoms1 = wfn1.get_atoms();
    const auto atoms2 = wfn2.get_atoms();
    bool atoms = compare_Atoms(atoms1, atoms2);

    const auto mos1 = wfn1.get_MOs_vec();
    const auto mos2 = wfn2.get_MOs_vec();
    bool mos = compare_MOs(mos1, mos2);
    return centers && types && exps && atoms && mos;
}
