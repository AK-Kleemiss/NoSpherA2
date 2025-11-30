//
// Created by lucas on 11/21/25.
//

#include "../Src/convenience.h"
#include <occ/qm/wavefunction.h>
#include "../Src/wfn_class.h"
#include "compareWfns.h"

using namespace occ::qm;
using namespace std;
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
