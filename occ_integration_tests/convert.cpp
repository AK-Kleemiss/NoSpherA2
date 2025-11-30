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
    auto results = conversion_benchmark(filepath, true);
    auto duration_occ = results.times.file2occ;
    auto duration_occ2nos = results.times.occ2nos;
    auto duration_wfn2nos = results.times.file2nos;
    fmt::println("Conversion from file to NoSpherA2 took {}µs.", duration_wfn2nos);
    fmt::println("Conversion from file to OCC took {}µs.", duration_occ);
    fmt::println("Conversion from OCC to NospherA2 took {}µs.", duration_occ2nos);
    fmt::println("Total, OCC and OCC2NOS took {}µs.", duration_occ2nos + duration_occ);
    fmt::println("---------------------------------------------");
    auto wfn_from_occ = results.occ2NOS;
    auto wfn_nos = results.file2NOS;

    wfnequal(wfn_from_occ, wfn_nos);
    return 0;
}
