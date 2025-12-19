//
// Created by lucas on 11/30/25.
//
#include <occ/main/occ_scf.h>
#include <occ/io/occ_input.h>
#include <occ/qm/wavefunction.h>
#include "compareWfns.h"
#include "../Src/wfn_class.h"

int main()
{
    occ::io::OccInput config = occ::io::read_occ_input_file("/home/lucas/CLionProjects/NoSpherA2/occ_integration_tests/alanine.toml");
    auto wfn_direct = occ::main::run_scf_external(config, true);
    auto wfn = occ::qm::Wavefunction::load("/home/lucas/CLionProjects/NoSpherA2/occ_integration_tests/alanine.owf.fchk");
    return 0;

}
