//
// Created by lucas on 11/17/25.
//

#include "WfnAdapter.h"
#include <fstream>
#include <occ/core/atom.h>
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <occ/qm/wavefunction.h>
#include <occ/core/util.h>
#include <fmt/ostream.h>
#include <occ/core/linear_algebra.h>
#include <occ/core/util.h>
#include <occ/io/fchkreader.h>
#include <occ/io/fchkwriter.h>
#include <occ/io/gaussian_input_file.h>
#include <occ/io/moldenreader.h>
#include <occ/io/wavefunction_json.h>
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <sstream>
#include <iostream>

using occ::format_matrix;
using occ::Mat;
using occ::qm::HartreeFock;
using occ::util::all_close;
using occ::io::JsonWavefunctionReader;
using occ::io::JsonWavefunctionWriter;
using occ::qm::Wavefunction;
int main()
{
    Wavefunction wfn = Wavefunction::load("/home/lucas/CLionProjects/NoSpherA2/tests_occ/water.owf.json");
    // WfnAdapter wfn2(wfn);
    // printf("NAtoms: %d", wfn2.get_ncen());
    WFN wfn2(wfn, WfnOrigin::UNKNOWN);
    return 0;
}
