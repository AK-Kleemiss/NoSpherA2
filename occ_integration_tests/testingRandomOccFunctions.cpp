//
// Created by lucas on 11/23/25.
//
#include <algorithm>
#include <cstdio>
#include <occ/qm/wavefunction.h>
#include <occ/io/conversion.h>
#include <Eigen/Eigen>
#include <iostream>
#include <ostream>
#include "../Src/convenience.h"
#include "../Src/wfn_class.h"
using namespace occ::qm;
using namespace std;


int main()
{
    const string filepath("/home/lucas/CLionProjects/NoSpherA2/occ_integration_tests/alanine.owf.fchk");
    auto wfn_nos = WFN(filepath);
    auto wfn = Wavefunction::load(filepath);
    volatile auto orb =  occ::io::conversion::orb::to_gaussian_order(wfn.basis, wfn.mo);
    return 0;
}
