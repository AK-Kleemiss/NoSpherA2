//
// Created by lucas on 11/23/25.
//
#include <algorithm>
#include <cstdio>
#include <occ/qm/wavefunction.h>
#include <occ/gto/gto.h>
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
    const string filepath("alanine.owf.fchk");
    auto wfn_nos = WFN(filepath);
    auto wfn = Wavefunction::load(filepath);

    auto test= occ::gto::transform_density_matrix_cartesian_to_spherical(wfn.basis, wfn.mo.D);

    wfn.save("alanine2.fchk");
    // for (auto c : wfn_nos.get_MO(0).get_coefficients())
    // {
    //     cout << c << ",\n";
    // }
    return 0;
}
