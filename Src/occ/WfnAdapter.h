//
// Created by lucas on 11/17/25.
//

#ifndef NOSPHERA2_WFN_ADDAPTER_H
#define NOSPHERA2_WFN_ADDAPTER_H
#include "../wfn_class.h"
#include <occ/qm/wavefunction.h>

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

class WfnAdapter : public WFN
{
    public:
        WfnAdapter(occ::qm::Wavefunction& wf)
            : m_originalWfn(wf),
            m_ncen(static_cast<int>(wf.atoms.size())),
            m_nfunc(static_cast<int>(wf.nbf)),
            m_nmo() {}

        const int& get_ncen() const override { return m_ncen; }

    private:
        occ::qm::Wavefunction& m_originalWfn;
        //Number of centers/atoms present in the wavefunction
        int m_ncen;
        //Number of basis functions/primitive Gaussians present in the wavefunction
        int m_nfunc;
        //Number of molecular orbitals present in the wavefunction
        int m_nmo;
        //Number of primitive exponents present in the wavefunction
        int m_nex;
        //Total charge of the wavefunction
        int m_charge;
        //Number of the ECPs modi, as specified in the ECPs_corrections.h file
        int m_ECP_m;
};


#endif //NOSPHERA2_WFN_ADDAPTER_H
