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

class wfn_addapter : public WFN
{
    public:
        wfn_addapter(const occ::qm::Wavefunction& wf) : m_originalWfn(wf) {};
    private:
        const occ::qm::Wavefunction& m_originalWfn;
};


#endif //NOSPHERA2_WFN_ADDAPTER_H
