#pragma once

#include "convenience.h"
// These values are made to fit the spherical density of the valence orbitals to the real radial nodal behaviour for the ECPs of the def2 basis sets
//  def2-SVP, def2-TZVP, def2-TZVPP, def2-QZVP, def2-QZVPP
namespace def_corrections
{
    extern const vec2 z;
    extern const vec2 c;

} // namespace def_corrections

namespace ptb_corrections
{
    extern const vec2 z;
    extern const vec2 c;
}

namespace xtb_corrections
{
    extern const vec2 z;
    extern const vec2 c;
}