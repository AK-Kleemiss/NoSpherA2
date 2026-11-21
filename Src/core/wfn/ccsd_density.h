#pragma once
#include <occ/io/occ_input.h>
#include <occ/qm/cc/density.h>
#include <occ/qm/wavefunction.h>
#include <iosfwd>
#include <string>

//the correlated densities (relaxed MP2, unrelaxed CCSD and CCSD(T) from the Lambda equations, local CCSD(T)) live in
//occ, occ/qm/cc/density.h; these names keep NoSpherA2's callers and tests
using CcsdDensityOptions = occ::qm::cc::CorrelatedDensityOptions;
using CcsdDensityResult = occ::qm::cc::CorrelatedDensity;

//occ's densities with NoSpherA2's log summary; occ's exceptions end the run through err_checkf
CcsdDensityResult ccsd_natural_orbitals(const occ::qm::Wavefunction &scf, const CcsdDensityOptions &opts, std::ostream &log);
CcsdDensityResult local_ccsd_natural_orbitals(const occ::qm::Wavefunction &scf, const CcsdDensityOptions &opts, std::ostream &log);
CcsdDensityResult mp2_natural_orbitals(const occ::qm::Wavefunction &scf, const CcsdDensityOptions &opts, std::ostream &log);
//unrelaxed MP3 (order 3) or MP4 (order 4, MP4(SDQ) unless opts.triples), as ORCA's AutoCI
CcsdDensityResult mpn_natural_orbitals(const occ::qm::Wavefunction &scf, int order, const CcsdDensityOptions &opts, std::ostream &log);

//-occ with method = "mp2", "mp3", "mp4(sdq)", "mp4", "ccsd" or "ccsd(t)": occ's RHF, then the density; ccsd_* keys of the toml's [scf] table
CcsdDensityResult ccsd_density_from_occ_input(const std::string &toml_path, occ::io::OccInput config, std::ostream &log);

//every -occ toml goes through here: MP2, MP3, MP4, CCSD and CCSD(T) give their natural orbitals, anything else is occ's SCF
occ::qm::Wavefunction wavefunction_from_occ_input(const std::string &toml_path, const occ::io::OccInput &config, std::ostream &log);
