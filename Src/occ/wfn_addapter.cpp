//
// Created by lucas on 11/17/25.
//

#include "wfn_addapter.h"
#include <occ/core/atom.h>
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <occ/qm/wavefunction.h>
#include <occ/core/util.h>
const char *json_contents = R"(
{
    "Molecule": {
        "Atoms": [
            {
                "BasisFunctions": [
                    {
                        "Coefficients": [
                            0.1543289707029839,
                            0.5353281424384732,
                            0.44463454202535485
                        ],
                        "Exponents": [
                            3.42525091,
                            0.62391373,
                            0.1688554
                        ],
                        "Shell": "s"
                    }
                ],
                "Coords": [
                    0.0,
                    0.0,
                    0.0
                ],
                "ElementLabel": "H",
                "ElementNumber": 1,
                "Idx": 0,
                "NuclearCharge": 1.0
            },
            {
                "BasisFunctions": [
                    {
                        "Coefficients": [
                            0.1543289707029839,
                            0.5353281424384732,
                            0.44463454202535485
                        ],
                        "Exponents": [
                            3.42525091,
                            0.62391373,
                            0.1688554
                        ],
                        "Shell": "s"
                    }
                ],
                "Coords": [
                    0.0,
                    0.0,
                    2.6456165874897524
                ],
                "ElementLabel": "H",
                "ElementNumber": 1,
                "Idx": 1,
                "NuclearCharge": 1.0
            }
        ],
        "BaseName": "h2",
        "Charge": 0,
        "CoordinateUnits": "Bohrs",
        "H-Matrix": [
            [
                -0.8422289108173067,
                -0.3785030189433606
            ],
            [
                -0.3785030189433606,
                -0.8422289108173067
            ]
        ],
        "HFTyp": "RHF",
        "MolecularOrbitals": {
            "EnergyUnit": "Eh",
            "MOs": [
                {
                    "MOCoefficients": [
                        0.621202118970972,
                        0.621202118970972
                    ],
                    "Occupancy": 2.0,
                    "OrbitalEnergy": -0.3773228237944341
                },
                {
                    "MOCoefficients": [
                        -0.8425697665943784,
                        0.8425697665943784
                    ],
                    "Occupancy": 0.0,
                    "OrbitalEnergy": 0.25890197203084375
                }
            ],
            "OrbitalLabels": [
                "0H   1s",
                "1H   1s"
            ]
        },
        "Multiplicity": 1,
        "S-Matrix": [
            [
                1.0,
                0.29569907102006404
            ],
            [
                0.29569907102006404,
                1.0
            ]
        ],
        "T-Matrix": [
            [
                0.7600318835666087,
                0.019745102188373418
            ],
            [
                0.019745102188373418,
                0.7600318835666087
            ]
        ]
    },
    "ORCA Header": {
        "Version": "5.0 - current"
    }
}
)";

using occ::format_matrix;
using occ::Mat;
using occ::qm::HartreeFock;
using occ::util::all_close;
using occ::io::JsonWavefunctionReader;
using occ::io::JsonWavefunctionWriter;
int main()
{

    std::vector<occ::core::Atom> atoms{{8, -1.32695761, -0.10593856, 0.01878821},
                                     {1, -1.93166418, 1.60017351, -0.02171049},
                                     {1, 0.48664409, 0.07959806, 0.00986248}};

    auto obs = occ::qm::AOBasis::load(atoms, "6-31G**");
    obs.set_pure(true);
    HartreeFock hf(obs);
    occ::qm::SCF<HartreeFock> scf(hf);
    scf.convergence_settings.energy_threshold = 1e-8;
    double e = scf.compute_scf_energy();

    occ::qm::Wavefunction wfn = scf.wavefunction();
    wfn_addapter wfn2(wfn);

    return 0;
}
