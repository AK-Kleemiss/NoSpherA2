#include "pch.h"
#include "TestRunner.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace NosTestFramework;

// ---------------------------------------------------------------------------
// Integration tests: each method mirrors one entry in tests/tests.toml.
// The test copies the test directory to a temp folder, runs NoSpherA2 in-process
// through NoSpherA2_DLL.dll, then compares the output log against the
// corresponding .good file using the same 1% relative numeric tolerance as the
// Python pytest harness.
//
// Prerequisites:
//   - NoSpherA2_DLL.dll must be built (it lives in the same output directory as
//     this test DLL: Windows/x64/<Config>/).
//   - MKLROOT / OCC_DATA_PATH are set automatically; OCC_DATA_PATH is
//     injected per-process from <repo>/Lib/occ/share/occ.
//   - Full tests (fchk_conversion, fourier_transform_full) are skipped
//     unless the environment variable RUN_FULL_TEST is set.
// ---------------------------------------------------------------------------

namespace NoSpherA2IntegrationTests
{
    TEST_CLASS(IntegrationTests)
    {
    public:

        TEST_METHOD(alanine_occ)
        {
            RunTest({
                "alanine_occ", "alanine_occ", "alanine_occ.good", "",
                {"-cif","alanine.cif", "-dmin","0.5", "-wfn","alanine.owf.fchk",
                 "-acc","1", "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(alanine_integrated_occ)
        {
            // subprocess=true: OCC's TBB parallel integrals cause heap corruption
            // when run in-process inside the managed vstest host. The subprocess
            // path gives process isolation and is exactly what CI (pytest) uses.
            RunTest({
                "alanine_integrated_occ", "alanine_integrated_occ",
                "alanine_integrated_occ.good", "",
                {"-occ","alanine.toml", "-cif","alanine.cif",
                 "-dmin","0.5", "-acc","1", "-all_charges", "-no_date"},
                false, true
            });
        }

        TEST_METHOD(disorder_THPP)
        {
            RunTest({
                "disorder_THPP", "disorder", "disorder_THPP.good", "",
                {"-cif","thpp.cif", "-hkl","thpp.hkl", "-acc","1",
                 "-mtc",
                     "olex2/Wfn_job/Part_1/thpp.wfx","0.1",
                     "olex2/Wfn_job/Part_2/thpp.wfx","0.2",
                 "-mtc_mult","1","1",
                 "-mtc_ECP","0","0",
                 "-mtc_charge","0","0",
                 "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(fourier_transform)
        {
            RunTest({
                "fourier_transform", "SALTED", "fourier_transform.good", "",
                {"-test_analytical", "-acc","1", "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(fractal)
        {
            // Actual output is a cube fractal plot, not the standard log file
            RunTest({
                "fractal", "sucrose_fchk_SF", "fractal.good",
                "sucrose_diff.cube_fractal_plot",
                {"-fractal","sucrose_diff.cube", "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(grown_water)
        {
            RunTest({
                "grown_water", "grown", "grown_water.good", "",
                {"-cif","water.cif", "-hkl","water.hkl", "-wfn","water.wfx",
                 "-acc","1", "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(Hybrid_mode)
        {
            RunTest({
                "Hybrid_mode", "Hybrid", "Hybrid_mode.good", "",
                {"-cif","ZP2.cif", "-dmin","0.9", "-acc","1",
                 "-mtc",
                     "ZP2_part1.gbw","0.1",
                     "ZP2_part2.gbw","0.2",
                 "-mtc_mult","1","1",
                 "-mtc_charge","0","0",
                 "-mtc_ECP","0","0",
                 "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(malbac_SF_ECP)
        {
            RunTest({
                "malbac_SF_ECP", "ECP_SF", "malbac_SF_ECP.good", "",
                {"-cif","malbac.cif", "-hkl","malbac.hkl", "-wfn","malbac.gbw",
                 "-acc","1", "-ECP","1", "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(openBLAS)
        {
            RunTest({
                "openBLAS", "OpenBLAS", "openBLAS.good", "",
                {"-blastest", "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(properties)
        {
            RunTest({
                "properties", "sucrose_fchk_SF", "properties.good", "",
                {"-wfn","olex2/Wfn_job/sucrose.wfx", "-cif","sucrose.cif",
                 "-lap", "-rdg", "-DEF", "-eli",
                 "-resolution","1.5",
                 "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(reading_SALTED)
        {
            RunTest({
                "reading_SALTED", "SALTED", "reading_SALTED.good", "",
                {"-test_reading_SALTED_binary", "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(ri_fit)
        {
            RunTest({
                "ri_fit", "epoxide_gbw", "ri_fit.good", "",
                {"-wfn","epoxide.gbw", "-cif","epoxide.cif",
                 "-dmin","0.4", "-ri_fit","combo_basis_fit",
                 "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(rubredoxin_cmtc)
        {
            RunTest({
                "rubredoxin_cmtc", "rubredoxin_cmtc", "rubredoxin_cmtc.good", "",
                {"-hkl","1yk4_h.hkl", "-acc","1",
                 "-cmtc",
                     "residues/1.gbw","residues/1.cif","0",
                     "residues/2.gbw","residues/2.cif","0.1",
                     "residues/3.gbw","residues/3.cif","0.2",
                     "residues/4.gbw","residues/4.cif","0",
                 "-mtc_mult","1","1","1","1",
                 "-mtc_charge","0","0","0","0",
                 "-mtc_ECP","0","0","0","0",
                 "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(SALTED)
        {
            RunTest({
                "SALTED", "SALTED", "SALTED.good", "",
                {"-SALTED","Model", "-cif","test_cysteine.cif",
                 "-wfn","test_cysteine.xyz", "-dmin","0.73",
                 "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(sucrose_IAM)
        {
            RunTest({
                "sucrose_IAM", "sucrose_IAM_SF", "sucrose_IAM.good", "",
                {"-cif","sucrose.cif", "-hkl","sucrose.hkl",
                 "-xyz","sucrose.xyz", "-IAM",
                 "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(sucrose_ptb)
        {
            RunTest({
                "sucrose_ptb", "sucrose_IAM_SF", "sucrose_ptb.good", "",
                {"-cif","sucrose.cif", "-dmin","0.8", "-wfn","wfn.xtb",
                 "-mult","0", "-charge","0", "-acc","1", "-ECP","3",
                 "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(sucrose_SF)
        {
            RunTest({
                "sucrose_SF", "sucrose_fchk_SF", "sucrose_SF.good", "",
                {"-cif","sucrose.cif",
                 "-hkl","olex2/Wfn_job/sucrose.hkl",
                 "-wfn","olex2/Wfn_job/sucrose.wfx",
                 "-acc","1", "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(sucrose_twin)
        {
            RunTest({
                "sucrose_twin", "sucrose_fchk_SF", "sucrose_twin.good", "",
                {"-cif","sucrose.cif",
                 "-hkl","olex2/Wfn_job/sucrose.hkl",
                 "-wfn","olex2/Wfn_job/sucrose.wfx",
                 "-acc","1",
                 "-twin","1","0","0","0","-1","0","0","1","-2",
                 "-all_charges", "-no_date"}
            });
        }

        TEST_METHOD(wfn_reading)
        {
            RunTest({
                "wfn_reading", "wfn_reading", "wfn_reading.good", "",
                {"-hkl","test.hkl", "-acc","1",
                 "-cif","test.cif", "-wfn","test.wfn",
                 "-all_charges", "-no_date"}
            });
        }

        // ---- Full tests (skipped unless RUN_FULL_TEST is set) ---------------

        TEST_METHOD(fchk_conversion)
        {
            RunTest({
                "fchk_conversion", "NiP3_fchk", "good.fchk", "",
                {"-b","dev2-TZVP", "-d","./", "-wfn","in.ffn",
                 "-all_charges", "-no_date"},
                /*full=*/true
            });
        }

        TEST_METHOD(fourier_transform_full)
        {
            RunTest({
                "fourier_transform_full", "SALTED", "fourier_transform_full.good", "",
                {"-test_analytical","full", "-acc","1",
                 "-all_charges", "-no_date"},
                /*full=*/true
            });
        }
    };
}
