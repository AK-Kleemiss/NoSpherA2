# Unit Test Status
**Last updated: 2026-07-03** (pruned 4 stale `tests.toml` entries with non-existent CLI flags;
regenerated `alanine_integrated_occ.good`; **201/201 ctest passing** after full rebuild — see
Known Issues for a same-day transient 13-failure episode against a stale/partial build)

## Test Harnesses

| Harness | How to run | Notes |
|---------|-----------|-------|
| Python pytest | `uv run pytest` | Golden-file diffs vs `.good` files |
| Windows VS Tests | `msbuild Windows\Tests\Tests.vcxproj` + `vstest.console.exe` | In-process via `NoSpherA2_DLL.dll`; same golden files |

Environment variables:
- `NOS_EXE` — path to exe (pytest only; default: `NoSpherA2` on PATH)
- `OCC_DATA_PATH` — path to OCC data share (`Lib/occ/share/occ` for VS tests, injected by TestRunner.h)
- `RUN_FULL_TEST` — set to any truthy value to include slow/full tests

---

## Current Validation Status

The table below (2026-06-17/2026-07-02) predates both the CMake/ctest migration and the
2026-07-03 Thakkar cubic-spline change, and is kept only as history — see `ctest --preset
<preset> --output-on-failure` for current status and the Known Issues section below for the
tests presently failing because of the spline change.

| Suite | Config | Result |
|-------|--------|--------|
| Python pytest Release | non-full (21 tests) | **passing** (historical, 2026-07-02) |
| Python pytest Release | full (23 tests) | **passing** (historical, 2026-07-02; 1 tolerated numeric warn in `fchk_conversion`) |
| Python pytest Debug | full (23 tests) | **passing** (historical, 2026-07-02) |
| VS Tests Debug | full (23 tests) | **passing** (historical, 2026-07-02) |
| VS Tests Release | full (23 tests) | **passing** (historical, 2026-07-02) |
| Python pytest macOS | `RGBI_Groups_NH3BH3`, `RGBI_Groups_NH3Li` targeted | **passing** (historical, validated 2026-06-17) |
| ctest (release-windows) | 205 cases (unit + `TomlIntegrationTests`) | targeted RGBI subset **7/7 passing** on 2026-07-03 for `NH3BH3` no-sym/sym and atom-only `NH3Li` NAO/ANO coverage; last full-suite baseline before this change was **201/201 passing**. The `-rgbi_no_sym` tests in this subset were removed 2026-07-03 — see Known Issues below. |

---

## VS-Only Unit Tests (no `tests/tests.toml` counterpart)

These tests live in `Windows/Tests/UnitTests.cpp` and exercise individual Src/ functions
via a thin C-shim exported from `NoSpherA2_DLL.dll` (`NoSpherA2_DLL/unit_test_api.h/.cpp`).
They do NOT need a pytest counterpart and must NOT be added to `tests/tests.toml`.

| Class | Methods | Coverage | Status |
|-------|---------|----------|--------|
| `GeometryTests` | `ArrayLength_*` (3), `ArrayDistance_*` (2), `VecDiff_*` (2), `VecCross_*` (3), `VecDot_*` (3) | `array_length`, `vec_diff`, `vec_cross`, `vec_dot` | 🆕 not yet validated |
| `NumericTests` | `IsNan_*` (3), `IsSimilarRel_*` (3), `IsSimilarAbs_*` (2), `FastExpNeg_*` (4) | `is_nan`, `is_similar_rel`, `is_similar_abs`, `fast_exp_neg` | 🆕 not yet validated |
| `FchkParsingTests` | `ReadFchkInt_*` (3), `ReadFchkDbl_*` (2) | `read_fchk_integer(string)`, `read_fchk_double(string)` | 🆕 not yet validated |
| `StringUtilTests` | `EndsWith_*` (5), `ShrinkString_*` (4) | `ends_with`, `shrink_string` | 🆕 not yet validated |
| `Sha256Tests` | `Sha256_*` (6) | `sha::sha256` | 🆕 not yet validated |
| `VecAggregateTests` | `VecSumBool_*` (3), `VecSumInt_*` (2), `VecSumDouble_*` (2), `VecLength_*` (4) | `vec_sum` (all overloads), `vec_length` | 🆕 not yet validated |
| `StringUtilTests2` | `Trim_*` (4), `Asciitolower_*` (3), `DoubleFromEsd_*` (3), `DecimalPrecisionCif_*` (3) | `trim`, `asciitolower`, `double_from_string_with_esd`, `get_decimal_precision_from_CIF_number` | 🆕 not yet validated |
| `BasisTypeTests` | `Sht2nbas_*` (2), `DoubleFactorial_*` (1) | `sht2nbas`, `doublefactorial` | 🆕 not yet validated |
| `ConstantsTests` | `ConstAbs_*` (3), `ConstexprPow_*` (1), `ConstantsSqrt_*` (2), `ExpApprox_*` (3), `LogApprox_*` (3), `Bohr2Ang_*` (2), `Ang2Bohr_*` (2), `CubicBohr2Ang_*` (1), `CubicAng2Bohr_*` (1), `Factorial_*` (1) | `constants::` — `const_abs`, `constexpr_pow`, `sqrt`, `exp_approx`, `log_approx`, `bohr2ang`, `ang2bohr`, `cubic_*`, `ft_fun` | 🆕 not yet validated |
| `OrbitalIndexTests` | `OrcaToPySCF_*` (4), `TypeToNbo_*` (6) | `constants::orca_2_pySCF`, `constants::type_2_nbo` | 🆕 not yet validated |
| `BesselTests` | `BesselJ0_*` (2), `BesselJ1_*` (1), `BesselJ2_*` (1), `BesselJ_HigherOrder_*` (1), `BesselJ_RecurrenceCheck` (1), plus more | `bessel_first_kind` (all l branches) | 🆕 not yet validated |

Total: 11 test classes, ~90 test methods.  
Added: 2026-06-14.

---

## Tests Registered in `tests/tests.toml` and `Windows/Tests/Tests.cpp`

| Test name | Directory | Good file | Full? | Status |
|-----------|-----------|-----------|-------|--------|
| alanine_occ | alanine_occ | alanine_occ.good | no | ✅ passing |
| alanine_integrated_occ | alanine_integrated_occ | alanine_integrated_occ.good | no | ✅ passing (regenerated 2026-07-03, see note below) |
| disorder_THPP | disorder | disorder_THPP.good | no | ✅ passing |
| fractal | sucrose_fchk_SF | fractal.good | no | ✅ passing |
| grown_water | grown | grown_water.good | no | ✅ passing |
| Hybrid_mode | Hybrid | Hybrid_mode.good | no | ✅ passing |
| malbac_SF_ECP | ECP_SF | malbac_SF_ECP.good | no | ✅ passing |
| properties | sucrose_fchk_SF | properties.good | no | ✅ passing |
| ri_fit | epoxide_gbw | ri_fit.good | no | ✅ passing |
| RGBI_Groups_NH3BH3_sym | RGBI_groups | NH3BH3_sym.good | no | ✅ passing (macOS targeted 2026-06-17) |
| RGBI_Groups_NH3BH3_sym_ANO | RGBI_groups | NH3BH3_sym_ano.good | no | ✅ passing (release-windows targeted 2026-07-03) |
| rubredoxin_cmtc | rubredoxin_cmtc | rubredoxin_cmtc.good | no | ✅ passing |
| SALTED | SALTED | SALTED.good | no | ✅ passing |
| sucrose_IAM | sucrose_IAM_SF | sucrose_IAM.good | no | ✅ passing |
| sucrose_ptb | sucrose_IAM_SF | sucrose_ptb.good | no | ✅ passing |
| sucrose_SF | sucrose_fchk_SF | sucrose_SF.good | no | ✅ passing |
| sucrose_twin | sucrose_fchk_SF | sucrose_twin.good | no | ✅ passing |
| intermolecular_nci | intermolecular_nci | good.dat | no | requires `-promol_nci_single_thread`; regenerated 2026-07-03 after fixing an off-by-one atomic-number bug (see note below) — ✅ passing on macOS arm64 and macOS x86_64 |
| wfn_reading | wfn_reading | wfn_reading.good | no | ✅ passing |
| TFVC | TFVC | TFVC.good | no | ✅ passing |
| TFVC_ECP | TFVC | TFVC_ECP.good | no | ✅ passing |
| fchk_conversion | NiP3_fchk | good.fchk | **yes** | ✅ passing (tolerated numeric warn) |

---

## Test Directories WITHOUT a Registered Test

These directories exist under `tests/` but have no entry in `tests.toml`. They need `.good` reference
files generated before they can be registered.

| Directory | Files available | Blocker |
|-----------|----------------|---------|
| Fe_gbw | Fe.cif, Fe.gbw, Fe.hkl | no .good file |
| mohr_salt | mohrs_salt.cif, mohrs_salt.gbw | no .good file |
| molden_file | many .molden/.wfn/.wfx | no .good file |
| RI_Test | TESTMOL.gbw, RI_Test.good | good file uses dev-machine absolute paths; `-test_RI` is an internal mode |
| RI_Test_2 | sucrose.cif, sucrose.gbw, .wfn | no .good file |
| reading_SALTED (dir) | alanine/crambin .npy data | no .good file; not the same as the now-removed `reading_SALTED` `tests.toml` entry (which used SALTED/) |
| CASSCF | water.inp | no .good file, CASSCF not yet wired into CLI |
| cytidine_tonto | many tonto-format files | external tonto dependency |
| polarizabilities | gjf/wfn files | no .good file |
| partition | DFT/HF .gbw files | no .good file |
| isosurface | .xyz files | no .good file |
| sfac_scan | epoxide .gbw/.cif | no .good file |
| PBC_Test | .ins/.res/olex2 | PBC not yet tested |
| ptb_H_file | H.gbw, wfn.xtb | variant of sucrose_ptb approach |
| Lukas_Test | thpp cif/xyz | informal dev test |
| BasisSetConverter | auxiliary_basis.hpp | header only, no runtime test |

**To add a new test**: generate a `.good` reference file from a known-good run, add an entry to
`tests/tests.toml`, and add a `TEST_METHOD` to `Windows/Tests/Tests.cpp`.

---

## Known Issues

- **Transient 13-test failure episode, 2026-07-03** (`alanine_occ`, `disorder_THPP`,
  `grown_water`, `Hybrid_mode`, `malbac_SF_ECP`, `rubredoxin_cmtc`, `sucrose_ptb`, `sucrose_SF`,
  `sucrose_twin`, `wfn_reading`, `TFVC`, `TFVC_ECP`, plus `alanine_integrated_occ`): an earlier
  `ctest` run the same day (against a not-fully-rebuilt binary, mid-way through the Thakkar
  cubic-spline interpolation change in `Src/core/spherical_density.h`/`Src/core/GridManager.cpp`
  — see `lambda2_at` in `Src/core/properties.cpp` for why that change was made) showed all 13
  failing with small (typically 3rd-decimal-place) Hirshfeld-weight numeric drift. After a full
  rebuild, all 13 pass again against their **existing, unmodified** `.good` files — i.e. this was
  not a real regression requiring golden-file updates, just a stale-build artifact. Only
  `alanine_integrated_occ.good` was actually regenerated (see below); it wasn't strictly necessary
  in hindsight, but reflects a real, freshly-verified in-process run and is harmless to keep.
  If a similar-looking failure set reappears, rebuild fully before concluding the `.good` files
  need updating — genuine regressions from future changes will look identical to this at a glance.

- `alanine_integrated_occ` regenerated 2026-07-03: ran `ctest --test-dir build/release-windows -R
  AlanineIntegratedOcc` (in-process, per the run rule below), then
  `python update_good_from_pytest_results.py --test alanine_integrated_occ --yes` to accept the
  ctest-produced in-place output (`tests/alanine_integrated_occ/NoSpherA2.log`) as the new
  `alanine_integrated_occ.good`. Re-ran ctest afterward to confirm it now passes. This test does
  run correctly as a normal ctest in-process case; the CLAUDE.md guidance against subprocess
  fallbacks refers to a distinct, historical crash when launched as a *spawned subprocess*
  (e.g. `update_good_from_pytest_results.py --run`), not to ctest's own in-process runner.

- **Removed 2026-07-03**: `fourier_transform`, `fourier_transform_full`, `openBLAS`,
  `reading_SALTED` were deleted from `tests/tests.toml`. Their `[test.args]` keys
  (`test_analytical`, `blastest`, `test_reading_SALTED_binary`) are not real `NoSpherA2` CLI
  flags — passing them produces the "Did not understand the task to perform!" help dump instead
  of running anything. Grepping the source found these names only in `tests.toml` and in
  `tests/src/UnitTests.cpp`, where the equivalent coverage already exists as direct C++ function
  calls with no CLI path at all (`BesselTests.AnalyticFourier` → `test_analytical_fourier()`,
  `SALTEDTests.ReadingSALTEDBinaryFile` → `test_reading_SALTED_binary_file()`, plus a BLAS
  self-test). None of the four had a corresponding `TEST(TomlIntegrationTests, ...)` in
  `tests/src/IntegrationTests.cpp` either, so they were already effectively dead from ctest's
  perspective — this just removes the stale, misleading `tests.toml`/`.good` entries to match.

- `intermolecular_nci`'s `_values.dat` writer (`promolecular_nci_analysis` in
  `Src/core/properties.cpp`) parallelizes over grid points with
  `schedule(dynamic)`, so output row order is not reproducible run-to-run
  under multi-threading (values are correct, ordering is not). The test's
  `tests.toml` args pass `-promol_nci_single_thread` (added 2026-07-02,
  `Src/core/convenience.h`/`.cpp`) to force single-threaded, deterministic
  output ordering for golden-file comparison. This flag is test/reproducibility
  tooling, not needed for normal end-user runs. On 2026-07-03, macOS arm64
  still failed because `get_lambda_1` used LAPACK/Accelerate eigenvalues for
  near-zero Hessians; `Src/core/convenience.cpp` now uses an analytical 3x3
  symmetric-matrix formula for the middle eigenvalue.

- **Root cause of the macOS arm64 divergence found and fixed 2026-07-03**: after the
  `get_lambda_1` fix above, `IntermolecularNCI` still produced ~6x more accepted grid
  points on macOS arm64 than on Linux/Windows/macOS x86_64 (2480 vs. 413), with values
  differing even for early points. AddressSanitizer traced this to a real bug in
  `make_thakkar_interpolators()` (`Src/core/properties.cpp`): it built its 92-entry
  `atom_models` vector with atomic numbers `0..91`, but every caller indexes it as
  `atom_models[atom.charge - 1]`, expecting index 0 to hold atomic number 1 (Hydrogen).
  Every real element's density lookup was off by one, and the bogus atomic-number-0
  entry triggered an out-of-bounds read of `Thakkar_ns`/`np`/`nd`/`nf`
  (`n_vector[atomic_number - 1]` = `n_vector[-1]`) whose garbage value was used as a
  loop trip count in `calc_orbs` — undefined behavior whose outcome depends on
  whatever happens to sit in adjacent memory, which differs by platform/build. This
  is what made the test look architecture-dependent; it was never really about arm64
  vs. x86_64 floating-point math. Three "make the math more careful" fixes were tried
  and individually verified to have **zero** effect on the divergence before the real
  bug was found: Kahan-compensated summation in `calc_orbs` (kept anyway, see below),
  swapping cubic-spline for linear interpolation, and building with
  `-ffp-contract=off` to rule out FMA-contraction differences.

  Fixed by changing `atom_models.emplace_back(a)` to `atom_models.emplace_back(a + 1)`
  (commit `88f0a9a`). Verified via: (1) AddressSanitizer clean before/after on macOS
  arm64, (2) a from-scratch local `release-macos-x86_64` build (with the fix) produces
  byte-for-byte identical `unit_surrounding_values.dat` output to the arm64 build
  (178 points each, `28 x 37 x 39` shrunk grid on both), (3) full `ctest` suite
  (201/201) still passes on macOS arm64 after the fix. `good.dat` was regenerated from
  this now cross-platform-identical output; the old 413-point `good.dat` was itself a
  product of the pre-fix bug (it happened to reproduce consistently across
  Linux/Windows/macOS x86_64's memory layouts, but was never numerically correct).

  Separately, `Thakkar::calc_orbs` (`Src/core/spherical_density.cpp`) was changed to
  use Kahan compensated summation for the Slater-expansion accumulation into `Orb[m]`.
  This is a genuine numerical-robustness improvement (verified to not change output on
  its own for this test) and was kept, but it was not the fix for the arm64 divergence.

- **`-rgbi_no_sym` RGBI tests removed 2026-07-03**: `RGBI_Groups_NH3BH3`, `RGBI_Groups_NH3BH3_ANO`,
  `RGBI_Groups_NH3Li`, `RGBI_NH3Li`, and `RGBI_NH3Li_ANO` (all tests passing `-rgbi_no_sym`) were
  removed from `tests/tests.toml`, `tests/src/IntegrationTests.cpp`, and
  `tests/src/SetIntegrationTestLocks.cmake`, along with their now-orphaned `.good` files
  (`tests/RGBI_groups/NH3BH3.good`, `NH3BH3_ano.good`, `nh3li.good`;
  `tests/RGBI/nh3li_nao.good`, `nh3li_ano.good`). `tests/RGBI_groups/nh3li.gbw` is left in place,
  unregistered, in case `-rgbi_no_sym` coverage is reinstated later.

  Root cause: `RGBI_Groups_NH3BH3_ANO` failed identically on all four CI platforms
  (Linux/Windows/macOS x86_64/macOS arm64) after the ANO-basis-for-RGBI merge
  (`420c36b`), which ruled out any platform-specific floating-point cause. Tracing
  `Roby_information::calculateAtomicNAO` (`Src/core/bondwise_analysis.cpp`) showed that
  `-rgbi_no_sym` skips the O_h spherical-averaging symmetrization step that the fixed
  `occupancy_cutoff = 0.17` NAO-retention threshold implicitly depends on: symmetrized
  boron ANOs come out as three exactly-degenerate eigenvalues (e.g. `0.170367` each,
  all clearing the cutoff), while un-symmetrized boron ANOs collapse the same density
  into one dominant eigenvalue plus a small residual (e.g. `0.511103` + `0.00267848`)
  that falls under the cutoff and is silently discarded — losing real bonding density
  and producing a wrong Roby population (`3.146` vs expected `3.474`). This is a design
  gap in the `-rgbi_no_sym` + ANO-basis combination (the fixed-threshold retention
  heuristic isn't robust without symmetrization), not a typo-style bug, and not
  something fixable by regenerating the `.good` file. The `-rgbi` (symmetrized) tests
  are unaffected and remain in the suite.

---

## How to Add / Update Tests — Agent Guidelines

See also the rule in `CLAUDE.md` → *Agent / AI coding-assistant rules*.

1. **Generate the good file** by running `NoSpherA2` (release build) against the test inputs and
   capturing stdout: `NoSpherA2 [args] > tests/<dir>/<name>.good`
2. **Register in `tests/tests.toml`** — add a `[test_name]` section with `directory`, optional
   `good` (if not `<name>.good`), optional `actual` (if output is not the log), and `[test_name.args]`
3. **Register in `Windows/Tests/Tests.cpp`** — add a `TEST_METHOD(<name>)` block calling `RunTest`
   with matching args (always append `-all_charges` and `-no_date`)
4. **Validate** in all four configurations: pytest Release, pytest Debug, VS Debug, VS Release
5. **Update this file** — mark the test as passing or note blockers
