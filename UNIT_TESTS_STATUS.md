# Unit Test Status
**Last updated: 2026-06-14**

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

## Current Validation Status (2026-06-14)

| Suite | Config | Result |
|-------|--------|--------|
| Python pytest Release | non-full (19 tests) | **passing** |
| Python pytest Release | full (21 tests) | **passing** (1 tolerated numeric warn in `fchk_conversion`) |
| Python pytest Debug | full (21 tests) | **passing** |
| VS Tests Debug | full (21 tests) | **passing** |
| VS Tests Release | full (21 tests) | **passing** |

---

## VS-Only Unit Tests (no `tests/tests.toml` counterpart)

These tests live in `Windows/Tests/UnitTests.cpp` and exercise individual Src/ functions
via a thin C-shim exported from `NoSpherA2_DLL.dll` (`NoSpherA2_DLL/unit_test_api.h/.cpp`).
They do NOT need a pytest counterpart and must NOT be added to `tests/tests.toml`.

| Class | Methods | Coverage | Status |
|-------|---------|----------|--------|
| `GeometryTests` | `ArrayLength_*` (3), `ArrayDistance_*` (2), `VecDiff_*` (2), `VecCross_*` (3), `VecDot_*` (3) | `array_length`, `vec_diff`, `vec_cross`, `vec_dot` in `convenience.cpp` | 🆕 not yet validated |
| `NumericTests` | `IsNan_*` (3), `IsSimilarRel_*` (3), `IsSimilarAbs_*` (2), `FastExpNeg_*` (4) | `is_nan`, `is_similar_rel`, `is_similar_abs`, `fast_exp_neg` | 🆕 not yet validated |
| `FchkParsingTests` | `ReadFchkInt_*` (3), `ReadFchkDbl_*` (2) | `read_fchk_integer(string)`, `read_fchk_double(string)` in `fchk.cpp` | 🆕 not yet validated |
| `StringUtilTests` | `EndsWith_*` (5), `ShrinkString_*` (4) | `ends_with`, `shrink_string` in `convenience.cpp` | 🆕 not yet validated |
| `Sha256Tests` | `Sha256_*` (6) | `sha::sha256` in `convenience.cpp` | 🆕 not yet validated |

Total: 5 test classes, 34 test methods.  
Added: 2026-06-14.

---

## Tests Registered in `tests/tests.toml` and `Windows/Tests/Tests.cpp`

| Test name | Directory | Good file | Full? | Status |
|-----------|-----------|-----------|-------|--------|
| alanine_occ | alanine_occ | alanine_occ.good | no | ✅ passing |
| alanine_integrated_occ | alanine_integrated_occ | alanine_integrated_occ.good | no | ⚠️ **crashing** — see INVESTIGATION_STATUS.md |
| disorder_THPP | disorder | disorder_THPP.good | no | ✅ passing |
| fourier_transform | SALTED | fourier_transform.good | no | ✅ passing |
| fractal | sucrose_fchk_SF | fractal.good | no | ✅ passing |
| grown_water | grown | grown_water.good | no | ✅ passing |
| Hybrid_mode | Hybrid | Hybrid_mode.good | no | ✅ passing |
| malbac_SF_ECP | ECP_SF | malbac_SF_ECP.good | no | ✅ passing |
| openBLAS | OpenBLAS | openBLAS.good | no | ✅ passing |
| properties | sucrose_fchk_SF | properties.good | no | ✅ passing |
| reading_SALTED | SALTED | reading_SALTED.good | no | ✅ passing |
| ri_fit | epoxide_gbw | ri_fit.good | no | ✅ passing |
| rubredoxin_cmtc | rubredoxin_cmtc | rubredoxin_cmtc.good | no | ✅ passing |
| SALTED | SALTED | SALTED.good | no | ✅ passing |
| sucrose_IAM | sucrose_IAM_SF | sucrose_IAM.good | no | ✅ passing |
| sucrose_ptb | sucrose_IAM_SF | sucrose_ptb.good | no | ✅ passing |
| sucrose_SF | sucrose_fchk_SF | sucrose_SF.good | no | ✅ passing |
| sucrose_twin | sucrose_fchk_SF | sucrose_twin.good | no | ✅ passing |
| wfn_reading | wfn_reading | wfn_reading.good | no | ✅ passing |
| **TFVC** | TFVC | TFVC.good | no | 🆕 **added 2026-06-14** — not yet validated |
| **TFVC_ECP** | TFVC | TFVC_ECP.good | no | 🆕 **added 2026-06-14** — not yet validated |
| fchk_conversion | NiP3_fchk | good.fchk | **yes** | ✅ passing (tolerated numeric warn) |
| fourier_transform_full | SALTED | fourier_transform_full.good | **yes** | ✅ passing |

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
| reading_SALTED (dir) | alanine/crambin .npy data | no .good file; different from `reading_SALTED` test (which uses SALTED/) |
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

### `alanine_integrated_occ` — heap corruption crash (all platforms)
- The standalone EXE crashes for **all** OCC density-fitting jobs (`-occ` flag with DF basis)
- The VS in-process DLL also crashes
- Root cause: heap corruption in OCC's integral engine (`~IntegralEngineDF`)
- Exit code: `0xC0000374` (STATUS_HEAP_CORRUPTION via `RtlFailFast` — uncatchable by SEH)
- The output fchk file is **never written** because the crash happens before destructors finish
- Full details: see `INVESTIGATION_STATUS.md`

### TFVC / TFVC_ECP — newly added, not yet validated
- Good files exist and were generated from prior runs, but these tests have never been run through
  the VS test harness
- Run against both Debug and Release before marking as passing

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
