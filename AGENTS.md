# AGENTS.md

This file provides guidance to Codex and other AI coding agents when working in this repository.

## What is NoSpherA2

NoSpherA2 is a C++ CLI tool that generates `.tsc` non-spherical scattering-factor files and property grids from quantum-chemistry wavefunction files. It is used with Olex2/olex2.refine for Hirshfeld Atom Refinement (HAR) in crystallography.

## Current Build System

The active build system is CMake/Ninja with presets. The top-level CMake project owns dependency configuration for `libcint`, `occ`, `mdspan`, `featomic`, Intel MKL or Apple Accelerate/OpenMP, the basis-set converter, the core library, the CLI app, and optional tests/DLL targets.

Bootstrap the local micromamba toolchain before configuring a fresh checkout:

```sh
cmake -P scripts/BootstrapMicromamba.cmake
```

Use the preset names that exist in `CMakePresets.json`; do not use old names such as `windows-msvc-release-full`, `linux-gcc`, or `macos-release-full-arm64`.

Common configure/build presets:

```sh
# Windows
cmake --preset release-windows
cmake --build --preset release-windows

cmake --preset debug-windows
cmake --build --preset debug-windows

# Linux
cmake --preset release-linux
cmake --build --preset release-linux

# macOS per architecture
cmake --preset release-macos-arm64
cmake --build --preset release-macos-arm64

cmake --preset release-macos-x86_64
cmake --build --preset release-macos-x86_64
```

The executable is written to the preset build tree under `bin`, for example `build/release-windows/bin/NoSpherA2.exe`.

Important CMake options:

- `NOSPHERA2_BUILD_TESTS=ON` builds the C++ GTest/CTest suite in `tests/src`.
- `NOSPHERA2_BUILD_DLL=ON` builds the optional Windows DLL target from `Src/dll`.
- `NOSPHERA2_DEPENDENCIES_ONLY=ON` configures/builds dependency targets and then returns before building NoSpherA2. CI uses this to cache dependency build trees.

CI configures with `NOSPHERA2_BUILD_TESTS=ON` and `NOSPHERA2_DEPENDENCIES_ONLY=OFF` after restoring or creating a dependency-only cache.

## Running Tests

Preferred CMake/CTest flow:

```sh
cmake --preset <preset> -DNOSPHERA2_BUILD_TESTS=ON
cmake --build --preset <preset>
ctest --preset <preset> --output-on-failure
```

Python golden-file tests are still available:

```sh
uv run pytest
```

Useful environment variables:

- `NOS_EXE` - path to the executable for Python tests.
- `OCC_DATA_PATH` - OCC runtime data directory. CTest defaults this to `occ/share`; local Windows Python runs commonly need `D:\git\NoSpherA2\Lib\occ\share\occ` or the current equivalent.
- `RUN_FULL_TEST` - truthy value enables slow/full Python integration tests.

Tests are golden-file comparisons: the harness runs NoSpherA2 with args from `tests/tests.toml`, captures stdout, and diffs against `.good` reference files in each test subdirectory.

### Windows Visual Studio Tests

The CMake flow is the default for new validation. The legacy Visual Studio solution and test project still exist under `Windows/` and remain useful for Windows/OCC integration debugging:

```ps
msbuild Windows\Tests\Tests.vcxproj /m /p:Configuration=Release /p:Platform=x64 /p:SolutionDir=D:\git\NoSpherA2\Windows\
vstest.console.exe Windows\Tests\Release_x64\Tests.exe /Platform:x64
```

Visual Studio tests must run in-process through `NoSpherA2_DLL.dll`; do not add subprocess fallbacks for failing VS tests, including `alanine_integrated_occ`. If all tests complete in under a second, the test binary or working directory is probably wrong.

The VS solution's `NoSpherA2`/`NoSpherA2_LIB` projects use a shared `OutDir`
(`Windows\Windows_utils\NoSpherA2_universal.props`) that writes into the same top-level `build/`
tree the CMake presets use, but under `build\$(Configuration)_$(Platform)\` (e.g.
`build\Release_x64\NoSpherA2.exe`) rather than the CMake preset's
`build\release-windows\bin\NoSpherA2.exe`. `Windows\Tests\Tests.vcxproj` was not repointed and
still builds to its own `Windows\Tests\$(Configuration)_$(Platform)\`.

## Adding or Changing Tests

Register Python/golden tests in `tests/tests.toml`:

```toml
[my_test_name]
directory = "my_test_dir"
# good = "custom.good"

[my_test_name.args]
cif = "structure.cif"
hkl = "data.hkl"
wfn = "wavefunction.wfx"
acc = 0
```

CLI args are always under `<testname>.args`; the test framework passes these as `--key value`.

When touching tests:

- Keep `tests/tests.toml`, `tests/src`, and `Windows/Tests/Tests.cpp` in sync where the same behavior is covered by more than one harness.
- Update `UNIT_TESTS_STATUS.md` whenever adding, removing, modifying, or investigating tests.
- Update the validation-status block in this file and `CLAUDE.md` if the overall validated pass count changes.

## Code Style

This codebase is terse, dense and comment-sparse. There is no formatter, so style
is enforced by imitation: **match the file you are editing**, not a global rule.
Code that is correct but laid out differently has been rejected by reviewers here.

Before committing any C++ change, re-read the staged diff for layout alone:

- **No blank lines inserted inside a function body to group statements.** Core
  files run 0-9% blank lines; `fchk.cpp` is 0%.
- **Comment lines should be 0-13% of added lines**, the baseline band for this
  repo. Rationale belongs in the commit message, not the source. Comments here
  explain mathematics, not decisions.
- **No markdown inside C++ comments** (`**bold**`, bullets, ASCII tables) and no
  benchmark results parked above a default value or function.
- **No defensive validation that has never fired.** The house guard style is one
  or two compact `err_checkf` statements at the top of the function.
- **Index-`for` with `int` and a short name** (`i`, `j`, `a`, `s`, `x`/`y`/`z`).
  2155 of 2620 loops. Do not introduce range-`for` or `std::accumulate` into
  index-`for` code.
- **Do not brace single-statement bodies** in files that leave them off; 40% of
  `if` and 24% of `for` bodies are braceless.
- **Use the `convenience.h` typedefs** (`vec`, `vec2`, `ivec`, `cvec`,
  `hkl_list`, `dMatrix*`), never `std::vector<double>`.
- **`#pragma omp` at column 0**, always, regardless of loop nesting depth.
- **Preserve whitespace convention per file.** `scattering_factors.cpp`,
  `XCW.cpp`, `basis_set.cpp`, `convenience.h`, `XCW.h` and `cell.cpp` are
  tab-indented; most others use 4 spaces.
- **Do not reformat untouched code**, delete commented-out code, or "fix" the
  repo-wide splits in brace style, `&` placement, `NULL` vs `nullptr` or function
  naming. Those are genuinely mixed and no migration is in progress.

Keep commit-message bodies short and prose-shaped; a one-line subject naming the
change is the norm here.

## Source Layout

```text
app/                 CLI executable target (`NoSpherA2`)
Src/core/            Core C++ library target (`NoSpherA2Core`, alias `NoSpherA2::Core`)
Src/dll/             Optional DLL target when `NOSPHERA2_BUILD_DLL=ON`
tests/src/           C++ GTest/CTest integration and unit tests
tests/               Python pytest harness and golden-file test cases
BasisSetGenerator/   Basis-set converter used to generate `Src/basis_data.cpp`
occ/                 OCC submodule, the main quantum-chemistry dependency
featomic/            Rust/C++ submodule for ML feature generation via metatensor
libcint/             C submodule for electron integral evaluation
mdspan/              C++ multidimensional-array reference implementation
cmake/               Project CMake modules
scripts/             Micromamba/bootstrap/Visual Studio helper scripts
Windows/             Legacy Visual Studio solution/projects and Windows test harness
Linux/               Linux packaging/build helper scripts
```

Key core modules live in `Src/core`:

- I/O: `fchk.cpp/h`, `cif.cpp/h`, `wfn_class.cpp/h`
- Molecular representation: `atoms.cpp/h`, `molecule.cpp/h`, `basis_set.cpp/h`
- Integration engine: `AtomGrid.cpp/h`, `GridManager.cpp/h`, `integrator.cpp/h`, `integration_params.cpp/h`
- Math/integrals: `nos_math.cpp/h`, `libCintMain.cpp/h`, `libCintKernels.cpp/h`
- Analysis/output: `scattering_factors.cpp/h`, `bondwise_analysis.cpp/h`, `isosurface.cpp/h`, `properties.cpp/h`, `cube.cpp/h`, `tsc_block.h`
- ML density: `SALTED_*.cpp/h`

## Dependency and Linking Notes

- The top-level project consumes bundled dependencies with `add_subdirectory`.
- `libcint` is patched parent-side for MSVC `ivdep` pragmas by generating rewritten sources under the build tree; do not edit the submodule for that workaround.
- `OCC` is the heaviest dependency. Check whether OCC already provides needed quantum-chemistry functionality before implementing new logic in `Src/core`.
- Most dependencies are linked statically.
- oneTBB is intentionally deployed dynamically. Expected runtime files beside packaged executables include `tbb12.dll`, `libtbb.so*`, or `libtbb.dylib`.
- Do not reintroduce `tbbmalloc_proxy` or `tbbmalloc` linking for NoSpherA2.

## GPU and CPU Must Agree Numerically

Every GPU path exists only as a faster way to compute what the CPU path already
computes. **A GPU path is not finished until its agreement with the CPU has been
measured on a real workload and written down.** Speed is not a result on its own; a
kernel that is fast and wrong is worse than no kernel, because it is harder to notice.

Rules:

- **Measure agreement, do not assume it.** Run the same input through the CPU and GPU
  paths and compare the physics the code exists to produce - scattering factors, the XCW
  convergence table, predicted coefficients - not an intermediate you happen to have.
- **Quote the wR2 shape and the maximum absolute difference.** The maximum *relative*
  difference is not a stable statistic: it is set by whatever the smallest near-zero
  component in the sample happens to be, and re-parsing the same data moved it three
  orders of magnitude once already.
- **Record the number beside the code and in the commit.** A tolerance nobody wrote down
  cannot be checked later, and the next reader cannot tell a regression from intent.
- **Check the physics before the clock.** The first I tensor kernel ran 11x faster and
  was completely wrong (GooF 69.1 against 4.665). Diff the science first, every time.
- **A path the hardware never selects is a path no test covers.** Give every branch a
  switch that forces it - `-gpu_fp64`, `-gpu_fp32`, `NOSPHERA2_GPU_BATCH`,
  `NOSPHERA2_GRID_LOCAL`, `NOSPHERA2_BLAS_GPU_MIN_FLOP` - or it rots unnoticed. Two paths
  had already rotted when this rule was written down: `sucrose_SF_gpu_grid` passed for a
  session while its kernel failed and fell back to the CPU, because a silent fallback
  reproduces the CPU reference exactly, and `blas_gpu_dgemm` shipped reached by no test at
  all because its size gate is above anything in the test data.

### Use `-gflops` before optimising anything

`-gflops` reports achieved GFLOP/s per stage for both the CPU and GPU paths, the host
thread count, and the slowest single call in each row. Each of those three exists because
something was got wrong without it:

- **The rate** stops a threshold calibrated on one card being trusted on another. A V100
  runs fp64 at half rate; this laptop's card at a sixty-fourth.
- **The thread count** because a CPU row that looks like slow hardware is usually a thread
  count of one. Every Linux build was single-threaded until 29 Aug 2026 and nobody noticed,
  because an ignored `#pragma omp` is legal and warns about nothing.
- **The slowest call** because creating the device context costs ~95 ms and lands inside
  whichever call touches the GPU first. That has twice been mistaken for a slow kernel -
  once reported as 63 ms of transfers that were really 2.4, and once as SALTED losing to
  the host 2.5x when per steady-state call it was winning 4.5x.

Measure before and after, in the same sitting, best of three. A speedup is a ratio, and the
denominator deserves as much suspicion as the numerator: the "20x faster XCW on the
cluster" that stood for several hours was a correct GPU measurement divided by a CPU
baseline that was accidentally serial. The honest figure was 1.68x.

Where a kernel deliberately trades precision for speed, the trade must be stated with its
measured cost, not left implicit. The current paths, all measured on real workloads:

| path | precision | agreement with the CPU |
|---|---|---|
| `calc_SF`, `-gpu_fp64` | fp64 throughout | wR2 1.23e-14 |
| `calc_SF`, default | fp64 phase, fp32 transcendental and sum | wR2 1.84e-8 |
| I tensor, `-gpu_itensor` | fp32 GEMM | 2.1e-11 worst over the lambda scan |
| SALTED, `-gpu_salted` | fp64 throughout | wR2 4.95e-13 |

The two fp32 paths are the deliberate trades. Both were measured before being kept, both
sit orders below the experimental uncertainty of the data being fitted, and the fp64
alternative is one flag away in each case.

## Known Pitfalls

- Use `cmake --list-presets=all` or inspect `CMakePresets.json` before assuming preset names. Host-conditional configure presets may hide non-host presets from a plain `cmake --list-presets`.
- On Windows, use an x64 Visual Studio Developer environment before configuring CMake so Ninja can find `cl.exe` and the MSVC standard headers.
- If MSVC builds fail on missing headers such as `stdlib.h`, `string.h`, or `stdio.h`, treat that as a toolchain-shell setup issue before suspecting CMake source rewriting.
- Golden-file float parsing must support scientific notation and multi-value lines without weakening existing thresholds.
- Clone with `--recursive`; if submodules are missing, run `git submodule update --init --recursive`.
- OCC source fixes must be committed inside the `occ` submodule first, then the parent repository must commit the updated submodule pointer.
- When local tests pass but CI fails, compare `.github/workflows/c-cpp_all.yml` step-by-step with local commands, especially `NOS_EXE`, `OCC_DATA_PATH`, and CMake options.
- New libcint parallel call sites should pre-allocate per-thread scratch buffers instead of passing `cache=nullptr` inside TBB parallel regions. Use the existing `three_center_max_cache_size<kind>` and `tbb::enumerable_thread_specific` pattern in OCC.

## Current Validation Notes

As of 2026-09-02, `ctest --preset release-windows` reports **275/275 passing, 0 failed**
(253 s; 6 not run: the four `full = true` XCW cases and two disabled `DeltaSeriesTests`).
This baseline covers the pTB cartesian-f fix in `WFN::read_ptb` and the new
`-no_date_but_gpu` flag. On a machine with a CUDA device, 16 golden-file cases had been
failing on GPU notes absent from the references; those notes now follow `-no-date`, with
`-no_date_but_gpu` for `sucrose_SF_gpu_grid`, whose reference must keep the note. See
`UNIT_TESTS_STATUS.md`.

As of 2026-08-25, `ctest --preset release-windows` reports **258/258 passing, 0 failed**
(569 s; 5 skipped, all pre-existing: the four `full = true` XCW cases and the optional
`Nbo47.EpoxideGennboMatchesReferenceWhenAvailable` fixture). This baseline includes the
new `-fukui` feature: 7 `FukuiTests` unit cases and the `TomlIntegrationTests.Fukui` /
`TomlIntegrationTests.FukuiPBC` golden-file cases. It also confirms that fixing
`WFN::read_fchk` to record `path` — which changes fchk-driven cube filenames from a
stem-less `_rho.cube` to `<stem>_rho.cube` — regressed nothing. See
`UNIT_TESTS_STATUS.md` for the Fukui traps and the grid-convergence validation.

As of 2026-07-18, `TomlIntegrationTests.P1_test_XCW` (in-process, `release-windows`,
`OMP_NUM_THREADS=20`) passes after fixing a real access violation: `GridManager::calculateMBISWeights`/
`calculateEMBISWeights` ignored `config_.no_density_eval` and forced an invalid WFN density
evaluation on XCW's MO-pruned dummy wavefunction. A separate `exit(0)` in
`XCW::run_XCW_fitting()` was also removed — it silently killed the in-process test binary before
GoogleTest's own pass/fail assertion ran, so earlier "passing" runs of this test were a false
positive rather than a real pass. See `UNIT_TESTS_STATUS.md` Known Issues for the full
root-cause writeup, including a still-open follow-up: any XCW orbital basis other than the
hardcoded default `def2-svp` (selectable via the newly wired `-b` flag) crashes with a separate,
unrelated access violation in OCC's SOAD initial guess. This was not re-run against the full
`ctest` suite this session; the last full-suite baseline remains the 2026-07-03 note below.

As of 2026-07-03, `ctest --preset release-windows` reports **201/201 passing** after a full
rebuild (following the Thakkar cubic-spline interpolation change in
`Src/core/spherical_density.h` / `Src/core/GridManager.cpp`, Hirshfeld-weight density
partitioning). A same-day run against a not-fully-rebuilt binary showed 13 `TomlIntegrationTests`
cases failing with small numeric drift; after rebuilding, all 13 pass against their existing
`.good` files unchanged — see `UNIT_TESTS_STATUS.md` Known Issues. Only
`alanine_integrated_occ.good` was actually regenerated. Treat this as historical status unless
you have rerun the current CMake presets in this checkout.

As of 2026-07-03, `ctest --preset release-macos-arm64` reports **202/202 passing**, including
`TomlIntegrationTests.IntermolecularNCI`, after fixing an off-by-one atomic-number indexing bug
in `make_thakkar_interpolators()` (`Src/core/properties.cpp`, commit `88f0a9a`) that had been
causing undefined-behavior-driven divergence between macOS arm64 and other platforms. Verified
identical (byte-for-byte) `unit_surrounding_values.dat` output between a from-scratch local
`release-macos-x86_64` build and the arm64 build. Separately, the 5 `-rgbi_no_sym` RGBI tests
(merged in via `420c36b`) were removed after `RGBI_Groups_NH3BH3_ANO` was found to fail
identically on all four CI platforms — a real design gap in the no-symmetrization + ANO-basis
combination, not a platform or golden-file issue. `RGBI_NH3Li`/`RGBI_NH3Li_ANO` were then
reinstated using `-rgbi` (symmetrized) instead. See `UNIT_TESTS_STATUS.md` Known Issues for
both root-cause writeups.

The June-2026 note below (21/21, VS/pytest harnesses) predates the CMake/ctest migration and is
kept only as history.

Known fixed areas from June 2026:

1. `alanine_integrated_occ` passes in-process in Release and Debug.
2. `IntegralEngineDF::~IntegralEngineDF()` clears Eigen buffers before libcint engines are destroyed.
3. oneTBB is deployed shared; `tbbmalloc_proxy` is not linked.
4. MSBuild instruction-set selection follows `NOS_AVX` to keep Eigen ABI consistent.
5. Parallel three-center libcint scratch buffers avoid heap corruption under vstest.
6. `fchk_conversion` writes and compares `log.fchk` and handles restricted/beta coefficients, CMO bounds, virtual orbitals, and G shells.
7. Debug CRT assertion dialogs are suppressed in the DLL test host so failures print instead of blocking invisible UI.

See `UNIT_TESTS_STATUS.md` for test inventory details, but prefer live CMake files and CI workflow files for current build-system facts.
