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

As of 2026-07-03, `ctest --preset release-windows` reports **201/201 passing** after a full
rebuild (following the Thakkar cubic-spline interpolation change in
`Src/core/spherical_density.h` / `Src/core/GridManager.cpp`, Hirshfeld-weight density
partitioning). A same-day run against a not-fully-rebuilt binary showed 13 `TomlIntegrationTests`
cases failing with small numeric drift; after rebuilding, all 13 pass against their existing
`.good` files unchanged — see `UNIT_TESTS_STATUS.md` Known Issues. Only
`alanine_integrated_occ.good` was actually regenerated. The June-2026 note below (21/21,
VS/pytest harnesses) predates the CMake/ctest migration and is kept only as history. Treat this
as historical status unless you have rerun the current CMake presets in this checkout.

`TomlIntegrationTests.IntermolecularNCI` was failing on macOS arm64 CI; fixed 2026-07-03 by
replacing the LAPACK-based `get_lambda_1` in `Src/core/convenience.cpp` with an analytical
trigonometric Cardano formula for the middle eigenvalue of a 3×3 real symmetric matrix. See
`UNIT_TESTS_STATUS.md` Known Issues for details.

Known fixed areas from June 2026:

1. `alanine_integrated_occ` passes in-process in Release and Debug.
2. `IntegralEngineDF::~IntegralEngineDF()` clears Eigen buffers before libcint engines are destroyed.
3. oneTBB is deployed shared; `tbbmalloc_proxy` is not linked.
4. MSBuild instruction-set selection follows `NOS_AVX` to keep Eigen ABI consistent.
5. Parallel three-center libcint scratch buffers avoid heap corruption under vstest.
6. `fchk_conversion` writes and compares `log.fchk` and handles restricted/beta coefficients, CMO bounds, virtual orbitals, and G shells.
7. Debug CRT assertion dialogs are suppressed in the DLL test host so failures print instead of blocking invisible UI.

See `UNIT_TESTS_STATUS.md` for test inventory details, but prefer live CMake files and CI workflow files for current build-system facts.
