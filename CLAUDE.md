# CLAUDE.md

This file provides guidance to Claude Code and other AI coding assistants when working in this repository.

## What is NoSpherA2

NoSpherA2 is a C++ CLI tool that generates `.tsc` non-spherical scattering-factor files and property grids from quantum-chemistry wavefunction files. It is used with Olex2/olex2.refine for Hirshfeld Atom Refinement (HAR) in crystallography.

## Current Build System

The active build system is CMake/Ninja with presets. The top-level CMake project owns dependency configuration for `libcint`, `occ`, `mdspan`, `featomic`, Intel MKL or Apple Accelerate/OpenMP, the basis-set converter, the core library, the CLI app, and optional tests/DLL targets.

Bootstrap the local micromamba toolchain before configuring a fresh checkout:

```sh
cmake -P scripts/BootstrapMicromamba.cmake
```

Use the preset names that exist in `CMakePresets.json`; do not use old names such as `windows-msvc-release-full`, `linux-gcc`, `macos-release-full-arm64`, or `windows-msvc-debug`.

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

### Windows Agent/CLI Notes

For CMake preset builds, run from an x64 Visual Studio Developer PowerShell or initialize the MSVC environment before invoking CMake/Ninja. If Ninja can find `cl.exe` but compilation fails on missing standard headers such as `stdlib.h`, the shell/toolchain environment is incomplete.

The legacy Visual Studio solution still exists for Windows/OCC integration debugging. If using MSBuild from a normal PowerShell session, pass the full `MSBuild.exe` path and preserve `/p:` switches:

```powershell
$msbuild = "C:\Program Files\Microsoft Visual Studio\18\Community\MSBuild\Current\Bin\amd64\MSBuild.exe"
$solDir = "D:\git\NoSpherA2\Windows\"

& $msbuild "D:\git\NoSpherA2\Windows\Tests\Tests.vcxproj" `
    /p:Configuration=Release /p:Platform=x64 "/p:SolutionDir=$solDir" /v:minimal
```

Do not invoke MSBuild through Bash for this project; Bash can strip leading `/` from `/p:` switches.

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

The CMake flow is the default for new validation. The Visual Studio test project remains useful when debugging Windows DLL/OCC behavior:

```powershell
msbuild Windows\Tests\Tests.vcxproj /m /p:Configuration=Release /p:Platform=x64 /p:SolutionDir=D:\git\NoSpherA2\Windows\
vstest.console.exe Windows\x64\Release\Tests.dll /Platform:x64
```

Visual Studio tests must run in-process through `NoSpherA2_DLL.dll`; do not add subprocess fallbacks for failing VS tests, including `alanine_integrated_occ`. If all tests complete in under a second, the test binary or working directory is probably wrong.

For coverage, keep `coverage.runsettings` aligned with the current DLL/test output layout before trusting the report.

## Agent Rules

### Unit-test documentation is mandatory

Whenever an agent adds, removes, modifies, or investigates a test in `tests/tests.toml`, `tests/src`, `Windows/Tests`, `tests/run_test.py`, `TestRunner.h`, or any `.good` reference file, it must:

1. Update `UNIT_TESTS_STATUS.md`.
2. Update the validation-status block in this file and `AGENTS.md` if the overall pass count changes.
3. Keep the file commitworthy; do not leave unit-test documentation stale after test work.

### Keep test harnesses aligned

`tests/tests.toml`, `tests/src`, and `Windows/Tests/Tests.cpp` should stay in sync where they cover the same workflows. Check all relevant harnesses when adding or removing tests.

### Respect submodule boundaries

Do not casually edit vendored submodules. For local parent-side workarounds, prefer generated-source or CMake target augmentation in the parent project. If an OCC source fix is required, commit it inside `occ` first, then update the parent submodule pointer.

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
- OCC is the heaviest dependency. Check whether OCC already provides needed quantum-chemistry functionality before implementing new logic in `Src/core`.
- Most dependencies are linked statically.
- oneTBB is intentionally deployed dynamically. Expected runtime files beside packaged executables include `tbb12.dll`, `libtbb.so*`, or `libtbb.dylib`.
- Do not reintroduce `tbbmalloc_proxy` or `tbbmalloc` linking for NoSpherA2.

## Known Pitfalls

- Use `cmake --list-presets=all` or inspect `CMakePresets.json` before assuming preset names. Host-conditional configure presets may hide non-host presets from a plain `cmake --list-presets`.
- On Windows, use an x64 Visual Studio Developer environment before configuring CMake so Ninja can find `cl.exe` and the MSVC standard headers.
- Golden-file float parsing must support scientific notation and multi-value lines without weakening existing thresholds.
- Clone with `--recursive`; if submodules are missing, run `git submodule update --init --recursive`.
- CI and local runs can differ because CI uses a dependency-only cache flow. Compare `.github/workflows/c-cpp_all.yml` step-by-step when a failure only appears in CI.
- VS tests are always in-process. Never reintroduce subprocess execution to work around a crash; fix the underlying OCC/libcint/DLL issue instead.
- New libcint parallel call sites should pre-allocate per-thread scratch buffers instead of passing `cache=nullptr` inside TBB parallel regions. Use the existing `three_center_max_cache_size<kind>` and `tbb::enumerable_thread_specific` pattern in OCC.

## Current Validation Notes

As of the latest documented June 2026 validation, the native test suites had passed in Release and Debug configurations with 21/21 tests. Treat this as historical status unless you have rerun the current CMake presets in this checkout.

Known fixed areas from June 2026:

1. `alanine_integrated_occ` passes in-process in Release and Debug.
2. `IntegralEngineDF::~IntegralEngineDF()` clears Eigen buffers before libcint engines are destroyed.
3. oneTBB is deployed shared; `tbbmalloc_proxy` is not linked.
4. MSBuild instruction-set selection follows `NOS_AVX` to keep Eigen ABI consistent.
5. Parallel three-center libcint scratch buffers avoid heap corruption under vstest.
6. `fchk_conversion` writes and compares `log.fchk` and handles restricted/beta coefficients, CMO bounds, virtual orbitals, and G shells.
7. `Hybrid_mode` Debug vstest passes after the restricted `read_gbw` coefficient-span fix and Debug CRT dialog suppression.

See `UNIT_TESTS_STATUS.md` for test inventory details, but prefer live CMake files and CI workflow files for current build-system facts.
