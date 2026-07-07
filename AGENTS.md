# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

## What is NoSpherA2

A C++ CLI tool that generates `.tsc` (non-spherical scattering factor) files and property grids from quantum chemistry wavefunction files. Used with Olex2/olex2.refine for Hirshfeld Atom Refinement (HAR) in crystallography.

## Build Commands

Preferred: use CMake presets.

**Windows (MSVC, full build):**
```ps
cmd /c "`"%VSINSTALLDIR%\Common7\Tools\VsDevCmd.bat`" -arch=x64 && cmake --preset release-windows && cmake --build --preset release-windows"
```

**Linux (GCC, full build):**
```sh
cmake --preset linux-gcc
cmake --build --preset linux-gcc
```

**macOS (full build, per-arch):**
```sh
cmake --preset macos-release-full-arm64
cmake --build --preset macos-release-full-arm64
```

Output: the executable is produced in the preset build folder (e.g. `build/release-windows/bin/NoSpherA2.exe`).

## Running Tests

```sh
# Run all tests (requires uv)
uv run pytest

# Skip slow/full tests (default)
NOS_EXE=./NoSpherA2 OCC_DATA_PATH=occ/share uv run pytest

# Enable full/slow tests
RUN_FULL_TEST=1 uv run pytest

# Or via CMake/CTest
cmake --build --preset <your-preset>
ctest --preset <your-preset>
```

Key environment variables:
- `NOS_EXE` — path to the executable (defaults to `NoSpherA2` in PATH)
- `OCC_DATA_PATH` — must point to OCC data files; on local Windows builds this is normally `D:\git\NoSpherA2\Lib\occ\share\occ`
- `RUN_FULL_TEST` — set to any truthy value to enable slow integration tests

Tests are golden-file comparisons: the harness runs `NoSpherA2` with args from `tests/tests.toml`, captures stdout, and diffs against `.good` reference files in each test subdirectory.

### Windows Visual Studio Tests

The repository also has a Visual Studio test project under `Windows/Tests`. Prefer this path when validating Windows/OCC integration issues:

```ps
msbuild Windows\Tests\Tests.vcxproj /m /p:Configuration=Release /p:Platform=x64 /p:SolutionDir=D:\git\NoSpherA2\Windows\
vstest.console.exe Windows\x64\Release\Tests.dll /Platform:x64
```

Visual Studio tests must run in-process through `NoSpherA2_DLL.dll`; do not add subprocess fallbacks for failing VS tests, including `alanine_integrated_occ`. If all tests complete in under a second, the test binary or working directory is probably wrong. The default suite should take tens of seconds; `RUN_FULL_TEST=1` enables the longer cases.

### Current Validation Status

As of 2026-06-13, the full suite passes in all actively validated native configurations:
- Python pytest Release full: `21 passed` (one tolerated `fchk_conversion` numeric warning)
- Python pytest Debug full: `21 passed`
- Visual Studio native Debug full: `21 passed`
- Visual Studio native Release full: `21 passed`

`Release+Copy` was intentionally not rerun in the latest validation pass.

### Adding a Test

Register in `tests/tests.toml`:
```toml
[my_test_name]
directory = "my_test_dir"         # subdirectory under tests/
# good = "custom.good"            # optional: override default good file name

[my_test_name.args]
cif = "structure.cif"
hkl = "data.hkl"
wfn = "wavefunction.wfx"
acc = 0
```

CLI args are always under `<testname>.args`. The test framework passes these as `--key value` to the executable.

## Architecture

### Source Layout

```
Src/           C++ source for the NoSpherA2 library and CLI
occ/           OCC submodule — quantum chemistry library (26+ modules)
featomic/      Rust submodule — ML feature generation via metatensor
libcint/       C submodule — electron integral evaluation
mdspan/        C++ submodule — multidimensional array reference impl
tests/         Python pytest harness + golden-file test cases
Windows/       VS project files and Windows-specific build scripts
Linux/         Linux build scripts
Mac/           macOS build scripts (universal binary: x86_64 + arm64)
```

### Key `Src/` Modules

- **Entry point**: `NoSpherA2.cpp` — CLI option parsing, dispatches to computation
- **I/O**: `fchk.cpp/h`, `cif.h`, `wfn_class.cpp/h` — parse Gaussian fchk/wfx, CIF, and wavefunction formats
- **Molecular representation**: `atoms.cpp/h`, `molecule.cpp/h`, `basis_set.cpp/h`
- **Integration engine**: `AtomGrid.cpp/h`, `GridManager.cpp/h`, `integrator.cpp/h`, `integration_params.cpp/h`
- **Math/integrals**: `nos_math.cpp/h`, `libCintMain.cpp/h`, `libCintKernels.cpp/h`
- **Analysis**: `scattering_factors.cpp/h`, `bondwise_analysis.cpp/h`, `isosurface.cpp/h`
- **ML density (SALTED)**: `SALTED_*.cpp/h` — machine-learning-based density predictions
- **Output**: `tsc_block.h`, `properties.cpp/h`, `cube.cpp/h` — TSC format, cube files, property grids

The `Src/CMakeLists.txt` builds both a static `libnosphera2` and the `NoSpherA2` executable.

### Dependency on OCC

The `occ/` submodule is the heaviest dependency — it provides HF/DFT, crystal analysis, isosurface computation, and more. When editing quantum chemistry logic, check whether OCC already provides the needed functionality before implementing it in `Src/`.

### Linking

Most dependencies (OCC, featomic, libcint, Intel MKL static libraries) are linked statically. oneTBB is intentionally built and deployed dynamically because static/proxy TBB caused Windows runtime instability in OCC integration tests.

Expected runtime files beside packaged executables:
- Windows: `tbb12.dll`
- Linux: `libtbb.so*`
- macOS: `libtbb.dylib`

Do not reintroduce `tbbmalloc_proxy` or `tbbmalloc` linking for NoSpherA2. The build scripts copy the core TBB runtime from `Lib/occ`.

## Known Pitfalls

- **Golden-file float parsing**: the diff parser in `tests/run_test.py` can fail when a line contains multiple values. When fixing comparison logic, support scientific notation and multi-value lines without breaking existing thresholds.
- **Windows toolchain**: before any Windows CMake configure, build, or test command, initialize the MSVC environment with the Visual Studio developer setup script. In this environment, run `VsDevCmd.bat -arch=x64` first, or chain it in the same `cmd /c` invocation as CMake. If standard headers such as `cstddef`, `vector`, `windows.h`, or libraries such as `kernel32.lib` are missing, treat that as an uninitialized MSVC/Windows SDK environment, not as a source patch failure.
- **Submodules**: always clone with `--recursive`. If submodules are missing, run `git submodule update --init --recursive`.
- **OCC submodule commits**: OCC source fixes must be committed in the `occ` repository first, then the parent NoSpherA2 repo must commit the updated submodule pointer.
- **CI vs local**: when a test passes locally but fails in CI, compare `.github/workflows/c-cpp_all.yml` step-by-step with your local commands, particularly the `NOS_EXE` and `OCC_DATA_PATH` values.

## OCC Integration and Windows Tests — Fixed (June 2026)

The `alanine_integrated_occ` test passes in both Release and Debug (21/21). Current fixes include:

1. **Heap corruption** — added `IntegralEngineDF::~IntegralEngineDF()` to clear Eigen buffers before libcint engines are destroyed.
2. **TBB instability** — oneTBB is now built shared (`TBB_BUILD_SHARED=ON`); `tbbmalloc_proxy` is not linked.
3. **AVX/ABI mismatch** — `Directory.Build.targets` selects the MSBuild instruction set from `NOS_AVX`, keeping NoSpherA2, the DLL, tests, and OCC on the same Eigen ABI.
4. **MSVC Debug dependency build** — `windows-msvc-debug` now has a workflow preset, matching the CI dependency action.
5. **fchk conversion** — `fchk_conversion` now writes and compares `log.fchk`; `Src/fchk.cpp` handles the restricted/beta coefficient path and G shells correctly.

All `[NOS-DBG]` diagnostic traces and Windows heap-check helpers have been removed from the source.

See `HANDOFF.md` and `INVESTIGATION_STATUS.md` for full details and validation commands.
