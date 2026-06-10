# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is NoSpherA2

A C++ CLI tool that generates `.tsc` (non-spherical scattering factor) files and property grids from quantum chemistry wavefunction files. Used with Olex2/olex2.refine for Hirshfeld Atom Refinement (HAR) in crystallography.

## Build Commands

**Windows** (in VS 2022 Developer PowerShell, with `make.exe` on PATH):
```ps
make.exe
```

**Linux / macOS:**
```sh
make
```

Output: `NoSpherA2` (or `NoSpherA2.exe` on Windows) at the repo root.

CMake presets also exist for targeted builds: `linux-occ-gcc`, `macos-release-*`, `windows-clang-cl`.

## Running Tests

```sh
# Run all tests (requires uv)
uv run pytest

# Skip slow/full tests (default)
NOS_EXE=./NoSpherA2 OCC_DATA_PATH=occ/share uv run pytest

# Enable full/slow tests
RUN_FULL_TEST=1 uv run pytest
```

Key environment variables:
- `NOS_EXE` — path to the executable (defaults to `NoSpherA2` in PATH)
- `OCC_DATA_PATH` — must point to `occ/share` for basis set and data files
- `RUN_FULL_TEST` — set to any truthy value to enable slow integration tests

Tests are golden-file comparisons: the harness runs `NoSpherA2` with args from `tests/tests.toml`, captures stdout, and diffs against `.good` reference files in each test subdirectory.

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

All dependencies (OCC, featomic, libcint, Intel MKL) are linked **statically** to produce a portable, self-contained executable.

## Known Pitfalls

- **Golden-file float parsing**: the diff parser in `tests/run_test.py` can fail when a line contains multiple values. When fixing comparison logic, support scientific notation and multi-value lines without breaking existing thresholds.
- **Windows toolchain**: use an x64 VS Developer shell for the `windows-clang-cl` CMake preset flows.
- **Submodules**: always clone with `--recursive`. If submodules are missing, run `git submodule update --init --recursive`.
- **CI vs local**: when a test passes locally but fails in CI, compare `.github/workflows/c-cpp_all.yml` step-by-step with your local commands, particularly the `NOS_EXE` and `OCC_DATA_PATH` values.
