# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is NoSpherA2

A C++ CLI tool that generates `.tsc` (non-spherical scattering factor) files and property grids from quantum chemistry wavefunction files. Used with Olex2/olex2.refine for Hirshfeld Atom Refinement (HAR) in crystallography.

## Build Commands

Preferred: use CMake presets.

**Windows (MSVC, full build):**
```ps
cmake --preset windows-msvc-release-full
cmake --build --preset windows-msvc-release-full
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

The CMake configuration bootstraps a micromamba environment by default (see README) to ensure dependencies are available.

### Building on Windows from an Agent/CLI session (no Developer Shell)

The VS Developer Shell (`VsDevCmd.bat`) cannot be sourced inside the Bash or PowerShell tools because environment variables set by a `call` inside `cmd /c` do not persist to the next tool invocation, and `msbuild` is not on the default PATH.

Use the full path to MSBuild and invoke it from **PowerShell** (not Bash — Bash strips leading `/` from `/p:` switches):

```powershell
$msbuild = "C:\Program Files\Microsoft Visual Studio\18\Community\MSBuild\Current\Bin\amd64\MSBuild.exe"
$solDir   = "D:\git\NoSpherA2\Windows\"

# Build the DLL (contains all Src/ production code)
& $msbuild "D:\git\NoSpherA2\NoSpherA2_DLL\NoSpherA2_DLL.vcxproj" `
    /p:Configuration=Release /p:Platform=x64 "/p:SolutionDir=$solDir" /v:minimal

# Build the VS unit-test project
& $msbuild "D:\git\NoSpherA2\Windows\Tests\Tests.vcxproj" `
    /p:Configuration=Release /p:Platform=x64 "/p:SolutionDir=$solDir" /v:minimal
```

Replace `Release` with `Debug` for a debug build. Supported configurations: `Debug`, `Release`, `Profile`, `GenCube`, `Release + Copy`.

Run the VS unit tests:

```powershell
$vstest = "C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\TestWindow\vstest.console.exe"
& $vstest "D:\git\NoSpherA2\Windows\x64\Release\Tests.dll" /Platform:x64
```

For code coverage (instruments `NoSpherA2_DLL.dll` + `Tests.dll`):

```powershell
& $vstest "D:\git\NoSpherA2\Windows\x64\Profile\Tests.dll" /Platform:x64 `
    /Settings:"D:\git\NoSpherA2\Windows\coverage.runsettings" /EnableCodeCoverage
```

**Key pitfalls:**
- `NoSpherA2_DLL.vcxproj` lives at `NoSpherA2_DLL\` (repo root), **not** under `Windows\`.
- Always pass `/p:SolutionDir=` so relative paths inside the vcxproj resolve correctly.
- Never use `Bash` for MSBuild invocations — the `/p:` flag prefix is silently stripped.

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
- `OCC_DATA_PATH` — must point to OCC data files; on local Windows builds this is normally `D:\git\NoSpherA2\Lib\occ\share\occ`
- `RUN_FULL_TEST` — set to any truthy value to enable slow integration tests

Tests are golden-file comparisons: the harness runs `NoSpherA2` with args from `tests/tests.toml`, captures stdout, and diffs against `.good` reference files in each test subdirectory.

### Windows Visual Studio Tests

The repository also has a Visual Studio test project under `Windows/Tests`. Prefer this path when validating Windows/OCC integration issues:

```powershell
# See "Building on Windows from an Agent/CLI session" above for the full invocation.
# Quick reference (run from a VS Developer PowerShell where msbuild is on PATH):
msbuild Windows\Tests\Tests.vcxproj /p:Configuration=Release /p:Platform=x64 /p:SolutionDir=D:\git\NoSpherA2\Windows\
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

**Also add a `TEST_METHOD` to `Windows/Tests/Tests.cpp`** — every entry in `tests.toml` must have a matching VS test method so both harnesses stay in sync.

After adding a test, update `UNIT_TESTS_STATUS.md` with the new row and its validation result.

## Agent / AI Coding-Assistant Rules

These rules apply to Claude Code, Codex, Copilot, and any other AI agent working in this repo.

### Unit-test documentation is mandatory

Whenever an agent **adds, removes, modifies, or investigates** a test — whether in `tests/tests.toml`, `Windows/Tests/Tests.cpp`, `tests/run_test.py`, `TestRunner.h`, or any `.good` reference file — it **must**:

1. Update `UNIT_TESTS_STATUS.md`:
   - If a new test is added: add a row to the registered-test table with status `🆕 not yet validated`.
   - If a test is validated as passing: change its status to `✅ passing` and update the validation-status table at the top.
   - If a test is found to be failing or crashing: change its status to `❌ failing` or `⚠️ crashing` and add a note to the Known Issues section.
   - If a test is removed: delete its row and note why in a comment commit message.
2. Update the "Current Validation Status" block in **this file (`CLAUDE.md`)** if the overall pass count changes.
3. Never leave `UNIT_TESTS_STATUS.md` stale after making test changes — the file is the single source of truth for test coverage and should be commitworthy on its own.

### Do not break the VS test / pytest parity

`tests/tests.toml` and `Windows/Tests/Tests.cpp` must stay in sync — every test registered in the TOML must have a matching `TEST_METHOD`, and vice versa. Check both files when adding or removing tests.

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
- **Windows toolchain**: use an x64 VS Developer shell for `windows-clang-cl` and `windows-msvc-debug` CMake preset flows.
- **Submodules**: always clone with `--recursive`. If submodules are missing, run `git submodule update --init --recursive`.
- **OCC submodule commits**: OCC source fixes must be committed in the `occ` repository first, then the parent NoSpherA2 repo must commit the updated submodule pointer.
- **CI vs local**: when a test passes locally but fails in CI, compare `.github/workflows/c-cpp_all.yml` step-by-step with your local commands, particularly the `NOS_EXE` and `OCC_DATA_PATH` values.
- **VS tests are always in-process**: the `subprocess` field has been removed from `TestDef`. Never reintroduce subprocess execution to work around an in-process failure — fix the underlying OCC/libcint issue instead. If a VS test crashes in-process, attach the debugger directly (the managed vstest host already has the debugger attached).
- **New libcint parallel call sites**: always pre-allocate per-thread scratch buffers instead of passing `cache=nullptr` inside TBB parallel regions. Use `three_center_max_cache_size<kind>` + `tbb::enumerable_thread_specific` (see `occ/src/qm/detail/df_kernels.h`).

## OCC Integration and Windows Tests — Fixed (June 2026)

The `alanine_integrated_occ` test passes in both Release and Debug (21/21) with all tests running in-process. Current fixes include:

1. **Heap corruption** — added `IntegralEngineDF::~IntegralEngineDF()` to clear Eigen buffers before libcint engines are destroyed.
2. **TBB instability** — oneTBB is now built shared (`TBB_BUILD_SHARED=ON`); `tbbmalloc_proxy` is not linked.
3. **AVX/ABI mismatch** — `Directory.Build.targets` selects the MSBuild instruction set from `NOS_AVX`, keeping NoSpherA2, the DLL, tests, and OCC on the same Eigen ABI.
4. **In-process heap corruption in TBB parallel integrals** — `compute_three_center_integrals_tbb` and `three_center_aux_kernel` now pre-allocate per-thread libcint scratch buffers so libcint never calls `malloc`/`free` inside the parallel region. This eliminates the `STATUS_HEAP_CORRUPTION` (0xC0000374) fast-fail that occurred under the busy .NET vstest heap.
5. **MSVC Debug dependency build** — `windows-msvc-debug` now has a workflow preset, matching the CI dependency action and fixing the missing-preset build failure.
6. **fchk conversion** — `fchk_conversion` now writes and compares `log.fchk`; `Src/fchk.cpp` handles restricted/beta coefficient emission, CMO bounds, virtual orbitals, and G shells correctly.

All `[NOS-DBG]` diagnostic traces and Windows heap-check helpers have been removed from the source.

See `HANDOFF.md` for full details and validation commands.

## Hybrid_mode Debug Fix — Fixed (June 2026)

`Hybrid_mode` Debug vstest now passes as part of the full 21/21 Debug suite. Two issues were resolved:

1. **OOB vector access in `read_gbw`** — `coefs_2D_s2_span` was unconditionally constructed with `coefficients[1].data()` even when `operators==1` (restricted wavefunction). With MSVC IDL=2 this raises a CRT assertion dialog blocking the vstest host (WaitReason=UserRequest, 0% CPU). Fix: use `coefficients[0].data()` as a fallback when `operators!=2`.
2. **Debug CRT dialog suppression** — `_CrtSetReportMode` + `_set_invalid_parameter_handler` added to `dllmain.cpp` (Debug-only) so future CRT assertions print to stderr rather than showing invisible modal dialogs in the test host.
