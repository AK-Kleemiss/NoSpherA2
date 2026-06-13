# Handoff: Current Windows/OCC/FCHK Status

## Current Status: COMPLETE

Last updated: 2026-06-13

The current branch has the Windows/OCC integration and fchk conversion work fixed
and validated. Full test results from the latest run:

- Python pytest Release full: `21 passed` in 81.88s, with one tolerated
  `fchk_conversion` numeric warning at line `261148`
- Python pytest Debug full: `21 passed` in 1103.46s
- Visual Studio native Debug full: `21 passed` in 24.4986 minutes
- Visual Studio native Release full: `21 passed` in 1.3985 minutes

`Release+Copy` was intentionally excluded from the latest requested validation pass.

---

## Latest Fixes

### 1. MSVC Debug dependency build preset

The Windows dependency action expected a `windows-msvc-debug` workflow preset, but
`CMakePresets.json` only had the configure/build presets. Added the missing workflow
preset so the dependency build no longer fails with:

```text
No such workflow preset in D:/a/NoSpherA2/NoSpherA2: "windows-msvc-debug"
```

### 2. fchk conversion test and writer

`fchk_conversion` was not actually comparing the generated fchk output. The test now
generates `log.fchk` and compares it against `good.fchk`.

The writer crash was fixed in `Src/fchk.cpp`:

- Bound Alpha MO coefficient output by the actual `CMO.size()` instead of
  `nao * wave.get_nmo()`
- Only write beta coefficient/energy blocks when the wavefunction is genuinely
  unrestricted (`get_is_unrestricted()` and multiplicity greater than 1)
- Include virtual orbitals in the restricted alpha CMO path
- Add consistent G-shell handling during fchk normalization and coefficient assembly

The stale basis typo `dev2-TZVP` was corrected to `def2-TZVP`.

### 3. wfn_reading golden output

The `wfn_reading.good` reference was updated for the current integration output:

- Grid points: `25738` -> `26030`
- Becke: `135.005819` -> `135.005930`
- TVFC: `134.990915` -> `134.991026`

---

## OCC / alanine_integrated_occ Status

`alanine_integrated_occ` passes in both Release and Debug with all Visual Studio tests
running in-process through `NoSpherA2_DLL.dll`.

Resolved OCC/Windows issues:

1. **Heap corruption** - `IntegralEngineDF::~IntegralEngineDF()` clears large Eigen
   buffers before libcint engines are destroyed.
2. **TBB instability** - oneTBB is built shared (`TBB_BUILD_SHARED=ON`); NoSpherA2
   links the core TBB runtime only.
3. **AVX/Eigen ABI mismatch** - `NOS_AVX` controls the MSBuild instruction set
   through `Directory.Build.targets`, so NoSpherA2, the DLL, tests, and OCC use the
   same Eigen ABI. Do not hardcode `EnableEnhancedInstructionSet` in individual
   `.vcxproj` files.
4. **Parallel libcint heap traffic** - OCC TBB integral kernels pre-allocate per-thread
   libcint scratch buffers instead of passing `cache=nullptr` inside parallel loops.

All temporary `[NOS-DBG]` traces and Windows heap-check diagnostics have been removed.

---

## Changed Files of Interest

### Parent repo

- `CMakePresets.json` - added `windows-msvc-debug` workflow preset
- `Directory.Build.targets` - central Windows instruction-set policy from `NOS_AVX`
- `Src/fchk.cpp` - fixed fchk coefficient output and G-shell handling
- `Src/test_functions.h` - corrected `def2-TZVP` basis spelling
- `tests/tests.toml` - `fchk_conversion` now compares generated `log.fchk`
- `tests/wfn_reading/wfn_reading.good` - updated current reference output
- `Windows/Tests/Tests.cpp` - native test registration updates
- `Windows/Tests/TestRunner.h` - in-process Visual Studio test execution
- `Windows/NoSpherA2.vcxproj` and `NoSpherA2_DLL/NoSpherA2_DLL.vcxproj` - no
  hardcoded AVX; instruction set comes from `Directory.Build.targets`
- `Windows/Build_Dependencies/build_deps.ps1` - Debug dependency preset selection
- `Windows/sync_occ_variant.ps1` - OCC variant freshness checks

### OCC submodule

- `CMakeLists.txt` - MSVC Debug flags including `_DEBUG` and `/bigobj`
- `src/qm/detail/df_kernels.h` - per-thread libcint scratch buffers in TBB paths
- `include/occ/core/element.h` and `src/core/element.cpp` - current OCC changes in this branch

OCC source changes must be committed in the `occ/` repository first, then the parent
NoSpherA2 repository must commit the updated submodule pointer.

---

## Rules Going Forward

- Do not reintroduce `tbbmalloc_proxy`, `tbbmalloc`, or
  `#include <tbb/tbbmalloc_proxy.h>`
- Keep OCC using shared oneTBB (`TBB_BUILD_SHARED=ON`)
- Never add subprocess fallbacks for Visual Studio tests. The VS test harness must
  stay in-process so failures are natively debuggable.
- When adding OCC/libcint parallel call sites, pre-allocate per-thread scratch buffers;
  do not pass `cache=nullptr` inside TBB parallel regions.
- When a test passes locally but fails in CI, compare `.github/workflows/c-cpp_all.yml`
  and `.github/actions/build-deps` with local commands, especially `NOS_EXE`,
  `OCC_DATA_PATH`, and the selected CMake preset.

---

## Validation Commands

```powershell
# Python full suite, Release
$env:RUN_FULL_TEST = '1'
$env:OCC_DATA_PATH = 'D:\git\NoSpherA2\Lib\occ\share\occ'
$env:NOS_EXE = 'D:\git\NoSpherA2\Windows\x64\Release\NoSpherA2.exe'
uv run pytest

# Python full suite, Debug
$env:NOS_EXE = 'D:\git\NoSpherA2\Windows\x64\Debug\NoSpherA2.exe'
uv run pytest
```

```powershell
# Build native Visual Studio tests
$msbuild = 'C:\Program Files\Microsoft Visual Studio\18\Community\MSBuild\Current\Bin\amd64\MSBuild.exe'
& $msbuild Windows\Tests\Tests.vcxproj /m /p:Configuration=Debug /p:Platform=x64 /p:SolutionDir=D:\git\NoSpherA2\Windows\ /v:minimal
& $msbuild Windows\Tests\Tests.vcxproj /m /p:Configuration=Release /p:Platform=x64 /p:SolutionDir=D:\git\NoSpherA2\Windows\ /v:minimal

# Run native Visual Studio tests
$env:RUN_FULL_TEST = '1'
$env:OCC_DATA_PATH = 'D:\git\NoSpherA2\Lib\occ\share\occ'
$vstest = 'C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\TestWindow\vstest.console.exe'
& $vstest D:\git\NoSpherA2\Windows\x64\Debug\Tests.dll /Platform:x64 /Logger:"console;Verbosity=normal"
& $vstest D:\git\NoSpherA2\Windows\x64\Release\Tests.dll /Platform:x64 /Logger:"console;Verbosity=normal"
```
