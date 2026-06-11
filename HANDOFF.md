# Handoff: OCC dynamic TBB and alanine_integrated_occ fix

## Current Status: COMPLETE

All crashes are fixed, all diagnostics removed. The full Visual Studio test suite
(`21/21`) passes in both Release and Debug modes.

---

## What Was Fixed

### 1. Windows heap corruption / `STATUS_HEAP_CORRUPTION` (Release EXE)
`IntegralEngineDF` was destroying its Eigen buffers in the wrong order relative to
the libcint-backed integral engines. Fixed by adding an explicit
`IntegralEngineDF::~IntegralEngineDF()` that clears `m_split_rij`, `m_integral_store`,
and `V_LLt` before the member engines go out of scope.

### 2. TBB runtime instability (all modes)
Linking `tbbmalloc_proxy` / static TBB caused allocator conflicts. Fixed by building
oneTBB shared (`TBB_BUILD_SHARED=ON` in all OCC CMake presets) and deploying
`tbb12.dll` / `libtbb.so*` / `libtbb.dylib` beside the executable.

### 3. AVX / Eigen ABI mismatch (Debug DLL)
`NoSpherA2_DLL` Debug|x64 lacked `/arch:AVX` while OCC static libs were built with
AVX. Fixed by adding `AdvancedVectorExtensions` to Debug|x64, Profile|x64, and
GenCube|x64 in `NoSpherA2_DLL/NoSpherA2_DLL.vcxproj`.

---

## Changed Files

### OCC submodule (`occ/`)
- `3rdparty/CMakeLists.txt` — oneTBB built as shared
- `libocc/CMakeLists.txt` — TBB runtime installed into `Lib/occ`
- `include/occ/qm/integral_engine_df.h` + `src/qm/integral_engine_df.cpp` — explicit destructor
- `include/occ/core/log.h` + `src/core/log.cpp` — `close_log_file()` for in-process tests
- `src/main/occ_scf.cpp` — calls `close_log_file()` from shutdown
- `include/occ/qm/scf_impl.h` — removed `nos_heapcheck` debug helper and all call sites
- `src/gto/shell.cpp` — removed `[NOS-DBG]` traces and per-shell HeapValidate loop
- `src/qm/hf.cpp` — removed `[NOS-DBG]` traces and HeapValidate block
- `src/driver/single_point.cpp` — removed `[NOS-DBG]` traces
- `src/main/occ_scf.cpp` — removed `[NOS-DBG]` traces

### NoSpherA2 parent repo
- `CMakePresets.json` — `TBB_BUILD_SHARED=ON` for all OCC presets (including `macos-base`)
- `Windows/NoSpherA2.vcxproj` + `NoSpherA2_DLL/NoSpherA2_DLL.vcxproj`
  - Link `tbb12.lib`, drop `tbbmalloc*` and allocator proxy
  - AVX (`AdvancedVectorExtensions`) enabled for all configs in DLL project
- `Windows/Build_Dependencies/build_deps.ps1` — no longer builds `tbbmalloc_proxy`
- `Linux/makefile` + `Mac/makefile` — link `Lib/occ/lib/libtbb`, deploy runtime
- `.github/workflows/c-cpp_all.yml` + `.github/actions/build-deps/action.yml` — CI uses TBB runtime artifacts
- `Src/NoSpherA2.cpp` — removed `_heapchk` debug call and `<malloc.h>` include
- `Windows/Tests/TestRunner.h` — unique temp dirs, CWD restore before cleanup
- `Windows/Tests/Tests.cpp` — tests run in-process through `NoSpherA2_DLL.dll`
- `NoSpherA2_DLL/dllmain.cpp` — cleaned stale heap-corruption wording

---

## Rules Going Forward

- Do **not** reintroduce `tbbmalloc_proxy`, `tbbmalloc`, or `#include <tbb/tbbmalloc_proxy.h>`
- Do **not** use static TBB inside OCC — keep `TBB_BUILD_SHARED=ON`
- OCC source changes must be committed inside `occ/` first, then update the submodule pointer in the parent repo
- When adding new OCC features that own large Eigen buffers, write an explicit destructor to control teardown order relative to libcint engines

---

## Validation Commands

```powershell
# Initialize VS environment
$vcvars = & cmd /c "`"C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvarsall.bat`" amd64 && set" 2>&1
$vcvars | Where-Object { $_ -match '=' } | ForEach-Object {
    $name, $value = $_ -split '=', 2
    [System.Environment]::SetEnvironmentVariable($name, $value, 'Process')
}

# Build OCC and install
$vscmake = "C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
& $vscmake --build build-windows-clang-cl --config Release --target occ_qm occ_driver occ_main
& $vscmake --install build-windows-clang-cl

# Build test suite
$msbuild = "C:\Program Files\Microsoft Visual Studio\18\Community\MSBuild\Current\Bin\MSBuild.exe"
& $msbuild Windows\Tests\Tests.vcxproj /m /p:Configuration=Release /p:Platform=x64 /p:SolutionDir=D:\git\NoSpherA2\Windows\ /v:minimal

# Run full suite
$env:RUN_FULL_TEST = '1'
$vstest = "C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\TestWindow\vstest.console.exe"
& $vstest D:\git\NoSpherA2\Windows\x64\Release\Tests.dll /Platform:x64 /Logger:"console;Verbosity=normal"
Remove-Item Env:\RUN_FULL_TEST -ErrorAction SilentlyContinue
```

Manual alanine smoke test:
```powershell
$env:OCC_DATA_PATH = "D:\git\NoSpherA2\Lib\occ\share\occ"
cd D:\git\NoSpherA2\tests\alanine_integrated_occ
& D:\git\NoSpherA2\Windows\x64\Release\NoSpherA2.exe -occ alanine.toml -cif alanine.cif -dmin 0.5 -acc 1 -all_charges -no_date
```
