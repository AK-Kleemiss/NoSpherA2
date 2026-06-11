# alanine_integrated_occ Crash Investigation — RESOLVED
**Last updated: 2026-06-11**

## Status: FIXED and VERIFIED

The `alanine_integrated_occ` test passes in both:
- In-process DLL mode (Windows VS test runner — Release and Debug)
- Subprocess EXE mode

Full `Tests.dll` suite: **21/21 passed** (including `alanine_integrated_occ`).

---

## Root Cause (Confirmed)

Two independent issues caused crashes:

### Issue 1: Windows heap corruption in `~IntegralEngineDF` (Release)
`IntegralEngineDF` owned Eigen matrices (`m_integral_store`, Split-RI-J buffers) whose aligned
allocations were freed *after* the libcint-backed integral engines (`m_ao_engine`, `m_aux_engine`)
were destroyed. The destruction order of member variables in a statically linked OCC caused
`handmade_aligned_free` to read a corrupted back-pointer → `RtlFailFast` →
`STATUS_HEAP_CORRUPTION (0xC0000374)`.

**Fix**: added an explicit `IntegralEngineDF::~IntegralEngineDF()` in
`occ/src/qm/integral_engine_df.cpp` that clears `m_split_rij`, `m_integral_store`, and `V_LLt`
in a controlled order before the member engines are destroyed.

### Issue 2: TBB runtime instability (all modes)
Linking `tbbmalloc_proxy` / static TBB caused allocator-interception conflicts across the
DLL/EXE boundary on Windows, leading to additional heap instability.

**Fix**: OCC's TBB is now built shared (`TBB_BUILD_SHARED=ON`). NoSpherA2 links against
`tbb12.lib` only (no allocator proxy). The `tbb12.dll` runtime is deployed beside the
executable.

### Issue 3: AVX/Eigen ABI mismatch (Debug DLL)
Debug|x64 `NoSpherA2_DLL` was compiled without `/arch:AVX` while OCC static libs were built
with AVX, giving different Eigen object sizes across the link boundary. This surfaced as an
access violation inside `occ::gto::Shell::Shell`.

**Fix**: added `<EnableEnhancedInstructionSet>AdvancedVectorExtensions</EnableEnhancedInstructionSet>`
to Debug|x64, Profile|x64, and GenCube|x64 in `NoSpherA2_DLL/NoSpherA2_DLL.vcxproj`.

---

## What Changed

### OCC submodule (`occ/`)
- `3rdparty/CMakeLists.txt` — builds oneTBB as shared
- `libocc/CMakeLists.txt` — installs TBB runtime artifacts into `Lib/occ`
- `include/occ/qm/integral_engine_df.h` + `src/qm/integral_engine_df.cpp` — adds explicit destructor
- `include/occ/core/log.h` + `src/core/log.cpp` — adds `close_log_file()` to release file handles after in-process tests
- `src/main/occ_scf.cpp` — calls `close_log_file()` from shutdown

### NoSpherA2 parent repo
- `CMakePresets.json` — `TBB_BUILD_SHARED=ON` for all OCC presets, including `macos-base`
- `Windows/NoSpherA2.vcxproj` + `NoSpherA2_DLL/NoSpherA2_DLL.vcxproj` — link `tbb12.lib`, drop `tbbmalloc*`; AVX enabled for all configs in DLL project
- `Windows/Build_Dependencies/build_deps.ps1` — no longer builds `tbbmalloc_proxy`
- `Linux/makefile` + `Mac/makefile` — link `Lib/occ/lib/libtbb`, deploy runtime beside executable
- `.github/workflows/c-cpp_all.yml` + `.github/actions/build-deps/action.yml` — CI artifacts include TBB runtime

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

# Run tests
$vstest = "C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\TestWindow\vstest.console.exe"
& $vstest D:\git\NoSpherA2\Windows\x64\Release\Tests.dll /Platform:x64 /Logger:"console;Verbosity=normal"
```

Manual smoke test:
```powershell
$env:OCC_DATA_PATH = "D:\git\NoSpherA2\Lib\occ\share\occ"
cd D:\git\NoSpherA2\tests\alanine_integrated_occ
& D:\git\NoSpherA2\Windows\x64\Release\NoSpherA2.exe -occ alanine.toml -cif alanine.cif -dmin 0.5 -acc 1 -all_charges -no_date
```

---

## Do Not Reintroduce
- `tbbmalloc_proxy` / `tbbmalloc` linking in any NoSpherA2 target
- `#include <tbb/tbbmalloc_proxy.h>` in `Src/pch.h` or anywhere else
- Static TBB linkage inside OCC (keep `TBB_BUILD_SHARED=ON`)
- Debug diagnostic `fprintf(stderr, "[NOS-DBG]...")` traces — these were removed in June 2026
