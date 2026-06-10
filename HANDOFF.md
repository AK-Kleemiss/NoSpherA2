# Handoff: OCC dynamic TBB and alanine_integrated_occ fix

## RESOLVED: Debug Test Explorer crash (latest)

The Debug crash in `occ::gto::Shell::Shell` was the same instruction-set/Eigen-ABI mismatch as the earlier Release crash: Debug|x64 `NoSpherA2_DLL` was built without `/arch:AVX` while the OCC static libs in `Lib/occ/lib` are built with AVX, giving different Eigen object layouts across the link boundary.

Fix: added `<EnableEnhancedInstructionSet>AdvancedVectorExtensions</EnableEnhancedInstructionSet>` to Debug|x64, Profile|x64, and GenCube|x64 in `NoSpherA2_DLL/NoSpherA2_DLL.vcxproj` (Release|x64 already had it).

Verified: `vstest.console.exe Windows\x64\Debug\Tests.dll /Platform:x64 /Tests:alanine_integrated_occ` → Passed (40 s, in-process).

## Previous handoff (historical): Debug Test Explorer crash

Repository: `D:\git\NoSpherA2`

Environment: Windows host, PowerShell commands, Visual Studio 18 tools.

### Release state

The Release Visual Studio test path was fixed and verified:

- `alanine_integrated_occ` runs in-process through `NoSpherA2_DLL.dll`.
- Full `Tests.dll` suite with `RUN_FULL_TEST=1` passed: `21/21`.
- The Release crash was traced to an instruction-set mismatch: `NoSpherA2_DLL` Release x64 was missing AVX while `NoSpherA2.exe` had AVX enabled.
- `NoSpherA2_DLL/NoSpherA2_DLL.vcxproj` now enables `AdvancedVectorExtensions` for Release x64.
- OCC now releases its file logger during SCF shutdown so `NoSpherA2_OCC.log` is not held open after in-process tests.
- `Windows/Tests/TestRunner.h` now restores the original current directory before cleanup and uses unique temp directories per test run.

### New remaining problem

Debug config still crashes in Visual Studio Test Explorer when debugging `alanine_integrated_occ`.

Exception:

```text
0xC0000005: Access violation reading location 0xFFFFFFFFFFFFFFFF
```

Observed stack:

```text
occ::gto::Shell::Shell(...)
occ::gto::AOBasis::load(...)
occ::driver::load_basis_set(...)
occ::driver::single_point_driver(...)
occ::driver::single_point(...)
occ::main::run_scf_external(...)
main(...) at Src/NoSpherA2.cpp:434
nos_cpp_dispatch(...) at NoSpherA2_DLL/dllmain.cpp:53
nos_seh_dispatch(...) at NoSpherA2_DLL/dllmain.cpp:90
NoSpherA2_run(...) at NoSpherA2_DLL/dllmain.cpp:155
NosTestFramework::RunTest(...) at Windows/Tests/TestRunner.h:292
IntegrationTests::alanine_integrated_occ(...) at Windows/Tests/Tests.cpp:39
```

This means Release is handled, but Debug in-process has a separate config/runtime/ABI issue.

### Changed files relevant to latest state

Parent repo:

- `NoSpherA2_DLL/NoSpherA2_DLL.vcxproj`
  - Added AVX for Release x64.
- `NoSpherA2_DLL/dllmain.cpp`
  - Cleaned stale OCC heap-corruption wording; SEH guard remains a narrow last-resort guard.
- `Windows/Tests/TestRunner.h`
  - Uses unique temp directories.
  - Restores original CWD before cleanup.
  - Keeps subprocess only as opt-in fallback.
- `Windows/Tests/Tests.cpp`
  - Documents that tests normally run in-process through `NoSpherA2_DLL.dll`.
- `Src/NoSpherA2.cpp`
  - Only final-newline churn from this latest pass.

OCC submodule:

- `occ/include/occ/core/log.h`
  - Declares `occ::log::close_log_file()`.
- `occ/src/core/log.cpp`
  - Implements `close_log_file()` by flushing, restoring non-file sinks, and dropping the previous file logger from spdlog if needed.
- `occ/src/main/occ_scf.cpp`
  - Calls `occ::log::close_log_file()` from `occ::main::shutdown()`.

### Important git/status notes

Root currently has:

```text
M NoSpherA2_DLL/NoSpherA2_DLL.vcxproj
M NoSpherA2_DLL/dllmain.cpp
M Src/NoSpherA2.cpp
M Windows/Tests/TestRunner.h
M Windows/Tests/Tests.cpp
m occ
?? Windows/x64/
```

OCC submodule content diffs:

```text
M include/occ/core/log.h
M src/core/log.cpp
M src/main/occ_scf.cpp
```

OCC may also show line-ending/status noise in:

```text
src/driver/single_point.cpp
src/gto/shell.cpp
src/io/json_basis.cpp
```

Check `git -C occ diff --name-only` and `git -C occ diff --stat` before touching those. Do not revert user changes.

`Windows/x64/` is generated build output.

### Likely next checks

1. Compare Debug x64 settings between:

```text
Windows/NoSpherA2.vcxproj
NoSpherA2_DLL/NoSpherA2_DLL.vcxproj
```

Focus on:

```text
EnableEnhancedInstructionSet
RuntimeLibrary
Optimization
DebugInformationFormat
PreprocessorDefinitions
OpenMPSupport / OpenMP flags
_ITERATOR_DEBUG_LEVEL
/MDd vs /MD
linked OCC libraries under Lib/occ/lib
```

2. Main suspicion:

Debug `NoSpherA2_DLL.dll` may be linking Release-built OCC static libraries from `Lib/occ/lib`. That can produce Debug/Release STL or Eigen ABI friction around `std::vector` and allocations, surfacing in `occ::gto::Shell::Shell`.

Possible fixes to investigate:

- Build/install OCC Debug libraries into a config-specific location and link Debug NoSpherA2 against those.
- Or force Debug NoSpherA2 DLL settings to match the Release OCC ABI where appropriate, especially iterator debug/runtime assumptions.

3. Reproduce outside Visual Studio Test Explorer:

```powershell
$vstest='C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\TestWindow\vstest.console.exe'
& $vstest 'D:\git\NoSpherA2\Windows\x64\Debug\Tests.dll' /Platform:x64 /Tests:alanine_integrated_occ /Logger:"console;Verbosity=normal"
```

4. If needed, add short-lived targeted traces or breakpoints in:

```text
occ/src/gto/shell.cpp
occ/src/io/json_basis.cpp
occ/src/driver/single_point.cpp
```

Remove any tracing before finishing.

### Known-good Release validation commands

```powershell
$vcvarsPath='C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvarsall.bat'
$vcvars=& cmd /c "`"$vcvarsPath`" amd64 && set" 2>&1
$vcvars | Where-Object { $_ -match '=' } | ForEach-Object {
  $name,$value=$_ -split '=',2
  [Environment]::SetEnvironmentVariable($name,$value,'Process')
}

$vscmake='C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe'
& $vscmake --build build-windows-clang-cl --config Release --target occ_core occ_main
& $vscmake --install build-windows-clang-cl

$msbuild='C:\Program Files\Microsoft Visual Studio\18\Community\MSBuild\Current\Bin\MSBuild.exe'
& $msbuild Windows\Tests\Tests.vcxproj /m /p:Configuration=Release /p:Platform=x64 /p:SolutionDir=D:\git\NoSpherA2\Windows\ /v:minimal

$env:RUN_FULL_TEST = '1'
$vstest='C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\TestWindow\vstest.console.exe'
& $vstest 'D:\git\NoSpherA2\Windows\x64\Release\Tests.dll' /Platform:x64 /Logger:"console;Verbosity=normal"
Remove-Item Env:\RUN_FULL_TEST -ErrorAction SilentlyContinue
```

## Current status

The `alanine_integrated_occ` crash has been fixed and the Visual Studio test suite passes.

Verified on Windows with Visual Studio 18 tools:

```powershell
MSBuild Windows\Tests\Tests.vcxproj /m /p:Configuration=Release /p:Platform=x64 /p:SolutionDir=D:\git\NoSpherA2\Windows\
vstest.console.exe D:\git\NoSpherA2\Windows\x64\Release\Tests.dll /Platform:x64
```

Result: `21/21 passed`, including `alanine_integrated_occ`. A full run with `RUN_FULL_TEST=1` also passed.

## What changed

### OCC submodule

- `occ/3rdparty/CMakeLists.txt`
  - Builds oneTBB as shared while keeping OCC itself static.
- `occ/libocc/CMakeLists.txt`
  - Installs the TBB runtime/import artifacts into `Lib/occ`.
  - On Windows this provides `Lib/occ/bin/tbb12.dll` and `Lib/occ/lib/tbb12.lib`.
  - On Linux/macOS this installs `libtbb.so*`/`libtbb.dylib` while excluding allocator proxy files.
- `occ/include/occ/qm/integral_engine_df.h`
- `occ/src/qm/integral_engine_df.cpp`
  - Adds `IntegralEngineDF::~IntegralEngineDF()`.
  - The destructor explicitly clears DF-owned Eigen buffers and Split-RI-J storage before libcint-backed integral engines are destroyed.
  - This is the actual `alanine_integrated_occ` crash fix.
- `occ/src/io/occ_input.cpp`
  - Adds TOML parsing for `use_split_ri_j`.
- `occ/include/occ/main/occ_scf.h`
- `occ/src/main/occ_scf.cpp`
  - Adds a small shutdown helper that flushes logs and clears timing state.

### NoSpherA2 parent repo

- `CMakePresets.json`
  - Sets `TBB_BUILD_SHARED=ON` consistently for OCC presets.
- `Windows/NoSpherA2.vcxproj`
- `NoSpherA2_DLL/NoSpherA2_DLL.vcxproj`
  - Link against `tbb12.lib`.
  - Do not link `tbbmalloc.lib` or `tbbmalloc_proxy.lib`.
  - Copy `Lib/occ/bin/tbb12.dll` beside the built EXE/DLL.
- `Windows/Build_Dependencies/build_deps.ps1`
  - No longer builds `tbbmalloc_proxy`.
  - Copies `tbb12.lib` and `tbb12.dll` from the OCC build/install.
- `Linux/makefile`
  - Links against `Lib/occ/lib/libtbb`.
  - Copies `libtbb.so*` beside the executable.
- `Mac/makefile`
  - Links against `Lib/occ_<arch>/lib/libtbb`.
  - Copies/lipos `libtbb.dylib` for single-arch and universal builds.
- `.github/workflows/c-cpp_all.yml`
- `.github/actions/build-deps/action.yml`
  - GitHub Actions artifacts now include `tbb12.dll`, `libtbb.so*`, or `libtbb.dylib` instead of allocator proxy files.

## Important notes

- Do not reintroduce `tbbmalloc_proxy` as a workaround. Static/proxy TBB was part of the instability.
- The old heap-corruption diagnostics and `[HEAP]`/`[DETECT]` tracing were removed. If they appear again, they are stale local edits.
- The OCC changes live inside the `occ` submodule. Commit them in the OCC repository first, then update the `occ` submodule pointer in NoSpherA2.

## Recommended submodule workflow

From the parent repo:

```powershell
cd D:\git\NoSpherA2\occ
git switch -c nosphera2-dynamic-tbb-df-teardown
git add 3rdparty/CMakeLists.txt libocc/CMakeLists.txt include/occ/main/occ_scf.h include/occ/qm/integral_engine_df.h src/io/occ_input.cpp src/main/occ_scf.cpp src/qm/integral_engine_df.cpp
git commit -m "Use shared TBB and fix DF integral teardown"
git push -u origin nosphera2-dynamic-tbb-df-teardown
```

After the OCC PR/branch is accepted, update NoSpherA2 to the desired OCC commit:

```powershell
cd D:\git\NoSpherA2\occ
git fetch origin
git switch main
git pull

cd D:\git\NoSpherA2
git add occ
git commit -m "Update OCC submodule for shared TBB fix"
```

## Validation commands

Windows:

```powershell
$vcvars = & cmd /c "`"C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvarsall.bat`" amd64 && set" 2>&1
$vcvars | Where-Object { $_ -match '=' } | ForEach-Object {
    $name, $value = $_ -split '=', 2
    [System.Environment]::SetEnvironmentVariable($name, $value, 'Process')
}

$vscmake = "C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
& $vscmake --build build-windows-clang-cl --config Release --target occ_qm occ_driver occ_main
& $vscmake --install build-windows-clang-cl

$msbuild = "C:\Program Files\Microsoft Visual Studio\18\Community\MSBuild\Current\Bin\MSBuild.exe"
& $msbuild Windows\Tests\Tests.vcxproj /m /p:Configuration=Release /p:Platform=x64 /p:SolutionDir=D:\git\NoSpherA2\Windows\ /v:minimal

$vstest = "C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\TestWindow\vstest.console.exe"
& $vstest D:\git\NoSpherA2\Windows\x64\Release\Tests.dll /Platform:x64 /Logger:"console;Verbosity=normal"
```

Manual alanine smoke test:

```powershell
$env:OCC_DATA_PATH = "D:\git\NoSpherA2\Lib\occ\share\occ"
cd D:\git\NoSpherA2\tests\alanine_integrated_occ
& D:\git\NoSpherA2\Windows\x64\Release\NoSpherA2.exe -occ alanine.toml -cif alanine.cif -dmin 0.5 -acc 1 -all_charges -no_date
```
