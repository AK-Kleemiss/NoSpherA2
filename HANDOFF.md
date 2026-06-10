# Handoff: OCC dynamic TBB and alanine_integrated_occ fix

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
