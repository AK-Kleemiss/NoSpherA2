---
name: nosphera2-build-system
description: Use when working in the NoSpherA2 repository on CMake presets, dependency bootstrap/cache flow, build targets, tests, CI build failures, runtime-library copying, or agent guidance updates after the build-system restructure.
---

# NoSpherA2 Build System

## Start Here

Treat `CMakePresets.json`, the top-level `CMakeLists.txt`, and `.github/workflows/c-cpp_all.yml` as the current sources of truth. Older notes may mention removed preset names such as `windows-msvc-release-full`, `linux-gcc`, `macos-release-full-arm64`, or `windows-msvc-debug`; verify against the live preset file before using them.

For a fresh checkout, bootstrap dependencies first:

```sh
cmake -P scripts/BootstrapMicromamba.cmake
```

Then configure/build with an active preset:

```sh
cmake --preset release-windows
cmake --build --preset release-windows
```

Common preset families are:

- Windows: `release-windows`, `debug-windows`
- Linux: `release-linux`, `debug-linux`
- macOS: `release-macos-arm64`, `release-macos-x86_64`, `debug-macos-arm64`, `debug-macos-x86_64`
- macOS universal build target: `cmake --build --preset release-macos-universal`

Use `cmake --list-presets=all` for a complete view. A plain `cmake --list-presets` can hide configure presets that do not match the host platform.

## Target Layout

The top-level CMake project configures dependencies and then builds project targets:

- `libcint`, `occ`, `mdspan`, and `featomic/featomic` are consumed with `add_subdirectory`.
- `BasisSetGenerator` builds `BasisSetConverter`, which generates `Src/basis_data.cpp`.
- `Src/core` builds the static `NoSpherA2Core` target and alias `NoSpherA2::Core`. Its sources sit in subsystem folders (`base`, `wfn`, `integration`, `tsc`, `xcw`, `salted`, `nbo`, `analysis`, `properties`, `gpu`), all on the PUBLIC include path (`NOSPHERA2_CORE_INCLUDE_DIRS`, also handed to the MSVC-host HIP compile), so a new folder goes into that list; the sources are ordered by file name, as in the flat folder.
- `Src/dll` builds the shared library (`NoSpherA2.dll`, `libNoSpherA2.so`, `libNoSpherA2.dylib`) into `bin`; it carries all the code.
- `app` builds the `NoSpherA2` executable, a launcher that loads that library from its own directory (`nosphera2_main`), so the exe never runs without it. Android links one static exe instead.
- `tests/src` builds `NoSpherA2_Tests` when `NOSPHERA2_BUILD_TESTS=ON`.

Runtime outputs are under the preset build tree's `bin` directory, for example `build/release-windows/bin/NoSpherA2.exe`.

## CMake Options

- Use `-DNOSPHERA2_BUILD_TESTS=ON` before running CTest.
- Use `-DNOSPHERA2_DEPENDENCIES_ONLY=ON` only for dependency-cache or dependency-prebuild workflows. The top-level `CMakeLists.txt` returns before adding NoSpherA2 targets in this mode.

CI restores or builds a dependency-only tree, saves that cache, then reconfigures the same build directory with `NOSPHERA2_BUILD_TESTS=ON` and `NOSPHERA2_DEPENDENCIES_ONLY=OFF`.

## Tests

Preferred CMake validation:

```sh
cmake --preset <preset> -DNOSPHERA2_BUILD_TESTS=ON
cmake --build --preset <preset>
ctest --preset <preset> --output-on-failure
```

Python golden-file tests remain available:

```sh
uv run pytest
```

Set `NOS_EXE` for Python tests when the executable is not on `PATH`. Set `OCC_DATA_PATH` when the runtime data path is not discoverable; CTest defaults this to `occ/share`.

When adding or changing tests, keep `tests/tests.toml`, `tests/src`, and `UNIT_TESTS_STATUS.md` aligned for the affected coverage.

## Windows Notes

For preset builds on Windows, use an x64 Visual Studio Developer shell or another environment where Ninja can find `clang-cl`, `lld-link` and the MSVC headers. `release-windows`/`debug-windows` compile with `clang-cl` from `.mambaenv/llvm` (`environment-llvm-win.yaml`); the ARM64 presets keep `cl`. Linux presets use `clang`/`clang++` with `CMAKE_LINKER_TYPE=LLD` from `.mambaenv/llvm` (`environment-llvm-linux.yaml`: clang 19.1.7, glibc 2.28 sysroot, GCC 14 libstdc++ headers; aarch64 in `.mambaenv/env`); macOS uses Apple clang. GCC is a test-only override (`-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_LINKER_TYPE=BFD`). A build tree configured before the switch keeps its old compiler in `CMakeCache.txt`; delete the cache. Missing standard headers such as `stdlib.h`, `string.h`, or `stdio.h` usually indicate shell/toolchain setup, not a broken source rewrite.

The Visual Studio solution was removed in October 2026; CMake/Ninja is the only build system. IDE work on every platform goes through `.vscode/` (CMake Tools with the presets; `cmake.useVsDeveloperEnvironment` supplies the VS environment, `clang-cl` included, on Windows). In-process tests (e.g. `alanine_integrated_occ`) stay in-process in the gtest binary; do not add subprocess fallbacks for crashing tests.

## Dependency Boundaries

Do not edit submodules casually. Prefer parent-project CMake or generated-source workarounds when fixing integration issues from the parent repository.

The parent project currently rewrites selected `libcint` source files under the build tree to make `#pragma GCC ivdep` portable to MSVC. Preserve this parent-side approach unless the user explicitly asks for a submodule fix.
