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
- `Src/core` builds the static `NoSpherA2Core` target and alias `NoSpherA2::Core`.
- `app` builds the `NoSpherA2` executable.
- `Src/dll` builds the optional DLL when `NOSPHERA2_BUILD_DLL=ON`.
- `tests/src` builds `NoSpherA2_Tests` when `NOSPHERA2_BUILD_TESTS=ON`.

Runtime outputs are under the preset build tree's `bin` directory, for example `build/release-windows/bin/NoSpherA2.exe`.

## CMake Options

- Use `-DNOSPHERA2_BUILD_TESTS=ON` before running CTest.
- Use `-DNOSPHERA2_BUILD_DLL=ON` only when the Windows DLL target is needed.
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

When adding or changing tests, keep `tests/tests.toml`, `tests/src`, `Windows/Tests/Tests.cpp`, and `UNIT_TESTS_STATUS.md` aligned for the affected coverage.

## Windows Notes

For preset builds on Windows, use an x64 Visual Studio Developer shell or another environment where Ninja can find `cl.exe` and the MSVC headers. Missing standard headers such as `stdlib.h`, `string.h`, or `stdio.h` usually indicate shell/toolchain setup, not a broken source rewrite.

The Visual Studio solution under `Windows/` is legacy but still useful for DLL/OCC integration debugging. VS tests must remain in-process through `NoSpherA2_DLL.dll`; do not add subprocess fallbacks for crashing tests.

The VS `NoSpherA2`/`NoSpherA2_LIB` projects share an `OutDir` (`Windows\Windows_utils\NoSpherA2_universal.props`) pointing at `build\$(Configuration)_$(Platform)\` (e.g. `build\Release_x64\NoSpherA2.exe`) — a sibling of, but distinct from, the CMake preset output `build\release-windows\bin\NoSpherA2.exe`. `Windows\Tests\Tests.vcxproj` was not repointed and still builds to `Windows\Tests\$(Configuration)_$(Platform)\`.

## Dependency Boundaries

Do not edit submodules casually. Prefer parent-project CMake or generated-source workarounds when fixing integration issues from the parent repository.

The parent project currently rewrites selected `libcint` source files under the build tree to make `#pragma GCC ivdep` portable to MSVC. Preserve this parent-side approach unless the user explicitly asks for a submodule fix.
