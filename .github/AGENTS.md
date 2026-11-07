# AGENTS.md – Project-wide AI assistant instructions

## Quick overview
- **Project**: **NoSpherA2** – a C++ CLI tool that generates `.tsc` files and property grids from quantum-chemistry wavefunctions.
- **Key submodules**: `occ/` (quantum-chemistry core), `featomic/` (Rust ML features), `libcint/`, `mdspan/`.
- **Build system**: CMake presets (Ninja) are the only build entrypoint across platforms. Dependencies are handled via a micromamba bootstrap step (`cmake -P scripts/BootstrapMicromamba.cmake`).

## Build & test commands (agents can run these automatically)
| Platform | Command | What it does |
|----------|---------|--------------|
| **Windows** | `cmd /c "\"%VSINSTALLDIR%\Common7\Tools\VsDevCmd.bat\" -arch=x64 && cmake --preset release-windows && cmake --build --preset release-windows"` | Initializes the MSVC x64 environment, then configures + builds NoSpherA2. |
| **Linux** | `cmake --preset release-linux && cmake --build --preset release-linux` | Configures + builds NoSpherA2. |
| **macOS** | `cmake --preset release-macos-arm64 && cmake --build --preset release-macos-arm64` | Configures + builds NoSpherA2 for a specific arch. |
| **All** | `ctest --preset <preset>` | Runs the configured CTest suite for the selected preset. |

Use `cmake --list-presets=all` for the full preset list; host-conditional presets hide non-host entries from a plain `--list-presets`.

## Important build notes for AI agents
- **Submodule handling**: run `git submodule update --init --recursive` before any build.
- **Static linking**: most third-party libraries (`featomic`, `libcint`, `occ`) are linked statically. OCC links against shared oneTBB, so the matching TBB runtime is copied beside the executable (`tbb12.dll`, `libtbb.so*`, or `libtbb.dylib`). Do not build or link `tbbmalloc_proxy` for NoSpherA2.
- **Windows environment**: initialize the VS Developer environment (`VsDevCmd.bat -arch=x64`) in the same shell before running CMake/CTest. Missing standard headers (`cstddef`, `windows.h`) or `kernel32.lib` mean the environment was not initialized.
- **IDE work**: the `.vscode/` setup (CMake Tools with the presets) is the supported IDE flow on Windows, Linux and macOS. The Visual Studio solution was removed in October 2026.
- **In-process tests**: tests documented as in-process (e.g. `alanine_integrated_occ`) stay in-process in the gtest binary; do not add subprocess fallbacks.

## Common pitfalls (agents should warn the user)
- Golden-file test parsing can fail on lines with multiple floating-point values or scientific notation. See `tests/run_test.py` for the parser.
- On macOS, build per-arch using `release-macos-arm64` / `release-macos-x86_64`.

## Helpful documentation links (link, don't duplicate)
- [Project overview & build instructions](../CLAUDE.md)
- [README – high-level description and usage](../README.md)
- [Test inventory and status](../UNIT_TESTS_STATUS.md)
- [CMake presets for targeted builds](../CMakePresets.json)

---
*These instructions are kept minimal; agents should refer to the linked files for detailed guidance.*
