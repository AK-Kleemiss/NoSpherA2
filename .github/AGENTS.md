# AGENTS.md – Project‑wide AI assistant instructions

## Quick overview
- **Project**: **NoSpherA2** – a C++ CLI tool that generates `.tsc` files and property grids from quantum‑chemistry wavefunctions.
- **Key submodules**: `occ/` (quantum‑chemistry core), `featomic/` (Rust ML features), `libcint/`, `mdspan/`.
- **Build system**: CMake presets are the primary build entrypoint across platforms. Dependencies are handled via a micromamba bootstrap step during configure.

## Build & test commands (agents can run these automatically)
| Platform | Command | What it does |
|----------|---------|--------------|
| **Windows** | `cmd /c "\"%VSINSTALLDIR%\Common7\Tools\VsDevCmd.bat\" -arch=x64 && cmake --preset release-windows && cmake --build --preset release-windows"` | Initializes the MSVC x64 environment, then configures + builds NoSpherA2. |
| **Linux** | `cmake --preset linux-gcc && cmake --build --preset linux-gcc` | Configures + builds full NoSpherA2. |
| **macOS** | `cmake --preset macos-release-full-arm64 && cmake --build --preset macos-release-full-arm64` | Configures + builds full NoSpherA2 for a specific arch. |
| **All** | `ctest --preset <preset>` | Runs the configured CTest suite for the selected preset. |
| **Windows native tests** | `vstest.console.exe Windows\x64\<Config>\Tests.dll /Platform:x64` | Use Visual Studio test tools for native/in-process debugging. Ensure the test binaries are built and the appropriate test adapter is installed.

## Important build notes for AI agents
- **OCC static library** (`Lib/occ/lib/libocc.a`) is built by the `occ` target. OCC links against shared oneTBB, so package the matching TBB runtime beside the executable (`tbb12.dll`, `libtbb.so*`, or `libtbb.dylib`).
- **TBB runtime**: Do not build or link `tbbmalloc_proxy` for NoSpherA2. The build scripts copy the core TBB runtime from `Lib/occ` into the executable/artifact layout.
- **Submodule handling**: The repository relies on several git submodules. Before any build, run `git submodule update --init --recursive`.
- **Static linking**: Most third‑party libraries (`featomic`, `libcint`, `occ`) are linked statically, with OpenMP/MKL and TBB runtime files deployed alongside the executable where needed.
- **Visual Studio version 18**: Windows builds must always initialize the VS Developer environment (VS version 18 in the current environment) before running CMake, CTest, MSBuild, or VSTest. Ensure the command invokes `VsDevCmd.bat -arch=x64` first, in the same shell/session as the build command.
- **Windows CMake presets**: use `release-windows` and `debug-windows`. If CMake or Ninja cannot find standard headers such as `cstddef`, `vector`, `windows.h`, or libraries such as `kernel32.lib`, the MSVC/Windows SDK environment was not initialized; rerun after `VsDevCmd.bat -arch=x64`.
- **Test execution on Windows**: Use Visual Studio test tools (`vstest.console.exe`) when debugging native/OCC behavior. Also keep the Python pytest harness passing because CI and golden-file coverage depend on it.
  - Visual Studio tests must run NoSpherA2 in-process through `NoSpherA2_DLL.dll`. Do not add subprocess fallbacks for failing VS tests, including `alanine_integrated_occ`; keep them natively debuggable in the VSTest host and fix the underlying in-process issue.
 - **OCC build process**: The OCC source resides in `occ/`. Build it via CMake presets (e.g. `linux-occ-gcc`) or as part of the full-build presets.

## Common pitfalls (agents should warn the user)
- Golden‑file test parsing can fail on lines with multiple floating‑point values or scientific notation. See `tests/run_test.py` for the parser.
- On macOS, build per-arch using `macos-release-full-arm64` / `macos-release-full-x86_64`.
- Windows builds require the VS 2022 Developer PowerShell (or the appropriate environment) so that `msbuild` and the C++ toolchain are on `PATH`.
- Current validated status as of 2026-06-13: Python pytest Release full `21 passed` with one tolerated `fchk_conversion` numeric warning; Python pytest Debug full `21 passed`; Visual Studio native Debug full `21 passed`; Visual Studio native Release full `21 passed`. `Release+Copy` was intentionally not rerun in the latest pass.

## Helpful documentation links (link, don’t duplicate)
- [Project overview & build instructions](CLAUDE.md)
- [README – high‑level description and usage](README.md)
- [Hand‑off notes for contributors](HANDOFF.md)
- [Investigation status & open issues](INVESTIGATION_STATUS.md)
- [CMake presets for targeted builds](CMakePresets.json)

---
*These instructions are kept minimal; agents should refer to the linked files for detailed guidance.*
