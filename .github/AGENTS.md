# AGENTS.md – Project‑wide AI assistant instructions

## Quick overview
- **Project**: **NoSpherA2** – a C++ CLI tool that generates `.tsc` files and property grids from quantum‑chemistry wavefunctions.
- **Key submodules**: `occ/` (quantum‑chemistry core), `featomic/` (Rust ML features), `libcint/`, `mdspan/`.
- **Build system**: top‑level `Makefile` dispatches to platform‑specific builds (Windows via MSBuild, Linux/macOS via make). Most dependencies are built statically and placed under `Lib/`; TBB is deployed as a runtime library beside the executable.

## Build & test commands (agents can run these automatically)
| Platform | Command | What it does |
|----------|---------|--------------|
| **Windows** | `make.exe` (or `make` from a VS 2022 Developer PowerShell) | Builds the solution (`Windows/NoSpherA2.sln`) and produces `NoSpherA2.exe` at the repo root. |
| **Linux / macOS** | `make` | Invokes the appropriate makefile in `Linux/` or `Mac/` and produces the executable. |
| **All** | `make test` or `make tests` | Builds the executable (if needed) and runs the Python pytest harness (`uv run pytest`). |
| **Windows native tests** | `vstest.console.exe Windows\x64\<Config>\Tests.dll /Platform:x64` | Use Visual Studio test tools for native/in-process debugging. Ensure the test binaries are built and the appropriate test adapter is installed.

## Important build notes for AI agents
- **OCC static library** (`Lib/occ/lib/libocc.a`) is built by the `occ` target. OCC links against shared oneTBB, so package the matching TBB runtime beside the executable (`tbb12.dll`, `libtbb.so*`, or `libtbb.dylib`).
- **TBB runtime**: Do not build or link `tbbmalloc_proxy` for NoSpherA2. The build scripts copy the core TBB runtime from `Lib/occ` into the executable/artifact layout.
- **Submodule handling**: The repository relies on several git submodules. Before any build, run `git submodule update --init --recursive`.
- **Static linking**: Most third‑party libraries (`featomic`, `libcint`, `occ`) are linked statically, with OpenMP/MKL and TBB runtime files deployed alongside the executable where needed.
- **Visual Studio version 18**: Windows builds must always start the VS Developer Command Prompt (VS version 18 in the current environment) to set up the compiler, linker, and environment variables. Ensure the build step invokes `VsDevCmd.bat` (e.g., `call "%VSINSTALLDIR%\Common7\Tools\VsDevCmd.bat"`) before running `msbuild`.
- **Windows CMake presets**: `windows-clang-cl` and `windows-msvc-debug` are both valid targeted flows. The `windows-msvc-debug` workflow preset is required by the Windows dependency build action.
- **Test execution on Windows**: Use Visual Studio test tools (`vstest.console.exe`) when debugging native/OCC behavior. Also keep the Python pytest harness passing because CI and golden-file coverage depend on it.
  - Visual Studio tests must run NoSpherA2 in-process through `NoSpherA2_DLL.dll`. Do not add subprocess fallbacks for failing VS tests, including `alanine_integrated_occ`; keep them natively debuggable in the VSTest host and fix the underlying in-process issue.
 - **OCC build process**: The OCC source resides in `occ/`. Building OCC is performed in a dedicated `build*` directory at the repository root (e.g., `build-macos-release-arm64`, `build-linux-occ-gcc`). After running CMake configuration and `make install`, the resulting libraries are installed into `Lib/occ/`. Whenever OCC source files are modified, the full configure‑build‑install cycle must be rerun to update the binaries in `Lib/occ/` before rebuilding NoSpherA2.

## Common pitfalls (agents should warn the user)
- Golden‑file test parsing can fail on lines with multiple floating‑point values or scientific notation. See `tests/run_test.py` for the parser.
- On macOS, the `NATIVE_ARCH` detection may pick the wrong architecture; ensure `make` is invoked with the correct `NATIVE_ARCH` if you need a specific binary.
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
