# AGENTS.md – Project‑wide AI assistant instructions

## Quick overview
- **Project**: **NoSpherA2** – a C++ CLI tool that generates `.tsc` files and property grids from quantum‑chemistry wavefunctions.
- **Key submodules**: `occ/` (quantum‑chemistry core), `featomic/` (Rust ML features), `libcint/`, `mdspan/`.
- **Build system**: top‑level `Makefile` dispatches to platform‑specific builds (Windows via MSBuild, Linux/macOS via make). All dependencies are built statically and placed under `Lib/`.

## Build & test commands (agents can run these automatically)
| Platform | Command | What it does |
|----------|---------|--------------|
| **Windows** | `make.exe` (or `make` from a VS 2022 Developer PowerShell) | Builds the solution (`Windows/NoSpherA2.sln`) and produces `NoSpherA2.exe` at the repo root. |
| **Linux / macOS** | `make` | Invokes the appropriate makefile in `Linux/` or `Mac/` and produces the executable. |
| **All** | `make test` or `make tests` | Builds the executable (if needed) and runs the Python pytest harness (`uv run pytest`). |
| **Windows (new version)** | `vstest.console.exe` (or `dotnet test` if using .NET) | Use Visual Studio test tools to run the test suite instead of the `uv` Python path. Ensure the test binaries are built and the appropriate test adapter is installed.

## Important build notes for AI agents
- **OCC static library** (`Lib/occ/lib/libocc.a`) is built by the `occ` target. By default it **includes** the TBB scalable allocator (`tbbmalloc`). The Windows Visual Studio project also links `tbbmalloc_proxy.lib` and `tbbmalloc.lib` explicitly – this is optional but harmless.
- **tbbmalloc**: The top‑level `Makefile` only defines a `tbbmalloc` target for macOS and Linux. On Windows the allocator is provided by the *Build Dependencies* step and linked via the `.vcxproj`. Agents should not attempt to build `tbbmalloc` on Windows; just ensure the `Lib/occ/` binaries are present.
- **Submodule handling**: The repository relies on several git submodules. Before any build, run `git submodule update --init --recursive`.
- **Static linking**: All third‑party libraries (`featomic`, `libcint`, `occ`, `tbbmalloc`) are linked statically, producing a portable executable.
- **Visual Studio version 18**: Windows builds must always start the VS 2022 Developer Command Prompt (VS version 18) to set up the compiler, linker, and environment variables. Ensure the build step invokes `VsDevCmd.bat` (e.g., `call "%VSINSTALLDIR%\Common7\Tools\VsDevCmd.bat"`) before running `msbuild`.
- **Test execution on Windows**: For this local version we prefer using Visual Studio test tools (e.g., `vstest.console.exe`) instead of the `uv` Python test runner. Adjust CI or local test scripts accordingly.
 - **OCC build process**: The OCC source resides in `occ/`. Building OCC is performed in a dedicated `build*` directory at the repository root (e.g., `build-macos-release-arm64`, `build-linux-occ-gcc`). After running CMake configuration and `make install`, the resulting libraries are installed into `Lib/occ/`. Whenever OCC source files are modified, the full configure‑build‑install cycle must be rerun to update the binaries in `Lib/occ/` before rebuilding NoSpherA2.

## Common pitfalls (agents should warn the user)
- Golden‑file test parsing can fail on lines with multiple floating‑point values or scientific notation. See `tests/run_test.py` for the parser.
- On macOS, the `NATIVE_ARCH` detection may pick the wrong architecture; ensure `make` is invoked with the correct `NATIVE_ARCH` if you need a specific binary.
- Windows builds require the VS 2022 Developer PowerShell (or the appropriate environment) so that `msbuild` and the C++ toolchain are on `PATH`.

## Helpful documentation links (link, don’t duplicate)
- [Project overview & build instructions](CLAUDE.md)
- [README – high‑level description and usage](README.md)
- [Hand‑off notes for contributors](HANDOFF.md)
- [Investigation status & open issues](INVESTIGATION_STATUS.md)
- [CMake presets for targeted builds](CMakePresets.json)

---
*These instructions are kept minimal; agents should refer to the linked files for detailed guidance.*