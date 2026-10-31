# NoSpherA2

### XCW per-element orbital bases

The XCW settings keyword `basis_overrides <file>` reads an all-electron BSE JSON
file, replacing the orbital shells for its elements while retaining `basis_set`
for all other elements. Relative paths are resolved beside the settings file;
paths containing spaces can be quoted. General and combined angular-momentum
contractions are split into individual spherical shells. The assembled basis
is shared by the SCF, atomic grids and output wavefunction. ECP overrides are
rejected. Olex2's XCW advanced-basis controls can fetch `SARC-DKH2` from BSE for
selected heavy elements; light elements retain the fallback basis.

This supplies SARC-DKH2 basis functions, not a DKH2 Hamiltonian: XCW's
Hartree-Fock Hamiltonian remains non-relativistic.

![Build](https://github.com/AK-Kleemiss/NoSpherA2/actions/workflows/c-cpp_all.yml/badge.svg)
[![DOI](https://img.shields.io/badge/DOI-10.1039/D0SC05526C-blue.svg)](https://doi.org/10.1039/D0SC05526C)
[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/10849/badge)](https://www.bestpractices.dev/projects/10849)

This repository exists for NoSpherA2, a software to calculate a .tsc file and property grid files from a given wavefunction. The .tsc file is a format to be used in olex2.refine for the use of non-spherical scattering- or form-factors during the refinement, enabling refinements like HAR (Hirshfeld Atom Refinement).

NoSpherA2 is published in Chem. Sci., 2021, 12, 1675-1692 (https://pubs.rsc.org/en/content/articlelanding/2021/sc/d0sc05526c).

Olex2 is provided free of charge by OlexSys Ltd. (https://www.olexsys.org)

The software is provided as-is under the BSD-2 licence. Please see [LICENSE](./LICENSE) for further details!

## Acknowledgments
This software makes extensive use of the following external packages:
- [libcint](https://github.com/sunqm/libcint)
- [mdspan](https://github.com/kokkos/mdspan/tree/stable)
- [featomic](https://github.com/metatensor/featomic)
- [occ](https://github.com/peterspackman/occ)
  - OCC was extensively used in XCW and is the engine used for calculating WFNs inside NoSpherA2. Please cite it if those functionalities were used: Spackman, P. R. (2026). Open Computational Chemistry (OCC) - A portable software library and program for quantum chemistry and crystallography. *J. Open Source Softw.* 11(117), 9609, [doi:10.21105/joss.09609](https://doi.org/10.21105/joss.09609).

## Building NoSpherA2
To clone this repository with all of its dependencies do:

```sh
git clone --recursive https://github.com/AK-Kleemiss/NoSpherA2.git
```

If you already cloned it without submodules, update them manually:

```sh
git submodule update --init --recursive
```

#### Prerequisites

- **CMake** (min 3.25)
- **C++20 compiler** (MSVC 19.35, GCC 12.2, Clang 15.0)

#### Step 1: Bootstrap the micromamba environment

Please run the following command to bootstrap the micromamba environment. This will download and install micromamba, and create a local environment with all dependencies needed to build NoSpherA2.

```bash
cmake -P scripts/BootstrapMicromamba.cmake
```

#### Step 2: Configure the project

NoSpherA2 uses CMake presets to simplify the build process.

```sh
cmake --preset <preset>
```

Common presets include:

- Windows: `release-windows`, `debug-windows`
- Linux: `release-linux`, `debug-linux`
- macOS (per-arch): `release-macos-arm64`, `release-macos-x86_64`, `debug-macos-arm64`, `debug-macos-x86_64`

##### Configuration options

Pass options while configuring a preset, for example:

```sh
cmake --preset release-windows -DNOSPHERA2_BUILD_TESTS=ON
```

- `NOSPHERA2_BUILD_TESTS` (default: `OFF`): build the CTest/GTest suite.
- `NOSPHERA2_BUILD_DLL` (default: `OFF`): build the optional Windows DLL.
- `NOSPHERA2_DEPENDENCIES_ONLY` (default: `OFF`): build dependency targets only; this is intended for dependency-cache preparation.
- `NOSPHERA2_GPU_AUTO` (default: `ON`): select CUDA for an NVIDIA driver or HIP for an AMD driver detected on the build host. Set it to `OFF` when choosing a backend explicitly.
- `NOSPHERA2_USE_CUDA` (default: `OFF`): compile the CUDA GPU paths. This can be enabled explicitly on a GPU-less build host when a CUDA compiler is available.
- `NOSPHERA2_USE_HIP` (default: `OFF`): compile the HIP GPU paths. Together with `NOSPHERA2_USE_CUDA` it produces one binary carrying both backends; at run time the backend with a device is used (`-tune NOSPHERA2_GPU_BACKEND=cuda|hip` overrides the choice).
- `NOSPHERA2_CUDA_PORTABLE` / `NOSPHERA2_HIP_PORTABLE` (default: `OFF`): compile for every supported NVIDIA / AMD architecture instead of only the build machine's GPU. Use this for a binary distributed to different GPUs.
- `NOSPHERA2_USE_CUTLASS` (default: `ON`): use CUTLASS headers for the CUDA single-precision I-tensor GEMM. It has no additional runtime dependency.

The CUDA runtime is linked statically and cuBLAS is opened by name on its first use rather than linked: it is on by default (`-no_gpu_cublas` pins the I-tensor GEMM to CUTLASS or the built-in kernel, which is what the reference tests do) and treated as absent when no matching library is installed. A CUDA-enabled executable therefore has no required CUDA DLL import and starts normally on a machine without an NVIDIA GPU; GPU requests fall back to the CPU when no usable device is present. HIP builds do not link the HIP runtime either: `libamdhip64.so` / `amdhip64_<major>.dll` is opened by name on the first GPU call and treated as "no device" when it is absent, so the same executable starts on a machine without ROCm. No CUDA, HIP, or cuBLAS runtime is copied into the executable directory. The CI artifacts `NoSpherA2-linux-x86_64-gpu` and `NoSpherA2-windows-x64-gpu` are such combined CUDA+HIP builds for every supported architecture.

At runtime, supported GPU paths are enabled by default: Fourier transforms, XCW I-tensor contractions, SALTED descriptor combinations, and Becke/TFVC atomic-grid weights. Each automatically falls back to the CPU if the device, memory budget, or input layout is unsuitable. Use `-no_gpu` to disable every GPU path, or `-no_gpu_grid`, `-no_gpu_itensor`, and `-no_gpu_salted` to pin an individual calculation to the CPU. `-gpu_blas` remains opt-in because its transfers only pay off for sufficiently large dense matrix products.

On Windows, configure from the x64 Visual Studio developer environment and use a CUDA/toolset combination supported by NVIDIA. If CMake cannot validate the selected CUDA compiler, configuration continues as a CPU-only build and reports the fallback; choose a compatible MSVC toolset or CUDA release before relying on GPU support.

On Linux aarch64 and Windows ARM64 the [Arm Performance Libraries](https://learn.arm.com/install-guides/armpl/) replace OpenBLAS whenever one is installed (`-DNOSPHERA2_USE_ARMPL=OFF` keeps OpenBLAS). CMake looks in `NOSPHERA2_ARMPL_DIR`, `.mambaenv/armpl`, `ARMPL_DIR`, `~/armpl`, `/opt/arm` and `Program Files\Arm Performance Libraries`, and links ArmPL statically, so nothing extra ships. Arm's site refuses scripted downloads: fetch the package by hand (on Windows the `_Windows.msi`, not the `_Windows_Arm64EC.msi`) and either name it in `NOSPHERA2_ARMPL_PACKAGE` or drop it into `.mambaenv/bootstrap`, and the bootstrap installs it into `.mambaenv/armpl`, which accepts Arm's licence.

#### Step 3: Build the project

```sh
cmake --build --preset <preset>
```

The final executable will be located in the preset build directory, e.g. `build/release-windows/bin/NoSpherA2.exe`.

---

### Developing with VS Code

The repository ships a VS Code setup (`.vscode/`) that works the same on Windows, Linux and macOS. Install the recommended extensions (C/C++ and CMake Tools; VS Code offers them on first open), then:

1. Pick a configure preset in the CMake Tools status bar (e.g. `debug-windows`, `debug-linux`, `debug-macos-arm64`) — run the micromamba bootstrap first, as described above.
2. Build with F7. On Windows, CMake Tools sets up the MSVC environment itself (`cmake.useVsDeveloperEnvironment`), so VS Code does not have to be started from a developer prompt; the integrated terminal profile "VS x64 Developer PowerShell" gives the same environment for running `cmake`/`ctest` by hand.
3. Pick a launch target (`NoSpherA2` or `NoSpherA2_Tests`) and debug with the "CMake target" launch configuration — `cppvsdbg` on Windows, gdb/lldb on Linux/macOS. Set the run's arguments and working directory in `.vscode/launch.json`; the "native NBO on epoxide" entry shows the pattern.
4. Tests appear in the Testing panel once configured with `-DNOSPHERA2_BUILD_TESTS=ON`.

IntelliSense comes from the CMake configure (`CMAKE_EXPORT_COMPILE_COMMANDS` is on in every preset), so there is no hand-maintained `c_cpp_properties.json`.

### Run Tests

After building the project using the `NOSPHERA2_BUILD_TESTS` flag, you can run the tests using `ctest`.

Note: the C++ test executables (and some Python tests) require `OCC_DATA_PATH` to point at the OCC runtime data directory (default in this repo: `occ/share`).

```bash
ctest --preset <preset> --output-on-failure
```

Example (Windows, Debug):

```powershell
ctest --preset debug-windows --output-on-failure
```

The GTest binary can also be run on its own from `tests/src`, e.g. `..\..\build\release-windows\bin\NoSpherA2_Tests.exe --gtest_filter=Qct*`.

### Building a macOS Universal Binary

To build a universal binary for macOS, build both per-arch release presets, then:

```bash
cmake --build --preset release-macos-universal
```

---

## Contributing: Pull Requests & Adding Tests

### Creating a Pull Request

1. **Fork the repository** and create a new branch for your feature or bugfix.
2. **Write clear, concise commit messages** describing your changes.
3. **Ensure your code builds and passes all tests** on your platform.
4. **Push your branch** to your fork and open a Pull Request (PR) against the `master` branch of this repository.
5. In your PR description, explain the purpose of your changes and reference any related issues.
6. Be responsive to review feedback and update your PR as needed.

### Adding a Test for New Features

1. **Locate the appropriate test file or directory** (commonly in the `tests/` folder or as specified in the codebase).
2. **Add your test** following the style of existing tests. For C++ code, this may be a new `.cpp` file or a new GTest file under `tests/src/`.
3. **Register your test in `tests/tests.toml`:**

   - Add a block for it; the `directory` names the subfolder of `tests/` holding the input files.
   - Add your test executable or call to the list of tests, as appropriate.
   - Example (a test calculating a tsc file for sucrose reading in an hkl file and wfn):

   ```toml
   [sucrose_SF]
   directory = "sucrose_fchk_SF"

   [sucrose_SF.args]
   cif = "sucrose.cif"
   hkl = "olex2/Wfn_job/sucrose.hkl"
   wfn = "olex2/Wfn_job/sucrose.wfx"
   acc = 1
   ```
4. If the reference output is not named after the test, add the `good` parameter:

   ```toml
   [disorder_THPP]
   directory = "disorder"
   good = "disorder_THPP.good"

   [disorder_THPP.args]
   cif = "thpp.cif"
   hkl = "thpp.hkl"
   acc = 1
   ```
5. Command line arguments are always passed in the block `<testname>.args`.
6. **Run `pytest` (or `ctest`)** to ensure your test runs.
7. **Document your test** if needed and mention it in your Pull Request.

---
