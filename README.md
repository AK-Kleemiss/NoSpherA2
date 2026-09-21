# NoSpherA2

![Build](https://github.com/AK-Kleemiss/NoSpherA2/actions/workflows/c-cpp_all.yml/badge.svg)
[![DOI](https://img.shields.io/badge/DOI-10.1039/D0SC05526C-blue.svg)](https://doi.org/10.1039/D0SC05526C)
[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/10849/badge)](https://www.bestpractices.dev/projects/10849)

This repository exists for NoSpherA2, a software to calculate a .tsc file and property grid files from a given wavefunction. The .tsc file is a format to be used in olex2.refine for the use of non-spherical scattering- or form-factors during the refinement, enabling refinements like HAR (Hirshfeld Atom Refinement).

NoSpherA2 is published in Chem. Sci., 2021, 12, 1675-1692 (https://pubs.rsc.org/en/content/articlelanding/2021/sc/d0sc05526c).

Olex2 is provided free of charge by OlexSys Ltd. (https://www.olexsys.org)

The software is provided as-is under the BSD-2 licence. Please see [LICENSE](./LICENSE) for further details!

## Building NoSpherA2

NoSpherA2 relies heavily on the submodule [featomic](https://github.com/metatensor/featomic) and the [mdSpan](https://github.com/kokkos/mdspan/tree/d34b447fbfdddfad63d2204923917e889ebe2e20) reference implementation. To clone this repository with all of its dependencies do:

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
- `NOSPHERA2_USE_HIP` (default: `OFF`): compile the HIP GPU paths. Together with `NOSPHERA2_USE_CUDA` it produces one binary carrying both backends; at run time the backend with a device is used (`NOSPHERA2_GPU_BACKEND=cuda|hip` in the environment overrides the choice).
- `NOSPHERA2_CUDA_PORTABLE` / `NOSPHERA2_HIP_PORTABLE` (default: `OFF`): compile for every supported NVIDIA / AMD architecture instead of only the build machine's GPU. Use this for a binary distributed to different GPUs.
- `NOSPHERA2_USE_CUTLASS` (default: `ON`): use CUTLASS headers for the CUDA single-precision I-tensor GEMM. It has no additional runtime dependency.

The Visual Studio solution under `Windows/` makes the same choices from the installed toolkits; see [GPU builds with the solution](#gpu-builds-with-the-solution).

The CUDA runtime is linked statically and cuBLAS is loaded only when `-gpu_cublas` is requested and a matching library is installed. A CUDA-enabled executable therefore has no required CUDA DLL import and starts normally on a machine without an NVIDIA GPU; GPU requests fall back to the CPU when no usable device is present. HIP builds do not link the HIP runtime either: `libamdhip64.so` / `amdhip64_<major>.dll` is opened by name on the first GPU call and treated as "no device" when it is absent, so the same executable starts on a machine without ROCm. No CUDA, HIP, or cuBLAS runtime is copied into the executable directory. The CI artifacts `NoSpherA2-linux-x86_64-gpu` and `NoSpherA2-windows-x64-gpu` are such combined CUDA+HIP builds for every supported architecture.

At runtime, supported GPU paths are enabled by default: Fourier transforms, XCW I-tensor contractions, SALTED descriptor combinations, and Becke/TFVC atomic-grid weights. Each automatically falls back to the CPU if the device, memory budget, or input layout is unsuitable. Use `-no_gpu` to disable every GPU path, or `-no_gpu_grid`, `-no_gpu_itensor`, and `-no_gpu_salted` to pin an individual calculation to the CPU. `-gpu_blas` remains opt-in because its transfers only pay off for sufficiently large dense matrix products.

On Windows, configure from the x64 Visual Studio developer environment and use a CUDA/toolset combination supported by NVIDIA. If CMake cannot validate the selected CUDA compiler, configuration continues as a CPU-only build and reports the fallback; choose a compatible MSVC toolset or CUDA release before relying on GPU support.

#### Step 3: Build the project

```sh
cmake --build --preset <preset>
```

The final executable will be located in the preset build directory, e.g. `build/release-windows/bin/NoSpherA2.exe`.

---

### Windows Build Instructions

If you want to develop NoSpherA2 on Windows using Visual Studio, initialize the micromamba environment as described above. Then you **have** to run the following command in a Developer PowerShell / terminal to install the dependencies:

#### Prerequisites

- **Visual Studio 2022** or newer with the C++ build tools

#### Command

```powershell
cmake -P scripts/SetupVSEnvironment.cmake
```

It configures `build/release-windows` and `build/debug-windows` and installs the dependencies (`libcint`, `occ`, `featomic`, MKL, ...) into `deps-install-release` / `deps-install-debug`, which the projects under `Windows/` read through the `NoSpherA2_release.props` / `NoSpherA2_debug.props` sheets in `Windows/Windows_utils/`. Run it again whenever a dependency changes; a solution that suddenly fails to compile against `occ` or `libcint` after a `git pull` usually just needs this refresh. The script reconfigures the CMake build directories, so do not run it while a CMake build is in progress.

Then open `Windows/NoSpherA2/NoSpherA2.sln` or build it from the same developer shell:

```powershell
msbuild Windows\NoSpherA2\NoSpherA2.sln /p:Configuration=Release /p:Platform=x64 /m
```

The solution holds four projects: `NoSpherA2_LIB` (everything under `Src/core`), `NoSpherA2` (the executable), `NoSpherA2_DLL` (the Olex2 DLL) and `Tests` (the GTest suite). The executable and the DLL land in `build\Release_x64\`, the tests in `Windows\Tests\Release_x64\Tests.exe`. The projects list their sources by hand while CMake globs them, so a `.cpp` or `.cu` added on the CMake side has to be added to the `.vcxproj` as well; `python scripts/check_vcxproj.py` reports the difference and runs as the first CI step.

#### GPU builds with the solution

The solution picks its GPU backends from the machine, the way the CMake build does with `NOSPHERA2_GPU_AUTO`. `Windows/Windows_utils/NoSpherA2_gpu.props` decides at build time:

- **CUDA** is used when `CUDA_PATH` names a toolkit whose Visual Studio integration is installed, i.e. `$(VCTargetsPath)\BuildCustomizations\CUDA <version>.props` exists. The CUDA installer adds it when the Visual Studio integration component is selected; a toolkit installed after Visual Studio, or into a different Visual Studio version, is not seen and the build silently stays on the CPU. The `.cu` files are compiled by `nvcc` with the static CUDA runtime, so the executable has no CUDA DLL import.
- **HIP** is used when `HIP_PATH` names a ROCm / HIP SDK for Windows with `bin\clang++.exe`. The `.cu` files are compiled a second time as HIP for `gfx90a`, `gfx942` and `gfx1100`; the HIP runtime is not linked but opened by name at run time (`Src/core/hip_runtime_shim.cpp`).
- **Both** installed give one binary with both backends, compiled into separate namespaces and joined by `Src/core/gpu_dispatch.cpp`; at run time the backend with a device is used.
- **Neither** gives a CPU-only build. Nothing has to be configured for it.

Environment variables that steer the choice, set in the shell that runs `msbuild` or in the Visual Studio session:

- `NOS_GPU=OFF` builds for the CPU only even when a toolkit is installed.
- `NOS_CUDA_ARCH` overrides the CUDA architectures with a `CodeGeneration` string such as `compute_86,sm_86;compute_89,sm_89`; `NOS_CUDA_ARCH=all` compiles the portable list used by the CI artifacts (`sm_70` to `sm_90` for CUDA 12, `sm_75` to `sm_120` for CUDA 13). Without it `nvidia-smi` names the cards in the machine (`Windows/Windows_utils/detect_cuda_arch.ps1`), and the portable list is the fallback when there is none.
- `CUTLASS_DIR` points at a CUTLASS checkout for the single-precision I-tensor GEMM; without it the property sheet looks for the copy the CMake configure fetched under `build\*\_deps\cutlass-src`, and without that the built-in GEMM is compiled.

The build log says what was chosen: `NoSpherA2 CUDA 13.4: compute_75,sm_75, CUTLASS <dir>` from the LIB project, `NoSpherA2 GPU: fat binary, CUDA and HIP (<path>)` or `NoSpherA2 GPU: none (CPU build)`. As with the CMake build, a GPU-enabled executable starts on a machine without a device and falls back to the CPU; `-no_gpu` and the per-path `-no_gpu_*` options apply unchanged.

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

The GTest binary can also be run on its own from `tests/src`, e.g. `..\..\build\release-windows\bin\NoSpherA2_Tests.exe --gtest_filter=Qct*`, or the solution's `Windows\Tests\Release_x64\Tests.exe` in the same way.

### Building a macOS Universal Binary

To build a universal binary for macOS, you can use the following command:

```bash
cmake --preset release-macos-universal
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
   acc = 0
   ```
4. If the reference output is not named after the test, add the `good` parameter:

   ```toml
   [disorder_THPP]
   directory = "disorder"
   good = "disorder_THPP.good"

   [disorder_THPP.args]
   cif = "thpp.cif"
   hkl = "thpp.hkl"
   acc = 0
   ```
5. Command line arguments are always passed in the block `<testname>.args`.
6. **Run `pytest` (or `ctest`)** to ensure your test runs.
7. **Document your test** if needed and mention it in your Pull Request.

---
