# NoSpherA2

![Build](https://github.com/AK-Kleemiss/NoSpherA2/actions/workflows/c-cpp_all.yml/badge.svg)
[![DOI](https://img.shields.io/badge/DOI-10.1039/D0SC05526C-blue.svg)](https://doi.org/10.1039/D0SC05526C)
[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/10849/badge)](https://www.bestpractices.dev/projects/10849)

This repository exists for NoSpherA2, a software to calculate a .tsc file and property grid files form a given wavefunction. The .tsc file is a format to be used in olex2.refine for the use of non-spherical scattering- or form-factors during the refinement, enabling refinemnts like HAR (Hirshfeld-Atom-Refinement).

NoSpherA2 is published in Chem. Sci., 2021,12, 1675-1692 (https://pubs.rsc.org/en/content/articlelanding/2021/sc/d0sc05526c).

Olex2 is free for use by OlexSys Ltd. (www.olxsys.org)

The software is provided as-is under the BSD-2 licence. Please see [LICENSE](./LICENSE) for further details!

## Building NoSpherA2

NoSpherA2 relies heavily on the submodule [featomic](https://github.com/metatensor/featomic) and the [mdSpan](https://github.com/kokkos/mdspan/tree/d34b447fbfdddfad63d2204923917e889ebe2e20) reference implementation. To clone this repository with all of it's dependencies do:

```sh
git clone --recursive https://github.com/AK-Kleemiss/NoSpherA2.git
```

If you already cloned it without submodules, update them manually:

```sh
git submodule update --init --recursive
```

#### **Prerequisites**

- **CMake** (min 3.25)
- ***C++20 compiler** (MSVC 19.35, GCC 12.2, Clang 15.0)

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

##### Configuration Options:

`NOSPHERA2_BUILD_TESTS` (default: `OFF`) - Build the test suite.

#### Step 3: Build the project

```sh
cmake --build --preset <preset>
```

The final executable will be located in the preset build directory, e.g. `build/release-windows/bin/NoSpherA2.exe`.

---

### Windows Build Instructions

If you want to develop NoSpherA2 on Windows using Visual Studio please also initialize the micromamba environment as described above. Then you **have** to run the following command in a Developer PowerShell / terminal to install the dependencies:

#### **Prerequisites**

- **Visual Studio 2022** with C++ build tools

### **Command**

```powershell
cmake -P scripts/SetupVSEnvironment.cmake
```

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

### Building a MacOS Universal Binary
To build a universal binary for macOS, you can use the following command:

```bash
cmake --preset release-macos-universal
```   

---

## Contributing: Pull Requests & Adding Tests

### Creating a Pull Request

1. **Fork the repository** and create a new branch for your feature or bugfix.ö
2. **Write clear, concise commit messages** describing your changes.
3. **Ensure your code builds and passes all tests** on your platform.
4. **Push your branch** to your fork and open a Pull Request (PR) against the `main` branch of this repository.
5. In your PR description, explain the purpose of your changes and reference any related issues.
6. Be responsive to review feedback and update your PR as needed.

### Adding a Test for New Features

1. **Locate the appropriate test file or directory** (commonly in the `tests/` folder or as specified in the codebase).
2. **Add your test** following the style of existing tests. For C++ code, this may be a new `.cpp` file or a new function in `Src/test_functions.h` file.
3. **Register your test in the tests.toml file:**

   - Open the `tests.toml` in the test subdirectory.
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
4. If your test requires a different output file from the test name, you can add the good parameter:

   ```toml
   [disorder_THPP]
   directory = "disorder"
   good = "disorder_THPP.good"

   [disorder_THPP.args]
   cif = "thpp.cif"
   hkl = "thpp.hkl"
   acc = 0
   ```
5. Command line arguments are always passed in the block `<testname>`.args.
6. **Run `pytest` (or `ctest`)** to ensure your test runs.
7. **Document your test** if needed and mention it in your Pull Request.

---
