# NoSpherA2

![Build](https://github.com/AK-Kleemiss/NoSpherA2/actions/workflows/c-cpp_all.yml/badge.svg)
[![DOI](https://img.shields.io/badge/DOI-10.1039/D0SC05526C-blue.svg)](https://doi.org/10.1039/D0SC05526C)
[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/10849/badge)](https://www.bestpractices.dev/projects/10849)

This repository exists for NoSpherA2, a software to calculate a .tsc file and property grid files form a given wavefunction. The .tsc file is a format to be used in olex2.refine for the use of non-spherical scattering- or form-factors during the refinement, enabling refinemnts like HAR (Hirshfeld-Atom-Refinement).

NoSpherA2 is published in Chem. Sci., 2021,12, 1675-1692 (https://pubs.rsc.org/en/content/articlelanding/2021/sc/d0sc05526c).

Olex2 is free for use by OlexSys Ltd. (www.olxsys.org)

The software is provided as-is under the BSD-2 licence. Please see [LICENSE](./LICENSE) for further details!

## Building NoSpherA2
NoSpherA2 relies heavily on the submodule [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS/tree/a64b75a2e00691b126a3c342a265f96fac98514f), [llvm](https://github.com/llvm/llvm-project/tree/2db262886f0c06c079e1b2808c4c14c16f8861b5), [featomic](https://github.com/metatensor/featomic) and the [mdSpan](https://github.com/kokkos/mdspan/tree/d34b447fbfdddfad63d2204923917e889ebe2e20) reference implementation. To clone this repository with all of it's dependencies do:

```sh
git clone --recursive https://github.com/AK-Kleemiss/NoSpherA2.git
```

If you already cloned it without submodules, update them manually:

```sh
git submodule update --init --recursive
```

---

### Windows Build Instructions
#### **Prerequisites**
- **Visual Studio 2022** with C++ build tools
- GNU Make for Windows
- cargo for compilation of rust (required for featomic): https://www.rust-lang.org/tools/install

#### **Building NoSpherA2**
Run inside a developers command prompt or power shell (amking sure make.exe is on your PATH):

```ps
make.exe
```

✅ **The executable will be located in the main folder:** `NoSpherA2.exe`.

---

### Linux Build Instructions

#### **Prerequisites**
Ensure `build-essential` and `cargo` are installed (select package manager according to your distro):

```sh
sudo <apt/dnf/yum> update && sudo <apt/dnf/yum> install -y build-essential cargo
```

#### **Building NoSpherA2**
Inside the NoSpherA2 directory, simply run:

```sh
make
```

✅ **The executable will be located in the main folder:** `NoSpherA2`.

---

### macOS Build Instructions

#### **Prerequisites**
- **Xcode Command Line Tools**
- **CMake** (install via Homebrew if needed: `brew install cmake`)
- **rustup** (install via Homebrew: `brew install rustup-init && rustup-init`; might requrie restarting the terminal to load properly)

#### **Building NoSpherA2**

```sh
make
```

✅ **The executable will be located in the main folder:** `NoSpherA2`.

---

## Contributing: Pull Requests & Adding Tests

### Creating a Pull Request

1. **Fork the repository** and create a new branch for your feature or bugfix.
2. **Write clear, concise commit messages** describing your changes.
3. **Ensure your code builds and passes all tests** on your platform.
4. **Push your branch** to your fork and open a Pull Request (PR) against the `main` branch of this repository.
5. In your PR description, explain the purpose of your changes and reference any related issues.
6. Be responsive to review feedback and update your PR as needed.

### Adding a Test for New Features

1. **Locate the appropriate test file or directory** (commonly in the `tests/` folder or as specified in the codebase).
2. **Add your test** following the style of existing tests. For C++ code, this may be a new `.cpp` file or a new function in `Src/test_functions.h` file.
3. **Register your test in the `makefile`:**
    - Open the `makefile` in the test subdirectory.
    - Add your test executable or call to the list of tests, as appropriate.
    - Example (a test calcualting a tsc file for sucrose reading in an hkl file wnd wfn):
      ```
        sucrose_SF:
            $(call RUN_TEST,$@,sucrose_fchk_SF, \
                -cif sucrose.cif \
                -hkl olex2/Wfn_job/sucrose.hkl \
                -wfn olex2/Wfn_job/sucrose.wfx \
                -acc 0 \
                -no-date)
      ```
4. **Run `make test`** to ensure your test compiles and runs.
5. **Document your test** if needed, and mention it in your Pull Request.

---


