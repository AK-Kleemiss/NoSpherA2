# NoSpherA2

![Build](https://github.com/AK-Kleemiss/NoSpherA2/actions/workflows/c-cpp_all.yml/badge.svg)
[![DOI](https://img.shields.io/badge/DOI-10.1039/D0SC05526C-blue.svg)](https://doi.org/10.1039/D0SC05526C)

This repository exists for NoSpherA2, a software to calculate a .tsc file and property grid files form a given wavefunction. The .tsc file is a format to be used in olex2.refine for the use of non-spherical scattering- or form-factors during the refinement, enabling refinemnts like HAR (Hirshfeld-Atom-Refinement).

NoSpherA2 is published in Chem. Sci., 2021,12, 1675-1692 (https://pubs.rsc.org/en/content/articlelanding/2021/sc/d0sc05526c).

Olex2 is free for use by OlexSys Ltd. (www.olxsys.org)

The software is provided as-is under the BSD-2 licence. Please see [LICENSE](./LICENSE) for further details!

## Building NoSpherA2
NoSpherA2 relies heavily on the submodule [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS/tree/a64b75a2e00691b126a3c342a265f96fac98514f), [llvm](https://github.com/llvm/llvm-project/tree/2db262886f0c06c079e1b2808c4c14c16f8861b5) and the [mdSpan](https://github.com/kokkos/mdspan/tree/d34b447fbfdddfad63d2204923917e889ebe2e20) reference implementation. To clone this repository with all of it's dependencies do:

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

#### **Building OpenBLAS**
Navigate to the OpenBLAS subdirectory, open a **Developer PowerShell**, and execute:

```ps
cd OpenBLAS
mkdir build ; cd build
cmake -G "Visual Studio 17 2022" -DCMAKE_BUILD_TYPE=Release -DNOFORTRAN=ON ..
msbuild -nologo OpenBLAS.sln -p:Configuration=Release -m
cmake -DBUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=install -P cmake_install.cmake
```

#### **Building NoSpherA2**
Navigate to the Windows subdirectory and run:

```ps
cd ..\..\Windows
msbuild NoSpherA2.sln -nologo -p:Configuration=Release
```

✅ **The executable will be located at:** `NoSpherA2\Windows\x64\Release`.

---

### Linux Build Instructions

#### **Prerequisites**
Ensure `build-essential` is installed:

```sh
sudo apt update && sudo apt install -y build-essential
```

#### **Building NoSpherA2**
Inside the NoSpherA2 directory, simply run:

```sh
make
```

✅ **The executable will be located inside the NoSpherA2 directory.**

---

### macOS Build Instructions

#### **Prerequisites**
- **Xcode Command Line Tools**
- **CMake** (install via Homebrew if needed: `brew install cmake`)

#### **Building OpenMP**

```sh
cd llvm-project/openmp
mkdir build && cd build
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++ -DLIBOMP_ARCH=aarch64 -DLIBOMP_ENABLE_SHARED=OFF -DCMAKE_OSX_ARCHITECTURES=arm64 -DCMAKE_OSX_DEPLOYMENT_TARGET=14.0 ..
make
```

#### **Building OpenBLAS**

```sh
cd ../../../../OpenBLAS
mkdir installation
make -j
make PREFIX=installation install
```

#### **Building NoSpherA2**

```sh
cd Mac
make
mv NoSpherA2 ../NoSpherA2_Mac
```

✅ **The executable will be located at:** `NoSpherA2_Mac`.


