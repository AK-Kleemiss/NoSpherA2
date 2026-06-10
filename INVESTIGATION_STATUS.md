# alanine_integrated_occ Crash Investigation Status
**Last updated: 2026-06-07**

## Goal
Fix the `alanine_integrated_occ` test so it passes in both:
- In-process DLL mode (Windows VS test runner)
- Subprocess EXE mode

## IMPORTANT: Standalone EXE also fails for ALL OCC jobs
The standalone `NoSpherA2.exe` (release build) crashes for **any** job that goes through the OCC
single-point driver with density fitting enabled. This is not limited to the test runner — running
`NoSpherA2.exe` directly on the command line with `alanine.toml` also crashes at
`~IntegralEngineDF`. The subprocess approach in the test runner therefore cannot work as a
workaround until the underlying OCC heap corruption is fixed in the exe itself.

---

## Crash Description

### EXE (Release build) — actual crash site
```
handmade_aligned_free(ptr)          ← Line 118, Eigen memory.h
aligned_free / conditional_aligned_free
DenseStorage::~DenseStorage()
Matrix::~Matrix()                   ← destroying m_integral_store
IntegralEngineDF::~IntegralEngineDF()
run_method<DFT,0>
single_point_driver
single_point
run_scf_external
main Line 419
```
Exit code: **0xC0000374** (STATUS_HEAP_CORRUPTION)

`handmade_aligned_free(ptr)` reads `original = *(double**)((ptr) - 1)` (the original malloc pointer stored 8 bytes before the aligned data). This value is **corrupted** → `free(corrupted_ptr)` → `RtlFailFast` → fast-fail crash.

### DLL (Debug build) — heap corruption detected earlier
```
occ::gto::Shell::Shell(...)
occ::gto::AOBasis::load(...)
occ::driver::load_basis_set(...)
occ::driver::single_point_driver(...)   ← LINE 351 (first call in this function)
occ::driver::single_point(...)
occ::main::run_scf_external(...)
main() Line 419
```

The debug heap validates on every `HeapAlloc`. The corruption was introduced somewhere **before** `single_point_driver line 351`, detected at the next alloc (`Shell::Shell`).

---

## Root Cause

**Windows heap is corrupted by OCC code.** The corruption happens during or before the DFT density-fitting computation (`IntegralEngineDF`).

The crash **cannot be caught** with SEH (`__try/__except`) because:
- Windows heap corruption on modern Windows uses `RtlFailFast` / `__fastfail` (hardware trap `int 0x29`)
- This bypasses ALL user-mode SEH handlers in the call stack
- `nos_seh_dispatch` in dllmain.cpp with `__except(GetExceptionCode() == 0xC0000374)` does NOT work

---

## What Was Tried and Why It Failed

### 1. SEH Handler in `nos_seh_dispatch`
```cpp
__try { result = nos_cpp_dispatch(...); }
__except (GetExceptionCode() == 0xC0000374 ? EXCEPTION_EXECUTE_HANDLER : EXCEPTION_CONTINUE_SEARCH)
{ ... }
```
**Failed**: `__fastfail` bypasses all user-mode SEH. Still crashes.

### 2. TBB Malloc Proxy
- Added `tbbmalloc_proxy.lib`, `tbbmalloc.lib` to linker
- Added `/INCLUDE:__TBB_malloc_proxy` flag
- Added PostBuildEvent to copy `tbbmalloc_proxy.dll` and `tbbmalloc.dll` to output dirs
- Tried `GetProcAddress(hProxy, "TBB_malloc_replacement_init")` in DllMain

**Failed** for multiple reasons:
1. `tbbmalloc_proxy.dll` only exports 2 symbols: `TBB_malloc_replacement_log` and `__TBB_malloc_proxy`. Does NOT export `TBB_malloc_replacement_init` → `GetProcAddress` returns NULL.
2. Loading order: `tbbmalloc_proxy.dll`'s DllMain runs while `NoSpherA2_DLL.dll` is still loading → our DLL's IAT never patched.
3. **Fundamental**: even if TBB intercepted our malloc/free, OCC's overflow corrupts the **Windows process heap** used by ucrtbase/ntdll internals → ANY subsequent `HeapAlloc`/`HeapFree` anywhere in the process fast-fails.

### 3. Subprocess Mode
`TestDef.subprocess = true` exists in TestRunner.h but doesn't help because the EXE also crashes (in release mode, crash is at `~IntegralEngineDF` — AFTER the computation but BEFORE the output file is written). So the subprocess exits with 0xC0000374 and the output file is never produced.

---

## Current State of Source Files

### `NoSpherA2_DLL/dllmain.cpp`
Contains (currently):
- Dead-code TBB proxy init attempt in DllMain (does nothing since export doesn't exist)
- `nos_cpp_dispatch`: C++ exception layer
- `nos_seh_dispatch`: SEH layer (ineffective against fast-fail)
- `NoSpherA2_run`: exported entry point
  - Saves/restores CWD
  - Syncs `OCC_DATA_PATH` via `GetEnvironmentVariableA` + `_putenv_s` + `occ::set_data_directory`
  - Saves/restores `std::cout`/`std::cerr` buffers

### `NoSpherA2_DLL/NoSpherA2_DLL.vcxproj`
Has:
- `tbbmalloc_proxy.lib;tbbmalloc.lib` in AdditionalDependencies (all x64 configs)
- `/INCLUDE:__TBB_malloc_proxy` in all x64 Link sections
- PostBuildEvent to copy TBB DLLs (added last session, may still be present)

### `Src/pch.h`
Has `#include <tbb/tbbmalloc_proxy.h>` (from prior sessions).

---

## Where the Overflow Is — Partial Investigation

### Test configuration (`tests/alanine_integrated_occ/alanine.toml`)
```toml
threads = 8
[scf]
input = "alanine.xyz"
method = "r2scan"
basis = "def2-svp"
multiplicity = 1
spherical = true
output = "fchk"
df-basis = "def2-universal-jkfit"
dft_grid_min_angular = 86
dft_grid_max_angular = 110
dft_grid_radial_precision = 1e-4
```

**Molecule**: alanine (C3H7NO2), ~13 atoms, elements H/C/N/O only.

**Basis sizes** (approx):
- `def2-svp` AO: nbf ≈ 65, max l=2 (d-function), max spherical shell size = 5
- `def2-universal-jkfit` aux: ndf ≈ 460, max l=4 (g-function for H/C/N/O), max spherical shell size = 9
  - Elements 1,6,7,8 all confirmed max l=4 in the JSON file
  - NO general contractions in either basis set (confirmed via JSON inspection)

**m_integral_store size**: `65 × 65 × 460 × 8 bytes ≈ 15.5 MB`

### Buffer size in `compute_three_center_integrals_tbb` (occ/src/qm/detail/df_kernels.h)
```cpp
size_t bufsize = aobasis.max_shell_size() * aobasis.max_shell_size() * auxbasis.max_shell_size();
// = 5 * 5 * 9 = 225 doubles
auto buffer = std::make_unique<double[]>(bufsize);
```

### libcint buffer analysis (`occ/src/libcint/src/cint3c2e.c`, `cart2sph.c`)
- `CINT3c2e_drv`: allocates `gctr` from `cache` (NOT from `buffer`/`out`)
- `c2s_sph_3c2e1`: allocates scratch (`buf1/2/3`) from `cache`, writes ONLY spherical output to `bufijk`
- For single-contracted shells (x_ctr=1) with no general contractions, output size = `dims[0]*dims[1]*dims[2]` = `(2*l_i+1)*(2*l_j+1)*(2*l_k+1)` ≤ `5*5*9 = 225`

**CONCLUSION**: The buffer in `compute_three_center_integrals_tbb` appears correctly sized for the no-general-contractions case. The overflow source was NOT conclusively identified.

### IMPORTANT NOTE: Two libcint builds
There is a `libcint/` directory at the repo root (used by NoSpherA2 directly), AND OCC has its own internal libcint (pulled as part of the OCC submodule). When OCC calls libcint, it uses ITS OWN libcint, NOT `D:\git\NoSpherA2\libcint\`. The crash may be in OCC's bundled libcint, not the one we inspected.

**Next step**: Find OCC's internal libcint location (somewhere inside `occ/` submodule).

---

## Call Chain for the Computation

```
main() line 419
  → run_scf_external(config)
      → print_header()
      → read_input_file("alanine.xyz", config)   [XyzFileReader]
      → single_point(config)
          → single_point_driver(config)
              → load_basis_set(m, "def2-svp", spherical=true)    ← Line 351
              → print_configuration(m, config)
              → set_density_fitting_basis(config.basis.df_name, ...)
                  → AOBasis::load(atoms, "def2-universal-jkfit")  [dfbasis]
                  → IntegralEngineDF(atoms, ao_shells, df_shells)
                      → m_aux_engine.one_electron_operator(coulomb)  [V = (P|Q)]
                      → V_LLt = LLT(V)
              → run_method<DFT, Restricted>(m, basis, config)
                  → DFT proc = ...
                  → scf.compute()                      [SCF iterations]
                      → proc.coulomb_and_exchange(mo)
                          → compute_stored_integrals()  [fills m_integral_store]
                          → stored_coulomb/exchange kernels
                  → Wavefunction wfn = scf.wavefunction()
              // ← run_method returns here, then destructor runs:
              → ~DFT → ~HartreeFock → ~unique_ptr<IntegralEngineDF>
                  → ~IntegralEngineDF
                      → ~m_integral_store       ← CRASH HERE (release)
```

---

## Key Technical Facts

1. **`handmade_aligned_free` crash**: `Eigen::aligned_free` on clang-cl uses `handmade_aligned_malloc/free` because `__clang__` defines `__GNUC__`, preventing `EIGEN_HAS_WIN32_MALLOC`. Stores original malloc ptr at `*(aligned_ptr - 8)`.

2. **Debug vs Release detection**:
   - Debug: detects any heap corruption at the NEXT `HeapAlloc`
   - Release: detects at the specific `HeapFree` of the corrupted block

3. **The multi-test scenario**: If tests share a process (DLL runner), a previous test's heap corruption could be detected in `alanine_integrated_occ` at its first allocation (`Shell::Shell` at line 351). Check whether other tests run BEFORE `alanine_integrated_occ` in the test class.

4. **OCC's `set_data_directory`**: Just sets a static string. Does NOT modify env vars. The `NoSpherA2_run` function now calls both `occ::set_data_directory(occBuf)` AND `_putenv_s("OCC_DATA_PATH", occBuf)` to keep everything in sync.

5. **`ensure_occ_data_path` stdout issue**: Fixed — the message goes to `NoSpherA2.log` (after cout redirect), not to test output.

---

## Possible Next Steps

### Option A: Find the actual OCC overflow (correct fix)
1. Find OCC's internal libcint: `ls occ/` — look for `libcint` or `extern/libcint` subdirectory
2. Check if OCC's libcint has different buffer sizing than standalone libcint
3. Use Windows Page Heap (`gflags /p /enable NoSpherA2.exe /full`) to find exact overflow address
4. Look at the DFT grid/XC code for overflows (not just integral code)
5. Check if `one_electron_operator` (2-center coulomb V matrix) has buffer issues

### Option B: Avoid freeing `m_integral_store` (workaround)
Modify `IntegralEngineDF` destructor to intentionally leak `m_integral_store` (call `new(&m_integral_store) Mat()` via placement-new, which won't work, OR use `m_integral_store.data()` pointer to mangle the stored metadata to a valid-enough value).  
**Problem**: This is fragile and undefined behavior.

### Option C: Run the test as subprocess + write output early
Modify `single_point_driver` or `run_method` to write the wavefunction/fchk output BEFORE `proc` (DFT object) goes out of scope and crashes.  
**Problem**: Output path and filename need to be accessible inside `run_method`.

### Option D: Use VirtualAlloc for `m_integral_store`
Override Eigen's allocator for this specific matrix to use `VirtualAlloc` with a guard page after the allocation. An overflow would then hit the guard page (access violation = catchable with SEH) instead of corrupting adjacent heap blocks.  
**Problem**: Significant Eigen internals changes needed.

### Option E: Mark the test as "expected failure" or skip it in the DLL runner
Use `subprocess = true` in the test definition and accept that the test may fail due to the crash. Once the exe crash is fixed, the subprocess mode would work.

---

## File Locations

- Test config: `tests/alanine_integrated_occ/alanine.toml`
- Test geometry: `tests/alanine_integrated_occ/alanine.xyz`  
- DLL entry: `NoSpherA2_DLL/dllmain.cpp`
- DLL project: `NoSpherA2_DLL/NoSpherA2_DLL.vcxproj`
- Test runner: `Windows/Tests/TestRunner.h`
- OCC DF kernels: `occ/src/qm/detail/df_kernels.h`
- OCC DF engine: `occ/src/qm/integral_engine_df.cpp`
- OCC cint interface: `occ/include/occ/qm/cint_interface.h`
- OCC single point: `occ/src/driver/single_point.cpp`
- Standalone libcint: `libcint/src/cint3c2e.c`, `libcint/src/cart2sph.c`
- OCC's own libcint: **NOT YET LOCATED** — search inside `occ/` submodule

## TBB DLL Locations
- Source: `D:\git\NoSpherA2\build-tbb-malloc\clang_20.1_cxx_64_md_release\tbbmalloc_proxy.dll`
- Source: `D:\git\NoSpherA2\build-tbb-malloc\clang_20.1_cxx_64_md_release\tbbmalloc.dll`
- Copied to: `Windows\x64\Debug\`, `Windows\x64\Release\`, `Windows\x64\Release + Copy\`
