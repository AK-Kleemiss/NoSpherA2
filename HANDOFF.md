# Handoff: Debugging heap corruption crash in alanine_integrated_occ test

## What we're trying to fix

The `alanine_integrated_occ` test crashes with `STATUS_HEAP_CORRUPTION` (0xC0000374 → process exit 116) every run.
The test runs r2scan/def2-svp DFT on alanine using OCC's quantum chemistry library.
The SCF converges successfully (12 iterations, energy −323.278 Hartree), then crashes during cleanup.

The OCC log (`tests/alanine_integrated_occ/NoSpherA2_OCC.log`) always ends at:
```
[DEBUG] D_diff rows=119 cols=119
```

## Root cause hypothesis

`m_esp_evaluator` (a `std::unique_ptr<ESPEvaluator<double>>` member of `IntegralEngine`) gets corrupted to `0x0000001400000026` = `{38, 20}` as two int32 values:
- 38 = number of DFT methods loaded (from `occ/share/methods/dft_methods.json`)
- 20 = `PTR_ENV_START` = `#define PTR_ENV_START 20` in `Lib/LibCint/include/cint.h`

The pointer is never legitimately set in this code path — it should stay null. Something writes 8 bytes `{38, 20}` to the address where `m_esp_evaluator` lives. Windows heap manager detects the corrupted block header when DFT is torn down.

## What has already been tried

- Replaced `tbb::enumerable_thread_specific<ThreadLocalData>` with `std::vector<ThreadLocalData>` in `compute_K_dft` (`occ/include/occ/dft/dft.h`) — OCC rebuilt, crash persists
- Single-thread test (`threads=1`) also crashes — not a race condition
- Guard bytes on libcint buffers in `compute_three_center_integrals_tbb` — not triggered
- All 12 SCF iterations complete before crash
- The pure-DFT path uses `stored_coulomb_kernel_r` (serial, no ETS) — exchange never called

## Current diagnostic instrumentation

### `Lib/occ/include/occ/qm/scf_impl.h` (NoSpherA2-side header — no OCC rebuild needed)

After each SCF iteration, calls `_heapchk()` and logs result:
```
[HEAP] iter=N _heapchk=ok        — heap intact
[HEAP] iter=N _heapchk=M CORRUPTED  — heap corrupt after iteration N
```

After the SCF loop exits, before/after each `matrix.resize(0,0)` call:
```
[HEAP] before FD_comm.resize: _heapchk=...
[HEAP] before D_last.resize: _heapchk=...
[HEAP] before D_diff.resize: _heapchk=...
[HEAP] after D_diff.resize: _heapchk=...
```

The `_heapchk()` call itself may raise 0xC0000374 if heap is corrupt — whichever `[HEAP]` message is LAST in the log before the crash tells you when the corruption first becomes detectable.

**`_HEAPOK` = 1, `_HEAPBADNODE` = -3, `_HEAPBADBEGIN` = -2, `_HEAPEMPTY` = -1**

### `occ/include/occ/qm/integral_engine.h` (OCC source — requires OCC rebuild)

- Explicit `~IntegralEngine()`: if `m_esp_evaluator` is non-null (= corrupted), logs `[DETECT]` and calls `.release()` to prevent crash
- `= default` move constructor/assignment (needed because explicit dtor suppresses implicit move)
- `esp_is_corrupted()` public method

### `occ/src/qm/integral_engine_df.cpp` (OCC source — requires OCC rebuild)

At entry of `IntegralEngineDF::coulomb()`: checks `esp_is_corrupted()` on both engines, logs `[DETECT]` if true.

### `Lib/occ/include/occ/qm/integral_engine.h` (NoSpherA2-side copy)

Same explicit destructor + `= default` move + `esp_is_corrupted()` as the occ/ version.

## Current build state

NoSpherA2 is being rebuilt right now (background task `bgqya0es8`). OCC was already rebuilt and installed with the `~IntegralEngine` destructor + move fixes.

After the build finishes, run the test and check `NoSpherA2_OCC.log` for `[HEAP]` and `[DETECT]` messages.

## Build instructions

### Initialize build environment (PowerShell — required for both OCC and NoSpherA2 builds)
```powershell
$vcvars = & cmd /c "`"C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvarsall.bat`" amd64 && set" 2>&1
$vcvars | Where-Object { $_ -match '=' } | ForEach-Object {
    $name, $value = $_ -split '=', 2
    [System.Environment]::SetEnvironmentVariable($name, $value, 'Process')
}
$vscmake = "C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
```

### Rebuild OCC (when `occ/` source or `occ/include/` headers change)
```powershell
Set-Location "D:\git\NoSpherA2"
& $vscmake --build build-windows-clang-cl --config Release
& $vscmake --install build-windows-clang-cl
```

### Rebuild NoSpherA2 (when `Lib/occ/include/` or `Src/` changes — no OCC rebuild needed)
```powershell
Set-Location "D:\git\NoSpherA2"
make.exe
```

### Run the failing test
```powershell
Set-Location "D:\git\NoSpherA2\tests\alanine_integrated_occ"
Remove-Item NoSpherA2_OCC.log -ErrorAction SilentlyContinue
& "D:\git\NoSpherA2\NoSpherA2.exe" -occ alanine.toml -cif alanine.cif -dmin 0.5 -acc 1
# Then read NoSpherA2_OCC.log for [HEAP] and [DETECT] messages
```

## Key files

| File | Role | Rebuild needed |
|------|------|---------------|
| `occ/include/occ/qm/integral_engine.h` | IntegralEngine with explicit dtor | OCC + NoSpherA2 |
| `occ/include/occ/dft/dft.h` | compute_K_dft with ETS→vector fix | OCC |
| `occ/src/qm/integral_engine_df.cpp` | coulomb() with [DETECT] checks | OCC |
| `occ/src/qm/detail/df_kernels.h` | stored_coulomb_kernel_r (serial J), compute_three_center_integrals_tbb | OCC |
| `Lib/occ/include/occ/qm/scf_impl.h` | compute_scf_energy with [HEAP] _heapchk calls | NoSpherA2 only |
| `Lib/occ/include/occ/qm/integral_engine.h` | NoSpherA2-side copy of integral_engine.h | NoSpherA2 only |
| `tests/alanine_integrated_occ/NoSpherA2_OCC.log` | OCC log — check for [HEAP]/[DETECT] | — |

## IntegralEngine memory layout (critical)
```
IntegralEngine members (in order):
  double m_precision{1e-12}
  AOBasis m_aobasis, m_auxbasis
  ShellPairList m_shellpairs
  mutable IntEnv m_env          // ~80 bytes; last field is vector<double> m_env_data
  mutable unique_ptr<ESPEvaluator<double>> m_esp_evaluator  // 8 bytes — GETS CORRUPTED to {38,20}
  mutable bool m_esp_initialized{false}
  bool m_have_ecp{false}
  int m_ecp_ao_max_l{0}, m_ecp_max_l{0}
```

`IntegralEngineDF` contains:
```
  double m_precision
  mutable IntegralEngine m_ao_engine   // ← m_esp_evaluator here gets corrupted
  mutable IntegralEngine m_aux_engine
  Eigen::LLT<Mat> V_LLt
  Mat m_integral_store                 // stored 3-center integrals (119×582)
  Policy m_policy
  size_t m_integral_store_memory_limit
  CoulombMethod m_coulomb_method
  mutable unique_ptr<SplitRIJ> m_split_rij
```

## What to do next

1. Wait for NoSpherA2 build to finish (or rebuild: `make.exe` with VS env initialized)
2. Run the test, read `NoSpherA2_OCC.log`
3. Look for `[HEAP]` messages:
   - If `[HEAP] iter=N _heapchk=M CORRUPTED` appears → corruption first detectable after iteration N's fock build
   - If no `[HEAP] iter=N` appears at all → `_heapchk()` itself crashed, meaning heap is corrupt very early
   - If all iter checks show ok, but one of the resize checks shows CORRUPTED → corruption happens during stack cleanup
4. Look for `[DETECT]` from `~IntegralEngine` or `coulomb()` entry
5. Once we know WHICH iteration (or which cleanup step) first shows corruption, add `_heapchk()` inside that iteration's `compute_fock` → `stored_coulomb_kernel_r` to narrow further

## Suspected area

The `m_integral_store` matrix (119×582 doubles = ~540 KB) is allocated in `IntegralEngineDF`. It is built in `compute_stored_integrals()`. The `stored_coulomb_kernel_r` uses it read-only. However `compute_three_center_integrals_tbb` fills it during initial setup. 

`m_esp_evaluator` lives inside `m_ao_engine` which is the first `IntegralEngine` inside `IntegralEngineDF`. If something writes past the end of `m_env.m_env_data` (the `vector<double>` that stores libcint's env array), it could overflow into `m_esp_evaluator`. The libcint env array contains atom/basis data including counts like 38 (number of DFT methods or shells) and offsets like 20 (PTR_ENV_START).

**Specifically suspicious**: `IntegralEnvironment` assignment (`m_env = IntEnv(...)`) in `set_auxiliary_basis()` called from `IntegralEngineDF` constructor. If the old env vector is not properly freed and the new one overlaps, or if the assignment operator has a bug, it could corrupt adjacent memory.
