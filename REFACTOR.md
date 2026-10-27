# Refactoring the XCW workflow
## Where we want to go
### Intended final Structure
In the end, the current XCW workflow should be running using three classes.
1. Structure Factors:

The structure_factor class already exists. This class (together with auxiliaries like the cell class) handles all computations regarding the calculation of structure factors. This includes, but is not limited to, Debye-Waller factors, the structure factors themselves, phase factors, anomalous dispersion correction, all kinds of crystallographic refinement indicators, scale factor etc.
Up until recently, this class was used to compute what is called the *I tensor*, a giant tensor holding integral values. Recent changes made storing this massive object redundant and it can now be calculated after every SCF iteration the SCF solver does (basically direct SCF approach). The structure_factor class is already close to completely reworked and won't need much more changes. The class should still be used to calculate the I tensor and there should be different functions for the calculation of F calc, using a density matrix (currently the way the XCW workflow does it), using a tsc_calculator and using tsc values.

2. SCF solver:

The SCF_wrapper should be renamed to SCF_solver. It is basically a class that does SCF calculations using the *occ submodule*. There should be a public function called run that does one complete solution of a system. I imagine the SCF_solver class to be an abstract class with different child classes like an XCW_solver that can handle the current XCW workflow. The basic SCF_solver class should be able to just solve a system on a defined theory level. This SCF_solver class is the heart of the procedure and should be easily extendable to other functions like energy decompositions etc.

3. XCW_wrapper:

The XCW_solver should be renamed to XCW_wrapper. This is the class that gets build for the XCW workflow from main(). It has minimal attributes, one structure_factors object, one SCF_solver object and the opt object with all the settings necessary for the procedure.
This class is the one that calls SCF_solver's run() and mediates between the two classes. XCW_wrapper then basically only has to do one thing: Tell the SCF_solver how it is supposed to run, it should to the inter-lambda extrapolations and gaussian halting criteria etc. This class decides for which lambda the SCF solver runs etc.

### Decisions
- **The SCF never sees the structure factors directly (decided 2026-10-26).** The SCF talks to the crystallographic side only through the abstract `SCF_coupling` interface (`XCW_solver.h`: `evaluate`, `perturbation`, `perturbation_response`, `quality`, `save_state`/`restore_state`, `on_device`). XCW_wrapper implements it on top of its `structure_factors` object and hands itself to the XCW child of SCF_solver. This keeps the SCF_solver base free of crystallography, so other children (energy decompositions etc.) need no structure factors at all.
- **How SCF_solver is extended (run 1).** `SCF_solver` is not pure abstract: on its own it solves the plain Hartree-Fock system, as "the basic SCF_solver" above should. A child adds a *model term* to the energy through protected virtual hooks, and the DIIS, rescue, TRAH and convergence machinery then works on E + that term without knowing what it is:
  - `evaluate_model(scf, timed)`: called at the current density before every Fock build
  - `model_energy()`: the term added to E (XCW: lambda chi^2)
  - `add_model_fock(F)`: its derivative, added to the Fock matrix (XCW: lambda times the perturbation matrix)
  - `add_model_response(scf, dD, dF)` and `model_active()`: the term's part of the TRAH Hessian-vector product
  - `save_model_state()` / `restore_model_state()`: around the finite-difference Hessian check
  - `begin_log()`, `log_iteration(scf)`, `end_log(scf)`, `functional_name()`: what the SCF log prints

  `XCW_solver` (the child) overrides them with the coupling and holds `lambda`; `XCW_wrapper::run_lambda` sets `scf_solver.lambda` and calls `scf_solver.run(hf, wfn, use_guess)`.
- **Code style (decided 2026-10-26).** CLAUDE.md's house style applies (terse, no blank lines inside functions, project typedefs, `#pragma omp` at column 0, tabs), together with what is asked for below:
  - Comment-sparse, but enough comments to see what is going on and which part of the code does what: mostly one-liners, sometimes two or three lines. Longer explanations only for theory-heavy parts (SOSCF, TRAH), in simple language.
  - Loop indices short, but descriptive where it helps: a nested loop over symmetry operations runs over `sym_op`, not `s`. Simple loops keep `i`, `j`, `k`, `l`.
  - Comments that cite a specific structure or run (Fe(phen)2(SCN)2, TEMPO, ...) are replaced by the general situation they describe.
  - GPU blocks are marked as such with a one-line note on what they do (each has a CPU equivalent).

## How we will get there
I will give you a time limit and you will refactor the code according to the criteria found in this file. You will work on a local clone of a GitHub repo, on a branch called XCW_refactor. After running until the deadline, you are tasked with preapring a pull request onto another branch called density_source. I will then look at the pull request, make changes if necessary and commit it. To keep track of what needs to be done and keep track of things we discuss outside of this file, you are expected to modify this file so you can read it and see what is next on the agenda. I do not expect the refactor to just happen in one go, so it is completely fine to just get things done partially during one run and keep some things for consecutive runs.

When the deadline comes, merge the then current density_source into XCW_refactor first. Changes made to density_source in the meantime are improvements: convert them into the new structure rather than dropping them because they do not fit it.

### Runs
| Run | Deadline | Goal |
|-----|----------|------|
| 1 | 2026-10-27 09:00 | Baseline; SCF_wrapper -> SCF_solver (abstract) + XCW child; XCW_solver -> XCW_wrapper; comment pass; PR onto density_source |

### Agenda
Done in run 1:
- [x] Baseline timings and memory of the XCW gtests and P1 goldens on the unchanged code
- [x] SCF_wrapper -> SCF_solver (`SCF_solver.h/.cpp`): base class with a public `run(hf, wfn, use_guess)` that solves one system; the XCW part lives in the child `XCW_solver` (`XCW_solver.h/.cpp`, about 100 lines) through the model hooks above. `SCF_coupling` and `fit_quality` moved to `XCW_solver.h`
- [x] XCW_solver -> XCW_wrapper (`XCW_wrapper.h/.cpp`): holds `opt`, `sf`, the SCF solver and the SCF log (plus the Gaussian halting history); drives the lambda scan and implements `SCF_coupling` by forwarding to `structure_factors`. `run()` is split into named steps: `setup_structure_factors`, `setup_scf_integrals`, `read_start_guess`, `run_lambda`, `approach_lambda` (the step halving after a failed lambda) and `extend_scan` (Gaussian halting)
- [x] structure_factors: F_calc entry points `calc_F_calc(const dMatrix2& D)` (density matrix, through the I tensor or on the grid as `XG.mode` says), `calc_F_calc(const itsc_block&)` (tsc values) and `calc_F_calc(std::vector<WFN>&, GridManager*, std::ostream&)` with `tsc_table` (the tsc calculator), next to the existing `calc_F_calc(const cvec2&)` (atomic scattering factors) and `calc_F_calc(I_tensor&, const dMatrix2&)` (the tensor walk). The perturbation, its response, `contract_I`/`contract_tensor`, the tscb table (`result_table`, with the grid route's layout state and the tsc grids) and the F_calc save/restore moved here from XCW_wrapper
- [x] Comment and naming pass over the three classes: case-specific comments generalised, GPU blocks marked, symmetry loops over `sym_op`, scatterer loops over `scatterer`
- [x] `options::SCF_settings` thresholds default to the `normal` preset (they were uninitialised, harmless only while the settings parser always filled them)
- [x] `ScfSolverTests` (new, 3): the plain SCF_solver against OCC's own SCF (water RHF, methyl radical UHF, def2-SVP, 1e-7 Eh) and a warm restart that goes straight to TRAH and converges in at most 3 iterations
- [x] Tests unchanged, timings and memory compared against the baseline (below)
- [x] Merged density_source at the deadline (`0efc6add`, seven commits: `-tshb` for XCW steps, `-xcw_partition ri`, fit caches). Their XCW parts went into the new structure: the RI partition's setup (`structure_factors::setup_ri(ao)`, AutoAux or `-ri_fit`) and the `NOSPHERA2_XCW_ADJOINT_CHECK` check (`structure_factors::check_adjoint`) into structure_factors, the RI branch of the tscb layout into `result_table`, the `XR.response` switch into `perturbation_response`, and the step's `.tshb` fit (`write_tshb`, its aux wavefunction and `PartitionCache`) plus the RI restraint targets into `XCW_wrapper::write_step_results`
- [x] Merged density_source again (`33fa9685`, three commits: density-following TFVC/MBIS/EMBIS shares on XCW's grid, RI restraints for every scheme, a SALTED fix, the typing net). The share machinery (`xcw_weights.cpp`, `XG.w0`/`XG.ws`) is structure_factors' and merged as is; the driver side went into XCW_wrapper: the refit between the solves of a step (`run_lambda`, `refit_moves`), the per-evaluation refit of `NOSPHERA2_XCW_WEIGHT_UPDATE` (`evaluate`), the lambda 0 fit and share report (`report_lambda`), the grid charges (`write_step_results`) and `NOSPHERA2_XCW_WEIGHT_CHECK` (`weight_check`). `result_table` takes whether to refit the shares to the table's density, and `check_adjoint` runs after the first step. The RI targets set from each step's density are gone upstream, since they now come from the shares

Next:
- [ ] SCF_solver for other theory levels: it is tied to `occ::qm::SCF<occ::qm::HartreeFock>` (`scf_t`). DFT (`occ::dft::DFT`) needs the procedure as a template parameter or behind an interface; the stored-ERI Fock build and the TRAH Hessian product are Hartree-Fock specific as they stand
- [ ] `SCF_log_writer` still has the XCW columns (criterion, GooF, R1); the base SCF writes its own simpler lines. A log per child would be cleaner
- [ ] `structure_factors::build_I` is about 1150 lines: split it into the screening, the grid setup and the integration, and decide whether the I tensor route stays now that the grid route is the default
- [ ] `write_step_results` in XCW_wrapper still writes seven kinds of output (tscb or the `.tshb` RI fit, wfn, Fcalc.txt, fchk, NBO .47, RGBI); an output writer of its own, holding the `.tshb` fit's cache, would leave XCW_wrapper with the scan alone
- [ ] The Gaussian halting functions still have blank-line-grouped bodies from before; reformat them when they are next touched

## What to look out for
When coming accross a comment, that references a specific case-fault or something of this nature, replace it by a more general description. I want short and concise comments instead of fully fledged explanations for the most part, for theoretically demanding parts like SOSCF or TRAH it is fine to be more elaborate but try sticking to simple language. The code you're rewriting should for the most part be self-commenting, no really abstract abbreviations or variable names, this is seldomly necessary. For for loops that just iterate trivially it is fine to use indices like i and j but if there are a lot of nested for loops it would be better to also describe what is being iterated through the assigned variable. I want the code to be as readable as possible so somebody who is new to the subject could start working on it after a short introduction. There are also code blocks pertaining GPU implementation, these should be marked as such and a short comment should be written explaining what is happening in the code block (there will always be a CPU equivalent to these code blocks so that should be easily done). If you find parts of the code where there is obvious unnecessary overhead like uncollapsed for loops (that have no reason to be uncollapsed) or redundant copying of variables, you are free to remove them and improve the code's efficiency.

## Tracking
The codebase has various tests implemented in it. Changes you make to the code should never break any tests that worked beforehand or slow down the code significantly. Additionally, keep track of memory usage and timings and compare it against each other and the current state of the codebase. These tracked quantities should be visualized in this file using, at least, tables.

### How it is measured
Machine: AMD Ryzen 7 PRO 8700G (8 cores / 16 threads), 31 GB RAM, Windows 11, MSVC 14.51, `release-windows-tests` preset, CPU only (no HIP compiler, so no GPU path is compiled). The benchmark runs `NoSpherA2.exe -do_XCW -no_gpu -cif P1_test.cif -hkl <hkl> -XCW_settings <settings> -no_date` in a scratch copy of `tests/P1_test`, 5 repeats per case, and reads wall time, CPU time and peak working set / peak commit from the process handle after it exits (PowerShell, `GetProcessMemoryInfo`). The cases:

| Case | Reflections | Settings |
|------|------------:|----------|
| P1_full_F | 3216 (`P1_test.hkl`) | `test_settings.txt`: sto-3g, F, lambda 0 to 0.05 in 0.01 |
| P1_full_F2 | 3216 | `test_settings_f2.txt`: the same against F^2 |
| P1_full_h2 | 3216 | `h2_test_settings.txt`: 1/\|H\|^2-weighted |
| P1_subset_F | 402 (`P1_subset.hkl`) | `test_settings.txt` |
| P1_full_svp | 3216 | `test_settings.txt` with `basis_set def2-svp` (a larger SCF) |

### Baseline (`a1448cf6`, 2026-10-26)
`ctest --preset release-windows-tests -j 6`: **1343 run, 0 failed, 23 skipped** (2 disabled), 120.8 s. The 94 XCW tests (`Xcw*`, `XCW*`, `P1_*`, `Extinction*`, `HessianVector*`) add up to 149.4 s of test time. With `RUN_FULL_TEST=1` the four `P1_*_full` goldens pass too (8.1, 4.0, 8.1, 3.9 s).

| Case | Wall (s) | CPU (s) | Peak working set (MB) | Peak commit (MB) |
|------|---------:|--------:|----------------------:|-----------------:|
| P1_full_F | 4.43 ± 0.04 | 59.6 ± 0.6 | 267.5 | 281.2 |
| P1_full_F2 | 4.43 ± 0.05 | 59.8 ± 0.7 | 267.4 | 281.2 |
| P1_full_h2 | 4.16 ± 0.17 | 53.8 ± 1.7 | 267.1 | 281.1 |
| P1_subset_F | 4.51 ± 0.09 | 63.7 ± 1.3 | 263.0 | 276.6 |
| P1_full_svp | 21.21 ± 0.32 | 306.2 ± 1.2 | 4008.7 | 4059.8 |

(mean ± standard deviation of 5 runs)

### Run 1 (2026-10-26)
`ctest --preset release-windows-tests -j 6`: **1346 run, 0 failed, 23 skipped** (2 disabled), 110.9 s (baseline: 1343 run, 120.8 s; the 3 new tests are `ScfSolverTests`). With `RUN_FULL_TEST=1` the four `P1_*_full` goldens pass too.

The same cases on the same machine; `baseline_b` is the baseline binary measured again in the same sitting as `run 1`:

| Case | Build | Wall (s) | CPU (s) | Peak working set (MB) | Peak commit (MB) |
|------|-------|---------:|--------:|----------------------:|-----------------:|
| P1_full_F | baseline | 4.43 ± 0.04 | 59.6 ± 0.6 | 267.5 | 281.2 |
| P1_full_F | baseline_b | 4.43 ± 0.09 | 62.5 ± 0.9 | 267.2 | 281.1 |
| P1_full_F | run 1 | 4.33 ± 0.05 | 61.4 ± 0.7 | 267.7 | 281.2 |
| P1_full_F2 | baseline | 4.43 ± 0.05 | 59.8 ± 0.7 | 267.4 | 281.2 |
| P1_full_F2 | baseline_b | 4.34 ± 0.07 | 62.2 ± 0.9 | 267.7 | 281.5 |
| P1_full_F2 | run 1 | 4.35 ± 0.02 | 62.1 ± 0.4 | 267.4 | 281.2 |
| P1_full_h2 | baseline | 4.16 ± 0.17 | 53.8 ± 1.7 | 267.1 | 281.1 |
| P1_full_h2 | baseline_b | 3.94 ± 0.05 | 56.2 ± 0.7 | 267.6 | 281.4 |
| P1_full_h2 | run 1 | 3.93 ± 0.10 | 56.0 ± 1.4 | 267.5 | 281.3 |
| P1_subset_F | baseline | 4.51 ± 0.09 | 63.7 ± 1.3 | 263.0 | 276.6 |
| P1_subset_F | baseline_b | 4.35 ± 0.04 | 62.1 ± 0.8 | 263.2 | 277.1 |
| P1_subset_F | run 1 | 4.41 ± 0.07 | 63.1 ± 1.1 | 262.9 | 276.6 |
| P1_full_svp | baseline | 21.21 ± 0.32 | 306.2 ± 1.2 | 4008.7 | 4059.8 |
| P1_full_svp | baseline_b | 20.89 ± 0.08 | 304.3 ± 1.5 | 4008.5 | 4059.8 |
| P1_full_svp | run 1 | 20.98 ± 0.05 | 305.1 ± 1.1 | 4005.9 | 4059.5 |

Every difference is within the run-to-run scatter: the refactor costs no time and no memory. The outputs (`NoSpherA2.log`, `SCF.log`, `NA2_*_Fcalc.txt`, `.tscb`) differ from the baseline's by no more than two baseline runs differ from each other, since the threaded sums are not bit-reproducible: lambda-table energies by up to about 2e-11 relative, F_calc in the last printed digit. Every case needs the same number of Fock builds as before. The step-halving retry and `load_wfn`, which no test covers, were checked by hand (a `max_iter 5` scan with 0.04 steps from a read-in lambda 0): the same retry sequence and lambda table as the baseline.

### Merge with density_source (2026-10-27)
`0efc6add` merged (seven commits, `-tshb` step fits and `-xcw_partition ri`, ported as listed in the agenda). `ctest --preset release-windows-tests -j 6`: **1330 passed, 0 failed, 23 skipped** (2 disabled) of 1355, 120.0 s; the four `P1_*_full` goldens pass. The benchmark cases use neither new option, so they measure that the merge leaves the default path alone; `baseline_c` is the baseline binary again in the same sitting:

| Case | Build | Wall (s) | CPU (s) | Peak working set (MB) | Peak commit (MB) |
|------|-------|---------:|--------:|----------------------:|-----------------:|
| P1_full_F | baseline_c | 4.61 ± 0.19 | 61.1 ± 1.0 | 267.2 | 281.2 |
| P1_full_F | merge | 4.37 ± 0.11 | 62.1 ± 0.8 | 267.7 | 281.5 |
| P1_full_F2 | baseline_c | 4.51 ± 0.15 | 61.0 ± 0.5 | 267.4 | 281.2 |
| P1_full_F2 | merge | 4.29 ± 0.04 | 61.3 ± 0.5 | 267.8 | 281.8 |
| P1_full_h2 | baseline_c | 4.00 ± 0.05 | 55.4 ± 0.2 | 267.1 | 281.0 |
| P1_full_h2 | merge | 3.92 ± 0.04 | 56.0 ± 0.5 | 267.6 | 281.3 |
| P1_subset_F | baseline_c | 4.43 ± 0.06 | 61.9 ± 0.9 | 263.0 | 277.0 |
| P1_subset_F | merge | 4.35 ± 0.07 | 62.1 ± 1.1 | 263.3 | 276.8 |
| P1_full_svp | baseline_c | 21.12 ± 0.06 | 302.7 ± 1.2 | 4009.5 | 4060.2 |
| P1_full_svp | merge | 21.06 ± 0.06 | 305.6 ± 1.0 | 4007.0 | 4060.5 |

Again within the scatter, with the same Fock build counts and outputs within the baseline's own noise.

### Second merge with density_source (2026-10-27)
`33fa9685` merged (three commits: density-following TFVC/MBIS/EMBIS shares and RI restraints for every scheme, ported as listed in the agenda, a SALTED fix and the typing net). `ctest --preset release-windows-tests -j 6`: **1333 passed, 0 failed, 23 skipped** (2 disabled) of 1358, 138.4 s; the three new `XcwPartitionTests.*SharesFollowTheDensity` pass unchanged, and so do the four `P1_*_full` goldens. The benchmark cases run on Hirshfeld weights, so the new share code is idle in them; `merge` is the binary of the first merge, the two binaries alternating run by run (4 runs each, the machine busier than in the sittings above):

| Case | Build | Wall (s) | CPU (s) | Peak working set (MB) | Peak commit (MB) |
|------|-------|---------:|--------:|----------------------:|-----------------:|
| P1_full_F | merge | 4.54 ± 0.05 | 62.1 ± 1.9 | 267.6 | 281.5 |
| P1_full_F | merge 2 | 4.65 ± 0.20 | 62.7 ± 1.2 | 267.3 | 281.2 |
| P1_full_F2 | merge | 4.56 ± 0.16 | 61.1 ± 1.0 | 267.9 | 281.6 |
| P1_full_F2 | merge 2 | 4.56 ± 0.15 | 63.2 ± 1.0 | 268.0 | 281.1 |
| P1_full_h2 | merge | 4.20 ± 0.25 | 56.3 ± 0.5 | 266.9 | 281.2 |
| P1_full_h2 | merge 2 | 4.09 ± 0.09 | 57.1 ± 0.7 | 267.4 | 281.2 |
| P1_subset_F | merge | 4.74 ± 0.23 | 63.2 ± 1.4 | 263.2 | 276.9 |
| P1_subset_F | merge 2 | 4.64 ± 0.18 | 64.0 ± 1.5 | 263.1 | 276.9 |
| P1_full_svp | merge | 21.99 ± 0.74 | 305.8 ± 1.8 | 4007.2 | 4060.1 |
| P1_full_svp | merge 2 | 22.61 ± 0.89 | 305.0 ± 3.2 | 4007.1 | 4060.4 |

Within the scatter, with the same Fock build counts (68, 68, 47, 74, 62) and the lambda tables and F_calc within thread noise of each other.
