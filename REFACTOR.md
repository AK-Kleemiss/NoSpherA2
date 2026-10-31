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
- **The SCF never sees the structure factors directly (decided 2026-10-26).** The SCF talks to the crystallographic side only through the abstract `SCF_coupling` interface (`XCW_solver.h`: `evaluate_cryst_crit`, `cryst_crit`, `perturbation_matrix`, `perturbation_matrix_response`, `save_F_calc`/`restore_F_calc`, `grid_on_device`; renamed in run 2 so that each says what it computes). XCW_wrapper implements it on top of its `structure_factors` object and hands itself to the XCW child of SCF_solver. This keeps the SCF_solver base free of crystallography, so other children (energy decompositions etc.) need no structure factors at all.
- **How SCF_solver is extended (run 1).** `SCF_solver` is not pure abstract: on its own it solves the plain Hartree-Fock system, as "the basic SCF_solver" above should. A child adds a *model term* to the energy through protected virtual hooks, and the DIIS, rescue, TRAH and convergence machinery then works on E + that term without knowing what it is:
  - `evaluate_model_term(scf, timed)`: called at the current density before every Fock build
  - `model_energy()`: the term added to E (XCW: the perturbation lambda * Criterion)
  - `add_model_fock(F)`: its derivative, added to the Fock matrix (XCW: lambda times the perturbation matrix)
  - `add_model_response(scf, dD, dF)` and `model_active()`: the term's part of the TRAH Hessian-vector product
  - `save_model_state()` / `restore_model_state()`: around the finite-difference Hessian check
  - `begin_log()`, `log_iteration(scf)`, `end_log(scf)`, `functional_name()`: what the SCF log prints

  `XCW_solver` (the child) overrides them with the coupling and holds `lambda`; `XCW_wrapper::run_lambda` sets `scf_solver.lambda` and calls `scf_solver.run(hf, wfn, use_guess)`.
- **The XCW quantities (decided 2026-10-27, run 3).** `calc_criteria` computes GooF(F) and GooF(F^2) (`GooF1`, `GooF2`). The GooF(F2) column always prints GooF(F^2). Criterion = GooF^2 against what is refined against (F or F^2, weighted under `h2`), Perturbation = lambda * Criterion, `quant` (Target quantity) = E + Perturbation. `quant_` is the functional the SCF minimises and tests convergence on; there is no separate `phi`.
- **SCF_solver names (decided 2026-10-27, run 3).** `fock_2e` for the two-electron Fock part (was `G`), `grad_alpha`/`grad_beta` for the virtual-occupied Fock blocks (the orbital gradient), `ah_shift` for the augmented-Hessian eigenvalue (was `theta`), `alpha_damp` for the density damping weight (also `SCF_settings::alpha_damp`), `trah_*` for everything of the second-order solver (was `soscf_*`; the on/off flag is `trah_active_`). `quant`, `kappa` and the augmented-Hessian scaling `alpha` in `trah_solve` keep their names. `stored_eri` keeps its name: it holds exactly the four-centre (ab|cd); occ holds the three-centre density-fitting integrals.
- **Only the grid route stays (decided 2026-10-27).** The I tensor is no longer built: F_calc and the perturbation come from the grid every evaluation, which every partition scheme (TFVC/MBIS/EMBIS shares, the RI partition) works on. This replaces "the class should still be used to calculate the I tensor" above.
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
| 1 | 2026-10-27 09:00 | Baseline; SCF_wrapper -> SCF_solver + XCW child; XCW_solver -> XCW_wrapper; comment pass; PR onto density_source (merged as #232) |
| 2 | 2026-10-27 | I tensor route removed, build_I split, self-documenting names, XCW log child, `-print_F` |
| 3 | open | Criterion = GooF^2, quant = E + lambda * Criterion (phi gone), SCF_solver names |

### Agenda
Done in run 1 (details in the git history of #232):
- [x] SCF_solver (plain Hartree-Fock, model hooks) with the XCW child `XCW_solver`; `SCF_coupling` and `fit_quality` in `XCW_solver.h`
- [x] XCW_wrapper drives the scan in named steps (`setup_structure_factors`, `setup_scf_integrals`, `read_start_guess`, `run_lambda`, `approach_lambda`, `extend_scan`) and implements `SCF_coupling` on `structure_factors`
- [x] structure_factors holds the F_calc entry points, the perturbation and its response, the tscb table (`result_table`) and the F_calc save/restore
- [x] Comment and naming pass, `SCF_settings` defaults, `ScfSolverTests` (3)
- [x] density_source merged twice (`0efc6add`, `33fa9685`): `-tshb` step fits, `-xcw_partition ri`, density-following TFVC/MBIS/EMBIS shares, all in the new structure

Run 2 (asked for 2026-10-27):
- [x] Remove the I tensor route (decided 2026-10-27): only the grid route stays, the one every partition scheme works on. Gone: the tensor build and its streaming, the `safe`/`read`/`stream`/`i_tensor_mb`/`i_float` keywords (now unknown keywords), `NOSPHERA2_XCW_GRID` and its comparison mode, `-itensor_hybrid`, `-gpu_itensor_tensor`, `i_tensor_stream.h` and `P1_test_XCW_gpu_itensor`. `i_double` now sets the grid's precision. The five subset goldens run on the grid (tensor to grid: energies within 6e-6 Eh, criteria in the fourth decimal); structure_factors.cpp is 1270 lines instead of 2450
- [x] Output: the lines about the tensor are gone ("XCW on the grid: the I tensor is not built", "Time taken for XCW integrals", "Skipped evaluation of ... grids")
- [x] `build_I` / `eval_I` are `setup_grid` with the auxiliaries `link_symmetry_equivalents`, `screen_ao_pairs`, `setup_atom_grids` and `build_grid_blocks`; the developer statistics dumps of the tensor (`NOSPHERA2_ITENSOR_AOSTATS`, `_SKIPSTATS`) went with it, its tuning knobs are `NOSPHERA2_XCW_GRID_CHUNK`, `_AO_TOL`, `_NO_MORTON`
- [x] SCF log child: `SCF_log_writer` writes the plain SCF's head and iteration lines (iteration, energy), `XCW_log_writer` adds the crystallographic columns and the console lambda table
- [x] `-print_F`: `NA2_<lambda>_Fcalc.txt` only when given
- [x] One-line comments on the calls that hide a lot of work (`XCW_wrapper::run`, `setup_structure_factors`, `setup_scf_integrals`, `run_lambda`, `write_step_results`)
- [x] `SCF_coupling` functions renamed (see Decisions); the model hook `evaluate_model` is `evaluate_model_term`
- [x] REFACTOR.md shortened, tracking in one table
- [ ] GPU/Metal I tensor kernels (`itensor_gpu.cu` build/hold/rows/cols, their `gpu_dispatch.cpp` entries, `itensor_gemm.cuh`, `itensor_metal.mm` and with it `NOSPHERA2_USE_METAL`) are unused since the route went; removed in a commit a GPU build can compile, since nothing here compiles them. `itensor_gpu_available`/`_integrated` and `eri_gpu_*` stay
- [ ] A GPU golden for the grid route (the tensor one, `P1_test_XCW_gpu_itensor`, is gone; it needs a device to record)
- [ ] `-gpu_itensor` / `-cpu_itensor_fp32` and the `i_double` keyword keep their names (they now set the grid); renaming them is a user decision

Run 3 (asked for 2026-10-27):
- [x] Criterion printed and weighed as GooF^2, Perturbation = lambda * Criterion, quant = E + Perturbation; `phi` merged into `quant_` (see Decisions). The goldens' Criterion, Perturbation, Target quantity and Crit(all) columns change; energies and the SCF path do not, the functional was already E + lambda GooF^2
- [x] SCF_solver renames (see Decisions); `stored_eri` keeps its name
- [x] `grid_contract` deals its blocks with `schedule(static, 4)` instead of `dynamic`: the per-thread sums of the perturbation no longer depend on load (squaring the criterion made the ninth decimal of the target quantity vary between identical runs under `-j 6`)

Later:
- [ ] SCF_solver for other theory levels: it is tied to `occ::qm::SCF<occ::qm::HartreeFock>` (`scf_t`); DFT needs the procedure as a template parameter or behind an interface
- [ ] `write_step_results` writes seven kinds of output; an output writer of its own would leave XCW_wrapper with the scan alone
- [ ] The Gaussian halting functions still have blank-line-grouped bodies; reformat them when next touched

## What to look out for
When coming accross a comment, that references a specific case-fault or something of this nature, replace it by a more general description. I want short and concise comments instead of fully fledged explanations for the most part, for theoretically demanding parts like SOSCF or TRAH it is fine to be more elaborate but try sticking to simple language. The code you're rewriting should for the most part be self-commenting, no really abstract abbreviations or variable names, this is seldomly necessary. For for loops that just iterate trivially it is fine to use indices like i and j but if there are a lot of nested for loops it would be better to also describe what is being iterated through the assigned variable. I want the code to be as readable as possible so somebody who is new to the subject could start working on it after a short introduction. There are also code blocks pertaining GPU implementation, these should be marked as such and a short comment should be written explaining what is happening in the code block (there will always be a CPU equivalent to these code blocks so that should be easily done). If you find parts of the code where there is obvious unnecessary overhead like uncollapsed for loops (that have no reason to be uncollapsed) or redundant copying of variables, you are free to remove them and improve the code's efficiency.

## Tracking
The codebase has various tests implemented in it. Changes you make to the code should never break any tests that worked beforehand or slow down the code significantly. Additionally, keep track of memory usage and timings and compare it against each other and the current state of the codebase. These tracked quantities should be visualized in this file using, at least, tables.

How it is measured: AMD Ryzen 7 PRO 8700G (8 cores / 16 threads), 31 GB RAM, Windows 11, MSVC 14.51, `release-windows-tests`, CPU only. `NoSpherA2.exe -do_XCW -no_gpu -cif P1_test.cif -hkl <hkl> -XCW_settings <settings> -no_date` in a scratch copy of `tests/P1_test`, 4-5 runs per case, wall and CPU time and peak working set from the process handle. Cases: `F`, `F2`, `h2` are the full 3216 reflections with `test_settings.txt`, `test_settings_f2.txt`, `h2_test_settings.txt` (sto-3g, lambda 0 to 0.05); `sub` the 402-reflection `P1_subset.hkl`; `svp` the full set with def2-svp. Each build is compared with the one before it in the same sitting, since the machine's speed drifts between sittings.

Wall time in s (mean, scatter ±0.1 s for the small cases and ±0.3-0.9 s for `svp`); peak working set about 267 MB (`sub` 263, `svp` 4008) in every build:

| Build | F | F2 | h2 | sub | svp | ctest (passed/run) |
|-------|--:|---:|---:|----:|----:|-------|
| baseline `a1448cf6` | 4.43 | 4.43 | 4.16 | 4.51 | 21.21 | 1320/1343 |
| run 1 | 4.33 | 4.35 | 3.93 | 4.41 | 20.98 | 1323/1346 |
| merge `0efc6add` | 4.37 | 4.29 | 3.92 | 4.35 | 21.06 | 1330/1355 |
| merge `33fa9685` | 4.65 | 4.56 | 4.09 | 4.64 | 22.61 | 1333/1358 |
| run 2 (no I tensor) | 4.29 | 4.26 | 3.84 | 4.30 | 20.77 | 1322/1346 |
| run 3 (GooF^2, static contraction) | 4.62 | 4.55 | 4.09 | 4.66 | 22.01 | 1329/1346 (`RUN_FULL_TEST=1`) |

Every step was within the run-to-run scatter of the build before it measured in the same sitting (`merge 33fa9685` ran in a busier sitting, its predecessor measured 4.54, 4.56, 4.20, 4.74, 21.99 there), needed the same number of Fock builds (68, 68, 47, 74, 62) and gave lambda tables and F_calc within thread noise of each other. Run 2's predecessor measured 4.37, 4.39, 4.03, 4.41, 21.05 in its sitting; the 12 fewer tests are the tensor's (file, streaming, device golden). Run 3 sat in a slower sitting: there its predecessor (the same code with the dynamic schedule) measured 4.43, 4.48, 4.06, 4.62, 21.62; the 0.03-0.4 s are the static schedule, `static, 1` (22.17 s svp) and `static, 16` (23.30 s) were slower. Its lambda tables match its predecessor's. With `RUN_FULL_TEST=1` the four `P1_*_full` goldens pass in every build.
