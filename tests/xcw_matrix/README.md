# XCW permutation-test fixture (Brief 3)

A tiny XCW case that runs end to end in well under a second, for the weighting × anharmonic × extinction ×
grown × F/F² × target matrix. This README covers the fixture, its timings, the stored-ERI decision and the matrix.
The matrix runs in ctest as `XcwMatrixTests` (tests/src/XcwMatrixTests.cpp) against `matrix.good`; see
[ctest](#ctest) below.

Measured 31 Oct 2026 on FLOWOFFICE with `spaces/3f316b01-r1sym/NoSpherA2.exe`: density_source plus the grown-mode R1
fix, in which grown F is summed over the asymmetric unit with all CIF operations.

## Fixture

| file | what |
|---|---|
| `nh3.cif` | NH3 in P2₁2₁2₁ (4 ops, non-centrosymmetric), harmonic; N1 Uani, H Uiso 0.053, λ 0.71073 Å |
| `nh3_C.cif` | same + Gram-Charlier C (10 terms) on N1 |
| `nh3_CD.cif` | same + Gram-Charlier C and D (10 + 15 terms) on N1 |
| `nh3_grown.xyz` | asymmetric unit + the neighbour molecule that H1 points at (H1···N 2.398 Å), 8 atoms, Å |
| `nh3.hkl` | 364 unique reflections, h,k,l ≥ 0 to sin θ/λ 0.80 Å⁻¹, HKLF-4 (I, σ(I)); 317 have I/σ ≥ 2 |
| `make_fixture.py` | regenerates all of the above (about 2 s) |
| `run_cases.py` | runs, times and parses cases (timing grid, ERI test, reference grid, matrix) |

**Structure.** Ammonia: the P2₁3 structure from `tests/cytidine_tonto/nh3_XCW`, a = 5.1305 Å, described in its
subgroup P2₁2₁2₁. This makes the whole molecule the asymmetric unit: 4 atoms and 8 sto-3g functions, with 4
operations that F_calc has to sum over.

The grown model keeps those 4 atoms and adds 4 grown atoms (SCF.log reports `4 asymmetric unit atoms, 4 grown, 8 in
the calculation`). Both the grown path and the symmetry path are therefore exercised.

**Intensities.**
- **Truth model:** one XCW at lambda 0 with def2-SVP, ungrown, using `nh3_CD.cif` and `extinction shelx fixed 0.2`
  (SHELX x on NoSpherA2's absolute |F_calc| scale).
- **True amplitudes:** the `-print_F` column scale·√y·|F_calc|, divided by that file's scale.
- **Observed intensities:** I_obs = (10·F)² + N(0, σ), where σ(I) = √(I + (0.01·I)² + 2²) and the seed is 31.
- **Check:** refitting with the truth basis (def2-SVP, lambda 0) gives R1 2.2 % and GooF 1.00, so the data are
  consistent.

With sto-3g the model is deliberately too small. That gives R1 4.3–5.2 % and GooF ≈ 2, with room for lambda to
improve the fit (R1 drops by about 0.5 pp over lambda 0 → 0.02).

The model misfit also explains why refined extinction lands at 0.35–0.41 (normal target) or 0.23–0.25 (H²
target) rather than the true 0.2. That is not a bug.

```
python make_fixture.py --exe D:/git/NoSpherA2/spaces/3f316b01-r1sym/NoSpherA2.exe   # scratch in E:/Temp/xcw_matrix_scratch/gen
python run_cases.py --exe <exe> --cases timing|timing_H2|eri|reference|matrix|matrix_sigma [--cpus 1] [--repeat 2] [--conv "sloppy normal_conv"]
```

### Settings line used everywhere

```
<conv> params 22 basis_set sto-3g F|F2 charge 0 mult 1 rhf start 0 step_size 0.01 end 0|0.02 [grown] [weighted] [weight ...] [extinction shelx [fixed 0.2]]
```

- Grown runs add `-xyz nh3_grown.xyz` to the command line.
- Every run passes `-no_gpu -no_date -print_F -cpus N`.
- `extinction` used to have to be the last keyword: a bare `extinction shelx` read the next token as its optional
  argument and stopped with `Could not read 'save_eri' after 'extinction'`. Fixed in ed2fd992: the value is
  optional, `fixed` and `weighted` come in any order, and the list ends at the next settings keyword
  (`ExtinctionTests.SettingsKeywordStopsAtTheNextKeyword`).

## Timings

The machine was at about 94 % load from a parallel build, so wall time is noisy. The table uses the child's CPU time
(user + sys, from psutil). It covers all 144 cells: target × cluster × F/F² × {none, C, CD} × {off, fixed, refined} ×
{1, 3} lambdas, best of 2, `-cpus 1`, `normal normal_conv`.

| factor | CPU median (s) | note |
|---|---|---|
| all 144 cells | 0.34 (sum 62.7) | wall median 0.59 under load |
| 1 lambda / 3 lambdas | 0.16 / 0.53 | max 0.33 s for one lambda |
| ungrown / grown | 0.30 / 0.38 | grid 10 986 → 21 972 points |
| F / F² | 0.34 / 0.30 | |
| anh none / C / CD | 0.36 / 0.34 / 0.27 | anharmonicity costs nothing measurable |
| ext off / fixed / refined | 0.31 / 0.38 / 0.35 | |
| target normal / H² | 0.30 / 0.39 | |
| max | 5.14 | `H2_u_F_CD_refined_3`, the TRAH stall below; `H2_u_F_none_refined_3` 4.73 |

With `-cpus 4` the wall time is about the same: the median is 0.56 s, because everything this small is serial
overhead. Results are identical between `-cpus 1` and `-cpus 4` to the printed digits for all 32 reference cells.

### TRAH ↔ DIIS ping-pong (fixed in ed2fd992)

**Fix:** the SCF solver now remembers the target and gradient at each TRAH → DIIS hand-off. The second hand-off that
improves neither ends that SCF ("N hand-offs to DIIS without a lower ... giving up this SCF"). The lambda step is then
accepted or half-stepped as usual. On the `H2_u_F_*_refined_3` cells:

| cell | Fock builds | wall |
|---|---|---|
| none, cpus 1 | 261 → 112 | 5.6 → 2.5 s |
| CD, cpus 1 | 261 → 97 | 5.6 → 2.0 s |
| CD, cpus 4 | 261 → 97 | 2.7 → 0.9 s |

The final R1, wR2, GooF, scale and extinction are unchanged. The original diagnosis follows.

`weighted` + F + refined extinction, ungrown, lambda 0.02, `normal` convergence:
- Every criterion converges except the density RMSD, which sits at 1.7e-6 against a threshold of 1e-6.
- TRAH predicts −2e-8 Eh but the target rises by 3–5e-8 Eh: rho < 0, the trust radius drops to 1e-4, and the run
  goes "back to DIIS".
- DIIS returns to the same point and TRAH restarts.

This repeats until the 100-iteration cap (170 Fock builds), and then `retrying lambda 0.015` rescues it. The final
numbers are fine; the cost is 5 s instead of 0.6 s.

The stall depends on round-off. `H2_u_F_none_refined_3` stalls with `-cpus 1` but not with `-cpus 4`, while
`H2_u_F_CD_refined_3` stalls with both.

A likely cause is that extinction and scale are re-fitted inside the target, which leaves noise of about 1e-8 Eh in
the target near convergence. The TRAH ↔ DIIS hand-off has no exit for "no progress in either direction".

The former workaround was `sloppy normal_conv`. Over the 72 sigma-weighted H² cells it agreed to the printed
digits. The matrix no longer uses it: with the stall fixed, `normal` fits the budget. `sloppy` also moves the
`chebyshev_opt` rows (m19 wR2 0.1433 sloppy vs 0.1719 normal; ON moves too), because their coefficients are
re-fitted on a less converged model.

## Stored ERIs: recommendation (c), nothing

| model | `save_eri` file | read back? |
|---|---|---|
| ungrown (8 functions, 36 pairs) | 6 180 B | yes, `0.00508 MB ... read from nh3.eri` in SCF.log |
| grown (16 functions) | 76 372 B | yes |

Reading the ERIs back gives bit-identical R1, wR2, scale and extinction. It saves nothing measurable: 0.47 s versus
0.53 s CPU for 3 lambdas, inside the noise. The integrals for 8–16 sto-3g functions cost less than process start-up.

- A committed `.eri` would be tiny, but it buys nothing. It would also go stale silently: `read_eri` rebuilds on a
  basis or geometry mismatch, so a stale file costs nothing but hides that it is stale.
- A FIXTURES_SETUP step would add a serial dependency for no gain.

## Matrix

The factors are:
- weighting: 8 levels;
- anharmonic: none, C, C+D;
- extinction: off, fixed (0.2), refined;
- cluster: ungrown (u), grown (g);
- F or F²;
- target: normal, or H² (`weighted`, XWR_type 2).

Rows 1–24 form a pairwise covering design: every pair of levels of any two factors appears at least once. 24 is the
minimum, because weighting × anharmonic alone has 24 pairs. `run_cases.check_pairwise()` asserts the coverage.

Row D is the all-default row and row ON is the all-on row. Every row runs 3 lambda steps (0, 0.01, 0.02), so the
`_opt` and `robust` schemes get updates between steps. Every row uses `normal normal_conv`.

| id | weighting | anh | ext | cluster | F/F² | target |
|---|---|---|---|---|---|---|
| D | (default) | none | off | u | F | normal |
| 1 | sigma | none | refined | u | F | H² |
| 2 | sigma | C | fixed | g | F² | normal |
| 3 | sigma | CD | off | g | F² | H² |
| 4 | unit | none | off | g | F | normal |
| 5 | unit | C | refined | u | F² | normal |
| 6 | unit | CD | fixed | u | F | H² |
| 7 | stl | none | fixed | u | F² | normal |
| 8 | stl | C | off | g | F | H² |
| 9 | stl | CD | refined | g | F | normal |
| 10 | shelx | none | off | u | F² | H² |
| 11 | shelx | C | fixed | g | F | normal |
| 12 | shelx | CD | refined | g | F | normal |
| 13 | shelx_opt | none | off | g | F | H² |
| 14 | shelx_opt | C | fixed | u | F² | normal |
| 15 | shelx_opt | CD | refined | g | F | H² |
| 16 | chebyshev | none | refined | g | F | normal |
| 17 | chebyshev | C | off | u | F² | H² |
| 18 | chebyshev | CD | fixed | g | F | H² |
| 19 | chebyshev_opt | none | refined | g | F | H² |
| 20 | chebyshev_opt | C | off | u | F² | normal |
| 21 | chebyshev_opt | CD | fixed | g | F² | H² |
| 22 | robust(shelx) | none | refined | u | F² | normal |
| 23 | robust(shelx) | C | fixed | g | F | H² |
| 24 | robust(shelx) | CD | off | g | F | normal |
| ON | chebyshev_opt | CD | refined | g | F² | H² |

The weight keywords (`run_cases.WEIGHT_KW`, the d7e0eb86 syntax) are:
- `weight unit`
- `weight stl 1`
- `weight shelx 0.1 0`
- `weight shelx_opt 0.1 0`
- `weight chebyshev 1 0.5 0.2`
- `weight chebyshev_opt 3`
- `weight robust shelx 0.1 0 6`

`sigma` rows pass no keyword.

**unit, stl and chebyshev are rescaled (a88ceae4).** These weights are absolute. On F², w = 1 against Fo² ~ 1e3–1e4
gave a criterion of order 1e3, and the scan crawled in half steps: m5 took 463 s and m17 298 s, and both stopped at
lambda ≈ 0.003.

The first weight fill now applies one CRYSTALS REWEIGHT-style factor k = 1/S², so that S = 1 under the scan's own
criterion, and holds it fixed. The knobs are `weight_rescale off` and `-xcw_weight_rescale 0`.

| row | before | after | k |
|---|---|---|---|
| m5 | 463 s, lambda 0.00305 | 0.38 s, lambda 0.02 | 7.40e-5 |
| m17 | 298 s, lambda 0.00328 | 0.69 s, lambda 0.02 | 1.58e-5 |

k comes from the initial guess. Scale and extinction are first settled on the guess under 1/σ² (as `eval_scale` does);
the scheme weights and k are then filled from that model, and the lambda 0 SCF starts with them. The criterion at
lambda 0 is therefore 0.664 (m5) and 0.963 (m17), not exactly 1. Settling extinction first moved the rows that combine
refined extinction with a non-sigma scheme (m5, m9, m16, m19 beyond tolerance; m12, m15, m22, mON in the last digits).

## ctest

`XcwMatrixTests.MatchesGolden/m<id>` runs each of the 26 rows in process with `-cpus 1`, in a scratch directory under
the system temp. It compares against `matrix.good`:

| quantity | tolerance |
|---|---|
| final lambda | 1e-6 |
| R1, wR2 (`XCW final model:` line) | 2e-4 abs |
| GooF, scale (`NA2_*_Fcalc.txt` header), extinction (last `XCW lambda ...: extinction` line in SCF.log) | 2e-3 rel |

A `-` in the extinction column means no extinction line is allowed. Nothing in `matrix.good` is a time.

**Measured 31 Oct 2026, FLOWOFFICE:** `ctest -R XcwMatrix -j4` takes 5.8 s for all 26 rows. The slowest row takes
1.7 s. With `-cpus 1` and `-cpus 4`, the numbers agree to the printed digits.

**Regenerate** after a deliberate change, and say in the commit which numbers moved and why:

```
python run_cases.py --exe D:/git/NoSpherA2/build/release-windows/bin/NoSpherA2.exe --cases matrix --cpus 1 --golden matrix.good --work E:/Temp/xcw_matrix_scratch/m_golden
```

`--golden` asserts every row exits 0 and that `check_pairwise()` holds.

## Reference values (sigma weighting, runnable today)

These tables predate d7e0eb86, which prints wR2 with the scheme weights. The wR2 below is `run_cases.py`'s own
recompute, from `nh3.hkl` and the last lambda's `NA2_*_Fcalc.txt`:
wR2 = √(Σw(Fo² − Fc²)² / ΣwFo⁴) over all 364 reflections, with w = 1/σ(Fo²)² and Fc = scale·√y·|F_calc|.

The other values come from:
- R1(gt), R1(all), GooF(F²) and Crit: the last row of the lambda table in NoSpherA2.log;
- scale: the Fcalc header;
- extinction x: the last `XCW lambda ...: extinction` line in SCF.log.

All values are at the final lambda 0.02. The grown rows depend on the 3f316b01 R1 fix.

**Matrix rows** (`--cases matrix_sigma`; `sloppy` and `normal` agree to these digits):

| id | R1(gt) | R1(all) | wR2 | GooF | Crit | scale | ext x |
|---|---|---|---|---|---|---|---|
| D | 0.04931 | 0.05341 | 0.09634 | 2.151 | 4.764 | 9.970 | – |
| 1 | 0.04448 | 0.04862 | 0.08908 | 1.987 | 5.062 | 9.904 | 0.2515 |
| 2 | 0.04622 | 0.05065 | 0.09167 | 2.044 | 4.178 | 9.982 | 0.2 (fixed) |
| 3 | 0.05026 | 0.05452 | 0.09806 | 2.192 | 6.210 | 9.725 | – |

**Reference grid** (`--cases reference`, `normal normal_conv`, 3 lambdas):

| target | cluster | F/F² | anh | ext | R1(gt) | R1(all) | wR2 | GooF | Crit | scale | ext x |
|---|---|---|---|---|---|---|---|---|---|---|---|
| normal | ungrown | F | none | off | 0.04931 | 0.05341 | 0.09634 | 2.151 | 4.764 | 9.97 | – |
| normal | ungrown | F | none | refined | 0.04392 | 0.04812 | 0.09067 | 2.023 | 4.125 | 10.08 | 0.4071 |
| normal | ungrown | F | CD | off | 0.04848 | 0.05290 | 0.09532 | 2.127 | 4.649 | 9.957 | – |
| normal | ungrown | F | CD | refined | 0.04323 | 0.04775 | 0.08985 | 2.004 | 4.038 | 10.06 | 0.3944 |
| normal | ungrown | F² | none | off | 0.04956 | 0.05363 | 0.09601 | 2.144 | 4.598 | 9.932 | – |
| normal | ungrown | F² | none | refined | 0.04433 | 0.04850 | 0.09021 | 2.013 | 4.054 | 10.03 | 0.3632 |
| normal | ungrown | F² | CD | off | 0.04865 | 0.05305 | 0.09499 | 2.120 | 4.495 | 9.92 | – |
| normal | ungrown | F² | CD | refined | 0.04361 | 0.04810 | 0.08939 | 1.994 | 3.975 | 10.02 | 0.3515 |
| normal | grown | F | none | off | 0.04947 | 0.05357 | 0.09647 | 2.154 | 4.782 | 9.964 | – |
| normal | grown | F | none | refined | 0.04404 | 0.04825 | 0.09075 | 2.025 | 4.136 | 10.07 | 0.4090 |
| normal | grown | F | CD | off | 0.04864 | 0.05305 | 0.09543 | 2.130 | 4.663 | 9.951 | – |
| normal | grown | F | CD | refined | 0.04336 | 0.04787 | 0.08992 | 2.005 | 4.047 | 10.05 | 0.3962 |
| normal | grown | F² | none | off | 0.04972 | 0.05379 | 0.09615 | 2.147 | 4.612 | 9.925 | – |
| normal | grown | F² | none | refined | 0.04446 | 0.04864 | 0.09029 | 2.015 | 4.061 | 10.03 | 0.3651 |
| normal | grown | F² | CD | off | 0.04879 | 0.05319 | 0.09510 | 2.123 | 4.506 | 9.913 | – |
| normal | grown | F² | CD | refined | 0.04373 | 0.04822 | 0.08945 | 1.995 | 3.981 | 10.01 | 0.3534 |
| H² | ungrown | F | none | off | 0.05090 | 0.05485 | 0.09859 | 2.204 | 6.651 | 9.752 | – |
| H² | ungrown | F | none | refined | 0.04448 | 0.04862 | 0.08908 | 1.987 | 5.062 | 9.904 | 0.2515 |
| H² | ungrown | F | CD | off | 0.04959 | 0.05387 | 0.09707 | 2.169 | 6.528 | 9.748 | – |
| H² | ungrown | F | CD | refined | 0.04344 | 0.04786 | 0.08775 | 1.956 | 4.968 | 9.897 | 0.2481 |
| H² | ungrown | F² | none | off | 0.05132 | 0.05525 | 0.09944 | 2.224 | 6.305 | 9.737 | – |
| H² | ungrown | F² | none | refined | 0.04494 | 0.04906 | 0.08996 | 2.008 | 4.888 | 9.881 | 0.2334 |
| H² | ungrown | F² | CD | off | 0.05000 | 0.05427 | 0.09788 | 2.188 | 6.202 | 9.733 | – |
| H² | ungrown | F² | CD | refined | 0.04385 | 0.04825 | 0.08857 | 1.976 | 4.806 | 9.875 | 0.2304 |
| H² | grown | F | none | off | 0.05118 | 0.05512 | 0.09879 | 2.209 | 6.674 | 9.744 | – |
| H² | grown | F | none | refined | 0.04459 | 0.04874 | 0.08919 | 1.990 | 5.091 | 9.896 | 0.2499 |
| H² | grown | F | CD | off | 0.04984 | 0.05411 | 0.09724 | 2.173 | 6.544 | 9.739 | – |
| H² | grown | F | CD | refined | 0.04353 | 0.04795 | 0.08782 | 1.958 | 4.990 | 9.889 | 0.2464 |
| H² | grown | F² | none | off | 0.05160 | 0.05553 | 0.09965 | 2.229 | 6.319 | 9.729 | – |
| H² | grown | F² | none | refined | 0.04506 | 0.04918 | 0.09010 | 2.011 | 4.908 | 9.873 | 0.2317 |
| H² | grown | F² | CD | off | 0.05026 | 0.05452 | 0.09806 | 2.192 | 6.210 | 9.725 | – |
| H² | grown | F² | CD | refined | 0.04397 | 0.04837 | 0.08867 | 1.978 | 4.820 | 9.867 | 0.2287 |

Grown and ungrown agree to within 3e-4 in R1, which is the expected result once grown F uses the asymmetric unit
with all operations. `NOSPHERA2_XCW_GROWN_COPYSUM=1` (the old copy-sum) is not covered here.

## Open items

1. wR2 is printed now (`XCW final model:`, with the scheme weights). `run_cases.py` keeps its sigma-weighted
   recompute as a cross-check: it asserts the recompute against `wR2(all)` on sigma rows.
   - The scale still comes from the `-print_F` header.
2. The rescale k comes from the initial guess (settled scale and extinction), not from a converged lambda-0 model, so
   S at lambda 0 is not exactly 1. This is deliberate: the weights are fixed before the lambda 0 SCF starts.
3. Goldens never compare timings. The `wall` and `cpu` columns of `run_cases.py` are for this README only.
