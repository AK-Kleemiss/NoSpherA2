# The held-out set: 16 molecules the NAO change was never tuned on, and their baseline column

`PREREGISTRATION.md` beside this file is the design, and it was committed before the
first stage-2 number of the set existed. This file is the record of what was actually
built and what the baseline binary measures on it. Read the pre-registration first; a
prediction written after the measurement is not a prediction, and the only reason the
predictions in that file count for anything is that they are in git ahead of these
numbers.

**Nothing in here agrees with NBO 7 and nothing in here is expected to.** The acceptance
gate is an *exact* basis-function match against gennbo's own output, and `-nbo_native`
fails it on all 16 as it does on the 22: `compare_all.py` reports **0 PASS, 16 FAIL**
over 15980 compared data points (8529 in the gated view, 7451 dropped). Every number
below is a measured delta against gennbo's printed values with its instrument floor
beside it.

## The 16, and why each one is in the set

| molecule | charge, mult | what it is in the set for |
|---|---|---|
| `sih4` | 0, 1 | second-row closed shell, the plain silicon reference for `si2h6` |
| `si2h6` | 0, 1 | ethane one row down - ethane is one of the eight the change was selected on |
| `cl2` | 0, 1 | third-row homonuclear, no polarity for the cascade to lean on |
| `clf3` | 0, 1 | **hypervalent beyond sf6/pf5**, and T-shaped rather than symmetric |
| `thiophene` | 0, 1 | aromatic ring with a third-row heteroatom; benzene's homologue |
| `bh3` | 0, 1 | electron-poor, an empty p on boron, a **diffuse and well-populated Rydberg set** |
| `mg_h2` | 0, 1 | s-block metal hydride, diffuse valence |
| `nacl` | 0, 1 | ionic closed shell, the cascade's classification under a near-full charge transfer |
| `nh4_plus` | +1, 1 | **cation** |
| `bf4_minus` | -1, 1 | **anion**, and four equivalent ligands |
| `zn_cl2` | 0, 1 | **transition metal beyond ni_co_4/ticl4**, d10, so a full d shell in the valence set |
| `hi` | 0, 1 | **ECP-bearing heavy atom** (def2-ECP on iodine, 28 core electrons replaced) |
| `hs` | 0, 2 | **open-shell radical**, third row |
| `cf3` | 0, 2 | **open-shell radical**, pyramidal, three fluorines |
| `c2h5` | 0, 2 | open-shell hydrocarbon large enough for the instrument to have range |
| `ch2` | 0, 3 | triplet, two unpaired electrons rather than one |

d/f basis functions with a non-identity component-major -> shell-major permutation are
present in every def2-TZVP member (the f sets on the second- and third-row atoms), and
`zn_cl2` carries the d shell in the *valence* set rather than as polarisation.

## The ECP case works, and its cost is visible rather than hidden

`hi` was the one member the briefing allowed to be dropped if gennbo or the drivers could
not handle it. It was not dropped: ORCA applies def2-ECP, the wavefunction has **26
electrons** (H 1 + I 53 - 28 replaced), `-convert_to_47`, gennbo, `-nbo_parse` and
`-nbo_native` all return 0, and both JSONs are `parser_version: 2`.

What it costs is recorded and not smoothed over:

- `nao_class_leak.py` prints `hi 3 rows classified into a different class by the two
  sides` before its row, and `d(Cor) -0.00831 / d(Val) +0.00823` - a core/valence
  reassignment of three rows, not a population error.
- gennbo prints iodine's first NAO as `Cor( 1s) 2.00000`, which under an ECP cannot be a
  1s: NBO labels a shell by its order within `l`, not by its true principal quantum
  number. Nothing in this directory may be read off a printed shell label.
- The per-class sine metric **saturates** on `hi` at 1.00 for Cor and Val (see below).
  A saturated metric is an upper bound, not a magnitude.

## What each molecule directory holds

Same set as `../data/`, same names, so one reader learns one layout:

| file | stage | what it is |
|---|---|---|
| `<mol>.inp` | 1 | the ORCA input that was actually run |
| `<mol>_start.xyz` | 1 | the idealised geometry `make_holdout.py` writes |
| `<mol>.out` | 1 | ORCA's output, which acceptance reads the basis count and energy off |
| `<mol>.xyz`, `<mol>_trj.xyz` | 1 | the optimised geometry and the trajectory |
| `provenance_orca.json` | 1 | hostname, threads, node cores, slurm job id, ORCA path |
| `<mol>.nbo` | 2 | **gennbo 7's own output**, kept as text so a parser fix is a local re-parse |
| `<mol>.gennbo.nbo.json` | 2 | that text parsed by `-nbo_parse` |
| `<mol>.native.nbo.json` | 2 | NoSpherA2's `-nbo_native` answer, **renat5 arm** (no env var) |
| `<mol>.baseline.native.nbo.json` | 2 | the same, **baseline arm** (`NAO_LEGACY_CASCADE=1`), job 595959 |
| `provenance_nbo.json` | 2 | hostname, threads, NoSpherA2 commit, gennbo path, keylist, flags, per-step seconds |

Both arms are committed so the before/after is checkable without the cluster:
`py -3.12 ../holdout_selfcheck.py` rebuilds the `d(Ryd)` table below from these files alone and
asserts that every molecule's two arms differ by more than its own floor. It reproduces the
table to the last printed digit; if it ever does not, the README is quoting numbers the tree
cannot produce.

`accepted_holdout.json` is the acceptance record for all 16 at once. All 16 are
`accepted: true`. It is a **weaker** criterion than `../data/`'s: none of these molecules
is in `tests/nbo_reference/index.json`, so there is no stored basis-function count to
match exactly. What was checked instead is a converged optimisation, a converged last
SCF, the charge/multiplicity/basis/keyword line as asked, and `<S**2>` within 4 % of
`S(S+1)`; the basis-function count and final energy are *recorded* so the next
regeneration has the exact test this one could not have.

## Deliberately not committed

On the cluster at `/work/akkleemiss/florian/nbo_ref_v2_holdout/<mol>/`:

- `<mol>.gbw` - ORCA's wavefunction.
- `<mol>.47` and `<mol>_native.47` - the NBO archives.
- the per-class sine probe's matrices, in `.../aopnao/<mol>/`: `<mol>.32` (gennbo's
  pre-NAOs in the AO basis), `<mol>.33` (its NAOs), `<mol>.naocpre.txt` and
  `<mol>.naoc.txt` (native's own two, dumped by `NAO_DUMP_CPRE` / `NAO_DUMP_C`). 11 MB
  for the 16, and every one of them is one command away.

Same rule as `../data/`: derived, one command each, so the repository keeps only what
cannot be rebuilt from the inputs. Whether the `.47` set belongs in git after all is
Florian's call, not this branch's.

## Reproducing all of it from one command

Everything is on the cluster (`ssh AKL007`), scripts one level up in
`tests/nbo_reference_v2/`. `$H = /work/akkleemiss/florian/nbo_ref_v2_holdout`.

```
python make_holdout.py $H                     # 16 x <mol>/<mol>.inp + <mol>_start.xyz
for m in ...; do sbatch holdout_stage.sh $m; done   # stage 1, one job per molecule
python accept_holdout.py $H                   # -> accepted_holdout.json
for m in ...; do sbatch nbo_stage.sh $m; done  # stage 2: .47, .nbo, both JSONs, provenance
python compare_all.py $H                      # the gate: 0 PASS, 16 FAIL
sbatch --export=ALL bin/twoarm.sh             # BOTH arms of one binary, 16 x 2, md5 per run
python bin/twoarm_table.py $H $H/twoarm       # the before/after table
python bin/storedchk.py  $H $H/twoarm         # stored column == renat5 arm, and the real floors
python nao_class_leak.py $H/twoarm/baseline   # worst L, baseline arm  (and $H for renat5)
LEGACY=1 SRC=$H V2=$H/bin BIN=$H/bin/NoSpherA2_519cb79 OUT=$H/aopnao_base \
    MOLS="bf4_minus bh3 cl2 clf3 hi mg_h2 nacl nh4_plus si2h6 sih4 thiophene zn_cl2" \
    sbatch --export=ALL bin/aopnao_stage.sh    # gennbo's lfn 32/33 + native's two dumps
python aonao_compare.py --space $H/aopnao_base # column 3, baseline arm (drop LEGACY for renat5)
```

`bin/NoSpherA2_519cb79` is the one binary both arms come from: `NAO_LEGACY_CASCADE=1` selects
the pre-change cascade, unset selects renat5. `bin/NoSpherA2_baseline` is a hard link to the
same inode kept only so older job scripts resolve - **its name is a lie**, and
`bin/DEPLOYED_COMMIT.txt` now says so instead of naming a commit.

Four traps, all paid for already:

- The login node's `python3` has **no numpy**, so `aonao_compare.py` dies at its import.
  Use `/work/akkleemiss/florian/fp_fdp/bin/python3` (numpy 2.0.2, Python 3.9.25).
- `aopnao_stage.sh` writes `<mol>.aopnao.nbo.json` and `aonao_compare.py --space` reads
  `<mol>.aonao.nbo.json`. Copy one to the other before the comparison.
- **Never take a binary's identity from a file beside it.** `md5sum` it inside the job. The
  first version of this README got its whole baseline column from a file whose hand-written
  label said `60055d4b` and whose bytes were the treatment build.
- `aopnao_stage.sh` does not echo its environment, so `LEGACY=1` was verified by effect:
  its dumps must match `twoarm/baseline` at 0.000e+00 and differ from `twoarm/renat5`.

## CORRECTION: the column first committed here as "baseline 60055d4b" is the renat5 arm

Commit `d259879e` labelled its column "binary 60055d4b, before the NAO change". It is the
**treatment** column. The file it used, `bin/NoSpherA2_baseline`, is byte-identical to the
renat5 build:

| file | md5 | bytes |
|---|---|---|
| `nbo_ref_v2_holdout/bin/NoSpherA2_baseline` (used for that column) | `1e45bc35f8d5e9617ee42dc62ca7d180` | 87952392 |
| `nao_renat5/bin/NoSpherA2_519cb79` (the change) | `1e45bc35f8d5e9617ee42dc62ca7d180` | 87952392 |
| `nao_renat5/bin/NoSpherA2.baseline` (the true pre-change binary) | `99a8aa340337a4d73d043084b96234a3` | 87952216 |

Same md5 as the treatment, 176 bytes larger than the baseline. The cause is a race, not an
inference error: the file was copied out of the shared tree's
`build/release-linux/bin/NoSpherA2` at 16:34, and the sibling session had replaced that path
with its renat5 build at 16:29 (restoring it at 16:30 - the mtimes still read 16:30 and
16:34). `bin/DEPLOYED_COMMIT.txt` said `60055d4b` because it was **typed by hand**, and
`naoarms.sh` echoed that file as if it were provenance. **A label beside a binary is not a
property of the binary.** Every arm below prints the `md5sum` of the file it invoked and the
value of `NAO_LEGACY_CASCADE` it saw, per molecule, into the job output.

The measurement itself was sound: re-taking the native side with the same file and no
environment variable reproduces the committed column at **0.00000 e on all 16 molecules**
(`bin/storedchk.py`). The numbers were mislabelled, not miscomputed.

## Both columns, from one binary and one environment variable

`bin/NoSpherA2_519cb79` run with `NAO_LEGACY_CASCADE=1` is the pre-change cascade and with
the variable unset is renat5, so the before/after comes from **one file** - which is what a
two-build comparison cannot give. Job 595959, 16 molecules x 2 arms, 32 of 32 `rc=0`.

**Discrimination check before the wave, not after** (`bin/armdisc.py` on the first molecule,
`bf4_minus`): 335 numeric fields compared, worst arm-to-arm difference **3.828** (an NAO
energy). The arms are wired up. Timings are excluded from that check on purpose - two runs of
one binary differ in `timings` and in nothing else, so a differ-anywhere test would pass on a
knob that does nothing. Had this check existed at 16:34 it would have caught the mixup free.

`d(Ryd)` = total Rydberg population, native minus gennbo. Its floor is **not** 1e-5: it is a
sum over N atoms of a five-decimal printed quantity, so the floor is `N x 0.5e-5`, from
1.0e-05 for `cl2` to 4.5e-05 for `thiophene`. `worst L` = largest per-atom intra-atomic
Rydberg leak; `max |dq|` = largest per-atom NPA **charge** deviation; both at 1e-5.

| molecule | floor | d(Ryd) base | d(Ryd) ren5 | factor | worst L base | worst L ren5 | abs dq base | abs dq ren5 |
|---|---|---|---|---|---|---|---|---|
| `bf4_minus` | 2.5e-05 | +0.00444 | -0.00022 | 20.0 | +0.00127 | -0.00005 | 0.01462 | 0.00275 |
| `bh3` | 2.0e-05 | +0.03707 | +0.00084 | 44.1 | +0.02113 | +0.00041 | 0.04251 | 0.00598 |
| `c2h5` | 3.5e-05 | +0.10884 | -0.00009 | 1199.6 | +0.04625 | -0.00070 | 0.01129 | 0.02550 |
| `cf3` | 2.0e-05 | +0.00295 | **-0.00614** | **0.48** | -0.00375 | -0.00697 | 0.03476 | 0.01442 |
| `ch2` | 1.5e-05 | +0.01770 | -0.00132 | 13.4 | +0.00805 | -0.00199 | 0.06055 | 0.05138 |
| `cl2` | 1.0e-05 | -0.00010 | **-0.00024** | **0.42** | -0.00006 | -0.00013 | 0.00000 | 0.00000 |
| `clf3` | 2.0e-05 | -0.00247 | -0.00239 | 1.03 | -0.00269 | -0.00071 | 0.00692 | 0.00019 |
| `hi` | 1.0e-05 | +0.00272 | +0.00008 | 34.4 | +0.00338 | +0.00011 | 0.00066 | 0.00003 |
| `hs` | 1.0e-05 | +0.00394 | +0.00042 | 9.5 | +0.00646 | +0.00048 | 0.01301 | 0.00895 |
| `mg_h2` | 1.5e-05 | -0.03591 | -0.00480 | 7.5 | -0.02652 | -0.00188 | 0.03772 | 0.00510 |
| `nacl` | 1.0e-05 | -0.00292 | +0.00010 | 28.7 | -0.00288 | +0.00007 | 0.00035 | 0.00013 |
| `nh4_plus` | 2.5e-05 | +0.05591 | +0.00348 | 16.1 | +0.04496 | +0.00342 | 0.01078 | 0.00271 |
| `si2h6` | 4.0e-05 | +0.04931 | +0.00103 | 48.1 | +0.01333 | +0.00040 | 0.03602 | 0.00863 |
| `sih4` | 2.5e-05 | +0.02473 | +0.00002 | 1454.6 | +0.01027 | +0.00171 | 0.06938 | 0.01297 |
| `thiophene` | 4.5e-05 | +0.19334 | +0.00160 | 121.0 | +0.04990 | -0.00155 | 0.01396 | 0.00553 |
| `zn_cl2` | 1.5e-05 | -0.00904 | +0.00099 | 9.1 | -0.00557 | +0.00070 | 0.00150 | 0.00031 |

Mean over the 68 atoms: intra-atomic leak **0.00841 e** baseline against **0.00051 e**
renat5; charge deviation **0.01264 e** against **0.00466 e**.

Two factors in that table are not measurements of a magnitude. `sih4`'s 1454.6 is
+0.00002 e against a 2.5e-05 floor - the renat5 value is **at its floor**, so the factor is
bounded below by the instrument, not by the cascade. `c2h5`'s 1199.6 is the same shape at
0.3 floors. Both say "indistinguishable from gennbo on this instrument", not "1200x better".

### The pre-registered claims, decided

- **H1 (transfer, worst `dpop` improves >= 5x): HOLDS.** Worst over the held-out set
  0.19334 e baseline against 0.00614 e renat5, a factor **31.5**. Median per-molecule factor
  **18.0**.
- **H2 (no member worse than 1.5x): FAILS, on two members.** `cf3` 0.00295 -> 0.00614
  (2.08x worse, both far above its 2.0e-05 floor) and `cl2` 0.00010 -> 0.00024 (2.4x worse,
  10 and 24 floors - small, but resolved). So "closer on 8 of 8, no counter-molecule" does
  **not** generalise: held-out data produces two counter-molecules.
- **The nominated counter-molecule was the wrong one.** `zn_cl2` was pre-registered as the
  member expected to get worse and it improves 9.1x. The prediction of a direction was wrong;
  the prediction that the set would find a counter-molecule was right.
- **H3 (the gate does not move): HOLDS as expected.** 0 PASS / 16 FAIL, both arms,
  15980 points. Nothing here agrees with NBO 7, and a 31.5x improvement is not an agreement.
- **The "fit rather than improvement" reading is not available.** It required a median
  held-out factor below 2x; the measured median is 18.0.

### The renat5 arm in detail

The same arm as the `ren5` columns above, with the field carrying each molecule's largest
deviation named - because on four members it is not the charge.

| molecule | d(Ryd) / e | worst L / e | max abs dq / e | worst other NPA field / e |
|---|---|---|---|---|
| `bf4_minus` | -0.00021 | -0.00005 | 0.00276 | 0.00277 (B val) |
| `bh3` | +0.00083 | +0.00041 | 0.00599 | 0.00641 (B val) |
| `c2h5` | -0.00013 | -0.00070 | 0.02552 | **0.18791 (C spin)** |
| `cf3` | -0.00607 | -0.00697 | 0.01439 | **0.11801 (C spin)** |
| `ch2` | -0.00131 | -0.00199 | 0.05140 | **0.29174 (C spin)** |
| `cl2` | -0.00026 | -0.00013 | 0.00001 | 0.00012 (Cl ryd) |
| `clf3` | -0.00232 | -0.00071 | 0.00020 | 0.00074 (F val) |
| `hi` | +0.00008 | +0.00011 | 0.00003 | 0.00832 (I val) |
| `hs` | +0.00044 | +0.00048 | 0.00896 | **0.06533 (S spin)** |
| `mg_h2` | -0.00477 | -0.00188 | 0.00510 | 0.00510 (Mg charge) |
| `nacl` | +0.00011 | +0.00007 | 0.00013 | 0.00020 (Na val) |
| `nh4_plus` | +0.00343 | +0.00342 | 0.00271 | 0.00615 (N val) |
| `si2h6` | +0.00107 | +0.00040 | 0.00863 | 0.00901 (Si val) |
| `sih4` | -0.00002 | +0.00171 | 0.01298 | 0.01467 (Si val) |
| `thiophene` | +0.00164 | -0.00155 | 0.00553 | 0.00553 (S charge) |
| `zn_cl2` | +0.00104 | +0.00070 | 0.00027 | 0.00097 (Zn val) |

Over the set's 68 atoms, renat5 arm: mean absolute intra-atomic Rydberg leak **0.00051 e**,
mean absolute charge deviation **0.00466 e** (baseline arm: 0.00841 e and 0.01264 e).

This table's `d(Ryd)` comes from `nao_class_leak.py`, which sums gennbo's printed per-class
totals; the two-arm table above uses `twoarm_table.py`, which sums the per-atom NPA Rydberg
column. The two agree to within **1.3 floors** on every member, and `sih4` is where that
matters: -0.00002 here against +0.00002 there, i.e. **a sign the instrument cannot resolve**,
2e-05 being 0.8 of `sih4`'s 2.5e-05 floor. The two-arm table uses one definition for both
arms, so its factors are internally consistent.

Two things this table says out loud:

1. **The NPA charge is a blind instrument for this defect, again.** `bf4_minus` leaks
   0.00005 e within an atom and shows 0.00276 e of charge deviation - a ratio of 54. The
   charge column does not see the cascade's error; it sees something else.
2. **The large deviations on the four open shells are in the spin field, not the
   charge.** 0.29 e on `ch2`, 0.19 e on `c2h5`, 0.12 e on `cf3`, 0.065 e on `hs`, while
   every charge deviation in the set stays at or below 0.05 e. An open-shell claim made
   on the spin-summed table would miss all four.

### Column 3, the per-class principal sine

From gennbo's own written matrices (`AOPNAO=W AONAO=W NRTSYM=OFF`, lfn 32 and 33) against
native's `NAO_DUMP_CPRE` / `NAO_DUMP_C` dumps, in the shared AO basis. Each molecule's
floor is `sqrt(2 max|C^T S C - 1|)` over both sides - a sine's floor is the **square
root** of the coefficient floor. A class counts as different only above 10x floor.

Both arms, one binary again (job 595960, `LEGACY=1` into `aopnao_stage.sh` for the baseline
dumps). That stage does **not** echo its environment, so the knob was verified by its effect
instead: the baseline dumps agree with job 595959's baseline arm at **0.000e+00** over 147 to
415 numeric fields and differ from the renat5 arm by 2.9 to 5.5 (an NAO energy) on `sih4`,
`thiophene` and `mg_h2`. `d_VR` is the number of valence dimensions the two sides disagree
about, the form comparable across molecules.

| molecule | floor | sin Cor base | sin Cor ren5 | sin Val base | sin Val ren5 | factor | d_VR base |
|---|---|---|---|---|---|---|---|
| `bf4_minus` | 6.9e-05 | 4.41e-05 | 4.41e-05 | 9.67e-02 | 1.10e-03 | 87.9 | 0.036 |
| `bh3` | 6.7e-05 | 0.00e+00 | 0.00e+00 | 1.65e-01 | 1.14e-02 | 14.5 | 0.085 |
| `cl2` | 6.5e-05 | 4.85e-05 | 4.85e-05 | 2.78e-02 | 2.42e-03 | 11.5 | 0.001 |
| `clf3` | 1.3e-04 | 4.38e-05 | 4.38e-05 | 4.63e-02 | 4.81e-03 | 9.6 | 0.006 |
| `hi` | 5.2e-05 | **1.00e+00** | **1.00e+00** | **1.00e+00** | **1.00e+00** | 1.0 | 0.005 |
| `mg_h2` | 5.8e-05 | 2.23e-05 | 2.23e-05 | 1.06e-01 | 1.26e-02 | 8.4 | 0.021 |
| `nacl` | 5.5e-05 | 4.16e-05 | 4.16e-05 | 1.98e-02 | 7.54e-04 | 26.3 | 0.001 |
| `nh4_plus` | 6.1e-05 | 0.00e+00 | 0.00e+00 | 1.33e-01 | 7.40e-03 | 18.0 | 0.079 |
| `si2h6` | 5.9e-05 | 3.88e-05 | 3.88e-05 | 1.97e-01 | 8.61e-03 | 22.9 | 0.194 |
| `sih4` | 5.8e-05 | 2.57e-05 | 2.57e-05 | 1.45e-01 | 9.77e-03 | 14.8 | 0.085 |
| `thiophene` | 6.9e-05 | 4.65e-05 | 4.65e-05 | 2.23e-01 | 2.38e-02 | 9.4 | 0.384 |
| `zn_cl2` | 6.4e-05 | 4.84e-05 | 4.84e-05 | 3.40e-02 | 2.90e-03 | 11.7 | 0.005 |

The core sine is **identical between the arms to every printed digit on 12 of 12** - the
change does not touch the core partition, which is what its own source comment claims. `sin
Val` improves on 11 of 12 (the twelfth is `hi`, pinned at the metric's ceiling in both arms),
factor 8.4 to 87.9, median 14.6, and **every renat5 value is still 10x its floor or more**, so
the valence/Rydberg spaces still differ after the change. Note that `cl2` improves 11.5x here
while its `d(Ryd)` got 2.4x worse: `d(Ryd)` is a signed sum that cancels across atoms and the
sine cannot cancel, so on `cl2` the two instruments disagree in sign and the subspace
instrument is the one that cannot hide a defect by cancellation.

12 of 16 measured, **4 VOID**: `c2h5`, `cf3`, `ch2` and `hs` are refused by the
comparison's own self-check, because native and gennbo assign different numbers of rows
per `(atom, l)` to each class, and a space called valence on one side is then not the same
space on the other. That is a real hole, not a small one: **the sine instrument as written
cannot measure an open shell**, so one third of the held-out set has no column 3 and the
four members chosen to stress spin are exactly the four it drops.

`sin Val` and `sin Ryd` are equal to 5 to 7 digits on every measured member. That is
forced by complementarity - Rydberg is valence's complement inside Val+Ryd - so it is
**one** measurement, not two. `hi`'s 1.00 is the saturation described above: identical
subspace dimensions holding different orbitals give a principal sine of exactly 1, which
is the metric's ceiling rather than a magnitude.

Classes called different at 10x floor: Cor 1/12, Val 12/12, Ryd 12/12.

## RETRACTED: "the 0.3161 e figure does not belong to this branch"

That section claimed the pre-registration's 0.3161 e before-value was a property of the
`density_source` share build `c8130055` and not of this branch, on the strength of benzene
reading 0.12055 with what was labelled the baseline binary. **The claim is withdrawn.** The
binary that produced 0.12055 was the renat5 build (md5 `1e45bc35...`, see the correction
above), so 0.12055 is the *treatment* value and it never contradicted anything: the sibling's
independently stored arms read legacy 0.43159, pre-change binary 0.43159, renat5 0.12055,
gennbo 0.11544. My number reproduces their treatment arm to five decimals.

What is withdrawn with it:

- The **cross-build mechanism hypothesis** - that the difference lives in `nbo.cpp`, `nrt.cpp`,
  `wfn_density.cpp`, `properties.cpp`, `isosurface.cpp` or `esp_gpu.cu` on the
  `density_source` side. There is no unexplained cross-build difference to explain; the whole
  gap is one `nao.cpp` difference, the change itself.
- The open hole "the mechanism is not named". It is named: renat5.
- "**No arm of 60055d4b reproduces 0.43157**" from job 594823. All four arms of that job ran
  the renat5 binary and none of them set `NAO_LEGACY_CASCADE=1`, the only switch that turns
  the new Rydberg step off, so the legacy number was not reachable in that job by
  construction. The job's positive control still stands for what it actually tested -
  `NAO_OWSO_OFF=1` and `NAO_CLASS_SPLIT=1` move benzene's total Rydberg (0.28604 and 1.12682
  against 0.12055) so environment plumbing arrives at the binary, and `NAO_CORE_POOLED=1` is
  flat at 0.12054 - but it does **not** bracket the baseline and must not be quoted as doing
  so.

What survives, because it was a git fact rather than an inference from the binary:
`Src/core/nao.cpp` is the byte-identical blob `022b9c21` at `c8130055` and at the two
branches' merge base `8d1b43b4`, and of the six nao.cpp commits since, only `2114ee18`
changes default execution. That is still true and still worth knowing. The inference drawn
from it was wrong.

**What went wrong, as a rule.** Two sessions shared one output path,
`nos_nboref/build/release-linux/bin/NoSpherA2`, and I copied it five minutes after the
sibling had written its own build there. Nothing in my instrumentation could see that, because
the only provenance I carried was a `DEPLOYED_COMMIT.txt` I had typed myself, and I printed
that file in job output as though it were a measurement. A commit label is not a checksum. The
cheap guard is the one now in `twoarm.sh`: `md5sum` the file from inside the job, and check
that the two arms actually disagree on one molecule before running the other fifteen.

## One deviation, one closed invariance, one hole

- **The thread count does not move a number. Measured, closed** (job 594824, `xcheck.sh` at
  `threads=8` against the stored `threads=1` files for `si2h6`, `thiophene` and `c2h5`). Each
  pair of `<mol>.native.nbo.json` differs in exactly **4 lines**, and all four are wall-clock:
  the `timings` object and `nrt.seconds`. Every NPA, NAO, NBO, E2 and NRT numeric field is
  identical to the last printed digit - the files are one line per record, so a line-level
  diff is exhaustive here rather than indicative. Stage 2 ran at `threads=8` for 13 of the 16
  and `threads=1` for the three added later, and that difference is now measured to be
  nothing. `holdout_stage.sh`'s header claim holds.
- The four open shells have no column 3 at all, as above: the sine comparator refuses them.
- **CLOSED, 25 Sep:** `cf3`'s regression is localised to **one spin**, and no new reference
  was needed to do it - see the next section. The sine instrument still refuses open shells;
  what was missing was not a reference but the right key in the file we already had.



## `cf3` is localised to the alpha spin, and `cl2` is not a robust counter-molecule

Two things were wrong in how the two counter-molecules were first reported, and both are
fixed by instruments that need no cluster and no new reference. Run them with
`py -3.12 spinsplit.py data_holdout` and `py -3.12 twofloors.py data_holdout <mols>`;
`spinsplit.py --selfcheck data_holdout` is the assert-only version.

### The per-spin reference was already in git

The open-shell gap was scoped as a missing reference. It is not: gennbo's own **per-spin**
NAO tables are already parsed into **`nao_alpha` and `nao_beta`** in every
`<mol>.gennbo.nbo.json`, and `-nbo_native` emits the same two keys, matched row for row
(`cf3` 124/124/124, `ch2` and `hs` 43, `c2h5` 92). The composite `nao` key is the spin-summed
table, which is why a comparison against *it* is meaningless for one spin - `alpha + beta`
reproduces the composite occupancy column to within **2.0 floors** on all four molecules.
Nothing had to be rebuilt or re-run; the comparator was reading the wrong key.

And a hypothesis about that gap is **refuted**: gennbo and native agree on the Cor/Val/Ryd
class of **every basis function** of all four open shells, in both arms - `misclass = 0`,
124/124 and 43/43 and 92/92. Whatever the per-`(atom, l)` row-count difference is, it is not
spin-dependent classification on this set.

### `cf3`: alpha is 1.9x worse, beta is 8.0x better, and the baseline's total was cancellation

Per-spin Rydberg population, native - gennbo, each side classed by itself. Floor here is
`n_Ryd x 0.5e-5` = 5.2e-04 for `cf3`, because this sum adds 104 printed values, not 4.

| molecule | spin | baseline | renat5 | factor | floors (base / ren5) |
|---|---|---|---|---|---|
| cf3 | alpha | -0.00285 | **-0.00541** | **0.53 worse** | 5.5 / 10.4 |
| cf3 | beta | +0.00580 | -0.00073 | 7.96 better | 11.1 / 1.4 |
| cf3 | total | +0.00302 | -0.00607 | **0.50 worse** | 5.8 / 11.7 |
| c2h5 | alpha | +0.05318 | +0.00029 | 186.4 | 138.1 / 0.7 |
| c2h5 | beta | +0.05571 | -0.00033 | 170.9 | 144.7 / 0.8 |
| ch2 | alpha | +0.00340 | -0.00042 | 8.1 | 18.9 / 2.3 |
| ch2 | beta | +0.01431 | -0.00089 | 16.0 | 79.5 / 5.0 |
| hs | alpha | +0.00194 | +0.00035 | 5.6 | 11.7 / 2.1 |
| hs | beta | +0.00203 | +0.00010 | 20.4 | 12.3 / 0.6 |

`cf3` is the only one of the four that gets worse, and only in **alpha**. The mechanism is
visible in the signs: the baseline's two spins carried errors of **opposite** sign
(-0.00285 and +0.00580) which partly cancelled in the total (+0.00302). `renat5` drives both
**negative** (-0.00541, -0.00073), so they add instead of cancelling (-0.00607). So part of
"2.08x worse on cf3" is the loss of a cancellation the baseline was benefiting from, not 2x
more error everywhere - but **alpha genuinely worsens 1.9x at 10.4 floors**, so this explains
the size of the regression without excusing it. `c2h5` is the control: same open-shell
structure, both spins improve by >170x, so the alpha channel is not generically harmed.

### `cl2` flips verdict between two legitimate summations

The same `d(Ryd)` can be summed over the per-atom NPA Rydberg column (k = N_atoms) or over the
per-NAO `Ryd` rows (k = n_Ryd). Both are the same physical quantity - they agree in value to
well within the looser floor on all 16, which `twofloors.py` asserts - but the floor is
`k x 0.5e-5`, so the per-NAO sum is the **looser** instrument.

| molecule | convention | k | floor | baseline | renat5 | verdict |
|---|---|---|---|---|---|---|
| cl2 | per-atom | 2 | 1.0e-05 | -0.00010 | -0.00024 | worse (10 -> 24 floors) |
| cl2 | per-NAO | 56 | 2.8e-04 | -0.00012 | -0.00026 | **both below floor: no signal** |
| cf3 | per-atom | 4 | 2.0e-05 | +0.00295 | -0.00614 | worse (148 -> 307) |
| cf3 | per-NAO | 104 | 5.2e-04 | +0.00302 | -0.00607 | worse (5.8 -> 11.7) |

`cl2` is the **only** molecule of the 16 whose verdict depends on the summation. So the honest
headline is **one robust counter-molecule, not two**: `cf3` fails on both instruments and by a
wide margin, while `cl2`'s regression exists on the tighter convention and vanishes on the
looser one, at a magnitude of 1.4e-04 e. Quote the convention with the number or do not quote
the number.

### Trap: `nao[i].energy` on an open shell is the ALPHA energy

The open-shell composite NAO table prints `Occupancy` and `Spin` - it has **no Energy column**.
The `energy` on a composite record in the reference JSON is therefore taken from elsewhere, and
it is measurably gennbo's **alpha** energy: equal on **124/124** rows of `cf3`, 43/43, 43/43,
92/92, and to the beta energy on **0** rows of any of them. `spinsplit.py --selfcheck` asserts
this so it cannot drift silently. It is the same defect family as the `[:packed]` slice in
`read_47` (`53f11d7f`) and the label beside the binary - an alpha quantity wearing a total
label - and `nboref_sidebyside.py` prints that field in a column headed `E g`, so read that
column as alpha on any open-shell molecule. No measured number in this README depends on it:
every column here is built from occupancies and charges.

### Trap: a blank line inside gennbo's NAO table separates atoms

It does not end the table. A reader that stops at the first blank line after the header returns
**the first atom only** - 31 of `ch3`'s 49 rows, which looks like a plausible table and is not
one. Break on a blank line only when the next non-blank line is not itself a row.
