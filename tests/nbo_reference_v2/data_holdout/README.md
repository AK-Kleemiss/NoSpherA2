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
| `<mol>.native.nbo.json` | 2 | NoSpherA2's `-nbo_native` answer on the same wavefunction |
| `provenance_nbo.json` | 2 | hostname, threads, NoSpherA2 commit, gennbo path, keylist, flags, per-step seconds |

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
python nao_class_leak.py $H                   # column 1 and 2 of the table below
SRC=$H V2=$H/bin BIN=$H/bin/NoSpherA2_baseline OUT=$H/aopnao \
    MOLS="bf4_minus bh3 c2h5 cf3 ch2 cl2 clf3 hi hs mg_h2 nacl nh4_plus si2h6 sih4 thiophene zn_cl2" \
    sbatch --export=ALL bin/aopnao_stage.sh    # gennbo's lfn 32/33 + native's two dumps
python aonao_compare.py --space $H/aopnao     # column 3
```

Two traps, both paid for once already:

- The login node's `python3` has **no numpy**, so `aonao_compare.py` dies at its import.
  Use `/work/akkleemiss/florian/fp_fdp/bin/python3` (numpy 2.0.2, Python 3.9.25).
- `aopnao_stage.sh` writes `<mol>.aopnao.nbo.json` and `aonao_compare.py --space` reads
  `<mol>.aonao.nbo.json`. Copy one to the other before the comparison.

## The baseline column: binary 60055d4b, before the NAO change

Every number here is `-nbo_native` against gennbo 7 on the same wavefunction, measured
with the binary whose `DEPLOYED_COMMIT.txt` reads **60055d4b** ("NAO: the step-4-operator
hypothesis is refuted"), the head of `nbo_external_reference` before the change landed.

`d(Ryd)` = total Rydberg population, native minus gennbo. `worst L` = the largest
per-atom intra-atomic Rydberg leak. `max |dq|` = the largest per-atom NPA **charge**
deviation. Both population columns are printed by gennbo to five decimals, so their floor
is **1e-5**; a figure at 1e-5 is at the floor and is not a measurement of a magnitude.

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

Over the set's 68 atoms: mean absolute intra-atomic Rydberg leak **0.00051 e**, mean
absolute charge deviation **0.00466 e**.

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

| molecule | floor | sin Cor | sin Val | sin Ryd |
|---|---|---|---|---|
| `bf4_minus` | 6.9e-05 | 4.41e-05 | 1.10e-03 | 1.10e-03 |
| `bh3` | 6.7e-05 | 0.00e+00 | 1.14e-02 | 1.14e-02 |
| `cl2` | 6.5e-05 | 4.85e-05 | 2.42e-03 | 2.42e-03 |
| `clf3` | 1.3e-04 | 4.38e-05 | 4.81e-03 | 4.81e-03 |
| `hi` | 5.2e-05 | **1.00e+00** | **1.00e+00** | 2.27e-02 |
| `mg_h2` | 5.8e-05 | 2.23e-05 | 1.26e-02 | 1.26e-02 |
| `nacl` | 5.5e-05 | 4.16e-05 | 7.54e-04 | 7.55e-04 |
| `nh4_plus` | 6.1e-05 | 0.00e+00 | 7.40e-03 | 7.40e-03 |
| `si2h6` | 5.9e-05 | 3.88e-05 | 8.61e-03 | 8.61e-03 |
| `sih4` | 5.8e-05 | 2.57e-05 | 9.77e-03 | 9.77e-03 |
| `thiophene` | 6.9e-05 | 4.65e-05 | 2.38e-02 | 2.38e-02 |
| `zn_cl2` | 6.4e-05 | 4.84e-05 | 2.90e-03 | 2.90e-03 |

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

## The 0.3161 e baseline figure does not belong to this branch

The pre-registration quotes "the shipped native's worst is 0.3161 e (benzene: 0.4316
against gennbo's printed 0.1155)" and a factor 63 against `renat5`. That figure comes
from the *stored* column in `../data/`, whose `provenance_nbo.json` names
`/work/akkleemiss/share/NoSpherA2_RGBI_NBO/NoSpherA2`, `c8130055 built 2026-09-25 02:20
on AKL007, branch density_source`.

Re-running **only the native side** of all 22 with the baseline binary 60055d4b, reusing
each molecule's stored `.gbw`, `_native.47`, flags and thread count so that the single
thing that changes is which NoSpherA2 ran (`bin/xcheck.sh`, job 594815):

| molecule | stored column (c8130055) | baseline binary (60055d4b) | gennbo |
|---|---|---|---|
| benzene, total Rydberg | 0.43157 | 0.12055 | 0.11544 |
| benzene, d(Ryd) | +0.31613 | **+0.00511** | - |
| ethane, d(Ryd) | +0.12985 | +0.00038 | - |
| lif, d(Ryd) | +0.00119 | +0.00001 | - |
| mean over the 22 | 0.01229 | 0.00030 | - |

`whichryd.py` localises it inside native's own output rather than in a parser: per-carbon
Rydberg 0.06826 (stored) against 0.01861 (baseline) against gennbo's 0.01787, with the
valence column moving the other way by the same amount, and both files
`parser_version: 2`.

What has been **excluded** by measurement, with a positive control (job 594823, four arms
of the one 60055d4b binary, each arm printing the environment the process actually saw -
because a silent knob reads exactly like a knob with no effect):

| arm | benzene total Ryd | ethane total Ryd |
|---|---|---|
| default | 0.12055 | 0.02524 |
| `NAO_CORE_POOLED=1` (restores the pre-`2114ee18` pooled block) | 0.12054 | 0.02523 |
| `NAO_OWSO_OFF=1` | 0.28604 | 0.07738 |
| `NAO_CLASS_SPLIT=1` | 1.12682 | 0.43210 |

The last two move the number, so the environment plumbing demonstrably works and
`NAO_CORE_POOLED` really is flat - the core-block commit is not the cause, which is what
its own source comment claims. **No arm of 60055d4b reproduces 0.43157.**

And `Src/core/nao.cpp` is the byte-identical blob `022b9c21` at `c8130055` and at the two
branches' merge base `8d1b43b4`; the six nao.cpp commits since then are on the
`nbo_external_reference` side only, and of those only `2114ee18` changes default
execution. So this is **not** a nao.cpp difference at all. The two builds diverge in
`nbo.cpp`, `nrt.cpp`, `wfn_density.cpp`, `properties.cpp`, `isosurface.cpp` and
`esp_gpu.cu` on the `density_source` side.

**Left open.** The mechanism is not named. What follows from the measurement and no more:
a 63x ratio taken against 0.3161 e is a ratio against that share build, and on the head of
this branch the same molecule's baseline is 0.00503 e. Any before/after on the held-out
set has to be taken with **one** binary and an environment switch, not with two builds.

## Two deviations and one untested invariance

- Stage 2 ran at `threads=8` for 13 of the 16 and at `threads=1` for `si2h6`,
  `thiophene` and `c2h5` (the three added later). Their `provenance_nbo.json` records it.
- `holdout_stage.sh`'s header claims a thread count must not move a number. That claim is
  **unverified** for this set. Until it is measured, the three `threads=1` members carry
  an uncontrolled difference from the other 13.
- The four open shells have no column 3 at all, as above.
