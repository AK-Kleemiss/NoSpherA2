# The 22 regenerated wavefunctions, and everything needed to build them again

This directory is the answer to the reason `nbo_reference_v2` exists at all. The
22 JSONs in `tests/nbo_reference/` are numbers with no inputs: the ORCA jobs and
the wavefunctions that produced them are gone, and the one same-named survivor
turned out to be a *different molecule* (its `water` is OHHHe, 73 basis
functions, where the reference records OHH with 43). Nobody can re-run them,
nobody can tell a regression from a fix, and nobody noticed the swapped molecule
because only the answers were kept.

So this time the inputs are kept, and this file says how to get from them back to
every number in here.

Nothing in this directory replaces `tests/nbo_reference/`. Those 22 JSONs stay
exactly where they are.

## What each molecule directory holds

| file | stage | what it is |
|---|---|---|
| `<mol>.inp` | 1 | the ORCA input that was actually run (byte-identical to `../inputs/<mol>/`, verified) |
| `<mol>_start.xyz` | 1 | the idealised geometry the optimisation started from |
| `<mol>.out` | 1 | ORCA's output, including the basis-function count and final energy acceptance reads |
| `<mol>.xyz`, `<mol>_trj.xyz` | 1 | the optimised geometry and the optimisation trajectory |
| `provenance_orca.json` | 1 | hostname, threads, node core count, slurm job id, ORCA path and module |
| `<mol>.nbo` | 2 | **gennbo 7's own output** - the external reference in its original text form |
| `<mol>.gennbo.nbo.json` | 2 | that output parsed by `-nbo_parse` |
| `<mol>.native.nbo.json` | 2 | NoSpherA2's own `-nbo_native` answer on the same wavefunction |
| `provenance_nbo.json` | 2 | hostname, threads, NoSpherA2 commit, gennbo path, the exact keylist and flags, per-step seconds |

`accepted.json` is the acceptance record for all 22 at once.

## Deliberately not committed

Two binary artefacts are left on the cluster at
`/work/akkleemiss/florian/nbo_ref_v2/<mol>/`:

- `<mol>.gbw` - 27 MB over the 22. ORCA's wavefunction.
- `<mol>.47` and `<mol>_native.47` - 42 MB over the 22. The NBO archive.

Both are *derived*, and the derivation is one command each (below), so what they
would add to the repository is 69 MB of things a reader can rebuild. What is kept
is everything that is **not** rebuildable without them: the input, the output,
gennbo's text, and both JSONs. That split is the judgement call in this
directory; if a future reader would rather have the `.47` in git than trust this
paragraph, it is 42 MB and the inputs are all here.

Whether the `.47` set belongs in the repository after all is Florian's call, not
this branch's.

## Regenerating any of it

Scripts are one level up, in `tests/nbo_reference_v2/`. All compute is on the
cluster (`ssh AKL007`); `env.sh` carries the module and binary paths, including
two traps worth reading before you run anything - the `gennbo` wrapper everyone
points at is broken on this cluster, and the ORCA module is not on the default
`MODULEPATH`.

```
python make_inputs.py <out_dir>        # 22 x <mol>/<mol>.inp + <mol>_start.xyz
sbatch orca_opt.sh <mol>              # stage 1 -> <mol>.out, <mol>.gbw, provenance_orca.json
python accept.py <root>                # -> accepted.json; refuses a molecule that drifted
sbatch nbo_stage.sh <mol>             # stage 2 -> .47, .nbo, both JSONs, provenance_nbo.json
python compare_all.py <root>          # native vs gennbo, all 22
```

`nbo_all_stage.sh` runs stage 2 over all 22 in one job. To rebuild only the two
excluded binaries for one molecule:

```
sbatch orca_opt.sh <mol>                                   # the .gbw
$NOSPHERA2 -convert_to_47 <mol>.gbw -nbo_keywords "$(python config.py --nbo-keywords <mol>)"
```

Stage 2 keeps `<mol>_native.47` as a copy taken *before* the keylist edit,
because `so2` gets a hand-written `$NRTSTR` appended for gennbo and the native
reader has no business seeing it. `config.py` is the single place that knows the
per-molecule deviations: `so2`'s `$NRTSTR`, `sf6`'s `NRTSYM=off`, and the two
molecules that need a raised `NRTE2`.

## Why these wavefunctions may stand in for the lost ones

A regenerated wavefunction is only a stand-in if it is the same molecule at the
same level of theory. The test is `accept.py`, and it is deliberately blunt:

1. **`basis_functions` must match `index.json` exactly.** Not within a
   tolerance - exactly. One basis function of difference means a different
   molecule, which is precisely how the `water` / OHHHe mixup would have been
   caught the first time.
2. **`final_energy_hartree` within 2.0e-5 Eh** of the recorded value.
3. Charge, multiplicity, basis and keyword line as recorded; optimisation
   converged.

All 22 passed all four. The worst energy deviation is ammonia at **9.03e-06 Eh**,
half the tolerance; the median is 1.69e-06 Eh.

### Where the 2.0e-5 Eh tolerance comes from

It is not a guess, and it is not fitted to the results. `ch3` is the calibration
point: `tests/nbo_reference/ch3_orca_reference.47` is the one archive that
survived, so `make_inputs.py` starts `ch3` from the *optimised* geometry of the
original run, read out of that archive's `$COORD` block. Its regenerated energy
therefore measures the pipeline against the original with the geometry held
fixed, and it comes out at **1.31e-08 Eh**. That is the floor: same ORCA version,
same basis, same grid defaults reproduce the recorded energy to eight decimals.
The other 21 add a re-optimisation on top of that floor, so the tolerance has to
cover how far a converged minimum can wander - three orders of magnitude above
the floor does, and ammonia's 9.03e-06 sits comfortably inside it while still
being far too tight for any of the failure modes worth catching (a different
molecule, a different functional, a different basis, an unconverged optimisation
all miss by 1e-3 or more).

`accepted.json` still records `tolerance_hartree` 1.0e-4 for `ni_co_4`,
`nitromethane` and `ticl4`, because it was written before those three
per-molecule loosenings were removed from `accept.py`. Their measured deviations
are 2.94e-08, 1.83e-06 and 2.89e-06 - all inside 2.0e-5 - so removing the
loosening does not change any verdict in this file. An unexercised tolerance
cannot be told apart from an unnecessary one, which is why it is gone.

## Acceptance record

| molecule | chg | mult | basis fns | recorded E / Eh | dev / Eh | host | wall / s |
|---|---|---|---|---|---|---|---|
| acetylene | 0 | 1 | 74 | -77.313388 | 6.55e-06 | AKL064 | 19 |
| ammonia | 0 | 1 | 49 | -56.549318 | 9.03e-06 | AKL064 | 16 |
| benzene | 0 | 1 | 222 | -232.183320 | 5.10e-06 | AKL064 | 119 |
| ch3 | 0 | 2 | 49 | -39.826588 | 1.31e-08 | AKL064 | 16 |
| ethane | 0 | 1 | 98 | -79.798965 | 2.58e-06 | AKL064 | 141 |
| ethene | 0 | 1 | 86 | -78.565717 | 1.39e-07 | AKL064 | 34 |
| formaldehyde | 0 | 1 | 74 | -114.495086 | 1.55e-06 | AKL064 | 33 |
| formate | -1 | 1 | 99 | -189.197679 | 4.96e-06 | AKL064 | 60 |
| hcn | 0 | 1 | 68 | -93.412610 | 3.72e-06 | AKL064 | 24 |
| lif | 0 | 1 | 45 | -107.430642 | 4.99e-07 | AKL064 | 26 |
| n2 | 0 | 1 | 62 | -109.521128 | 6.23e-07 | AKL064 | 21 |
| ni_co_4 | 0 | 1 | 293 | -1961.652333 | 2.94e-08 | AKL064 | 375 |
| nitromethane | 0 | 1 | 142 | -245.001629 | 1.83e-06 | AKL064 | 129 |
| no | 0 | 2 | 62 | -129.891252 | 2.32e-06 | AKL064 | 24 |
| o2 | 0 | 3 | 62 | -150.329914 | 2.00e-06 | AKL064 | 23 |
| ozone | 0 | 1 | 93 | -225.421928 | 5.94e-07 | AKL064 | 41 |
| pf5 | 0 | 1 | 192 | -840.760146 | 7.11e-08 | AKL064 | 118 |
| pyridine | 0 | 1 | 216 | -248.224443 | 4.76e-06 | AKL064 | 244 |
| sf6 | 0 | 1 | 223 | -997.213510 | 6.98e-08 | AKL064 | 123 |
| so2 | 0 | 1 | 99 | -548.600329 | 1.14e-06 | AKL064 | 37 |
| ticl4 | 0 | 1 | 193 | -2690.326773 | 2.89e-06 | AKL006 | 37 |
| water | 0 | 1 | 43 | -76.426059 | 3.90e-07 | AKL064 | 30 |

`recorded E` is `index.json`'s value, i.e. the lost run's; `dev` is against the
regenerated one. Every row is `! B3LYP def2-TZVP TightSCF Opt`, ORCA 6.1.1
(`orca/avx2-6.1.1-xqm3jnz`), 4 MPI ranks, and every row optimised to convergence.
`ticl4` ran on AKL006 rather than AKL064 - both 96-core nodes, and both jobs got
4 cores; that difference is recorded rather than smoothed over because a sibling
branch once chased a phantom 6.5x regression that was a 4-thread job on a 96-core
node with nothing in the log saying so.

## Stage 2 provenance, in one place

All 22 ran in slurm job **581339** on **AKL028** (96-core node, 8 cores given,
`-nbo_threads 8`), 198 s of wall time over the whole set, no failures.
NoSpherA2 was `/work/akkleemiss/share/NoSpherA2_RGBI_NBO/NoSpherA2`, commit
**c8130055**, built 2026-09-25 02:20 on AKL007 from `density_source`; gennbo 7
from `/work/software/bin/NBO/bin`. The keylist is `NRT E2PERT=<value> NRTLST=0.1
NRTDTL` plus the per-molecule additions above, and the native side gets the
matching `-nbo_e2min` - `E2PERT` **must** carry its value, because `E2PERT` alone
is accepted and silently ignored, exit 0, no warning, which is how the surviving
keylist came to set nothing at all.

The 22 `.nbo` files in here were parsed by commit c8130055, which does **not**
contain the valency-spin fix on this branch (`bf73d301`: NBO prints three NRT
valency tables for an open shell and the composite one was stored as a second
beta table). The `.nbo` text is gennbo's own and is unaffected; only the
`.gennbo.nbo.json` files of `ch3`, `no` and `o2` carry the duplicate. Re-running
`-nbo_parse` on the stored `.nbo` files with a rebuilt binary fixes those three
JSONs without touching the cluster or ORCA - which is the point of having kept
the `.nbo` text.
