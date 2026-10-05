# The held-out set: 16 molecules the NAO change was never tuned on

The set covers chemistry that the 22 of `../data/` do not: cations, anions, an ECP, d10
metals, radicals and a triplet. None of these molecules is in
`tests/nbo_reference/index.json`.

## The 16

| molecule | charge, mult | what it is in the set for |
|---|---|---|
| `sih4` | 0, 1 | second-row closed shell, the plain silicon reference for `si2h6` |
| `si2h6` | 0, 1 | ethane one row down |
| `cl2` | 0, 1 | third-row homonuclear, no polarity |
| `clf3` | 0, 1 | **hypervalent beyond sf6/pf5**, and T-shaped rather than symmetric |
| `thiophene` | 0, 1 | aromatic ring with a third-row heteroatom; benzene's homologue |
| `bh3` | 0, 1 | electron-poor, an empty p on boron, a **diffuse and well-populated Rydberg set** |
| `mg_h2` | 0, 1 | s-block metal hydride, diffuse valence |
| `nacl` | 0, 1 | ionic closed shell, near-full charge transfer |
| `nh4_plus` | +1, 1 | **cation** |
| `bf4_minus` | -1, 1 | **anion**, and four equivalent ligands |
| `zn_cl2` | 0, 1 | **transition metal beyond ni_co_4/ticl4**, d10, a full d shell in the valence set |
| `hi` | 0, 1 | **ECP-bearing heavy atom** (def2-ECP on iodine, 28 core electrons replaced) |
| `hs` | 0, 2 | **open-shell radical**, third row |
| `cf3` | 0, 2 | **open-shell radical**, pyramidal, three fluorines |
| `c2h5` | 0, 2 | open-shell hydrocarbon |
| `ch2` | 0, 3 | triplet, two unpaired electrons rather than one |

Every def2-TZVP member has f sets on its second- and third-row atoms, so each one exercises
the component-major to shell-major permutation of d/f functions.

`hi` runs end to end: 26 electrons, and `-convert_to_47`, gennbo, `-nbo_parse` and
`-nbo_native` all return 0. Under an ECP, gennbo labels iodine's first NAO `Cor( 1s)`,
because NBO names a shell by its order within `l` and not by its principal quantum number.
Do not read anything off a printed shell label.

## What each molecule directory holds

The layout is the same as `../data/`:

| file | stage | what it is |
|---|---|---|
| `<mol>.inp` | 1 | the ORCA input that was actually run |
| `<mol>_start.xyz` | 1 | the idealised geometry `make_holdout.py` writes |
| `<mol>.xyz` | 1 | the optimised geometry |
| `provenance_orca.json` | 1 | hostname, threads, node cores, slurm job id, ORCA path |
| `<mol>.gennbo.nbo.json` | 2 | gennbo 7's output parsed by `-nbo_parse` |
| `provenance_nbo.json` | 2 | hostname, threads, NoSpherA2 commit, gennbo path, keylist, flags, per-step seconds |

`accepted_holdout.json` marks all 16 `accepted: true`. Its criterion is weaker than
`../data/`'s, because no stored basis-function count exists to match against. The checks
are:

- a converged optimisation and a converged last SCF;
- the charge, multiplicity, basis and keyword line as requested;
- `<S**2>` within 4 % of `S(S+1)`.

The basis-function count and the final energy are recorded, so the next regeneration can
apply the exact test.

`.gbw`, `.47`, `.out` and `.nbo` stay on the cluster at
`/work/akkleemiss/florian/nbo_ref_v2_holdout/<mol>/`.

## Reproducing it

The scripts are one level up. Run them on the cluster (`ssh AKL007`), with
`$H = /work/akkleemiss/florian/nbo_ref_v2_holdout`.

```
python make_holdout.py $H                     # 16 x <mol>/<mol>.inp + <mol>_start.xyz
sbatch holdout_stage.sh <mol>                 # both stages, one job per molecule
python accept_holdout.py $H                   # -> accepted_holdout.json
python compare_all.py $H                      # native vs gennbo
```
