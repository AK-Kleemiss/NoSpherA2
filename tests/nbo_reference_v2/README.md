# nbo_reference_v2 — an NBO reference that can actually arbitrate

`tests/nbo_reference/` holds 22 gennbo 7 result JSONs and nothing that produced them. This
directory holds the missing half: the geometries, the ORCA inputs and outputs, the `.47`
archives and gennbo 7's own output, so the comparison can be re-run and re-audited rather
than believed.

The 22 old JSONs stay exactly where they are. They are still a useful regression gate on the
parser and on our own past output; they are only useless as an *external* reference.

## Why the old set cannot arbitrate

1. **Its wavefunctions and ORCA inputs are gone.** The directory ships answers and no
   geometries. The only two molecules whose same-named wavefunction still exists anywhere in
   the test tree are *different molecules*: `water` there is OHHHe with 73 NAOs, against the
   reference's OHH with 43.
2. **Its NRT weights are wrong.** The parser read two columns by position: NRT's `RS` is a
   *rank*, not a structure number (weight sums came out 112 / 145 / 200 %), and an open-shell
   NAO table prints `Spin` where the parser read `Energy` (42 of ch3's 49 energies landed at
   zero). Both bugs are fixed in the code as of `b4584267`, but the 22 stored JSONs still
   carry the broken numbers. Do not treat a stored NRT weight or open-shell NAO energy as
   truth.
3. **`compare_nbo.py --all .` run inside `tests/nbo_reference` compares the dataset with
   itself** and reports 22/22 at 0.00000. That is not a gate. A comparison needs two sides
   with independent provenance.

What survived, and what makes this fixable, is `tests/nbo_reference/index.json`: it records
per molecule the ORCA `keyword_line` (uniformly `! B3LYP def2-TZVP TightSCF Opt`), basis,
charge, multiplicity, **`basis_functions`**, **`final_energy_hartree`** and whether the
optimisation converged — read back from what the programs printed. That is enough to rebuild
a wavefunction and *prove* it is the same molecule and state.

## The acceptance gate (`accept.py`) — the crux

A regenerated wavefunction is accepted only if all of these hold against `index.json`:

- the optimisation converged;
- charge, multiplicity, HF type, basis and ORCA version match;
- **`basis_functions` matches exactly** — this is the identity test. A count that differs by
  even one means a different molecule was built. It is exactly how the `water` / OHHHe mixup
  would have been caught;
- `|ΔE|` against the recorded `final_energy_hartree` is inside tolerance.

Tolerance is `2.0e-5` Ha by default — four times ORCA's own `5e-6` Eh per-step convergence
criterion, so it accepts the same minimum reached by a different path and rejects a different
minimum. Three molecules are given `1.0e-4` Ha and the reason is recorded in `accept.py`:
`ni_co_4` and `ticl4` for soft metal–ligand bending modes, `nitromethane` for an almost free
methyl rotor. **`accept.py` prints both basis-function counts and the actual deviation for
every molecule, pass or fail** — a column of "PASS" cannot be re-audited.

`ch3` is the calibration point. `tests/nbo_reference/ch3_orca_reference.47` is the one
surviving archive and its `$COORD` is the *optimised* geometry, so ch3 is the only molecule
that starts from the original structure. It reproduces the recorded energy to 1.3e-8 Ha,
which is what says the pipeline itself is faithful and the 1e-6-scale deviations elsewhere
are starting geometries relaxing into the same minimum.

## What the comparison actually compares

Stage 2 runs both analyses **on the same wavefunction and the same `.47` archive**:

- the **external** side: `NoSpherA2 -convert_to_47` writes the archive, **gennbo 7** analyses
  it, and `NoSpherA2 -nbo_parse` turns gennbo's output into JSON with the *fixed* parser;
- the **native** side: `NoSpherA2 -nbo_native` on the same archive.

`compare_all.py` then reports NPA charges, NAO occupancies and energies, NBO occupancies and
hybrid composition, E2 entries, and NRT bond orders, valencies and weights. gennbo is the
reference side, so a quantity gennbo printed and we did not is a `missing`, not a silent
pass. **NRT weights are read off the by-rank comparison, never by label**: starting the
search from a different Lewis structure renumbers the structure list and swaps the
descriptions of equivalent structures without moving a single weight.

Both sides are put on the thresholds NBO itself reported (`index.json`
`thresholds_reported_by_nbo`): `config.py` translates them into `-nbo_e2min` and `-nrt_e2` for
the native run, so the comparison measures the method and not a printing cut-off.

## Running it

Everything runs on the cluster (`ssh AKL007`); a local box is capped at ~30 % CPU.

```
python make_inputs.py <out_dir> --nprocs 4     # 22 <mol>/<mol>.inp + <mol>_start.xyz
# upload <out_dir> and bin/ to $ROOT on the cluster, then per molecule:
sbatch bin/orca_opt.sh  <mol>                  # stage 1: optimise, write <mol>.gbw
python bin/accept.py <root> --json accepted.json
sbatch bin/nbo_stage.sh <mol>                  # stage 2: .47, gennbo 7, -nbo_parse, -nbo_native
python bin/compare_all.py <root> --json comparison.json
```

`env.sh` carries the environment and is sourced by both batch scripts. Each stage writes a
`provenance_*.json` next to the molecule with the **hostname, the threads the job was given
and the core count of the node it landed on** — a number without its thread count and node is
not a measurement anyone can reuse, and an audit here once chased a 6.5x "regression" that
was a 4-thread job landing on a 96-core node.

## Traps, all of them paid for

- **`module load orca/...` reports "unknown"** until `module use
  /work/software/spack/share/spack/lmod/linux-rocky9-x86_64/Core` — the ORCA modulefile is in
  the spack lmod tree, which is not on the default MODULEPATH. And do not pipe `module load`
  into `head`: it runs in a subshell and its `setenv` is lost.
- **`/work/akkleemiss/share/gennbo` is broken on this cluster.** It sets
  `NBOBIN=/opt/nbo7/bin`, which exists on no node; the working install is
  `/work/software/bin/NBO/bin`. The shipped wrapper also runs `rm -i` (which blocks a batch
  job on a prompt) and hardcodes `NBOMEM=100gb`. `gennbo7` here is a bash replacement without
  either.
- **`BASH_SOURCE` does not locate a batch script's directory**: slurm copies it into
  `/var/spool/slurmd/job<id>/`. Both scripts take the path from `$NBOREF_BIN`.
- **Cancel a running job before resubmitting the same name**, or the submit helper hands back
  the old job's id — that has cost this project a wave already.
- **`--time=12:00:00` sat pending behind the Active queue; `--time=02:00:00` backfills
  immediately.**
- **`sf6` needs `NRTSYM=off`** or NBO aborts the whole analysis after 0.08 s and prints
  nothing.
- **`so2` cannot be started from NBO's own Lewis structure** and needs a hand-written
  `$NRTSTR`; the four 9-valence-pair structures are in `config.py` as `NRTSTR_SO2`. It is
  appended only to the gennbo copy of the archive — `<mol>_native.47` is taken *before* the
  append, because the native reader has no business seeing it.
- **`%maxcore` in a Python f-string**: a bare `%%` reaches ORCA literally and it dies with
  "expected an identifier after '%%'".

## Layout

```
inputs/<mol>/<mol>.inp        ORCA input at the recorded keyword line
inputs/<mol>/<mol>_start.xyz  the starting geometry, built by make_inputs.py
bin/ (this directory)         make_inputs.py, accept.py, config.py, compare_all.py,
                              env.sh, gennbo7, orca_opt.sh, nbo_stage.sh
compare_nbo.py                a copy of tests/nbo_reference/compare_nbo.py, the comparison
                              engine, imported unchanged by compare_all.py
results/<mol>/                <mol>.out, <mol>.47, <mol>_native.47, <mol>.nbo,
                              <mol>.gennbo.nbo.json, <mol>.native.nbo.json,
                              provenance_orca.json, provenance_nbo.json
accepted.json, comparison.json
```
