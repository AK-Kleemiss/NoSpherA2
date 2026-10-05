# nbo_reference_v2: the NBO 7 reference, with its inputs

`tests/nbo_reference/` holds 22 gennbo 7 result JSONs and nothing that produced them. This
directory holds what is needed to rebuild them, so that the numbers can be re-run instead of
taken on trust. Kept here are the inputs, the geometries, the acceptance records, gennbo 7's
parsed output as the golden files, and the scripts that regenerate all of it. External
program output (ORCA `.out`, gennbo `.nbo`) is not published.

The old 22 JSONs stay where they are, as a regression gate on the parser. Do not treat them
as an external reference: they predate two parser fixes (NRT `RS` is a rank and was read as
a structure number; an open-shell NAO table prints `Spin` where `Energy` was read), so their
NRT weights and open-shell NAO energies are wrong.

## Acceptance (`accept.py`)

`index.json` records, per molecule, the ORCA keyword line, the charge, the multiplicity,
`basis_functions` and `final_energy_hartree`. A regenerated wavefunction is accepted only if:

- the optimisation converged;
- the charge, multiplicity, basis and keyword line match;
- **`basis_functions` matches exactly**, which is the identity test. A difference of one
  function means a different molecule. This is how the old set's `water`, which was really
  OHHHe, should have been caught;
- `|ΔE| <= 2.0e-5` Eh, four times ORCA's own per-step criterion.

`ch3` starts from the optimised `$COORD` in `../nbo_reference/ch3_orca_reference.47`, the one
archive that survived. It reproduces the recorded energy to 1.3e-8 Eh. That figure is the
pipeline's floor.

`make_holdout.py` and `make_radicals.py` add molecules that were never in `index.json`.
`accept_holdout.py` judges those on the criteria that can still be checked: convergence,
charge and multiplicity.

## Running it

Everything runs on the cluster (`ssh AKL007`).

```
python make_inputs.py <out_dir> --nprocs 4     # <mol>/<mol>.inp + <mol>_start.xyz
sbatch orca_opt.sh  <mol>                      # stage 1: optimise, write <mol>.gbw
python accept.py <root> --json accepted.json
sbatch nbo_stage.sh <mol>                      # stage 2: .47, gennbo 7, -nbo_parse, -nbo_native
python compare_all.py <root> --json comparison.json
```

`nbo_all_stage.sh` runs stage 2 over the whole set. `holdout_stage.sh` runs both stages for
one hold-out molecule. `config.py` is the only place that holds the per-molecule keylist
deviations. `env.sh` carries the environment. Each stage writes a `provenance_*.json` with
the host, the threads given and the node's core count.

## Traps

- **`module load orca/...` reports "unknown"** until you run `module use
  /work/software/spack/share/spack/lmod/linux-rocky9-x86_64/Core`. Do not pipe `module load`
  into `head`: the `setenv` is then lost in a subshell.
- **`/work/akkleemiss/share/gennbo` is broken.** Its `NBOBIN=/opt/nbo7/bin` exists on no
  node, and its `rm -i` blocks a batch job. `gennbo7` here is the replacement, pointing at
  `/work/software/bin/NBO/bin`.
- **`BASH_SOURCE` does not locate a batch script**, because slurm copies the script into
  `/var/spool/slurmd/job<id>/`. The scripts take their directory from `$NBOREF_BIN`.
- **`sf6` needs `NRTSYM=off`**, or NBO aborts after 0.08 s and prints nothing.
- **`so2` needs a hand-written `$NRTSTR`** (`NRTSTR_SO2` in `config.py`). It is appended
  only to gennbo's copy. `<mol>_native.47` is taken before the append.
- **`E2PERT` without `=<value>` is accepted and silently ignored**: exit code 0 and no
  warning. Always spell the value out.
- **An open shell prints three NRT tables**: alpha, beta, then the composite. A parser that
  takes the spin from the last `Beta spin orbitals` header files the composite table as a
  second beta table.

## Layout

```
inputs/<mol>/                  ORCA input and starting geometry, as make_inputs.py writes them
data/<mol>/                    the 22 + 3 radicals: <mol>.inp, _start.xyz, .xyz,
                               <mol>.gennbo.nbo.json (golden), provenance_orca/nbo.json
data/accepted*.json            acceptance verdicts
data_holdout/<mol>/            the 16 hold-out molecules, same layout
../nbo_reference/compare_nbo.py  the comparison engine, imported by compare_all.py
```

`.gbw`, `.47`, `.out` and `.nbo` stay on the cluster under
`/work/akkleemiss/florian/nbo_ref_v2/<mol>/`. Any of them can be rebuilt with one stage
command.
