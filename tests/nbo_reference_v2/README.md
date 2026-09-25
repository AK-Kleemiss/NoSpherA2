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

`compare_all.py` takes the two JSONs from the same molecule directory and has **no default
candidate directory** — the defect that started this whole thread was `compare_nbo.py`'s
candidate path defaulting to the reference directory, so that `--all .` compared the dataset
with itself. Here both paths are named per molecule and printed, a molecule missing either
JSON is a failure rather than a silent skip, and the summary states how many molecules and
how many data points the comparison rests on.

It reports NPA charges, NAO occupancies and energies, NBO occupancies and
hybrid composition, E2 entries, and NRT bond orders, valencies and weights. gennbo is the
reference side, so a quantity gennbo printed and we did not is a `missing`, not a silent
pass. **NRT weights are read off the by-rank comparison, never by label**: starting the
search from a different Lewis structure renumbers the structure list and swaps the
descriptions of equivalent structures without moving a single weight.

Both sides are put on the thresholds NBO itself reported (`index.json`
`thresholds_reported_by_nbo`): `config.py` translates them into `-nbo_e2min` and `-nrt_e2` for
the native run, so the comparison measures the method and not a printing cut-off.

### Two views, both printed

`compare_all.py` reports every molecule twice. `ungated` compares every printed orbital;
`gated` drops the empty Rydberg NBOs and the E2 rows that accept into one. The engine already
does this for Rydberg *NAOs* and its reason applies verbatim to NBOs: any unitary rotation
inside the empty space is an equally valid NBO set, so its labels, ordering and one-electron
energies are that code's arbitrary choice. ch3's first five gennbo RY energies (3.50, 26.26,
1.35, 12.33, 1.99) against native's (2.95, 1.56, 5.68, 1.40, 2.21) are a permutation-like set,
not a disagreement — while their occupancies agree to 9.8e-4 across all 36.

**The gated column never replaces the ungated one, and the ungated verdict is the verdict.** A
gate that is only ever shown after it has removed most of the failures cannot be told apart
from a gate that was tuned until the test passed, so both columns are printed with the number
of points each dropped. On water that is 429 points / 206 failed against 217 / 72, and the
gated view still fails — which is the point: the gate says *which* comparisons failed, it does
not grant permission to stop counting them.

### A surplus is a disagreement too

`compare_keyed` in the engine walked the reference side only, so an entry the *candidate* printed
and gennbo did not was neither `compared` nor `missing` — ch3's twelve native-only E2 rows left no
trace at all in the very check that exists to find disagreements. The counter now lives in
`tests/nbo_reference/compare_nbo.py` itself, one line at the end of `compare_keyed`, and
`Quantity.ok()` fails on a surplus exactly as it fails on a missing entry. It belongs in the
engine rather than in a wrapper because the old 22-molecule gate imports the same function and
had the same blind spot; the byte-identical second copy of the engine that used to sit in this
directory is gone for the same reason. `compare_nbo.py --selftest` covers both directions, and
`--all <dir>` now **refuses** the reference directory itself — run that way it compared every
molecule with its own file and reported 22/22 at max dev 0.00000.

### Is a disagreement a factor of two? `ratio_probe.py`

`ratio_probe.py <dir> <mol> --expect 2.0` prints `native/gennbo` for every matched pair, per
quantity, plus the spread. "Exactly twice" as a summary cannot distinguish twelve pairs at 2.00
from eleven at 2.00 and one at 1.6, and an *empty* gennbo table against twelve native entries is
merely *consistent* with a factor of two — which is why `e2_threshold_probe.sh` lowers the print
threshold on **both** sides to 0.02 kcal/mol until gennbo prints its own numbers.

The test is on the residual `native/2 - gennbo` against half of gennbo's own last printed digit,
not on a ratio band. gennbo's E2 table carries two decimals, so a printed 0.20 is anything in
[0.195, 0.205) and the ratio against native's 0.39517 is 1.976 whatever the truth is; calling
that a 1.2 % failure would be measuring NBO's report format. Covalent/ionic *shares* are printed
separately because a share is a ratio of two quantities that are both doubled — it survives the
factor, so anything left in it is a second, independent finding.

What it found on ch3 (slurm 580673, AKL064, 4 threads):

- **NRT: a factor of two, exactly.** 16 of 16 matched pairs at ratio 2.0000, residual 0.00000 —
  valency 1.50000 → 3.00000, electron count 4 → 8 and 3 → 6, bond-order total 0.50000 → 1.00000.
  Native doubles *each spin's own* value rather than printing a spin-summed total: gennbo's two
  electron counts differ (4 and 3) and native's differ correspondingly, where a sum would be 7
  in both blocks.
- **E2: not a factor of two, and the missing alpha table is not a missing table.** At 0.02
  kcal/mol gennbo printed 21 entries and native 44, with beta at gennbo 6 × 0.0600 against native
  6 × 0.3913–0.3952 (ratio 6.52–6.59) and **no** alpha rows from native at all. Rerun at 0.001
  kcal/mol on both sides (slurm 581532, AKL064, 4 threads) native prints them: its alpha
  C–H → C–H\* entries are **0.0148 kcal against gennbo's 0.2000**, a ratio of 0.071. So the two
  spins are wrong in *opposite* directions, 13.5× too small in alpha and 6.5× too large in beta.
  Over all 54 matched pairs the ratio runs 0.0714 to 27.57, mean 6.14. It is not a spin swap
  either: a swap would have printed gennbo's own 0.0600 and 0.2000.
- **The covalent/ionic partition says the same thing.** `cov + ion = val` holds exactly on both
  sides and the totals are a clean 2×, yet the ionic share is 19.20 % gennbo against 24.79 %
  native in alpha and 13.61 % against 5.15 % in beta: **+5.59 points in one spin, −8.46 in the
  other**, which no single factor can produce. (An earlier draft said 16.40 % for beta. That was
  NBO's composite alpha+beta table read as a beta table — see the composite trap below.)
- **The lead is the spin-resolved NBO/NAO construction, not the E2 code.** Three independent
  quantities put alpha and beta delocalisation the wrong way round: the E2 energies above, the
  antibond occupancies (alpha 0.00002 e native vs 0.00056 gennbo, beta 0.00174 vs 0.00015) and
  the ionic share. `E2 = q F(i,j)²/Δe` and `q* ≈ 2q(F/Δe)²` share the same off-diagonal Fock
  element, so one wrong spin-resolved Fock or density matrix produces all three; an E2 bug
  produces only the first. An absolute occupancy tolerance cannot see any of it: every value is
  within 1.6e-3 of 1 or of 0, so a 28× relative difference in an antibond population passes a
  2e-3 tolerance untouched.
- **NPA is the one block that nearly survives**, ratios 0.9099 to 1.0126 — carbon −0.49127 gennbo
  against −0.44722 native, 0.044 e. Being spin-summed it partly cancels the inversion, which is
  consistent with the same cause rather than a fourth one.

None of this changes a convention in the native code. A factor of two applied in the wrong place
is how this thread started.

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
- **`E2PERT` without a value is accepted and silently ignored.** Only `E2PERT=<value>` sets the
  threshold; `E2PERT` and `E2PERT <value>` run to completion, exit 0, print no warning and leave
  NBO's default (0.5 kcal closed shell, 0.25 open shell). The surviving ch3 archive carried the
  bare form, so the keylist this project inherited never set anything. The two sides happened to
  match anyway — because the value `index.json` records *is* what NBO printed under the ignored
  keyword — which is agreement by luck, and `config.py` now spells the value out.
- **NBO prints THREE NRT tables for an open shell**, alpha, beta and "(composite alpha+beta)",
  and the composite ones come *after* the last `Beta spin orbitals` header. A parser that takes
  the spin from that header files the composite table as a second beta table:
  `-nbo_parse` did exactly that for valencies (not for bond orders, which read their own title),
  so ch3's JSON has 12 valency rows, 8 of them tagged `beta`, and any consumer keying by
  `(spin, atom)` kept whichever came last — carbon's real beta row `1.5000 1.2959 0.2041 3.0000`
  replaced by the composite `3.0000 2.5079 0.4921 7.0000`. Fixed in `Src/core/nbo_run.cpp`
  (`valency_spin`, read from the table's own title), but a *stored* JSON written before that fix
  still carries the duplicate, so `nrt_spin_test.py` and `ratio_probe.py` keep the **first** row
  for a key — always the right one, since the composite table comes last — and print a note
  saying they dropped a duplicate rather than averaging it away. The check that the exclusion is
  right: the valency and bond-order ionic shares then agree at 13.60/13.61 %, where before they
  differed by 2.8 points.

## Layout

```
inputs/<mol>/<mol>.inp        ORCA input at the recorded keyword line
inputs/<mol>/<mol>_start.xyz  the starting geometry, built by make_inputs.py
bin/ (this directory)         make_inputs.py, accept.py, config.py, compare_all.py,
                              env.sh, gennbo7, orca_opt.sh, nbo_stage.sh
../nbo_reference/compare_nbo.py  the one comparison engine, imported by compare_all.py,
                              nrt_spin_test.py and the old 22-molecule gate. There is no copy
                              in this directory any more: it was byte-identical and tracked,
                              so "make the engine count a surplus" could not be trusted in two
                              places at once.
ratio_probe.py                native/gennbo per matched pair, resolution-aware
nrt_spin_test.py              each spin's NRT value against its own gennbo counterpart
e2_threshold_probe.sh         both sides at one low E2 print threshold
results/<mol>/                <mol>.out, <mol>.47, <mol>_native.47, <mol>.nbo,
                              <mol>.gennbo.nbo.json, <mol>.native.nbo.json,
                              provenance_orca.json, provenance_nbo.json
accepted.json, comparison.json
```
