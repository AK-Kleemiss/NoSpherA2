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

## The NAO class partition, and why NPA cannot referee it

`nao_class_leak.py <root>` decomposes each side's NAO table into Core/Valence/Rydberg per atom
and pairs the rows **by rank within an (atom, l) block**, never by the printed shell label: NBO 7
numbers ethane's C 1s block `1s 2s 3s 5s 4s` where native numbers it `1s 2s 3s 4s 5s`, and pairing
on that label reported 11 of 22 molecules "not comparable" until the rank pairing replaced it. The
printed principal quantum number is each side's private bookkeeping. **Rank, never label** - the
same trap costs the E2 table matches through its `BD ( 2)` ordinal, and `label_pairing_audit.py`
measures how many.

What it found, over the 22: the core is exact everywhere (|d| <= 3e-5 e), and `d(Val) = -d(Ryd)`
to five decimals in every single molecule. The error is an intra-atomic valence -> Rydberg leak,
and it is 47x larger than the NPA charge comparison lets you see - benzene's valence set is
0.31611 e short and its Rydberg set 0.31605 e long while its worst NPA charge deviation is
0.00673 e. It is not a per-Rydberg-function normalisation error: `ryd_per_function.py` shows d per
function spanning -0.000694 to +0.001939 e *including sign changes*, with the top 3 of 26 functions
carrying 38-88 % of the magnitude where flat would be 11.5 %. It is not spin-resolved either -
ch3's leak (0.03418 e) is smaller than closed-shell ethane's (0.12989 e).

### The two-arm experiment (`nao_split_arms.sh`, `NAO_CLASS_SPLIT`)

`nao.cpp` step 4 re-diagonalises the m-averaged density inside each (atom, l) block with valence
and Rydberg **together**, and its comment claimed that separating them "inflates the Rydberg
occupancies tenfold". That was an assertion, so `NAO_CLASS_SPLIT=1` now runs the separated variant
and `nao_split_arms.sh` runs both arms over the same 22 wavefunctions and the same gennbo JSONs.

Pre-registered before the run, and written into the script: split will be **worse**, because the
top eigenvalue of the m-averaged block is an upper bound on any single shell's own diagonal, so the
valence shell can only come out with more population from one block than from a split one.

Measured: mean |intra-atomic leak| 0.01229 e one-block vs **0.06897 e** split, over 111 atoms.
Benzene's Rydberg set 0.43157 -> 1.12682 e. Split also loses the one structural feature the
baseline has: pf5, so2 and sf6 are the only three molecules whose leak has the *opposite* sign in
the one-block form, and splitting drives all 22 the same way. One block is the better of the two,
the comment is accurate, and step 4 is exonerated - what is left of the leak is upstream, in step
3's class orthogonalisation or in which (n, l) count as valence per element.

### And it is not step 4 at all: 92 % of the excess is inter-l

`ryd_excess_by_l.py <root>` splits the Rydberg deviation by (atom, l) block and asks whether the
block holds a valence shell. Over the 22, **+0.10250 e lands in valence-bearing l blocks and
+1.15335 e in blocks with no valence shell** - the d and f polarisation shells of first-row atoms,
where the worst single block is ethane's C d at +0.03977 e. Step 4 is unitary *within* a block and
those blocks contain no valence population to take, so it cannot be the source: the leak crosses l,
and inside this construction only step 3's Schmidt projection of one class out of another can move
population between l values. That is where to look next, and the valence-shell definition is
exonerated too - `nao_class_leak.py` reports zero rows classified differently by the two sides once
they are paired by rank.

The hypervalent three behave the other way round in this decomposition as well: pf5 -0.01703 e,
sf6 -0.01802 e and so2 -0.00322 e in the no-valence blocks, against +0.25135 e for benzene and
+0.10830 e for ethane. Their baseline Rydberg totals are 0.93-0.98x gennbo's, i.e. essentially
right, while benzene's is 3.74x and ethane's 6.22x. Any candidate fix has to keep them right -
that is what makes them the acceptance test rather than three more data points.

For the record on the rejected arm: split takes benzene's Rydberg population to 9.75x gennbo's
0.11552 e and ethane's to 17.38x. The old comment's "tenfold" was not a figure of speech.

**And the metric that cannot see any of it.** The two arms differ by 0.695 e in benzene's Rydberg
population and by **1e-10 in every one of 111 NPA charges**. Step 4's mixing is intra-atomic and
unitary, so it moves no charge between atoms - which means NPA agreement is not evidence about the
class partition, in either direction. Any future claim about the NAO construction has to be made
on the NAO table, not on charges.

## Three open shells that were never in the lost run

`make_radicals.py` adds allyl, hco and no2 - doublets with an alpha/beta asymmetry large enough to
tell "doubles each spin's own value" apart from "sums the two spins". They have no
`thresholds_reported_by_nbo` record because they were never in the original run, so `config.py`
gives them the same thresholds as the open shells that do (`NO_RECORD`, written out rather than
defaulted silently) and `accept.py` takes `<S**2>` near 0.75 plus a converged optimisation in place
of a stored basis-function count. Accepted: allyl 123 basis functions, `<S**2>` 0.778384,
E -117.227034018565; hco 68, 0.753460, -113.848612257692; no2 93, 0.753692, -205.080839945861.

### They refute "just divide by two"

`nrt_ratio_stats.py` and `nrt_halving_test.py` apply `ratio_probe.py`'s resolution-aware residual
test (`|native/2 - gennbo| <= 0.0005 + 0.01|gennbo|`, gennbo's NRT tables carry three decimals) to
**every** spin-resolved NRT number rather than the subset one probe matched:

| molecule | numbers | halving repairs | halving fails | worst residual |
|---|---|---|---|---|
| o2 | 34 | 34 | 0 | - |
| no | 34 | 22 | 12 | 0.0365 |
| ch3 | 74 | 46 | 28 | 0.1268 |
| hco | 54 | 29 | 25 | 0.0565 |
| no2 | 60 | 20 | 40 | 0.1211 |
| allyl | 154 | 55 | 99 | 0.1838 |

o2 is the only molecule where the error *is* a factor of two. On ch3 the two worst failures are
the C's covalency (+0.1268) and its electrovalency (-0.1268) - equal and opposite, so halving fixes
the total valency and leaves the covalent/ionic partition wrong. no2 additionally breaks a symmetry
gennbo keeps: its two symmetry-equivalent oxygens come out at 1.5152 and 1.7753 alpha valency
against gennbo's 0.8355 for both.

`nrt_spin_test.py` had been skipping all of this. Its `discriminates()` guard bailed when a molecule
did not have exactly two spin labels, and the fixed parser now correctly stores gennbo's composite
alpha+beta valency table under a third label, `composite` - so a guard meant to stop a
non-discriminating molecule passing a wrong fix was instead stopping every open shell from being
tested at all. A spin is alpha or beta; the composite block's absence on the native side is now one
reported line instead of one failure per row.

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
data/<mol>/                   the committed record of the run: <mol>.inp, _start.xyz, .xyz,
                              _trj.xyz, .out, .nbo, both .nbo.json, both provenance files.
                              NOT .gbw (27 MB) or .47 (42 MB) - those stay on the cluster and
                              are one command to rebuild. data/README.md carries the recipe and
                              the acceptance table, which is what tests/nbo_reference/ lacks.
data/accepted.json            per-molecule acceptance verdict, keyed by molecule name
$ROOT/<mol>/ on the cluster   the live run directory both stages write, by default
                              /work/akkleemiss/florian/nbo_ref_v2; also holds .gbw and .47
```
