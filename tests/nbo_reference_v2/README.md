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
and the NPA charge comparison is **blind** to it rather than merely 47x coarser: benzene's valence
set is 0.31611 e short and its Rydberg set 0.31605 e long while its worst NPA charge deviation is
0.00673 e, and the two class arms below move 0.695 e of benzene's Rydberg population while all 111
NPA charges agree to 1e-10.  A charge is a sum over classes, so a leak between two classes cancels
in it exactly.  The distinction decides what to do next: a coarse metric would be worth tightening,
a blind one has to be replaced - which is what the per-orbital AO -> NAO comparison does. It is not a per-Rydberg-function normalisation error: `ryd_per_function.py` shows d per
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

## One truth per molecule: `parser_version`

Three parser defects were fixed in this lane (gennbo's NRT `RS` column is a rank, not a structure
number; its composite alpha+beta valency table is its own table and must not overwrite `beta`; an
open-shell NAO table prints `Spin` where the old code read `Energy`).  That made every stored JSON
written before the fix *wrong*, not merely old — and for a while there were two copies of the
reference on the cluster, `nbo_ref_v2/` with the old parse and `nbo_ref_v2_fixed/` with the new one.
A whole halving measurement was made against the old copy and had to be thrown away.

A warning in a README does not stop that happening again, so:

* `NBO_JSON_PARSER_VERSION` in `Src/core/nbo_run.h` is stamped into every JSON the writer produces
  as `"parser_version"`.  Bump it whenever a change alters what a stored reference *means*.
* `config.load_nbo(path)` is the only way the scripts here read a stored reference, and it raises
  `SystemExit` on a file whose stamp is not the current one.  Eleven readers were converted;
  `index.json` and `provenance_nbo.json` are not parser output and are left alone.
* Every `.nbo` kept in `data/<mol>/` and on the cluster was re-parsed **in place** by
  `regen_refs.sh`, and the native side was re-run in place with the same per-molecule flags, so
  both files in a directory carry version 2 and nothing in the tree predates the fix.
  `verify_one_truth.py` then compared the re-parse with the retired `nbo_ref_v2_fixed/` copy field
  by field: 23 of 23 molecules identical outside `parser_version`, `source` and `timings`.
* The retired copy keeps its files **unstamped on purpose**, with a `README-SUPERSEDED.txt`.  A
  second copy of a reference is how the wrong one gets measured, so it now fails loudly instead of
  answering.

Re-parsing does not weaken the arbitration: `-nbo_parse` re-reads gennbo 7's own output text, which
is kept and unchanged.  The reference is still NBO 7's numbers; only our reading of them is new.

One file did not survive that sweep and could not: `ch3_original/ch3_original.native.nbo.json` was
computed from `ch3.gbw` — a *different molecule* that happens to share the name, which is exactly
the pairing that made the old dataset unusable.  The archive's own wavefunction is gone, so it is
renamed `ch3_original.native-from-ch3-gbw.NOT-A-PAIR.json` and cannot be picked up as a pair again.

## Step 3 was not where the leak has to be

`nao_class_leak.py` and the `NAO_CLASS_SPLIT` arm between them exonerate step 4: the split variant
is 5.6x worse, and `ryd_excess_by_l.py` puts +1.15335 e of the +1.256 e Rydberg excess in `(atom, l)`
blocks that hold no valence shell at all — blocks a unitary intra-block rotation cannot feed.  Step
1 is per-`(atom, l)` as well.  Only step 3's Schmidt projection crosses `l`.

`NAO_DUMP_STEP3=1` now prints, for every shell, both the occupancy step 4 inherits and the pre-NAO
occupancy it started from (the dump goes to `NoSpherA2.log`, not stdout).  `step3_stages.py` reads
those logs and puts four stages side by side — pre-NAO, after step 3, final, gennbo — and
`step3_stages.sh` produces them for all 25 molecules in place.  The measurement on ethane and
benzene that motivated it:

| molecule | after step 3 | final | gennbo | step 4 recovered |
|----------|-------------|-------|--------|------------------|
| ethane   | 0.43210 | 0.15471 | 0.02486 | 68 % of the excess |
| benzene  | 1.12682 | 0.43157 | 0.11552 | 69 % |
| sf6      | 0.71526 | 0.34846 | 0.35735 | 102 % |
| pf5      | 0.50271 | 0.24899 | 0.26845 | 108 % |

So step 4 is not a cosmetic step: it repairs most of what step 3 leaves, and it repairs the
hypervalent cases *completely*.  What it cannot repair is population parked in a block with no
valence shell — which is precisely the residual.

The pre-registration and the acceptance gate for any later change to step 3 are in
`step3_stages.py`'s docstring, fixed before the run so they cannot be relaxed afterwards: pf5, so2
and sf6 must keep a final Rydberg total within 0.90–1.05 of gennbo's, no molecule may get worse by
more than 0.01 e, and the no-valence-block excess must fall.  A change that improves the mean while
flattening those three has found a fudge factor, not a fix.

### What the arms then measured (job 586936, `nao_step1_arms.sh`, all 25 molecules)

The localisation above is wrong, and the gate caught the first candidate that would have
exploited it.

`NAO_OWSO_OFF=1` drops the occupancy weighting inside each class and lets the plain Löwdin do the
whole within-class orthogonalisation.  It is not a candidate — it throws the OWSO away — it is the
test of one claim: after step 3 the Rydberg set is the S-orthogonal complement of the span of the
natural minimal pre-NAOs, so a class **total** at that stage cannot depend on anything inside step 3.
Confirmed: every molecule's `s3_Ryd` is unchanged to all five printed decimals while every `fin_Ryd`
moves.

But the *distribution* over `(atom, l)` blocks is not invariant, and that is what the final answer
sees — the no-valence-block excess went 1.39828 → 1.14700 e and the final excess 1.25585 → 1.12478 e.
So "only step 3's Schmidt projection crosses `l`" was the wrong localisation: the weighting inside
step 3 moves the final numbers too.

And it is still not a candidate, because it is what the gate was written to catch.  The mean improves
while pf5 goes 0.93 → 1.12, so2 0.94 → 1.05 and sf6 0.98 → 1.19 × gennbo, sf6 worse by 0.077 e on its
own.  That is the fudge factor: flattening the three molecules that were already right in order to
improve the average.

`NAO_PRENAO_NET=1` builds the pre-NAOs from the atom's own block of `P` with `S^A` as the metric —
the net atomic population, the way the 1985 paper reads — instead of from `(S P S)^A`.  It was
rejected years ago on epoxide's NPA charges, which the class-split arm proved blind to exactly this
kind of error, so it needed re-running against the Rydberg metric.  It is not a near miss: the
Rydberg total after step 3 goes 10.40557 → 194.78163 e and the final one 4.00571 → 56.45473 against
gennbo's 2.74986, every molecule between 8× and 72× worse.  Gross is confirmed, and no longer on a
charge argument.

Together those two prune the search rather than narrow it.  The class totals step 4 inherits are
fixed by the span of the natural minimal pre-NAOs alone, and both inputs to that span are already at
their better setting — so there is nothing upstream of step 4 left to change.  What is left is which
vectors step 3 hands to each block and what step 4 does with them.  NBO 7 prints no intermediate
*occupancy* table - but it does emit the transformation itself, so the vectors are arbitrable even
where the intermediate occupancies are not.


## The arbitrated intermediate (`aonao_probe.sh`, `aonao_stage.sh`, `aonao_compare.py`)

`$NBO AONAO=W $END` makes NBO 7.0.9 write **its own AO -> NAO transformation matrix** to logical
file 33 at nine decimals, and `AONAO` alone prints the same matrix into the `.nbo`.  Probe job
588007 established that on this install, against the binary's own `$NBO HELP $END` text rather than
a remembered manual - and it also killed the guessed unit number: `AONAO=W48` is refused, because
lfn 48 is reserved.  The whole matrix-output family exists on the same terms (AOPAO, AOPNAO, AONHO,
AONBO, AONLMO, AORNBO, AOMO, PAOPNAO, NAONHO, NAONBO), so AO -> NBO is available later without a
second capability question.  The native side needed no new code beyond a dump: `NAOResult::C`
already holds our AO -> NAO matrix, and `NAO_DUMP_C=1` prints one `NAOC` line per NAO with its
(atom, l, m, shell, class, occupancy) and its AO coefficients.

What makes this a fair comparison rather than a convention argument: **both sides read the same
`<mol>_native.47`**.  The reader in `nbo.cpp` builds its AO map from that file's own CENTER/LABEL
arrays and takes S and P from `$OVERLAP`/`$DENSITY`, so the AO order and both matrices are literally
the same arrays on both sides and nothing is converted.

The metric is an m-averaged mixing matrix between the two sides' shells inside each (atom, l) block,
`M[a][b] = (1/(2l+1)) sum_{m,m'} (c_native(a,m)^T S c_gennbo(b,m'))^2`, paired **by rank** within the
block.  Squaring makes it sign-blind and the m sum makes it blind to component order, so neither
convention can masquerade as an error.  It is reported as two separate numbers: the **in-block leak**
`1 - M[a][a]`, which an intra-atomic step could produce, and the **out-of-block remainder**, which it
could not.

Three gates run before any number is quoted, because three metrics have already been retired on this
branch for being incommensurable:

1. the lfn 33 layout is decided by which reshape satisfies `C^T S C = 1`, and the check asserts the
   loser **fails** - with `S = 1` a matrix and its transpose would both pass, so the demo uses a
   non-trivial metric;
2. the occupancies printed in the same run must come back out of the matrix as
   `diag(C_g^T S P S C_g)` (that is the NAO-basis density because `C^-1 = C^T S`; the wrong form
   `C_g^T P C_g` is computed too and the report names which one matched);
3. the AONAO runs use a reduced keylist (`AONAO=W` alone, since the NAO stage is upstream of NBO, E2
   and NRT), so that run's own NAO table must equal the stamped reference JSON's to 1e-5 - otherwise
   the molecule is declared **VOID** instead of compared.

The number of VOID molecules is printed next to the number compared, so a shrinking denominator
cannot flatter a mean - the failure that has just cost another lane a corpus headline.

Three more gates test the **metric** rather than `nao.cpp`.  The exact one is completeness: each native
shell's m-averaged mixing summed over *all* gennbo NAOs must be 1, because both sides are
S-orthonormal sets spanning the same AO space, so any deviation is an error in the matrix read, the
pairing or the block bookkeeping and nothing may be quoted until it passes.  The soft one is the
contrast: mean in-block leak on pf5/so2/sf6, whose final table is already right, against ethane and
benzene, which are worst - if the metric does not separate those, it is measuring something other
than the leak it exists to explain.  The third closes completeness's blind spot: a row sum is
invariant under any relabelling of gennbo's shells, so it cannot see a pairing error - the defect
that once discarded 11 of 22 molecules - and the pairing diagnostic therefore prints, per molecule,
how many shells have their argmax on the rank partner and the distribution of (argmax - rank)
offsets.  It is deliberately diagnostic rather than pass/fail, because a genuine leak *will* move an
argmax off the diagonal and that is the finding; a systematic off-by-one instead appears as most
shells sharing one non-zero offset, and the report says so in those words.  When completeness does
fire, the report distinguishes a **shortfall** - missing weight, so lfn 33 is not the full n x n set
and the fault is on the read side - from an excess, which would have to be block bookkeeping: 0.93
instead of 1 reads like a small error and means a missing column.

A third gate was proposed and **does not survive the algebra**, so it is not implemented: that the
out-of-block remainder should be exactly 0 on LiF or N2 because the pipeline is intra-atomic.  NAOs
are orthonormal over the whole molecule on each side, not per atom, so a native NAO on one atom
generically has amplitude on another atom's and on other l when expanded in gennbo's basis - and that
is the inter-l leak this metric exists to measure (92 % of the final-table excess crosses l), not an
artefact of it.  A diatomic does not change it: two atoms' NAO sets are mutually orthogonal within
one side, never between sides.

### What it says (jobs 589060 + 589634, AKL012, 2 threads, 8 molecules, 0 VOID)

The staging job ran in seconds once it was shaped like the probe, and every gate passed: lfn 33 is
column-major on all eight molecules with `|C^T S C - 1| <= 2.8e-09`, the printed occupancies come back
out of it as `diag(C^T S P S C)` to `5.0e-06` while the wrong form does not, all eight agree with the
stamped reference JSON at `0.0e+00`, completeness holds to `2.74e-09`, and the contrast gate separates
the molecules whose final table is already right (mean in-block leak 0.17423 over 182 shells) from the
two worst (0.43950 over 136).

sf6 was **VOID on the first pass** with a 0-byte lfn 33 and `rc=0`, and the cause is worth keeping: NBO
stopped at `SYMOPS: generated 48 symmetry operator(s) for Th but expected 24` after 0.02 CPU seconds,
which is *before* the NAO stage, so it produced a valid-looking `.nbo` of 1190 bytes and no matrix at
all.  `NRTSYM=OFF` suppresses it and job 589634 recovered the full 1.15 MB matrix.  A job that exits 0
and writes an empty file is exactly what the VOID denominator exists for.

That left one keyword in sf6's keylist and not in the other seven's - a per-molecule input difference
inside the very gate that rests on sf6 - so it was measured rather than argued about.  Job 589969 reran
lif and water with `AONAO=W NRTSYM=OFF`, on the same node, against the same `.47`, and both lfn 33 files
came back **byte-identical** to the committed ones (`md5 4b7206f8...` and `a2d5f3e1...`; the only
difference anywhere in the input is the keylist line itself).  `NRTSYM` does not reach the NAO stage, at
`0.0e+00` and not at a tolerance, so sf6's matrix is comparable with the other seven's and the contrast
gate is not mixing two inputs.  The record is `data_nao_split/aonao/nrtsym_probe_589969.txt`.

Three arms, because a rank pairing is not physics and a unit-normalised shell is not a population:

| by class | mean `1-M[a][a]` (rank-paired) | mean `1-max_b M[a][b]` | mean weight outside the (atom,l) block | shells |
|---|---|---|---|---|
| Core    | 0.00000 | 0.00000 | 0.00000 | 34 |
| Valence | 0.00783 | 0.00783 | 0.00367 | 72 |
| Rydberg | 0.37253 | 0.27308 | 0.17676 | 273 |

Read per shell, the Rydberg set is the whole problem - 48x the valence error - and the cores are
*exactly* right, which exonerates the AO read, S, P and the core partition in one number.  Of the
Rydberg 0.373, 0.099 is only an ordering difference (LiF's 2p Rydberg pair scores 0.00007 on its rank
partner and 0.98861 one rank over: right shape, different order), 0.096 is genuine rotation inside the
(atom,l) block, and 0.177 leaves the block altogether.  The `+-1` offsets in the pairing diagnostic are
symmetric (benzene `-1:24 +0:42 +1:24`), so they are mutual swaps, not the shared non-zero offset that
would mean a systematic pairing error; sf6's are one-sided (`+1:7` of 79) but too few to be that either.

Weighted by occupancy the ranking **inverts**, and that is the finding:

| electrons, sum `(2l+1)*occ*leak` | rank-paired | ordering-insensitive | out-of-block |
|---|---|---|---|
| Core    | 0.00026 | 0.00026 | 0.00001 |
| Valence | 0.94266 | 0.94266 | 0.51486 |
| Rydberg | 0.39467 | 0.31255 | 0.19472 |
| total   | 1.33759 | 1.25547 | 0.70960 |

**70 % of the mis-shaped density is in the valence shells**, and 0.51486 e of it lands outside its own
(atom,l) block - against 0.19472 e for the whole Rydberg set, whose per-shell error is 48x larger.
None of the valence number is an ordering artefact: the two arms agree to five decimals, because
valence shells are not near-degenerate enough to swap.  So the per-shell picture ("the Rydberg
construction is wrong") and the charge picture ("the valence shells carry the error that moves
populations") are both true, and only the second is commensurable with the NPA failure this exists to
explain.  It is also what `d(Val) = -d(Ryd)` looks like one level down: the valence shells themselves
are built with a small shape error that carries a lot of charge, and it surfaces as Rydberg population.

This is the fourth metric on this branch to rank the wrong thing until it was weighted the way the
failure is, so the charge-scale arm is printed next to the mean and neither is quoted alone.

### Where that weight goes, and whether it reaches the final table

The out-of-block weight can only be two things, and they mean different things, so they are split:

| electrons | in-block, other shell | same atom, other l | another atom |
|---|---|---|---|
| Core    | 0.00025 | 0.00000 | 0.00001 |
| Valence | 0.42780 | 0.14081 | 0.37405 |
| Rydberg | 0.19995 | 0.05639 | 0.13834 |
| total   | 0.62800 | 0.19719 | 0.51240 |

Of the valence set's 0.94266 e, 0.42780 e mixes another shell of the **same** (atom, l) - which is the
valence/Rydberg swap the final tables report, seen on the donor side - 0.14081 e crosses l on the same
atom, and 0.37405 e sits on **another atom entirely**.  So the intra-atomic part, the only part that can
move population between classes of one atom, is 0.56861 e, and within it the same-l channel is three
times the cross-l one.  That is *not* the same decomposition as the final table's "92 % of the excess is
inter-l" (that one is a receiver-side statement about (atom, l) blocks with no valence shell), and the
two must not be quoted as the same number: the cross-l channel is real here but it is the smaller half
of the intra-atomic weight.

The inter-atomic 0.37405 e has nowhere to go in a class table, and it does not show up as charge either:
the worst atomic charge deviation across the eight is 0.02802 e, a factor of 13 smaller.  Inter-atomic
mis-shape is largely reciprocal between bonded partners, so it cancels - the same cancellation that made
the NPA charge agreement meaningless, now measured one level upstream.

Per molecule, against the final table's own class error `d(Val)` (native minus gennbo, from the stamped
reference JSONs, not copied from an earlier printout):

| molecule | d(Val) | intra Val | intra Ryd | ratio | inter Val | worst dq |
|---|---|---|---|---|---|---|
| ammonia | -0.02959 | 0.02990 | 0.00352 | 1.13 | 0.01592 | 0.00996 |
| benzene | -0.31611 | 0.36191 | 0.13222 | 1.56 | 0.22158 | 0.00673 |
| ethane  | -0.12989 | 0.13181 | 0.07932 | 1.63 | 0.05841 | 0.00903 |
| lif     | -0.00116 | 0.00037 | 0.00060 | 0.83 | 0.00010 | 0.00197 |
| pf5     | +0.01945 | 0.00884 | 0.01435 | 1.19 | 0.02523 | 0.01308 |
| sf6     | +0.00885 | 0.01162 | 0.01179 | 2.65 | 0.03157 | 0.01632 |
| so2     | +0.01701 | 0.01424 | 0.01213 | 1.55 | 0.01357 | 0.02802 |
| water   | -0.01268 | 0.00992 | 0.00241 | 0.97 | 0.00768 | 0.00165 |

The intermediate and the final table are the same size on every molecule, and on the two that matter -
benzene at 0.31611 e and ethane at 0.12989 e, the largest class errors in the set - the intra-atomic
valence mis-shape alone is 0.36191 e and 0.13181 e, i.e. 114 % and 101 % of what the final table shows.
The ratio never exceeds 2.65 and never falls below 0.83.  That is the bridge: the AO -> NAO matrix
carries enough mis-shape, on the right atoms and in the right class, to produce the population error
that `-nbo_native` fails on, and nothing downstream has to be invoked to explain it.

What it is **not** is proof that nothing downstream also contributes.  The occupancy weight is an
estimate - the population a mis-shaped shell actually moves is `w*(occ_a - occ_b)` and this uses
`w*occ_a`, and M is m-averaged, so a shell's components are assumed to share its occupancy.  On the four
molecules whose `d(Val)` is 0.001-0.03 e the ratio scatters 0.83-1.19 either side of 1, which is the
estimate's own resolution; nothing here can check the **sign**, since a mixing weight is positive while
pf5/so2/sf6 have native's valence set too large and the other five have it too small.

Two read defects were caught by the gates on the way, both of the family this branch keeps meeting:

- `unpack_upper` filled the `.47`'s packed triangle row by row, but UPPER is packed **column by
  column**, so S and P were scrambled and the layout gate refused all five molecules.  It was right
  and the fault was mine one function earlier.  Its doctest had used `n = 2`, where the two orders
  coincide - a test that could not fail; it uses `n = 3` now.
- gennbo prints an (atom, l) block **component-major** (water's oxygen p set is `2px 3px 4px 5px`,
  then `2py ...`), so slicing it in groups of `2l+1` grouped four shells of one component together.
  Every p shell then scored `M[a][a] ~ 1/3` - one of three m pairs matching - and the run passed all
  three gates while producing nonsense.  The grouping now reads the self-identifying `lang` field, and
  the demo carries the case that scores 1/3 under the old slicing and 1 under the new.

The kept inputs are in `data_nao_split/aonao/<mol>/`: gennbo's own lfn 33 (`<mol>.33`, the external
reference in its original form, so a reader fix is a local re-parse), the native `NAOC` dump, the
AONAO run's parsed JSON, both jobs' stdout and the full comparison table.  The `.47` is not duplicated
here - it is the same `<mol>_native.47` listed under "Deliberately not committed" in `data/README.md`.

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
