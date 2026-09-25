# The held-out set: what is expected of it, written before any number exists

`renat5` was selected on eight molecules - lif, water, ammonia, ethane, benzene, pf5, so2,
sf6 - out of an enumerated space of 120 definitions in which **nothing closed a single
molecule**. It was nominated by name before it was measured and the 120-member sweep then
failed to displace it on either headline figure, which is a stronger claim than "best of
144". But every number attached to it comes from those eight, and eight molecules that are
seven second-row closed shells plus two third-row hypervalents is a narrow view of a cascade
whose defect is an intra-atomic valence -> Rydberg leak crossing angular momentum.

This file is committed **before the first job of the set finished**, so that the predictions in
it cannot have been written around the answers. A prediction written after the measurement is
not a prediction.

## What is being tested, and against what

Nothing here is expected to agree with NBO 7. `-nbo_native` fails the acceptance gate - an
**exact** basis-function match against gennbo - on 22 of 22 molecules in every arm measured so
far, and the held-out set is not expected to change that. What is being tested is whether an
**improvement** selected on eight molecules is still an improvement on molecules it was never
shown. Every number in the result table is a measured delta against gennbo's own printed
values with its instrument floor beside it, never an agreement.

## The three instruments, pinned before use

1. **Worst Rydberg population error, `dpop`.** Per molecule, the absolute difference between
   the total Rydberg population `-nbo_native` reports and the total Rydberg occupancy gennbo
   prints, both read out of the stored `parser_version: 2` JSONs. On the eight tuned
   molecules the shipped native's worst is **0.3161 e** (benzene: 0.4316 against gennbo's
   printed 0.1155) and `renat5`'s replica is **0.0050 e**, a factor 63.
2. **NPA charge delta.** Per molecule, the largest per-atom difference between the two sides'
   printed NPA charges. Its floor is gennbo's own print precision, five decimals.
3. **Per-class principal sine, Cor / Val / Ryd**, measured on gennbo's own basis from its
   `AONAO=W` / `AOPNAO=W` dumps against native's `NAO_DUMP_C=1` columns, each against that
   molecule's own **measured** floor (re-quantise every printed matrix by +/-0.5e-09 over four
   draws, take 3x the spread: 1.8e-08 for lif to 3.95e-05 for benzene). A sine's floor is the
   **square root** of the coefficient floor; a flat tolerance on a quantity that goes through
   a printed inverse is wrong by construction, and was already refuted once in this lane.

A threshold below the printing floor is not a gate. Wherever a number lands at its floor, the
floor is printed beside it and the reading is "at the floor", not "zero".

## The thirteen members, and what each is for

Every member is new: none of the 25 directories in `data/` contains it. Every member is 2-5
atoms on purpose - the point is coverage of the class structure, not of molecular size.

| molecule | chg | mult | what it covers that the eight do not |
|---|---|---|---|
| `sih4` | 0 | 1 | Si: third-row main group whose def2-TZVP d set is a genuinely populated diffuse shell |
| `cl2` | 0 | 1 | two equally polarisable third-row centres - the **inter-atomic** Val<->Ryd channel, where 39-66 % of gennbo's own transformation was found to sit |
| `clf3` | 0 | 1 | hypervalent beyond sf6/pf5: T-shaped, and the central atom **keeps two lone pairs** |
| `bh3` | 0 | 1 | electron-deficient: an empty valence p, so the Val/Ryd boundary has no occupancy to guide it |
| `nh4_plus` | +1 | 1 | a **cation** - the density is contracted and the Rydberg set squeezed |
| `bf4_minus` | -1 | 1 | an **anion** - the excess charge pushes population into the diffuse set |
| `nacl` | 0 | 1 | the third-row **homologue of lif**, which is one of the eight: a row transfer on otherwise identical bonding |
| `mg_h2` | 0 | 1 | group 2, absent from the eight: Mg 3s valence with 3p/3d above it |
| `hs` | 0 | 2 | an open shell with the spin on a **third-row** atom |
| `cf3` | 0 | 2 | an open shell whose spin sits on carbon with three fluorines pulling on it |
| `ch2` | 0 | 3 | a **triplet** carbene, the only non-O2 triplet in this tree |
| `zn_cl2` | 0 | 1 | a transition metal beyond ni_co_4/ticl4: **d10**, so 3d is filled valence and 4d is Rydberg, and def2-TZVP gives Zn several d shells plus an f - the component-major -> shell-major permutation is certainly not the identity here, and that permutation has already produced one defect in this lane |
| `hi` | 0 | 1 | an **ECP-bearing** heavy atom (def2-ECP, 28 core electrons absent from the wavefunction). Its first job is to find out whether gennbo and the `.47` writer handle an ECP at all - the answer is the measurement, and if it cannot be built that is recorded, not dropped |

## The per-member expectation

The mechanism being extrapolated is specific: `renat5` adds the published cascade's Rydberg
re-naturalisation, so it can only move a molecule that **has** Rydberg population to move. The
expectation therefore tracks how diffuse and how populated each member's Rydberg set is, and
it is stated as a factor on `dpop`, not as a direction.

| molecule | expected | why |
|---|---|---|
| `sih4` | **HELP**, >2x | populated Si d set, the class the change acts on |
| `cl2` | **HELP**, >2x | two polarisable centres, large diffuse population |
| `clf3` | **HELP**, >2x | pf5 and sf6 are the two nearest tuned members and both improved 35x and 15x on the sine |
| `nacl` | **HELP**, >2x | lif improved 12x on the sine; this is the same bond one row down. **The cleanest transfer test in the set: if nacl does not improve, the change is second-row-specific.** |
| `bf4_minus` | **HELP**, >2x, and the **largest baseline dpop** of the main-group members | an anion pushes population into exactly the set the change re-naturalises |
| `cf3` | **HELP**, >2x | three F give a large diffuse population; the spin path has no reason to behave differently, and if alpha and beta disagree that is itself a finding |
| `hs` | **HELP**, >2x | third-row open shell |
| `mg_h2` | **HELP**, but the **second most likely member to break it** | Mg's 3s valence is barely occupied, so an occupancy-weighted step has little to weight |
| `nh4_plus` | **FLAT to modest**, <2x | a contracted cation has little Rydberg population to redistribute |
| `ch2` | **FLAT**, <2x | small basis, small Rydberg population |
| `bh3` | **FLAT**, <2x - and this is the **control** | almost no Rydberg population. If bh3's `dpop` moves a lot, the instrument is not measuring what this file claims it measures |
| `zn_cl2` | **DIRECTION UNKNOWN - nominated to break it** | the change's Rydberg step was only ever measured where the Rydberg set is s/p/d polarisation. On a d10 metal the m-averaging runs over five d components of a **filled valence** 3d with 4d above it, which has no counterpart in the eight. This is the member expected to get **worse**, and it is in the set for that reason |
| `hi` | **UNKNOWN; feasibility first** | if it builds at all: the core is nearly untouched by this change in every arm measured, so the core sine should be flat and only Ryd should move |

## The falsifiable claims

- **H1 (transfer).** The worst `dpop` over the held-out members improves by at least **5x**
  between the baseline binary and the changed one. Deliberately weaker than the 63x claimed on
  the eight, so that a partial transfer still counts as a transfer.
- **H2 (no counter-molecule).** No held-out member gets worse than **1.5x** its baseline
  `dpop`. The claim on the eight is "closer on 8 of 8, no counter-molecule"; one held-out
  counter-molecule falsifies that as a general statement, and `zn_cl2` is the nominated
  candidate.
- **H3 (the gate does not move).** No held-out member passes the exact basis-function
  acceptance gate against gennbo, before or after the change. This is **expected**, stated
  here so that nobody reads a failure as news.
- **What would make this a fit rather than an improvement:** a median `dpop` improvement factor
  on the held-out set below **2x** while the eight show 63x. That reading is available and it
  is written down before the measurement, so it cannot be argued away afterwards.

## What the set cannot do

`index.json` is the record of the lost original run and none of these thirteen is in it, so
`accept.py`'s exact basis-function identity test **has nothing to compare against** - the same
weakening allyl, hco and no2 carry. `accept_holdout.py` checks what is available (converged
optimisation, converged last SCF, charge / multiplicity / basis / keyword line as asked,
`<S**2>` within 4 % of S(S+1)) and **records** the basis-function count so the next
regeneration of this set does have the exact test. It can prove the wavefunction is the one
this tree asked for. It cannot prove it reproduces somebody else's earlier wavefunction.

## Addendum, 25 Sep 2026 - written AFTER the baseline was measured, and marked as such

Nothing above is edited: the predictions stand as registered. What this addendum records is
that one number the predictions were *calibrated against* turned out to belong to a different
binary, which changes the scale H1 and H2 are read on and nothing else.

"The shipped native's worst `dpop` is 0.3161 e (benzene 0.4316 against gennbo's 0.1155)" is a
property of the share build `/work/akkleemiss/share/NoSpherA2_RGBI_NBO/NoSpherA2`,
`c8130055 built 2026-09-25 02:20 on AKL007, branch density_source`, which is what produced the
stored native column in `../data/`. Re-running only the native side of the same 22 with the
baseline binary 60055d4b - same stored `.gbw`, same `_native.47`, same flags, same thread count,
so the only variable is which NoSpherA2 ran - gives benzene 0.12055 against gennbo's 0.11544,
a `dpop` of **0.00503 e**, and a mean over the 22 of 0.00030 e rather than 0.01229 e.

Consequences, and only these:

- The factor 63 is a ratio against that share build, not against the head of this branch.
- H1's "5x improvement" and H2's "1.5x worse" are ratios, so they survive unchanged - but they
  must be taken on one binary with an environment switch, never on two builds, and the
  before-number must come from 60055d4b.
- The held-out set's own baseline sits at 0.00002 to 0.00607 e, which is one to two orders
  below the 0.3161 e the set was sized against. The "fit rather than improvement" reading above
  is therefore harder to reach than intended on this set: a 5x improvement on 0.006 e is close
  to the 1e-5 print floor. That is a weakness of the set, recorded rather than repaired.

See `README.md`, section "The 0.3161 e baseline figure does not belong to this branch", for the
excluded candidates and the positive control that makes the exclusions readable.

## Addendum 2, 25 Sep 2026 - the first addendum was wrong, and the claims are now decided

Addendum 1 above said the 0.3161 e calibration figure belonged to the share build `c8130055`
and not to this branch. **That is withdrawn.** The binary it rested on is byte-identical to the
renat5 treatment build (md5 `1e45bc35f8d5e9617ee42dc62ca7d180`, checked by `md5sum` against the
preserved pre-change binary `99a8aa34...`, which is 176 bytes smaller), so benzene 0.12055 was
the *treatment* value all along and agreed with the sibling's renat5 arm rather than
contradicting the 0.3161 e before-value. Both addenda stay in place: an addendum that is itself
corrected is part of the record.

The 0.3161 e figure therefore stands, and so does the scale the predictions were written on.
Held-out verdicts, both arms from one binary and one environment variable (job 595959):

- **H1 HOLDS.** Worst `d(Ryd)` 0.19334 e baseline -> 0.00614 e renat5, factor **31.5** against
  the 5x claimed. Median per-molecule factor 18.0.
- **H2 FAILS on two members**: `cf3` (2.08x worse, both values far above its 2.0e-05 floor) and
  `cl2` (2.4x worse, 10 floors against 24). "Closer on 8 of 8, no counter-molecule" does not
  survive held-out data.
- **The nominated counter-molecule was the wrong one.** `zn_cl2` was predicted to get worse and
  improves 9.1x. `bh3` was pre-registered as the control - "if bh3's dpop moves a lot, the
  instrument is not measuring what this file claims" - and it moved 44x. The reading taken is
  the narrower one: `bh3` does have Rydberg population to move (0.03707 e baseline), so the
  control was mis-specified rather than the instrument being wrong. That is a weakened control
  and it is recorded as such, not repaired after the fact.
- **H3 HOLDS**: 0 PASS / 16 FAIL in both arms. Nothing here agrees with NBO 7, and a 31.5x
  improvement is not an agreement.
- The **fit-rather-than-improvement** reading required a median below 2x. It is 18.0, so that
  reading is not available - a result of the set, not a concession to it.
