# Troubleshooting and gotchas

Everything here is a guard the code actually raises, or a trap that has cost
somebody a day. Ordered roughly by how often it bites.

## The run printed nothing

It printed into `NoSpherA2.log` in the working directory. `std::cout` is
redirected there before the options are even parsed.

```powershell
Get-Content .\NoSpherA2.log -Tail 40
```

`-out <file>` renames the log. Four commands deliberately restore the console
first, because their entire output is meant to be read: `-h`, `-topology`,
`-eli_family`, `-fukui_analysis`.

## "Did not understand the task to perform!"

The parser accepted every flag but nothing on the line selected a job, so
NoSpherA2 fell through to printing the help text. The two common causes:

* The job flag was consumed as the **argument** of an earlier flag that takes
  one more value than you thought.
* The binary is **older than the flag**. A flag added after your build was
  linked is just an unrecognised token. Check with the binary's own help:

```powershell
.\NoSpherA2.exe -h > help.txt ; Select-String -Path help.txt -Pattern "-topology"
```

If the flag is not in that output, rebuild before you debug the command line.

## Flag order matters

Most flags only set a field and can go anywhere. But a number of them **do their
work inside the option parser** and then stop parsing, so anything written after
them is read by nobody. There is no warning.

**One-shot options** (they finish the run): `-topology`, `-eli_family`,
`-qtaim_eli`, `-merge`, `-merge_nocheck`, `-tscb`, `-tsc_labels <table> <cif>`,
`-cube_convert`, `-rho_cube`, `-rho_at_points`, `-calc_dens_1D`,
`-atom_dens`, `-atom_dens_diff`, `-combine_mos`, `-dipole_moments`,
`-laplacian_bonds`, `-ewal_sum`, `-draw_orbits`, `-spherical_harmonic`,
`-spherical_atoms`, `-write_ri_coefs`, `-SALTED_COEFS`, `-SALTED_Training`,
`-interaction_energy`, `-convert_XCW`, `-sfac_diffuse`.

The rules that follow from that:

| Put this first | Before this | Otherwise |
| --- | --- | --- |
| `-wfn <file>` | `-rho_at_points`, `-calc_dens_1D`, `-SALTED_COEFS`, `-SALTED_Training` | hard error, "No wavefunction specified" |
| `-ri_fit <basis>` | `-calc_dens_1D`, `-SALTED_Training`, `-interaction_energy` | hard error or an unfitted density |
| `-acc <n>` | `-eli_analysis`, `-qtaim_eli` | the default accuracy 2 is used |
| `-radius`, `-resolution` | `-rho_cube` | the defaults 2.0 / 0.1 are used |
| `-no_gpu_salted` | `-SALTED_COEFS` | the GPU path runs anyway (the globals are set later) |
| `-old_tsc`, `-tsc_block` | `-merge`, `-tscb` | silently ignored |
| `-repulsion_exchange`, `-repulsion_overlap`, `-no_gpu_density` | `-interaction_energy` | never read |
| `-nbo_native` / `-nbo` / `-nbo_parse` / `-convert_to_47` | every `-nbo_*` and `-nrt*` sub-option | the sub-options are silently ignored |

The NBO case is the one that looks most wrong: the sub-options are collected by a
loop that starts at the position of the parent flag, so `-nrt -nbo_native x.gbw`
runs an NBO analysis with no NRT and says nothing about it. Always write the
parent flag first.

`-mtc` and `-cmtc` are the mirror image: they swallow tokens *forwards* until one
begins with `-`, in pairs (`wfn parts`) or triples (`wfn cif parts`). Two
consequences. `-mtc a.gbw -group 1 2 3` does not do what it reads like — there is
no `-group` in this syntax, `-group` terminates the list, and the atom indices
after it are then taken as further wavefunction names. And `-mtc_mult`,
`-mtc_charge`, `-mtc_ECP` must come **after** the complete `-mtc` list, in
fragment order.

## "Option X needs more values than were given after it"

Exactly what it says: the flag takes more arguments than are left on the line, or
the next token starts with `-` where a value was expected. The wrapper prints the
offending flag by name. `-twin` wants nine numbers, `-hkl_min_max` six,
`-sfac_diffuse` six values and two file names, `-polarizabilities` seven files.

## `-dmin` alone generates a full sphere

`-dmin` asks for every reflection in a sphere to that resolution. On a
**complete** dataset that is about twice what you need; on a **partial** dataset
it is catastrophic — one Fe test generated 286,998 rows where the combination
with `-hkl_min_max` gives 308.

Give both, and let the intersection do the work. The sphere carries a deliberate
1e-3 margin, so it is always a superset of the cctbx list and never short.

Fine resolutions are expensive, not broken: `-dmin 0.05` on sucrose is 12.06
million reflections and an 8.8 GB `.tscb`, with 54 s spent just building the
list.

With `-ED`, `-dmin` becomes `dmin/2 - 0.001` and the index box is ignored — every
beam is needed.

## The `.tsc` streaming block size does not buy speed

`-tsc_block` only trades memory. The table hash is identical at every block size
tested and there is no time trend. Peak memory is roughly
`queue_depth x scatterers x block x 16 B`. `-tsc_block 0` means "hold the whole
table": 3.5 GB for a medium protein, 40 GB for 8566 atoms against 293k
reflections.

## Disorder tables are keyed by atom ID

Rows in a multi-fragment (`-mtc`) table are keyed by atom ID in hex, not by
label, because atom labels are **not unique across disorder parts**. A
label-keyed table silently loses or misattributes rows, and the result is the
wrong answer at very nearly the right size — one case produced 1025 scatterers
where the structure has 1026. Check the scatterer count in the log against the
CIF; only a byte comparison against a non-streamed run catches the subtler
version.

Related: `-tsc_labels` on an existing table **refuses** any table without a
`SCATTERER_IDS` block, and older builds wrote a hardcoded `TITLE test` for
`-tscb`. Do not infer what produced a table from its header.

## Molden files

The `[GTO]` normalisation convention is assumed to be the **bare** one that ORCA
and `orca_2mkl` write, unconditionally and with no warning. A genuine
`normprim` molden is misread by orders of magnitude. Auto-detection is not
possible: no file in any corpus carries a `program=` keyword and `[Title]` is not
evidence, since ORCA leaves it empty.

**If the molden was not written by ORCA, feed a `.gbw` or `.fchk` instead.**

For the same reason, do not treat DGrid as ground truth for a molden file: DGrid
guesses the convention from `[Title]`, so an ORCA file with an empty title is
read as `normprim` and disagrees by a factor of several hundred at a fluorine
nucleus — with DGrid wrong. DGrid also silently converts an *unrestricted* ORCA
molden into a restricted wavefunction, dropping the beta orbitals, while its
header's per-spin electron counts still look right.

Spherical molden `p` shells were read in the wrong order (z,x,y instead of
x,y,z) in builds from mid-September 2026 until the fix, which gives a 1-5 %
density error near nuclei and mixes in-plane with out-of-plane `p`. If you are on
an old binary, cross-check against the `.gbw`.

## Cube grids

**`-resolution 0.1` is too coarse for QTAIM.** It invents non-nuclear attractors
at C-C bond critical points. 0.05 A is the real minimum.

**`-radius` hard-zeroes every point further than that from any atom.** A diffuse
density therefore loses electrons to the edge of the box, and it is not a bug: an
Yb/vDZP case integrates to 19.58 of 24 electrons at a 10 A radius and needs
about 30 A to reach 23.999. **Check the integrated electron count in the log
before you believe a cube.**

Box-mode (bondwise, basin) cube geometry is not trustworthy in detail: the
requested spacing is not reproduced exactly (`incr = size/np`), with `bohr` the
origin is converted but the increment stays in Angstrom, and mode 1 centres the
box on the second atom where the documentation says the first.

`-def` alone used to write plain rho as the "static deformation density" unless
`-hirsh` or `-HDEF` was also given, and `-s_rho` alone crashed. Both are fixed;
on an older binary add `-hirsh` and avoid a bare `-s_rho`.

Off-centre property and Hirshfeld grids in builds before 18 Sep 2026 scaled the
box by 1.89 about the origin for any molecule whose shortest bond exceeds
1.058 A — i.e. anything without an O-H or N-H — producing clipped, shifted
surfaces. Update the binary rather than fiddling with the box.

`-hdef` in lowercase used to be rejected while only `-HDEF` worked, which is what
older help texts disagreed about. Current builds accept both.

## Surfaces are not reproducible file-by-file

The `.obj` and `Hirshfeld_surface.dat` face order changes between runs because
marching cubes collects triangles in parallel. Compare sorted tables, not files.
Also, after the `Time to calculate Values` log line the stream is stuck at zero
decimals, so an isosurface value of `0.002` prints as `0`.

## GPU and CPU do not agree bitwise

Never expect bit-identity; expect these magnitudes (measured, on `wR2` or the
equivalent derived quantity):

| Path | Agreement with the CPU |
| --- | --- |
| `calc_SF` with `-gpu_fp64` | 1.2e-14 |
| `calc_SF` default (fp32 transcendentals and sum) | 1.8e-8 |
| I tensor with `-gpu_itensor` (fp32 GEMM) | 2.1e-11 |
| `-gpu_salted` | 5.0e-13 |

All of those are orders of magnitude below the ~1 % experimental uncertainty, and
the double-precision alternative is one flag away in each case. Judge agreement
on a **sensitive derived** quantity, not on a total: in one XCW run the total
energy agreed to 1.5e-12 while the halting statistic moved by 5.3e-7.

A silent CPU fallback is indistinguishable from a correct GPU result in the
numbers, so if it matters whether the device ran, use `-no_date_but_gpu` (plain
`-no_date` suppresses the "GPU in use" line as well as the timings) and check the
log.

## A fitted or predicted density has no orbitals

On the `-ri_fit` and `-SALTED` paths:

* **ELF exits** with "ELF needs orbitals".
* `-MO`, `-rdg` and `-elf` stay zero.
* The **Laplacian is exact** for the fitted density.
* **ELI-D is an orbital-free PC07 surrogate**: basin *positions* are right for
  rho > 0.01, but peaks are capped near 3.4 and the tail is lost.
* **QTAIM charges** from a fitted density are good to about 0.01 e.

## SALTED: pass the charge constraint

Without `-salted_charge_constraint` the prediction is short by about 0.42 e,
which shifts the whole ESP by roughly 0.08 au — larger than the ESP range on the
surface itself. There is no reason not to pass it.

## RI fitting

`-ri_fit auto_aux <elements...> <basis...>`: element symbols placed directly
after `auto_aux` restrict which elements the set is generated for, and the basis
names that follow fill the rest, first match winning. `auto_aux H
def2-universal-jkfit` is what fixes hydroxyl-hydrogen populations (0.13 e down to
under 0.013 e) at about twice the wall time.

`-multipole_centre` with `lmax >= 2` breaks the fit. Use partitioned moments —
i.e. the default `-multipole_partition`.

## ELI family: two of the four members carry no new information

* **ELI-q is exactly ELI-D^(-8/3)**: identical critical points, inverted order.
  Its gradient ascents shatter into 60-100 tail basins, so do not run basin
  analysis on it.
* **ELIA is identically `1 - zeta^2`** for any single determinant, hence
  identically 1 for a closed shell. It is deliberately excluded from basin
  analysis.
* **Triplet ELI-D** is informative only for genuinely open-shell wavefunctions,
  and its prefactor follows DGrid rather than Kohout's Eq. 5, so a comparison
  against the paper shows a constant factor.

If you generate DGrid reference fields: `ELI-D triplet` is a user-input error
(bare `triplet` means the density); the correct spellings are
`ELI-D triplet-pair` and `ELI-D beta-beta`. Never ask DGrid for `beta-beta` on a
restricted wavefunction — it writes a plausible constant artefact (about
7.5e-13 in every voxel) instead of refusing.

## NBO and NRT

* **`.wfn` and `.wfx` are rejected** (exit 255): they are primitive-only and the
  analysis needs the contracted basis. Use `.gbw`, `.fchk` or `.molden`. Shells
  are supported to `g`.
* Open-shell archives written by older builds were analysed as closed-shell —
  plausible and silently wrong. Check that the output has `Alpha spin` and
  `Beta spin` sections.
* With **external NBO 7.0.9**: `NRTSUB=` does not exist (use `NRT <atom list>`
  or `$NRTSTR`); `$CHOOSE` and `$NRTSTR` must go at the **end** of the `.47`, or
  you get `'END' is not an acceptable orbital type`; symmetry detection on
  SF6-like cases can abort the whole analysis silently, which `NRTSYM=off` fixes;
  transition-metal complexes need `NRTE2=5...20` to finish at all — and that
  changes the answer, not just the precision.
* Compare NRT weights **by rank, never by structure description**, and always
  quote the threshold alongside the number.
* **`-cpus` does not reach this code.** The native search reads `-nbo_threads`,
  and falls back to `omp_get_max_threads()` when it is absent — so `-cpus 1`
  still runs on every core, and a timing taken that way is not a serial timing.

## The tsc NaN guard was dead in Release builds

`x != x` was folded away under `/fp:fast`, so a table containing NaNs could be
written without complaint before the `std::isnan` fix. If you are looking at an
old table that refines pathologically, check for NaNs directly.

## Adding `-cif` makes cubes much slower

`-cif` switches cube evaluation onto the periodic path — the one Olex2 always
drives. Roughly 27x the cost of the molecular path on the same grid. That is
expected, not a regression.
