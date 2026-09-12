# VMD visualisation helpers

## Fukui functions and the dual descriptor

Run NoSpherA2 with `-fukui` to get the cubes, then generate a VMD script for them:

```sh
NoSpherA2 -wfn mol.fchk -fukui -resolution 0.2 -radius 3.5
python scripts/vmd/make_fukui_vmd.py .
vmd -e mol_fukui.vmd
```

`-fukui` writes four cubes and a summary beside the wavefunction:

| File | Quantity |
|---|---|
| `<stem>_fukui_plus.cube` | f+(r) = \|psi_LUMO\|^2 |
| `<stem>_fukui_minus.cube` | f-(r) = \|psi_HOMO\|^2 |
| `<stem>_fukui_zero.cube` | f0(r) = (f+ + f-) / 2 |
| `<stem>_dual_descriptor.cube` | df(r) = f+ - f- |
| `<stem>_fukui.dat` | frontier orbital indices/energies and the integrated norms |

## If VMD says "Cannot read file of type ..."

```
ERROR) Cannot read file of type r
ERROR) Loading of startup molecule files aborted.
```

That is a **mistyped flag**, not a problem with the file. VMD treats any option it
does not recognise as a *file-type specifier* — the `vmd -pdb foo.pdb` /
`vmd -cube foo.cube` syntax — so `-r` becomes "a file of type r" rather than an
"unknown option" error. The letter in the message is the flag you typed. The script
flag is `-e`.

To check a generated script without opening a window:

```sh
vmd -dispdev text -eofexit -e mol_fukui.vmd
```

A working run prints `Analyzing Volume...` once per cube (four times for a full set).
If you see no `Analyzing Volume` lines, VMD never executed the script.

As a fallback that bypasses command-line parsing entirely, start VMD normally and in
the Tk Console (Extensions > Tk Console) run:

```tcl
source D:/path/to/mol_fukui.vmd
```

## Reading the pictures

The naming is the single easiest thing to get backwards:

- **f+** is large where the molecule **accepts** density, i.e. where it is attacked by a
  **nucleophile**. Those are the **electrophilic** sites.
- **f-** is large where the molecule **gives density up**, i.e. where it is attacked by an
  **electrophile**. Those are the **nucleophilic** sites.
- The **dual descriptor** folds both into one signed field and is usually the more
  readable picture: positive (blue) electrophilic, negative (red) nucleophilic.

## Isovalues

The generator reads the actual data range out of each cube and picks a default isovalue
from it, because an isovalue above the data range renders an empty screen and that is by
far the most common way one of these plots appears "not to work". Carrying a value over
from a `rho` or `ELI` plot will reliably show nothing: these frontier densities peak
around 0.01-0.1.

The default is 30% of the **smallest** peak among f+/f-/f0, and it is shared by all three
surfaces. That is deliberate on both counts. f+ and f- routinely differ several-fold in
height (the LUMO is usually the more diffuse orbital and so the flatter one), so scaling
to the largest would put the isovalue above the f+ maximum and hide it entirely; and f+
and f- are only visually comparable when they are drawn at the same level, so per-cube
autoscaling would make a diffuse orbital look as concentrated as a tight one.

Override with `--iso` and `--dual-iso`. The script warns if a value you pass lies outside
a cube's range.

## Just analysing a wavefunction — no cubes

If all you want is *where does this molecule react*, skip the grid entirely:

```sh
NoSpherA2 -fukui_analysis mol.gbw
```

No resolution, no radius, no CIF. It reports the frontier orbitals and the HOMO-LUMO
gap, prints the condensed Fukui table for all five partitions, names the most
electrophilic and most nucleophilic atom, and writes `mol_fukui.dat`. Takes a fraction of
a second for a small molecule — the cost is one Becke-grid integration, not a cube.

```
Conceptual-DFT reactivity analysis of epoxide.gbw
Read 7 atoms and 62 molecular orbitals (12 occupied).
HOMO = MO 11  energy -0.230142  occupation 2.000
LUMO = MO 12  energy 0.075122  occupation 0.000
...
Summary (Hirshfeld partition):
  Most ELECTROPHILIC site (a nucleophile attacks here): C   df = 0.24873
  Most NUCLEOPHILIC site (an electrophile attacks here): O   df = -0.54337
```

The wavefunction may also be given with `-wfn`. If it has no virtual orbitals the command
says so and names the formats that do, rather than emitting zeros.

Use `-fukui` instead when you want the cubes as well.

## Condensed (atom-summed) Fukui functions

`-fukui` also prints a per-atom table to `NoSpherA2_cube.log` and writes it into
`<stem>_fukui.dat`:

```
f+_A = integral of w_A(r) |psi_LUMO(r)|^2      A is attacked by a nucleophile
f-_A = integral of w_A(r) |psi_HOMO(r)|^2      A is attacked by an electrophile
df_A = f+_A - f-_A                             > 0 electrophilic, < 0 nucleophilic
```

**Prefer these numbers to the isosurfaces when deciding where a molecule reacts.**
The point-wise Fukui function is dominated by near-nuclear density, so its maxima track
the heaviest atom rather than the reactive site — on ethylene oxide both f+ and f- peak
on the oxygen, which says nothing useful. The condensed values give the right answer
immediately: f- sits on O (0.76), f+ on the two carbons (0.28 + 0.32), so the dual
descriptor is strongly negative on O and positive on C. That is textbook epoxide
reactivity — nucleophilic ring-opening at carbon, oxygen as the basic site.

Five columns are printed, one per atomic partition: **Hirshfeld, Becke, TFVC, MBIS,
EMBIS**. They come from one integration pass, because NoSpherA2 already carries all five
weight sets on the same grid. Hirshfeld is the scheme the conceptual-DFT literature uses
(De Proft et al., *J. Comput. Chem.* **2002**, 23, 1198); the others are not otherwise
available for Fukui functions anywhere.

They do not agree, and the disagreement is the interesting part. On ethylene oxide the
carbon f- is 0.069 under Hirshfeld but 0.024 under TFVC — nearly a factor of three. The
condensed Fukui function is only defined relative to a choice of w_A, and which choice
behaves best is an open question that these columns let you measure rather than argue.

The **sum over atoms of each column is 1** for both f+ and f-, and 0 for the dual
descriptor. That is an exact result, so it is a free correctness check on every run.

## Sanity check before believing a picture

`<stem>_fukui.dat` reports the integral of f+ and f- over the evaluated grid. Each equals
**1** for the exact functions. A markedly smaller value means the grid or `-radius` is
clipping the orbital, not that the chemistry is odd — increase `-radius` first, then
refine `-resolution`. On sucrose the integrals go 1.0045 / 0.9327 at 0.8 A, to
0.9966 / 0.9994 at 0.2 A.

## Caveats

- The frontier-orbital approximation assumes a well-separated, non-degenerate HOMO and
  LUMO. It is unreliable for degenerate or near-degenerate frontier orbitals, which is
  common for transition-metal complexes; orbital-weighted variants exist for that case.
- Wavefunction files that store only occupied orbitals — plain `.wfn`, and some `.wfx` —
  have no LUMO, so no Fukui function can be formed. NoSpherA2 reports this and skips
  rather than emitting a zero cube. Use an `.fchk`, `.gbw` or `.molden` instead.
- **Unrestricted wavefunctions.** The frontier pair is the globally highest occupied and
  globally lowest virtual *spin orbital*, ignoring spin — i.e. the electron that is
  easiest to remove and the level an added electron would actually enter. That is the
  correct **spin-unresolved** Fukui function, but it is not a spin-polarised
  (f_N / f_S) treatment in the sense of Galván/Vela/Gázquez.

  In practice both members of the pair usually come from the *same* manifold, which makes
  the result effectively a one-spin Fukui function. NoSpherA2 now prints which manifold
  won and each manifold's own frontier, so this is visible rather than implicit:

  ```
  Chosen HOMO is beta, chosen LUMO is beta.
  alpha: HOMO = MO 65 (-0.129327), LUMO = MO 66 (0.021202)
   beta: HOMO = MO 553 (-0.103838), LUMO = MO 554 (-0.060801)
  Both come from the beta manifold, so this is effectively a beta-only Fukui function.
  ```

  One pathology to know: for a **high-spin** system the lowest virtual is often the
  spatial partner of a singly-occupied orbital, so f+ and f- describe nearly the same
  region and the dual descriptor collapses toward zero. That is an artefact of the
  frozen-orbital approximation, not a statement about the molecule.
