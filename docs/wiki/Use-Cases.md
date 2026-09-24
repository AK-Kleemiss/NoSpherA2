# Use cases

One PowerShell block per job, each runnable as it stands once the file names are
yours. `.\NoSpherA2.exe` assumes the binary is in the working directory; put it
on `PATH` and drop the `.\`.

Remember that **the output goes to `NoSpherA2.log`**, not to the console. Each
recipe below ends with the `Get-Content` line that actually shows you the
result.

The invocations are modelled on the regression suite in `tests/tests.toml`,
which is the authoritative set of command lines that are known to work.

---

## 1. Aspherical form factors for a refinement

The core job: an ORCA `.gbw` (or a Gaussian `.fchk`), a CIF and a reflection
list in, a binary `.tscb` out. Default Hirshfeld partitioning, accuracy 2.

```powershell
.\NoSpherA2.exe -wfn sucrose.gbw -cif sucrose.cif -hkl sucrose.hkl -acc 2 -mult 1 -charge 0 -all_charges
```

```powershell
Get-Content .\NoSpherA2.log -Tail 30
```

The table is written next to the CIF. Point `olex2.refine` at it, or set the
NoSpherA2 tab in Olex2 to use it.

## 2. The spherical reference for the same structure

Run the identical job with IAM form factors, so that any improvement you claim
is measured against the same reflection list and the same code path.

```powershell
.\NoSpherA2.exe -wfn sucrose.gbw -cif sucrose.cif -hkl sucrose.hkl -acc 2 -mult 1 -charge 0 -IAM
```

## 3. The same table with a different partitioning

Becke, TFVC, MBIS and EMBIS are one flag each. This one also asks for a
generated reflection sphere to 1.5 A instead of reading the `.hkl`, writes text
`.tsc`, and caps the streaming block size.

```powershell
.\NoSpherA2.exe -wfn epoxide.gbw -cif epoxide.cif -dmin 1.5 -acc 1 -mult 1 -charge 0 -tfvc -old_tsc -tsc_block 50
```

Swap `-tfvc` for `-becke`, `-mbis`, `-embis` or `-hirsh`. If you generate the
reflection list, give **both** `-dmin` and `-hkl_min_max` unless your dataset is
a complete sphere — see
[Troubleshooting](Troubleshooting.md#dmin-alone-generates-a-full-sphere).

## 4. A heavy-atom structure with ECPs

`-ECP 1` selects the ORCA def2-ECP treatment, which is what you want when the
wavefunction was computed with those ECPs.

```powershell
.\NoSpherA2.exe -wfn malbac.gbw -cif malbac.cif -hkl malbac.hkl -acc 2 -mult 1 -charge 0 -ECP 1
```

## 5. A disordered structure: one wavefunction per part

`-mtc` takes alternating pairs of its own — a wavefunction, then the disorder
parts it covers as one token. It is **not** the separate `-group` flag, and the
part token is a list: `0.1` means parts 0 and 1, part 0 being the ordered atoms
that carry no `_atom_site_disorder_group` in the CIF. So both fragments include
the backbone, and each adds its own part:

```powershell
.\NoSpherA2.exe -cif thpp.cif -hkl thpp.hkl -acc 1 -mtc part_1/thpp.wfx 0.1 part_2/thpp.wfx 0.2 -mtc_mult 1 1 -mtc_charge 0 0 -mtc_ECP 0 0
```

`-cmtc` is the same with a CIF per fragment — wavefunction, CIF, parts — for
fragments that come with their own asymmetric unit:

```powershell
.\NoSpherA2.exe -hkl 1yk4_h.hkl -acc 1 -cmtc residues/1.gbw residues/1.cif 0 residues/2.gbw residues/2.cif 0.1 residues/3.gbw residues/3.cif 0.2
```

Both keep consuming tokens until one begins with `-`, so put `-mtc_mult`,
`-mtc_charge` and `-mtc_ECP` after the whole list, in fragment order. Check the
scatterer count in the log against the number of atoms in the CIF before you
refine against the table.

## 6. Convert a wavefunction between formats

ORCA `.gbw` to `.wfn`:

```powershell
.\NoSpherA2.exe -gbw2wfn mol.gbw
```

A `.ffn` plus a named basis set to a Gaussian checkpoint:

```powershell
.\NoSpherA2.exe -wfn in.ffn -b def2-TZVP -d ./ -fchk log.fchk
```

Wavefunction geometry and basis into a CIF, for archiving alongside the
structure:

```powershell
.\NoSpherA2.exe -gbw2wfn mol.gbw -wfn_cif mol_wfn.cif
```

## 7. Convert or merge `.tsc` tables

Binary `.tscb` to text `.tsc` and back — the direction follows the input's
extension:

```powershell
.\NoSpherA2.exe -tscb experimental.tscb
```

Merge the tables of several fragments into one, checking that the reflection
sets agree:

```powershell
.\NoSpherA2.exe -merge part1.tsc part2.tsc
```

Rewrite an existing table's atom IDs into atom labels using the CIF:

```powershell
.\NoSpherA2.exe -tsc_labels experimental.tscb mol.cif labelled.tsc
```

## 8. Density, Laplacian, ELI-D and RDG cubes in one run

Every property is an independent flag, so ask for all of them at once and pay
for the grid once. Keep the resolution coarse while you are finding your feet: a
0.1 A grid on a large molecule is a large file.

```powershell
.\NoSpherA2.exe -wfn epoxide.gbw -resolution 0.5 -radius 2.0 -rho -lap -eli -rdg -DEF -hirsh
```

```powershell
Get-Content .\NoSpherA2.log -Tail 20
```

Add `-cubeb` to write the compact binary `.cubeb` format instead of text, and
convert back when you need to look at one:

```powershell
.\NoSpherA2.exe -cube_convert epoxide_rho.cubeb epoxide_rho.cube
```

## 9. Critical points of the density

Console output, no cube, no grid to choose: every nuclear, bond, ring and cage
critical point from the analytic Hessian, with the Poincare-Hopf relation as the
completeness check.

```powershell
.\NoSpherA2.exe -topology epoxide.gbw
```

## 10. QTAIM and ELI-D basins, charges and delocalisation indices

The integration job. Arguments are wavefunction, resolution and radius, in that
order. 0.05 A is the practical minimum resolution — a coarser grid invents
non-nuclear attractors at C-C bond critical points.

```powershell
.\NoSpherA2.exe -acc 4 -eli_analysis nh3li.gbw 0.05 3.0
```

```powershell
Get-Content .\NoSpherA2.log -Tail 60
```

`-acc` has to come **before** `-eli_analysis`. For a structure with diffuse
density, raise the radius and check the integrated electron count in the log
before you trust the charges.

The rest of Kohout's ELI family as point values, printed to the console:

```powershell
.\NoSpherA2.exe -eli_family nh3li.gbw
```

And an ELI-D field masked to the QTAIM basins of selected atoms (0-based
indices):

```powershell
.\NoSpherA2.exe -qtaim_eli nh3li.gbw 0,1,2,3 0.05 3.0
```

## 11. Intermolecular NCI between two fragments

Promolecular densities, so no wavefunction is needed — just the two geometries.
The trailing numbers are the optional cut-offs in the order `rcut1 rcut2
rho_abs_max rdg_max colour_max`.

```powershell
.\NoSpherA2.exe -promol_nci unit.xyz surrounding.xyz
```

## 12. ESP-coloured Hirshfeld surface for Olex2

The surface of the asymmetric unit against its packing environment. Olex2 builds
the two `.xyz` files for you; standalone, export them from the structure.

```powershell
.\NoSpherA2.exe -hirshfeld_surface asu.xyz pack.xyz -wfn mol.gbw -resolution 0.2 -esp
```

An ESP mapped onto a density isosurface instead, at the default 0.002 e/bohr^3:

```powershell
.\NoSpherA2.exe -wfn epoxide.gbw -resolution 0.5 -radius 2.0 -esp_isosurface
```

Both write an `.obj` mesh plus a `.dat` table. The triangle order in those files
is not reproducible between runs (marching cubes collects in parallel), so
compare sorted tables rather than diffing files.

## 13. NBO, E2 and NRT on a small molecule

Run entirely inside NoSpherA2 — no `gennbo` installation. Every sub-option must
come **after** `-nbo_native`.

```powershell
.\NoSpherA2.exe -nbo_native nico4.gbw -nrt -nrt_arrows -nbo_e2min 0.5 -nbo_threads 12 -nbo_json nico4_nbo.json
```

```powershell
Get-Content .\NoSpherA2.log -Tail 80
```

To use an external NBO 7 instead:

```powershell
.\NoSpherA2.exe -nbo nico4.gbw -nbo_exe C:\NBO7\bin\gennbo.exe -nbo_keywords "NRT NRTSYM=off" -nbo_keep47
```

Or just write the archive and stop:

```powershell
.\NoSpherA2.exe -convert_to_47 nico4.gbw -nbo_47 nico4.47
```

The input must carry a contracted basis: `.gbw`, `.fchk` or `.molden`. A `.wfn`
or `.wfx` is rejected.

## 14. Natural population analysis and resonance-group bond indices

```powershell
.\NoSpherA2.exe -wfn nh3bh3.gbw -npa
```

```powershell
.\NoSpherA2.exe -wfn nh3bh3.gbw -rgbi -rgbi_groups "0,4,5,7" "1,2,3,6" -rgbi_basis ano
```

## 15. RI-fitted density, with multipole restraints

Fit the density in an auxiliary basis and use the fit for the partitioning. The
restraints hold the fit's atomic multipole moments up to `lmax` at the values a
real partitioning gives.

```powershell
.\NoSpherA2.exe -wfn epoxide.gbw -cif epoxide.cif -hkl epoxide.hkl -acc 2 -ri_fit def2-universal-jkfit -multipole_moments Hirshfeld 2 -multipole_strength 0.5
```

Let the auxiliary set be generated, but pin hydrogen to a proper fitting basis —
this is what fixes hydroxyl-H populations, at roughly twice the wall time:

```powershell
.\NoSpherA2.exe -wfn sucrose.gbw -cif sucrose.cif -hkl sucrose.hkl -acc 2 -ri_fit auto_aux H def2-universal-jkfit
```

## 16. A machine-learned (SALTED) density

```powershell
.\NoSpherA2.exe -wfn mol.xyz -SALTED .\salted_model -salted_charge_constraint -esp -resolution 0.2
```

Always pass `-salted_charge_constraint`: without it the prediction is short by
about 0.4 e, which shifts the whole ESP by more than its range on the surface.

Coefficients only, no field. `-SALTED_COEFS` runs inside the parser, so
everything it needs — including `-no_gpu_salted` if you want the CPU path — has
to come first:

```powershell
.\NoSpherA2.exe -wfn mol.xyz -SALTED .\salted_model -salted_charge_constraint -SALTED_COEFS
```

## 17. Electrostatic interaction energy between two molecules

Both densities are RI-fitted (or SALTED-predicted, if `-SALTED` is given), then
their electrostatic interaction is evaluated. `-interaction_energy` runs inside
the parser and stops it, so it goes **last**.

```powershell
.\NoSpherA2.exe -ri_fit def2-universal-jkfit -repulsion_exchange pbe -repulsion_overlap 0.1 -interaction_energy A.gbw B.gbw
```

Every contacting pair in a crystal, from a job file:

```powershell
.\NoSpherA2.exe -ri_fit def2-universal-jkfit -interaction_energies pairs.job
```

## 18. X-ray constrained wavefunction fitting

The lambda scan. `-no_gpu` is not required, but the regression tests pin it
because the device I tensor differs from the host one in the tenth significant
figure.

```powershell
.\NoSpherA2.exe -wfn mol.gbw -cif mol.cif -hkl mol.hkl -acc 2 -mult 1 -charge 0 -do_XCW -XCW_settings xcw.txt -anom_disp mol.disp -no_gpu
```

```powershell
Get-Content .\NoSpherA2.log -Tail 60
```

The settings file controls where the two-electron integral tensor lives
(`stream`, `i_tensor_mb`, `save`, `read`, `df_basis`, `i_double`, `i_float`),
which is the whole memory question for a fit of any size.

## 19. Force the GPU, and prove it ran

The GPU paths are on by default and fall back to the CPU silently, so a run that
looks slow tells you nothing. `-no_date_but_gpu` keeps the "GPU in use" notes
while still suppressing timings, and `-gflops` reports the achieved rate.

```powershell
.\NoSpherA2.exe -wfn sucrose.gbw -cif sucrose.cif -hkl sucrose.hkl -acc 2 -mult 1 -charge 0 -gpu_grid -gpu_fp64 -gflops -no_date_but_gpu
```

```powershell
Select-String -Path .\NoSpherA2.log -Pattern "GPU","GFLOP"
```

Use `-gpu_fp64` when you need the CPU and GPU results to agree to round-off;
leave it off for production speed.

## 20. Run entirely on the CPU

For a reference number, a machine with no usable device, or to rule the GPU out
of a discrepancy.

```powershell
.\NoSpherA2.exe -wfn sucrose.gbw -cif sucrose.cif -hkl sucrose.hkl -acc 2 -mult 1 -charge 0 -no_gpu -cpus 8
```

## 21. Density at a list of points

For feeding another program. The points file holds the count on the first line
then `x y z` per line, **in bohr**; the result lands in `<points>.rho`.

```powershell
.\NoSpherA2.exe -wfn mol.gbw -rho_at_points grid.txt
```

`-wfn` has to precede `-rho_at_points`, which runs inside the parser.

## 22. Dipole moments, polarisabilities, Fukui functions

```powershell
.\NoSpherA2.exe -wfn mol.gbw -dipole_moments
```

Polarisability from seven finite-field wavefunctions, in the order the flag
expects:

```powershell
.\NoSpherA2.exe -e_field 0.005 -polarizabilities f0.gbw fxp.gbw fxm.gbw fyp.gbw fym.gbw fzp.gbw fzm.gbw
```

Condensed Fukui indices — this one prints to the console, and needs a
wavefunction that still has its virtual orbitals (a `.fchk`, not a `.wfx`
trimmed to the occupied set):

```powershell
.\NoSpherA2.exe -fukui_analysis mol.fchk
```
