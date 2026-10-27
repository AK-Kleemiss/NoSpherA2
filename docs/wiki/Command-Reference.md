# Command reference

Every flag below was read out of the option parser in
`Src/core/convenience.cpp` (`options::digest_io_options`,
`digest_run_options`, `digest_partition_options`, `digest_property_options`,
`digest_ri_options`, `digest_xcw_options`, `digest_dev_options` and
`options::look_for_debug`). Defaults come from `struct options` and
`struct properties_options` in `Src/core/convenience.h`.

Notation: `<x>` is required, `[x]` optional, `x|y` alternatives. Angstrom is
the unit unless stated otherwise.

Either separator works in every flag: `-rgbi_groups` and `-rgbi-groups` are the
same option, and the spelling in the tables is the canonical one. Case is not
normalised, so where a flag has an upper-case form it is listed separately.

**Before you compose a long command line, read
[flag order matters](Troubleshooting.md#flag-order-matters).** Several options
run inside the parser and stop it, so options placed after them are silently
ignored.

## Contents

1. [Input, readers and format conversion](#1-input-readers-and-format-conversion)
2. [Scattering factors and `.tsc`/`.tscb` tables](#2-scattering-factors-and-tsctscb-tables)
3. [Partitioning, grids and basis](#3-partitioning-grids-and-basis)
4. [Multi-fragment and disorder workflows](#4-multi-fragment-and-disorder-workflows)
5. [Density, properties and cubes](#5-density-properties-and-cubes)
6. [Topology and basin analysis](#6-topology-and-basin-analysis)
7. [Bonding analysis](#7-bonding-analysis)
8. [RI density fitting, SALTED and interaction energies](#8-ri-density-fitting-salted-and-interaction-energies)
9. [X-ray constrained wavefunction fitting](#9-x-ray-constrained-wavefunction-fitting)
10. [Performance and hardware](#10-performance-and-hardware)
11. [Diagnostics, testing and developer flags](#11-diagnostics-testing-and-developer-flags)

---

## 1. Input, readers and format conversion

Almost every job starts by naming a wavefunction. The reader is chosen from the
file extension, and the formats are not interchangeable in capability: `.gbw`,
`.fchk` and `.molden` carry a contracted basis and are the only ones the
NBO-family tools can use, while `.wfn`/`.wfx` are primitive-only. Unoccupied
MOs are dropped on read for most jobs, which matters for Fukui functions (they
need the virtuals).

| Flag | Arguments | Meaning |
| --- | --- | --- |
| `-wfn <file>` | 1 | The wavefunction. `.gbw`, `.fchk`, `.molden`, `.wfn`, `.wfx`, `.xtb` (pTB), `.ffn`, `.xyz`; ORCA's natural orbitals `base.mdci.nat` (CCSD) and `base.mp2nat` (relaxed MP2), and an `orca_2json` export `base.json` of an AutoCI density (CCSD(T), MP2-MP4). See [Correlated densities](#correlated-densities-mp2-ccsd-ccsdt). |
| `-xyz <file>` | 1 | Geometry only. With `-SALTED` this doubles as the wavefunction source. |
| `-cif <file>` | 1 | Crystal structure: cell, symmetry, atom labels. Switches cube evaluation onto the periodic path. |
| `-hkl <file>` | 1 | Reflection list to compute form factors for. |
| `-b <basis>` | 1 | Basis set name, needed when the input file does not carry one (`.ffn`, some `.fchk` paths). |
| `-d <dir>` | 1 | Directory the basis-set files are looked up in. |
| `-coef <file>` | 1 | Pre-computed density-fitting coefficient file. |
| `-fchk <file>` | 1 | Write a Gaussian-format checkpoint file. |
| `-wfn_cif <file>` | 1 | Write the wavefunction geometry and basis into a CIF. |
| `-gbw2wfn <file>` | 1 | Convert an ORCA `.gbw` into a `.wfn`. |
| `-convert_to_47 <wfn>` | 1 + sub-options | Write an NBO `.47` archive from a contracted-basis wavefunction. |
| `-cube_convert <in> <out>` | 2 | Text `.cube` <-> binary `.cubeb`. The output format follows the extension of `<out>`. One-shot. |
| `-combine_mos <a> <b>` | 2 | Combine the MO sets of two wavefunctions; select with `-cmos1`/`-cmos2`. One-shot. |
| `-cmos1 <n> [n ...]` | 1+ | MO indices taken from the first wavefunction. |
| `-cmos2 <n> [n ...]` | 1+ | MO indices taken from the second. |

## 2. Scattering factors and `.tsc`/`.tscb` tables

The core job. Give `-wfn`, `-cif` and `-hkl` and NoSpherA2 partitions the
density, transforms each atomic piece over the reflection list and writes a
table. Output is binary `.tscb` by default; `-old_tsc` writes text `.tsc`.
The table is streamed in blocks of reflections so that peak memory stays bounded
— the block size only trades memory against nothing, so leave it alone unless
you are memory-constrained.

| Flag | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `-acc <0..4>` | 1 | `2` | Integration-grid accuracy. Higher is denser and slower. |
| `-mult <n>` | 1 | `0` | Spin multiplicity of the wavefunction. |
| `-charge <n>` | 1 | `0` | Total charge. |
| `-dmin <d>` | 1 | `99.0` (off) | Generate a full reflection sphere down to this resolution instead of reading `-hkl`. |
| `-hkl_min_max <h1 h2 k1 k2 l1 l2>` | 6 | `-100..100` each | Bound the generated index range. |
| `-old_tsc` | 0 | off | Write text `.tsc` instead of binary `.tscb`. |
| `-tshb` | 0 | off | With `-ri_fit` or `-SALTED`, write the fitted densities as `experimental.tshb` (Gaussians times real spherical harmonics, keyed like the table's scatterers, version 4) instead of the table; olex2.refine evaluates them analytically. One file per run: every `-mtc`/`-cmtc` part goes into it (an atom from the first part that has it, as in the table) Thakkar-filled atoms and the atoms of a part without orbitals (an `.xyz`) as Slater records in their fill charge, with the relativistic (heavier than Kr) or beyond-+1 ion differences as fitted Gaussians. With `-do_XCW` each lambda step writes `NA2_<lambda>.tshb` instead of its table: an RI fit of the converged density in the `-ri_fit` basis (AutoAux without it), restrained by `-multipole_moments` if given; it is another partition than the XCW table, R1 1.6e-3 to 3.4e-3 between the two on the P1 test. The table is written instead, with a `-tshb:` line giving the reason, for `-IAM`, a grid partition and `-cif`. |
| `-tsc_block <n>` | 1 | derived from `-mem` | Reflections per streamed block. `0` means hold the whole table in memory. |
| `-mem <MB>` | 1 | `1000.0` | Memory budget the block size is derived from. When given, half of it also bounds occ's stored density-fitting integrals for `-occ` (see below). |
| `-tsc_labels` | 0 | off | Key the table rows by atom label instead of atom ID. |
| `-tsc_labels <table> <cif> [out]` | 2-3 | `<table>.labels.tsc` | Rewrite an existing table's IDs into labels using the CIF. One-shot. |
| `-merge <a.tsc> <b.tsc> ...` | 2+ | | Merge tables, checking that the reflection sets agree. Given `.tshb` files instead, merges them into `combined.tshb`, the first file winning for an atom in several, as the table merge does; `-merge a.tshb b.tshb -cif x.cif [y.cif ...]` checks the result against every atom of the CIFs, and the header says `COVERAGE: complete` only when every input does and no CIF atom is missing, `partial` otherwise (always without `-cif`). One-shot. |
| `-merge_nocheck <a> <b> ...` | 2+ | | Merge without the consistency check. One-shot. |
| `-tscb <file>` | 1 | | Convert between `.tsc` and `.tscb` (the direction follows the input extension). Rejects any other extension with "Wrong file ending!". One-shot. |
| `-IAM` | 0 | off | Write spherical (independent-atom-model) form factors — the reference case. |
| `-ED` | 0 | off | Electron-diffraction scattering factors. Changes `-dmin` to `dmin/2 - 0.001` and ignores the index box. |
| `-twin <9 numbers>` | 9 | | A twin law as a 3x3 matrix, row-major. Repeat the flag for several laws. |
| `-pbc <n>` | 1 | `0` | Number of periodic images included in the density. |
| `-group <n> [n ...]` | 1+ | | Atom group for the current fragment; `+n` marks a negated index. Olex2 emits this. |
| `-all_charges` | 0 | off | Print every atomic charge, not just a summary. |
| `-Anion <"Z Z ...">` | 1+ | | Elements to treat as anions when picking spherical references. |
| `-Cation <"Z Z ...">` | 1+ | | Elements to treat as cations. |
| `-ECP <0..3>` | 1 | `0` | ECP treatment. `1` = ORCA-style def2-ECP, `3` = pTB. Also spelled `-ecp`, `-Ecp`. |
| `-method <name>` | 1 | | Record the QM method used (goes into the table header). |
| `-refine [acc]` | 0-1 | `0.1` | Numerical integral accuracy for the refinement path. |
| `-occ <file>` | 1 | | An occ input (`.toml`) that occ runs, or an occ wavefunction file. With `method = "mp2"`, `"mp3"`, `"mp4(sdq)"`, `"mp4"`, `"ccsd"` or `"ccsd(t)"` the density is the correlated one (relaxed MP2, unrelaxed MP3 / MP4 / CCSD / CCSD(T)), not the SCF one; see [Correlated densities](#correlated-densities-mp2-ccsd-ccsdt) for methods, backends and keys. |

With density fitting (`df-basis` in `[scf]`), occ stores the three-centre integrals over the screened shell pairs,
one triangle, packed, so the SCF does not recompute them every Fock build. They get half of `-mem` when it is given, else
half of the memory the job can have; `df_store_memory_mb` in `[scf]` sets the budget explicitly. Pairs that do not fit
are recomputed every build for J, and exchange stays integral-direct unless every pair fits. `use_direct_df_kernels = true`
recomputes everything. Sucrose r2SCAN/def2-SVP: SCF 54.8 s direct, 15.8 s stored (1.16 GB), 28.3 s at
`df_store_memory_mb = 500`, the same energy to 1e-12 Eh.

### Correlated densities: MP2, CCSD, CCSD(T)

`-occ` with a correlated `method` runs occ's RHF and then the density of that method; the natural orbitals go on to the
scattering factors like any other wavefunction. Closed shell (RHF) only. The densities are occ's
(`occ/qm/cc/density.h`: `ccsd_density`, `local_ccsd_density`, `mp2_relaxed_density`, `mpn_density`); NoSpherA2 builds that code
against MKL (`OCC_DENSITY_CBLAS`), and reads the `ccsd_*` keys and writes the log summary.

| `method` | density |
| --- | --- |
| `"mp2"` | relaxed MP2: MP2 Lagrangian with the Z-vector orbital response (PySCF `grad/mp2.py`), the frozen core relaxed too. SCS/SOS-MP2 (`mp2_spin_scaling`) is refused. |
| `"mp3"`, `"mp4(sdq)"`, `"mp4"` | unrelaxed MP3 / MP4(SDQ) / MP4(SDTQ), ORCA AutoCI's density: oo and vv blocks from the MPn Lagrangian, and for MP4 the ov block of the second-order wavefunction's singles-doubles and doubles-triples couplings (MP3 has none). No orbital response; see the note on unrelaxed MP-n densities below. |
| `"ccsd"` | unrelaxed CCSD from the Lambda equations (PySCF `ccsd_lambda` / `ccsd_rdm`). |
| `"ccsd(t)"` | unrelaxed CCSD(T): the (T) terms in Lambda and in the density, the response form (PySCF `ccsd_t_lambda_slow` / `ccsd_t_rdm_slow`, `for_grad=True`). |

**Integral backends.** `exact` (default) builds the MO integrals from exact four-index integrals and runs the SCF with
exact J/K, as ORCA's `NoRI`. Density fitting is `ccsd_backend = "df"`, a `df-` or `ri-` method prefix (`"df-ccsd(t)"`,
`"ri-mp2"`) or, for MP2, occ's `mp2_backend = "df"`; the fitting basis is `ri-basis` in `[scf]`, else occ's correlation
fitting basis for the orbital basis. THC is not available for densities.

- DF never stores the O(v⁴) integrals (benzene cc-pVDZ 0.6 GB, naphthalene 3.6 GB in the exact backend), but the Lambda
  and (T) steps still hold dense o·v³ blocks, and the SCF becomes DF-JK too. At small sizes it is not faster.
- Benzene cc-pVDZ CCSD(T), DF against exact: 69 s vs 54 s (8 cores), E_corr 4.3e-4 Eh, E(T) 5.6e-5 Eh, scattering
  factors 7e-5 e (C) — about 25× below the (T) effect, about 40× above the occ/ORCA agreement.
- occ's `cc-pvdz-rifit` has no oxygen: with cc-pVDZ and O in the molecule set `ri-basis = "def2-svp-rifit"` (or
  `"cc-pvtz-rifit"`).

**`[scf]` keys** (all `method`s unless noted):

| Key | Default | Meaning |
| --- | --- | --- |
| `ccsd_frozen_core` | `-1` | `-1` chemical core (ORCA's default), `0` all electrons, `n` the n lowest occupied. |
| `ccsd_backend` | `exact` | `exact` or `df` (see above). |
| `ccsd_max_memory_gb` | `1.0` | Integral and (T) buffer budget. The (T) step needs about 8·(4·o·v³ + 8·v³·threads) bytes (naphthalene cc-pVDZ ~4 GB). |
| `ccsd_tol` | `1e-5` | CCSD and Lambda convergence: max \|D·Δt\| of the amplitudes, ORCA's MDCI `STol` with TightSCF. The MP2 Z-vector is solved to 1e-10. |
| `ccsd_local` | `false` | CCSD / CCSD(T) only: cluster-in-molecule densities (one cluster per Pipek-Mezey LMO) for larger systems, with `ccsd_cluster_radius` (8 Å), `ccsd_domain_threshold` (0.01), `ccsd_pair_threshold` (1e-5 Eh) and `ccsd_lno_threshold` (1e-6). Clusters that span the molecule give the canonical density back. |

**Accuracy and cost.** Against ORCA 6.1 (water and benzene cc-pVDZ, ORCA's tolerances): scattering factors within
~3e-6 e for MP2 (relaxed), CCSD (MDCI) and CCSD(T) (AutoCI); energies within ORCA's convergence (~1e-6 Eh). MP3 and
MP4 against AutoCI (water): natural occupations within 7e-8, energies within 1e-7 Eh, the same charge table. Canonical
CCSD(T) scales as o³v⁴: on 8 cores, water seconds, benzene ~1 min, naphthalene ~12 min (ORCA MDCI 14 min for its
CCSD(T) energy and CCSD density, AutoCI 13 min for benzene's CCSD(T) density).

**Densities from ORCA** (all read with `-wfn`):

- CCSD: `! CCSD NoRI` with `%mdci Density Unrelaxed NatOrbs true end` writes `base.mdci.nat`. MDCI's CCSD(T) run writes
  the *CCSD* density (no (T)); the log says so when `base.property.txt` lies beside the file.
- relaxed MP2: `! MP2 NoRI` with `%mp2 Density Relaxed NatOrbs true end` writes `base.mp2nat`.
- CCSD(T) (and MP2, MP3, MP4, unrelaxed): AutoCI, `! AUTOCI-CCSD(T) NoRI KeepDens` with `%autoci density unrelaxed end`
  (on Windows AutoCI is a separate download, `Orca6.1.0.Win64_autoci*.zip`, unpacked into the ORCA folder), then
  `orca_2json base.gbw` with `base.json.conf` = `{"MOCoefficients": true, "Densities": ["autocipur"]}` and `-wfn base.json`.
- Unrelaxed MP-n densities (ORCA's or occ's, which are the same) are poor for scattering factors: on water they close only a third to a half of the HF to
  CCSD(T) difference (max |Δf| on O against CCSD(T): HF 2.2e-2, MP2 1.5e-2, MP4 1.2e-2, CCSD 1.4e-3 e). Use the
  relaxed MP2 density, CCSD or CCSD(T).

## 3. Partitioning, grids and basis

How the molecular density is divided among atoms. This is the single choice that
most changes the resulting form factors. `PartitionType` is one of `Becke`,
`TFVC`, `Hirshfeld`, `RI`, `MBIS`, `EMBIS`; **Hirshfeld is the default** and is
what the NoSpherA2 paper describes.

| Flag | Arguments | Meaning |
| --- | --- | --- |
| `-hirsh` | 0 | Hirshfeld partitioning (the default; the flag is for being explicit). |
| `-becke` | 0 | Becke fuzzy-cell partitioning. Also `-Becke`, `-BECKE`. |
| `-tfvc` | 0 | Topological fuzzy Voronoi cells. Also `-TFVC`. |
| `-mbis` | 0 | Minimal-basis iterative stockholder. Also `-MBIS`. |
| `-embis` | 0 | Extended MBIS. Also `-EMBIS`. MBIS and EMBIS take their exponentials four at a time (`exp_kernel`, within 1 ulp of the CRT; sucrose charges and tsc identical to print precision); `-tune NOS_EXP_CRT` uses the CRT `exp` instead. |
| `-ri_fit <basis...>` | 1+ | RI partitioning: fit the density in an auxiliary basis and partition the fit. See [section 8](#8-ri-density-fitting-salted-and-interaction-energies). Also `-RI_FIT`. |
| `-sfac_diffuse <a> <b> <c> <cif> <wfn> <dmin>` | 6 | Diffuse-scattering form factors. One-shot. |
| `-atom_sfac <wfn1> <wfn2>` | 2 | Compare the form factors of one atom between two wavefunctions. |
| `-spherical_atoms` | 0 | Write the tabulated spherical atomic densities. One-shot. |

## 4. Multi-fragment and disorder workflows

For a structure that needs more than one wavefunction: disorder parts,
independent molecules, or a mixed treatment where some fragments are aspherical
and others IAM.

`-mtc` and `-cmtc` do **not** use the separate `-group` flag. Each takes
alternating tokens of its own and keeps reading until the next token that begins
with `-`:

* `-mtc <wfn> <parts> <wfn> <parts> ...`
* `-cmtc <wfn> <cif> <parts> <wfn> <cif> <parts> ...`

`<parts>` is a **single token** listing the CIF disorder parts that wavefunction
covers, separated by a comma or a dot — so `0.1` means parts 0 **and** 1, not
"part 0.1". Part 0 is every atom whose `_atom_site_disorder_group` is absent or
unreadable, i.e. the ordered backbone, which is why a disorder job is normally
`0.1` and `0.2` rather than `1` and `2`. The numbers are matched against that
CIF field directly.

| Flag | Arguments | Meaning |
| --- | --- | --- |
| `-mtc <wfn> <parts> ...` | 2 per fragment | Multi-tsc mode: several wavefunctions, each with the parts it covers. |
| `-cmtc <wfn> <cif> <parts> ...` | 3 per fragment | As `-mtc`, each fragment bringing its own CIF. |
| `-mtc_mult <n> [n ...]` | 1+ | Per-fragment multiplicities, in `-mtc` order. |
| `-mtc_charge <n> [n ...]` | 1+ | Per-fragment charges. |
| `-mtc_ECP <n> [n ...]` | 1+ | Per-fragment ECP modes. |

Disorder tables are keyed by atom ID, not by label — see
[Troubleshooting](Troubleshooting.md#disorder-tables-are-keyed-by-atom-id).

## 5. Density, properties and cubes

Point-wise fields on a grid. The grid is a box around the molecule: `-radius`
sets how far past the outermost atom it extends and every point further than
that from *any* atom is set to zero; `-resolution` is the spacing. Each property
is a separate flag and several can be written in one run.

| Flag | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `-resolution <d>` | 1 | `0.1` | Grid spacing. |
| `-radius <r>` | 1 | `2.0` | Grid extent past the atoms; also the hard cut-off radius. |
| `-rho` | 0 | off | Electron density cube. |
| `-lap` | 0 | off | Laplacian of the density. |
| `-eli` | 0 | off | ELI-D cube. |
| `-eli_spin` | 0 | off | Spin-resolved ELI-D cubes from one pass: `<stem>_eli_aa.cube`, `_eli_bb.cube`, `_eli_triplet.cube` (aa only for a closed shell, which equals `-eli`); 0 where a spin channel has no density. `scripts/vmd/eli_spin.py` turns them into a VMD scene. |
| `-elf` | 0 | off | Electron localisation function. |
| `-esp` | 0 | off | Electrostatic potential. |
| `-rdg` | 0 | off | Reduced density gradient (NCI). |
| `-def` | 0 | off | Static deformation density. Also `-DEF`. |
| `-HDEF` | 0 | off | Hirshfeld deformation density. Also `-hdef`. |
| `-MO <n>\|all` | 1 | | One molecular orbital, or every one. Repeat for several. |
| `-ibo_cube <list>` | 1 | | Runs `-ibo` and writes `<wfn>_IBO_<n>.cube` for the IBOs in `<list>`: the numbers `-ibo` prints, `all`, or `C2:C3` (every bond-like IBO between atoms 2 and 3, labels 1-based, element optional), comma-separated, e.g. `6,C2:C3`. Same grid as `-MO`. |
| `-s_rho` | 0 | off | Spin density. |
| `-fukui` | 0 | off | Fukui function. Also `-Fukui`. Needs the virtual orbitals. |
| `-esp_isosurface [v]` | 0-1 | `0.002` when bare, else off | ESP mapped onto a density isosurface at value `v`. |
| `-hirshfeld_surface <asu.xyz> <pack.xyz>` | 2 | | Hirshfeld surface of the asymmetric unit against its packing environment. Writes `.obj` and `Hirshfeld_surface.dat`. |
| `-cubeb` | 0 | off | Write cubes in the binary `.cubeb` format. |
| `-cube <file>` | 1 | | Read a density cube instead of computing one. Also `-cube_density`. |
| `-rho_cube <wfn>` | 1 | | Density cube for a wavefunction and nothing else. Set `-radius`/`-resolution` **before** this flag. One-shot. |
| `-density_difference <wfn2>` | 1 | | Subtract the density of a second wavefunction from the first. |
| `-atom_dens <wfn> [alpha_MOs] [beta_MOs]` | 1-3 | | Spherically averaged atomic density; MO lists are comma-separated. One-shot. |
| `-atom_dens_diff <wfn1> <wfn2>` | 2 | | Difference of two spherically averaged densities. One-shot. |
| `-spherical_aver_fukui <wfn1> <wfn2>` | 2 | | Radial Fukui profile to `fukui_averaged_density_wfn.dat`. |
| `-spherical_aver_hirsh <...>` | | | Radial Hirshfeld-weight profile. |
| `-rho_at_points <points>` | 1 | | Density at explicit points; the file holds a count then `x y z` per line **in bohr**, output goes to `<points>.rho`. Needs `-wfn` first. One-shot. |
| `-calc_dens_1D [a1] [a2] [n] [pad]` | 0-4 | `0 1 1000 2.0` | Density along the line between two atoms. Needs `-wfn` and `-ri_fit` first. One-shot. |
| `-dipole_moments` | 0 | | Dipole moments of the wavefunction. One-shot, output goes to the log. |
| `-polarizabilities <7 wfn files>` | 7 | | Polarisability from seven finite-field wavefunctions. |
| `-e_field <F>` | 1 | `0.005` | Field strength used by the finite-field jobs. |
| `-laplacian_bonds <wfn>` | 1 | | Laplacian profiles along each bond. One-shot. |
| `-fukui_analysis [wfn]` | 0-1 | | Condensed Fukui indices, printed to the console. |
| `-promol_nci <f1.xyz> <f2.xyz> [f3.xyz ...] [rcut1] [rcut2] [rho_max] [rdg_max] [colour_max]` | 2+ | `0.95 0.75 0.5 1.0 0.015` | Intermolecular NCI from promolecular densities of two or more fragments. |
| `-promol_nci_single_thread` | 0 | off | Serialise the promolecular NCI loop. |
| `-fractal <cube>` | 1 | | Fractal-dimension analysis of a cube. |
| `-get_g` | 0 | off | Also write the kinetic-energy density G. |
| `-ewal_sum <cube> [k_max] [alpha]` | 1-3 | | Ewald sum over a residual-density cube. One-shot. |
| `-draw_orbits <l,m[,res[,radius]]>` | 1 | `res 0.025`, `radius 3.5` | Cube of a single real spherical harmonic. One-shot. |
| `-spherical_harmonic` | 0 | | Self-test of the spherical-harmonic code. One-shot. |

## 6. Topology and basin analysis

Integration rather than plotting. `-topology` finds every critical point of the
density from the analytic Hessian and checks the set with the Poincare-Hopf
relation. `-eli_analysis` does the real work: it locates QTAIM and ELI-D basins,
integrates them, and reports charges, basin populations and localisation /
delocalisation indices (LI/DI).

| Flag | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `-topology <wfn>` | 1 | | Every critical point of rho — nuclear, bond, ring, cage — with the Poincare-Hopf check. Console output. One-shot. |
| `-eli_analysis <wfn> <resolution> <radius>` | 3 | | QTAIM + ELI-D basins, populations, LI/DI. |
| `-eli_family <wfn> [points]` | 1-2 | | ELI-D alpha/beta/triplet and ELI-q as point values. Console output. One-shot. |
| `-qtaim_eli <rho.cube> <eli.cube> <atoms> [bg]` | 3-4 | `bg 0` | Mask an ELI-D cube to the QTAIM basins of the listed atoms. `<atoms>` is comma-separated, **0-based**. One-shot. |
| `-qtaim_eli <wfn> <atoms> [res] [radius] [bg]` | 2-5 | | The same, computing both fields from a wavefunction. One-shot. |
| `-basin_grid <n>` | 1 | `1` | Basin-integration grid level. |
| `-basin_cube` | 0 | off | Integrate basins on a cube instead of walking the analytic gradient. |
| `-qct` | 0 | off | Interactive QCT (quantum chemical topology) menu. Also `-QCT`. |

Basins walk the analytic gradient by default; `-basin_cube` restores the older
cube-based path. ELI-D and the fitted densities still use a cube on purpose.
Both paths stop at the electron-density isosurface of 0.0001 e/bohr³; density
outside it is reported separately from the finite basin populations and volumes.

## 7. Bonding analysis

Localised-orbital bonding descriptors. Two independent routes:

* **In house**: `-nbo_native` runs the NBO search, the hybrid analysis, second-
  order perturbative (E2) donor-acceptor energies and natural resonance theory
  (NRT) inside NoSpherA2, with no external program. NRT is solved as a convex QP
  over the simplex and parallelised over candidate structures. Each NRT bond
  order comes with its range over the weights whose D(w) lies within 0.1 % of the
  fit's: the fit determines a bond order only that far (benzene C-C 1.44-1.57).
* **External**: `-nbo` writes a `.47` archive, calls `gennbo`, and parses the
  result. `-nbo_parse` parses an existing output.

`-rgbi` is a separate scheme (resonance-group bond indices) and `-npa` gives
natural population analysis on its own. `-ibo` gives IAO partial charges and
intrinsic bond orbitals (Knizia, JCTC 9 (2013) 4834) against the MINAO
reference basis, closed shell only.

All of these need a contracted basis, so `.gbw`, `.fchk` or `.molden` — not
`.wfn`/`.wfx`.

| Flag | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `-nbo_native <wfn>` | 1 | | NBO search, hybrids, E2 and NRT in house. |
| `-nbo <wfn>` | 1 | | Write a `.47`, run external `gennbo`, parse the output. |
| `-nbo_parse <file>` | 1 | | Parse an existing `gennbo` output. |
| `-convert_to_47 <wfn>` | 1 | | Write the `.47` archive and stop. |
| `-nbo_keywords <"...">` | 1 | | Extra `$NBO` keylist entries. |
| `-nbo_exe <path>` | 1 | | The `gennbo` executable. On Windows, `wsl:<path>` names a binary inside WSL. |
| `-nbo_dir <dir>` | 1 | | Working directory for the external run. |
| `-nbo_json <file>` | 1 | | Write the analysis as JSON. |
| `-nbo_47 <file>` | 1 | | Name for the written `.47`. |
| `-nbo_e2min <kcal>` | 1 | | Drop E2 interactions below this threshold. |
| `-nbo_threads <n>` | 1 | | Threads for the native NBO/NRT search. |
| `-nbo_keep47` | 0 | off | Do not delete the `.47` afterwards. |
| `-nrt` | 0 | off | Run NRT (native) / request it in the keylist (external). |
| `-nrt_e2 <n>` | 1 | | NRT E2 threshold. |
| `-nrt_arrows` | 0 | off | Report NRT arrow (charge-transfer) contributions. |
| `-nrt_bond_scale <f>` | 1 | | Scale factor on the NRT bond orders. |
| `-nrt_max <n>` | 1 | | Maximum number of resonance structures kept. |
| `-nrt_atoms <list>` | 1 | | Restrict the resonance search to these atoms. |
| `-nrt_exhaustive` | 0 | off | Enumerate candidates exhaustively. |
| `-nrt_no_symmetry` | 0 | off | Skip symmetry detection. |
| `-nrt_no_components` | 0 | off | Do not split into components. |
| `-nrt_no_ion` | 0 | off | Exclude ionic structures. |
| `-nrt_rotate` | 0 | off | Native: a structure that is a symmetry image of one whose orbital sweeps converged starts its own sweeps from that structure's orbitals, rotated into place in the NAO basis, instead of the greedy fill. Sweeps on single runs: SF6 59124 -> 31900, TiCl4 9355 -> 2921, benzene 6157 -> 3424; bond orders over the 41 v2 molecules move by up to 2.2e-4. |
| `-nrt_target <D(w)>` | 1 | one search | Native only, implies `-nrt`: walk arrows again from every structure the fit weights and refit, until D(w) reaches the value, no new structure takes weight or the `-nrt_max` pool is full; `0` = until one of the latter. |
| `-nrt_screen <r\|off>` | 1 | off | Native, with `-nrt_target`: a walk that yields 200 or more new structures builds in full only those whose reduced cost, from orbitals taken from the walk root's untouched slots plus three sweeps, is below `r` - 16 at a time with their symmetry images, refitting in between. At `r` = 0.2 grown benzene, thiophene and PF5 run 2.4-2.9 times faster, but a structure whose orbitals settle late is missed: SF6 builds none of its walk and ends with D(w) 0.23 % above the unscreened fit, while thiophene, pyridine and Ni(CO)4 end below it through a different walk. |
| `-nrt_no_mirror` | 0 | mirror on | Native: switches off what is otherwise default since 2026-10-14 (formerly opt-in `-nrt_mirror`, still accepted): a single move that adds no formal-charge separation (sum over atoms of \|valence electrons - Lewis electron count\|), such as the allyl cation's bond shift, counts as a whole resonance arrow instead of a dropped intermediate; symmetry images complete the candidate set even where the walk passed them as intermediates; the walk from a symmetry reference of the parent explores the structures the parent's walk already reached. With it allyl's and the allyl cation's mirror bonds are equal (without: 0.744/0.721 in the minority spin, 1.463/1.445 in the cation; NBO 7 weighs the two mirror structures equally), a cation distorted by 0.006 A stays within the symmetric fit's 0.1 % range, and grown Ni(CO)4 keeps Td. Mean \|dBO\| against NBO 7 over the 41 v2 molecules: 0.0356 with, 0.0384 without (LiF 0.997 against 0.000, NaCl 0.000 against 1.000, CF3 beta 0.031 against 0.167); r2SCAN 0.0451 against 0.0477. Grown thiophene finds 50 % more structures with it. |
| `-nrt_stall <n>` | 1 | `0` = off | Native: a structure's orbital sweeps stop once Anderson mixing has not lowered the largest vector change by 0.1 % for `n` sweeps, instead of running to the cap of 200; the log counts the structures stopped this way. With `n` = 20, CPU time of single runs: Ni(CO)4 -41 %, SF6 -23 %, PF5 -20 %; 29 of the 41 v2 molecules give identical output, ten move by at most 1.6e-6 in D(w), BF4- by 5.5e-5, and Ni(CO)4 ends 0.05 % higher with bond orders moved by up to 0.0017; agreement with NBO 7 (mean \|dBO\| 0.0384 -> 0.0383) unchanged. Grown Ni(CO)4 (`-nrt_target 0`) ends 0.12 % higher with bond orders moved by up to 0.006. |
| `-nrt_parent <file.json>` | 1 | search's own | Native only, closed shell, implies `-nrt`: the NBO search places this Lewis structure (NBO's `$CHOOSE`) and NRT starts from it. JSON `{"format": "NoSpherA2 Lewis structure", "version": 1, "atoms": N, "lone_pairs": {"O1": 2}, "bonds": [["O1", "H2", 1]]}`: atoms as element + 1-based index or a bare index, valence lone pairs only, bond orders 1-3; a `.native.nbo.json` is read through its `"parent_lewis"`. Refused with a one-line reason on an unknown label, wrong element, duplicate bond, order outside 1-3, an electron count that does not add up, `three_centre` or `alpha`/`beta` entries, or a bond beyond `-nrt_bond_scale`. The log compares its rho(NL) with the search's own best structure. |
| `-npa` | 0 | off | Natural population analysis with the per-orbital table. |
| `-npa_summary` | 0 | off | NPA, totals only. |
| `-ibo` | 0 | off | IAO charges and IBOs (Pipek-Mezey-like, exponent 4, 2x2 Jacobi), listed as core / lone pair / bond-like / delocalised with IAO and Mulliken populations per atom. Closed-shell `.gbw` or spherical `.molden`; refuses unrestricted. `-ibo_cube` writes them as cubes. |
| `-rgbi` | 0 | off | Resonance-group bond indices. |
| `-rgbi_groups <"i,j,k"> [...]` | 1+ | | Explicit atom groups (comma-separated, per group). |
| `-rgbi_basis nao\|ano` | 1 | `ano` | Orbital basis used for RGBI. |
| `-rgbi_no_sym` | 0 | off | Disable symmetry use in RGBI. |
| `-rgbi_theta` | 0 | off | Report the RGBI mixing angles. |
| `-rgbi_EVs` | 0 | off | Report the RGBI eigenvalues. |

**Every `-nbo_*` and `-nrt*` sub-option must appear *after* the
`-nbo_native`/`-nbo`/`-nbo_parse`/`-convert_to_47` flag it belongs to** — they
are read by a loop starting at that flag's position, so anything earlier on the
line is ignored without a warning.

### EQC energy decomposition

Rahm and Hoffmann's E = nX̄ + (Vnn − Eee) split of a bond-forming reaction,
fragments → parent, printed as X-analysis `-m 2` prints it (Q, covalency and
ionicity). Two modes:

* **Mode 1** (`-eqc <gbw> -eqc_frag ...`): the ORCA `.gbw` seeds occ, which
  re-converges the parent and then each fragment from the parent's density
  block, so every term comes from one program. `E(ORCA)` is read from the `.out`
  next to the gbw.
* **Mode 2** (`-eqc_wfn <parent> <frag> <frag> ...`): every term read from the
  files, no SCF. `.wfx` and `.fchk` carry their own energy; `.wfn`, `.gbw` and
  `.molden` need the ORCA `.out` of the same name beside them.

| Flag | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `-eqc [gbw]` | 0-1 | | Mode 1 on this ORCA gbw (or `-wfn`). |
| `-eqc_frag <atoms> <charge> <mult>` | 3+ | | One fragment per triple, atoms 0-based as `0,2-4`; every atom in exactly one fragment. |
| `-eqc_method hf\|<functional>` | 1 | `hf` | occ method for Mode 1; must be the gbw's. |
| `-eqc_basis <name>` | 1 | from the gbw | Load the basis by name (needed with ECPs); checked against the gbw. |
| `-eqc_cold` | 0 | off | Also converge every SCF from occ's own guess and report its iteration count. |
| `-eqc_wfn <files>` | 3+ | | Mode 2: parent first, then the fragments. |

Example, ethane into two methyl radicals:
`NoSpherA2 -eqc ethane.gbw -eqc_frag 0,2-4 0 2 1,5-7 0 2`.
Theory, every output line and the caveats: [EQC energy decomposition](EQC-Energy-Decomposition.md).

## 8. RI density fitting, SALTED and interaction energies

Two ways to get a density that is cheap to evaluate: fit the real one in an
auxiliary basis (`-ri_fit`), or predict the fitting coefficients with a trained
SALTED model (`-SALTED`). Both give the same object — a set of auxiliary-basis
coefficients — and both can then be partitioned, put on a grid, or used for
electrostatic interaction energies. Multipole restraints keep the fit's moments
honest.

A fitted or predicted density has **no orbitals**, so ELF, MO and RDG are not
available from it.

| Flag | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `-ri_fit <basis> [basis ...]` | 1+ | | RI-fit the density in these auxiliary basis sets, and use RI partitioning. Also `-RI_FIT`. Needs a wavefunction with its basis set (`.gbw`, `.fchk`, `.molden`; a `.wfx` or `.wfn` stops), and stops when the fit holds the wavefunction's electrons to worse than 0.1 % (1 % with `-multipole_moments`). |
| `-ri_fit auto_aux [elements...] [basis...]` | 1+ | | Generate the auxiliary set automatically; element symbols right after `auto_aux` restrict which elements it is generated for, the basis names that follow fill the rest (first match wins). |
| `-multipole_moments <scheme> [lmax]` | 1-2 | off; omitted `lmax` is `0` | Restrain atomic moments up to `lmax` (0-8). The fit restrains atom-centred auxiliary coefficients for charges and dipoles; orders 2 and above use grid-partitioned rows for stability. It does not constrain the total electron count exactly. `Hirshfeld`, `TFVC`, `MBIS` and `EMBIS` supply the targets. `Nuclear`, `Mulliken` and `Sanderson` use atom-centred charge-only order 0; any supplied order is ignored. |
| `-multipole_strength <w>` | 1 | `1.0` | Weight of the multipole restraints. Must be positive. |
| `-write_ri_coefs` | 0 | | Fit and write `RI_COEFS.npy`. One-shot. |
| `-ri_cube <...>` | | | Cube of the fitted density. Also `-RI_CUBE`. |
| `-SALTED <model>` | 1 | | Predict the density with a SALTED model — a directory or a `.salted` file. Also `-salted`. |
| `-SALTED_COEFS` | 0 | | Predict and write the coefficients only. Needs `-wfn` (or `-xyz`) first. One-shot. Also `-salted_coefs`. |
| `-SALTED_Training` | 0 | | Generate SALTED training data. Needs `-wfn` and `-ri_fit` first. One-shot. |
| `-interaction_energy <A> <B>` | 2 | | Electrostatic interaction between two fitted densities. With `-SALTED` both are predicted, otherwise both are RI-fitted. One-shot. |
| `-interaction_energy <A> <A.npy> <B> <B.npy>` | 4 | | The same from pre-computed coefficient files. One-shot. |
| `-interaction_energies <job>` | 1 | | Every contacting molecule pair in the crystal, from a job file. |
| `-repulsion_overlap <f>` | 1 | `0.0` | Overlap-repulsion weight in the interaction energy. Must be >= 0. |
| `-repulsion_exchange dirac\|pbe\|b88\|r2scan` | 1 | `dirac` | Exchange functional for the repulsion term. |
| `-classify_atoms <...>` | | | Geometry-based atom-type classification. |
| `-classify_atoms_list <...>` | | | The same over a list of structures. |
| `-calc_featomic_descriptor <...>` | | | featomic descriptor for one structure. |
| `-calc_featomic_descriptors <list>` | 1 | | One descriptor `.npy` per structure named in the list file. |
| `-geometry_aid_cutoff <r>` | 1 | `3.5` | Neighbour cut-off for the geometry-aid model. |
| `-geometry_aid_center_weight <w>` | 1 | `1.0` | Weight of the central atom in the descriptor. |
| `-geometry_aid_metals` | 0 | off | Include metals in the classification. |

## 9. X-ray constrained wavefunction fitting

Fit the wavefunction against the measured structure factors: a lambda scan in
which the X-ray agreement statistic is traded against the energy. The dominant
cost is the two-electron integral ("I") tensor, so much of the `-XCW_settings`
file is about where that tensor lives.

Everything about the crystal comes from `-cif`: cell, symmetry, positions,
U_iso, U_ij, the Gram-Charlier C and D, the wavelength and the anomalous
dispersion (`_atom_type_scat_dispersion_real/imag`, or per site
`_atom_site_dispersion_*`). For every lambda the run writes
`NA2_<lambda>.tscb`, `.wfn`, `.fchk` and `_Fcalc.txt`, where `<lambda>` loses
its decimal point and is padded to seven digits (0.01 gives `NA2_0010000`).
The SCF details go to `SCF.log`.

| Flag | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `-do_XCW` | 0 | off | Run the XCW lambda scan. Needs `-cif`, `-hkl` and `-XCW_settings`. |
| `-XCW_settings <file>` | 1 | | Settings file; see [settings-file keywords](#settings-file-keywords). |
| `-xcw_partition grid\|ri` | 1 | `grid` | How XCW shares the density among the atoms for F_calc. `grid`: the grid weights (`-hirsh` default). `ri`: an RI fit of the density in the `-ri_fit` basis (AutoAux without one), (mu nu\|P) from occ held in memory; with `-multipole_moments` its atomic moments are restrained to the grid partition of that scheme (`hirshfeld`, `tfvc`, `mbis`, `embis`, exact derivative at fixed shares), to Mulliken charges (linear in D, exact) or to fixed Nuclear / Sanderson charges. `-tfvc`, `-mbis` and `-embis` work on the `grid` route too: their shares are refit to XCW's own density (stand-alone `make_chi` / MBIS / EMBIS fits on the XCW grid, MBIS and EMBIS to 1e-6) inside an SCF solve while the density moves by more than 1e-3 (MBIS/EMBIS loose there, 1e-2 of that move), held once it moves less, and between solves of a step until the refit moves F by less than 1e-5 of max \|F\|; they need the grid route (not with `read`, `safe` or `NOSPHERA2_XCW_GRID`). |
| `-xyz <file>` | 1 | | The grown molecule, for a settings file with `grown`. |
| `-twin <9 numbers>` | 9 | | Twin laws applied to the reflections, as in [section 2](#2-scattering-factors-and-tsctscb-tables). |
| `-mem <MB>` | 1 | | When given, caps the I tensor, the density-fitting integrals and the stored two-electron integrals together; without it they are sized from the memory the job can have. |
| `-acc <0..4>` | 1 | `2` | Integration grid of the I tensor and of each `.tscb`. |
| `-xcw_int_precision <eps>` | 1 | `1e-10` | Integral screening threshold of the Fock build. |
| `-xcw_extrapolate` / `-no_xcw_extrapolate` | 0 | on | Seed each lambda step from the density extrapolated through the two previous steps (RHF only). |
| `-xcw_incremental` / `-no_xcw_incremental` | 0 | off | Build the two-electron Fock part from the change of the density each iteration. |
| `-rgbi` | 0 | off | Roby-Gould bond indices of each refined wavefunction, written to `NA2_<lambda>_RGBI.txt`. |
| `-convert_XCW <stdout> <step>` | 2 | | Import lambda steps from a Tonto XCW log. One-shot. |

The I tensor runs on the GPU by default; `-gpu_itensor`, `-gpu_itensor_tensor`,
`-itensor_hybrid` and `-cpu_itensor_fp32` in [section 10](#10-performance-and-hardware)
choose how. The partitioning flags of [section 3](#3-partitioning-grids-and-basis)
and `-all_charges` apply to each `.tscb` written. Developer knobs, through
`-tune NAME=value`: `NOSPHERA2_XCW_I_FLOAT` (single-precision tensor),
`NOSPHERA2_NO_INVERSION_SCREENING` (integrate both h and -h),
`NOSPHERA2_XCW_ERI_NO_BUDGET` (store the two-electron integrals outside `-mem`) and
`NOSPHERA2_DEBUG_ASYMFACT` (print the symmetry factors). For density-following shares: `NOSPHERA2_XCW_WEIGHT_UPDATE=k` (refit at every k-th evaluation inside the SCF instead of between solves; `0` freezes them per step, the reference), `NOSPHERA2_XCW_WEIGHT_FTOL` (default `1e-5`), `NOSPHERA2_XCW_WEIGHT_TOL` (MBIS/EMBIS fit, default `1e-6`), `NOSPHERA2_XCW_WEIGHT_DTOL` (density move that refits inside a solve, default `1e-3`, `0` between solves only), `NOSPHERA2_XCW_WEIGHT_MIX` (Anderson depth of the MBIS/EMBIS fits, default `6`, `0` the plain fixed point), `NOSPHERA2_XCW_TFVC_OLD` (the stand-alone Becke+TFVC weight pass instead of the TFVC-only one), `NOSPHERA2_XCW_WEIGHT_CHECK` (how far the shares move F and R1 over a step) and `NOSPHERA2_XCW_ADJOINT_CHECK` (now after the first step).

Example, run in `tests/P1_test`:
`NoSpherA2 -do_XCW -no_gpu -cif P1_test.cif -hkl P1_subset.hkl -XCW_settings test_settings_save.txt -no_date`,
with `test_settings_save.txt` holding
`normal normal_conv params 177 basis_set sto-3g max_iter 100 F charge 0 mult 1 rhf start 0 step_size 0.01 end 0.05 safe p1_i.bin save_eri p1_eri.bin`.

### Settings-file keywords

The file is one stream of whitespace-separated keywords, on one line or many,
and case does not matter. A keyword that takes a value is followed by it. An
unknown keyword stops the run, and so does a file without `basis_set`.

**Model and scan**

| Keyword | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `basis_set <name>` | 1 | required | Orbital basis from the built-in library, e.g. `sto-3g`, `def2-svp`. |
| `basis_overrides <file>` | 1 | | Per-element replacement bases from a BSE JSON file, relative to the settings file; other elements keep `basis_set`. |
| `params <n>` | 1 | `1` | Number of refined parameters, subtracted from the reflection count in the GooF and criterion denominators. Without it the run warns and uses 1. |
| `charge <n>` | 1 | `0` | Total charge of the model. |
| `mult <n>` | 1 | `1` | Spin multiplicity. |
| `rhf` \| `uhf` | 0 | `rhf` | Restricted or unrestricted Hartree-Fock. |
| `grown` | 0 | off | The model is a grown molecule rather than the asymmetric unit; needs `-xyz`. Symmetry operations it applies explicitly leave the sum and the atom weights are solved to count every position once. |
| `wavelength <lambda>` | 1 | from the CIF | Wavelength in Angstrom; overrides `_diffrn_radiation_wavelength`. |
| `start <x>` | 1 | `0` | First lambda. |
| `end <x>` | 1 | `1` | Last lambda, inclusive. |
| `step_size <x>` | 1 | `0.01` | Lambda increment. |

**Fit target**

| Keyword | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `f` \| `f2` | 0 | `f` | Fit against F or F^2. |
| `weighted` | 0 | off | Fit the 1/\|H\|^2-weighted residual self-energy criterion instead of the classical goodness of fit. |
| `i_sigma <x>` | 1 | `2` | Only reflections with I/sigma(I) >= x enter the fit and the scale, i.e. F/sigma(F) >= 2x. R1(all) and Crit(all) still use all. |
| `extinction <model> [iso\|aniso] [fixed] [start]` | 1-4 | | Extinction correction; `<model>` is `shelx` (`zachariasen`), `bc_gaussian` or `bc_lorentzian`. Isotropic, refined and starting at `1e-4` unless told otherwise. Needs a wavelength. |

**SCF convergence**

One convergence preset and one damping preset; the individual keywords after
them override single values.

| Preset | Target change | DIIS error | Orbital gradient | Max density change | RMS density change |
| --- | --- | --- | --- | --- | --- |
| `sloppy` | `3e-5` | `1e-4` | `5e-4` | `1e-4` | `1e-5` |
| `normal` (default) | `1e-6` | `1e-5` | `7e-5` | `1e-5` | `1e-6` |
| `tight` | `5e-7` | `5e-6` | `3e-5` | `1e-6` | `1e-7` |
| `very_tight` | `1e-7` | `1e-6` | `1e-5` | `1e-7` | `1e-8` |

| Preset | Damping | Level shift | Switched off below DIIS error |
| --- | --- | --- | --- |
| `slow_conv` | `0.8` | `1.0` | `1e-5` (both) |
| `normal_conv` (default) | `0.5` | `0.5` | `1e-3` (damping), `1e-2` (shift) |
| `fast_conv` | none | none | |

| Keyword | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `conv <x>` | 1 | preset | Change of the target quantity, energy + lambda x penalty. |
| `conv_diis <x>` | 1 | preset | DIIS error. |
| `gradient <x>` | 1 | preset | Orbital gradient. |
| `maxp_diff <x>` | 1 | preset | Largest change of the density matrix. |
| `rmsp_diff <x>` | 1 | preset | RMS change of the density matrix. |
| `max_iter <n>` | 1 | `100` | SCF iterations per lambda step. |
| `damp <alpha>` | 1 | preset | Density damping factor. |
| `shift <x>` | 1 | preset | Level shift. |
| `diis_damping <x>` | 1 | preset | DIIS error below which damping is switched off. |
| `diis_shift <x>` | 1 | preset | DIIS error below which the level shift is switched off. |
| `soscf` | 0 | | Accepted and ignored. A cold start runs DIIS, a warm start (the previous lambda step or a `.fchk` restart) TRAH from its first iteration. |
| `check_hessian` | 0 | off | Finite-difference check of the Hessian-vector product on the first second-order step, reported in `SCF.log`; enters TRAH once the DIIS error is below `1e-2`. |
| `df_basis <name>` | 1 | | Density-fit the Fock build with this auxiliary basis, e.g. `def2-universal-jkfit`. |
| `hessian_df_basis <name>` | 1 | automatic | Density-fit only the TRAH Hessian products; iterations keep exact Fock builds. Unset, `def2-universal-jkfit` is used when the ERIs do not fit in memory; `none` keeps the products exact. |
| `guess_basis <name>` | 1 | | Start the first lambda from a Hartree-Fock converged in this smaller basis, projected into the orbital basis. |
| `load_wfn` | 0 | off | Start the first lambda from `NA2_<start>.fchk` of an earlier run. |

**I tensor and stored integrals**

The I tensor depends on the structure, the basis and the reflections only, not
on any refinement setting, so a saved one serves every later scan of the same
data. Without `stream` or `i_tensor_mb` it is held in memory when it fits four
fifths of the memory the job can have (or of `-mem`) and streamed when it does
not.

| Keyword | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `stream` | 0 | | Keep the tensor on disk and read a window of reflections at a time, with a 2048 MB budget. |
| `i_tensor_mb <n>` | 1 | | The same with an n MB budget. |
| `safe <path>` | 1 | | Write the tensor to `<path>` while the refinement runs, for a later `read <path>`. |
| `read [path]` | 0-1 | `I_tensor_stream.bin` | Reuse the tensor at `path` instead of building it; built when the file is missing or does not match. A tensor that is read must not be written again with `safe`. |
| `i_float` | 0 | off | Hold the tensor in single precision: half the memory. |
| `i_double` | 0 | off | Keep it in double precision although it was built in single. |
| `save_eri <path>` | 1 | | Write the stored two-electron integrals for a later run. |
| `read_eri <path>` | 1 | | Reuse them; a file of another basis or geometry is refused and the integrals are built. |

**Analysis and output**

| Keyword | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `gaussian_halt` | 0 | off | Test the residuals of the strong reflections for normality at each lambda, add the A^2 column and recommend a halting lambda. |
| `strong_cutoff <x>` | 1 | `3.0` | \|F\|/sigma above which a reflection counts as strong for `gaussian_halt`. |
| `nbo` | 0 | off | Also write `NA2_<lambda>.47`, an NBO input, for every lambda. |

## 10. Performance and hardware

GPU acceleration is **on by default** for the Fourier transform, the XCW I
tensor, the SALTED equicomb step, the Becke/TFVC grid weights and the fitted
density, with a silent CPU fallback when no device or runtime is found. Single
precision is used where it is safe; the double-precision alternative is one flag
away in each case. See
[GPU and CPU do not agree bitwise](Troubleshooting.md#gpu-and-cpu-do-not-agree-bitwise)
for the sizes of the differences.

| Flag | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `-cpus <n>` | 1 | all | OpenMP threads. **Does not reach the native NBO/NRT search** — that reads `-nbo_threads` only, and without it takes every core regardless of `-cpus`. |
| `-mem <MB>` | 1 | `1000.0` | Memory budget (drives the tsc block size). |
| `-no_gpu` | 0 | | Disable every GPU path. |
| `-gpu_fp64` | 0 | off | Double precision on the GPU throughout. |
| `-gpu_fp32` | 0 | off | Force single precision. |
| `-gpu_grid` / `-no_gpu_grid` | 0 | on | Grid weights on the GPU. |
| `-gpu_density` / `-no_gpu_density` | 0 | on | Fitted-density evaluation on the GPU. |
| `-gpu_salted` / `-no_gpu_salted` | 0 | on | SALTED equicomb on the GPU. |
| `-gpu_itensor` / `-no_gpu_itensor` | 0 | on | XCW I tensor on the GPU. |
| `-gpu_itensor_tensor` / `-no_gpu_itensor_tensor` | 0 | off | Tensor-core path for the I tensor. |
| `-gpu_cublas` / `-no_gpu_cublas` | 0 | on | Use cuBLAS when the runtime is present. |
| `-gpu_blas` | 0 | off | Route the dense linear algebra through the GPU BLAS. |
| `-cpu_itensor_fp32` / `-no_cpu_itensor_fp32` | 0 | on | Single precision for the CPU I tensor. |
| `-itensor_hybrid` / `-no_itensor_hybrid` | 0 | off | Split the I tensor between CPU and GPU. |
| `-gflops` | 0 | off | Report achieved GFLOP/s. |
| `-basin_grid <n>` | 1 | `1` | Basin-integration grid level (also a cost knob). |

## 11. Diagnostics, testing and developer flags

These exist for development and regression testing. They are listed for
completeness; their arguments and behaviour are not part of any stable
interface and can change without notice.

| Flag | Meaning |
| --- | --- |
| `-test` | Run the built-in test suite. |
| `-profiling [root]` | Profiling run; `root` defaults to `tests`. Also `-profile`. |
| `-partitioning_test` | Partitioning regression test. |
| `-test_RI` | RI-fitting regression test. |
| `-RI_WFN_DIFF` | Compare a fitted density against the wavefunction it came from. |
| `-NNLS_TEST` | Non-negative-least-squares self-test. |
| `-lahvatest`, `-lukas_test` | Ad-hoc development tests. |
| `-rkpts`, `-skpts` | Read / save the k-point set. |
| `-out <file>` | Rename the log file (read before the parser, so position-independent). |
