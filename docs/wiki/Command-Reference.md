# Command reference

Every flag below was read out of the option parser in
`Src/core/convenience.cpp` (`options::digest_io_options`,
`digest_run_options`, `digest_partition_options`, `digest_property_options`,
`digest_ri_options`, `digest_xcw_options`, `digest_dev_options` and
`options::look_for_debug`). Defaults come from `struct options` and
`struct properties_options` in `Src/core/convenience.h`.

Notation: `<x>` is required, `[x]` optional, `x|y` alternatives. Angstrom is
the unit unless stated otherwise.

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
| `-wfn <file>` | 1 | The wavefunction. `.gbw`, `.fchk`, `.molden`, `.wfn`, `.wfx`, `.xtb` (pTB), `.ffn`, `.xyz`. |
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
| `-tsc_block <n>` | 1 | derived from `-mem` | Reflections per streamed block. `0` means hold the whole table in memory. |
| `-mem <MB>` | 1 | `1000.0` | Memory budget the block size is derived from. |
| `-tsc_labels` | 0 | off | Key the table rows by atom label instead of atom ID. |
| `-tsc_labels <table> <cif> [out]` | 2-3 | `<table>.labels.tsc` | Rewrite an existing table's IDs into labels using the CIF. One-shot. |
| `-merge <a.tsc> <b.tsc> ...` | 2+ | | Merge tables, checking that the reflection sets agree. One-shot. |
| `-merge_nocheck <a> <b> ...` | 2+ | | Merge without the consistency check. One-shot. |
| `-tscb <file>` | 1 | | Convert between `.tsc` and `.tscb` (the direction follows the input extension). Rejects any other extension with "Wrong file ending!". One-shot. |
| `-IAM` | 0 | off | Write spherical (independent-atom-model) form factors — the reference case. |
| `-ED` | 0 | off | Electron-diffraction scattering factors. Changes `-dmin` to `dmin/2 - 0.001` and ignores the index box. |
| `-anom_disp <file>` | 1 | | Add anomalous dispersion corrections from a file. |
| `-twin <9 numbers>` | 9 | | A twin law as a 3x3 matrix, row-major. Repeat the flag for several laws. |
| `-pbc <n>` | 1 | `0` | Number of periodic images included in the density. |
| `-group <n> [n ...]` | 1+ | | Atom group for the current fragment; `+n` marks a negated index. Olex2 emits this. |
| `-all_charges` | 0 | off | Print every atomic charge, not just a summary. |
| `-Anion <"Z Z ...">` | 1+ | | Elements to treat as anions when picking spherical references. |
| `-Cation <"Z Z ...">` | 1+ | | Elements to treat as cations. |
| `-ECP <0..3>` | 1 | `0` | ECP treatment. `1` = ORCA-style def2-ECP, `3` = pTB. Also spelled `-ecp`, `-Ecp`. |
| `-method <name>` | 1 | | Record the QM method used (goes into the table header). |
| `-refine [acc]` | 0-1 | `0.1` | Numerical integral accuracy for the refinement path. |
| `-occ <file>` | 1 | | Occupation-number input file. |

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
| `-embis` | 0 | Extended MBIS. Also `-EMBIS`. |
| `-ri_fit <basis...>` | 1+ | RI partitioning: fit the density in an auxiliary basis and partition the fit. See [section 8](#8-ri-density-fitting-salted-and-interaction-energies). Also `-RI_FIT`. |
| `-sfac_diffuse <a> <b> <c> <cif> <wfn> <dmin>` | 6 | Diffuse-scattering form factors. One-shot. |
| `-atom_sfac <wfn1> <wfn2>` | 2 | Compare the form factors of one atom between two wavefunctions. |
| `-spherical_atoms` | 0 | Write the tabulated spherical atomic densities. One-shot. |

## 4. Multi-fragment and disorder workflows

For a structure that needs more than one wavefunction: disorder parts,
independent molecules, or a mixed treatment where some fragments are aspherical
and others IAM. Each fragment gets its own wavefunction and its own `-group`
list; `-mtc` collects them.

| Flag | Arguments | Meaning |
| --- | --- | --- |
| `-mtc <wfn> ... -group <...> ...` | varies | Multi-tsc mode: several wavefunctions, each with its own atom group. |
| `-cmtc <wfn> ... -group <...> ...` | varies | As `-mtc`, for fragments related by crystallographic symmetry. |
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
| `-elf` | 0 | off | Electron localisation function. |
| `-esp` | 0 | off | Electrostatic potential. |
| `-rdg` | 0 | off | Reduced density gradient (NCI). |
| `-def` | 0 | off | Static deformation density. Also `-DEF`. |
| `-HDEF` | 0 | off | Hirshfeld deformation density. (Lowercase `-hdef` is **not** accepted.) |
| `-MO <n>\|all` | 1 | | One molecular orbital, or every one. Repeat for several. |
| `-s_rho` | 0 | off | Spin density. |
| `-fukui` | 0 | off | Fukui function. Also `-Fukui`. Needs the virtual orbitals. |
| `-esp_isosurface [v]` | 0-1 | `0.002` when bare, else off | ESP mapped onto a density isosurface at value `v`. |
| `-hirshfeld_surface <asu.xyz> <pack.xyz>` | 2 | | Hirshfeld surface of the asymmetric unit against its packing environment. Writes `.obj` and `Hirshfeld_surface.dat`. |
| `-cubeb` | 0 | off | Write cubes in the binary `.cubeb` format. |
| `-cube <file>` | 1 | | Read a density cube instead of computing one. Also `-cube_density`. |
| `-rho_cube <wfn>` | 1 | | Density cube for a wavefunction and nothing else. Set `-radius`/`-resolution` **before** this flag. One-shot. |
| `-density_difference <wfn2>` | 1 | | Subtract the density of a second wavefunction from the first. Also `-density-difference`. |
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

## 7. Bonding analysis

Localised-orbital bonding descriptors. Two independent routes:

* **In house**: `-nbo_native` runs the NBO search, the hybrid analysis, second-
  order perturbative (E2) donor-acceptor energies and natural resonance theory
  (NRT) inside NoSpherA2, with no external program. NRT is solved as a convex QP
  over the simplex and parallelised over candidate structures.
* **External**: `-nbo` writes a `.47` archive, calls `gennbo`, and parses the
  result. `-nbo_parse` parses an existing output.

`-rgbi` is a separate scheme (resonance-group bond indices) and `-npa` gives
natural population analysis on its own.

All of these need a contracted basis, so `.gbw`, `.fchk` or `.molden` — not
`.wfn`/`.wfx`.

| Flag | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `-nbo_native <wfn>` | 1 | | NBO search, hybrids, E2 and NRT in house. |
| `-nbo <wfn>` | 1 | | Write a `.47`, run external `gennbo`, parse the output. |
| `-nbo_parse <file>` | 1 | | Parse an existing `gennbo` output. |
| `-convert_to_47 <wfn>` | 1 | | Write the `.47` archive and stop. |
| `-nbo_keywords <"...">` | 1 | | Extra `$NBO` keylist entries. |
| `-nbo_exe <path>` | 1 | | The `gennbo` executable. |
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
| `-npa` | 0 | off | Natural population analysis with the per-orbital table. |
| `-npa_summary` | 0 | off | NPA, totals only. |
| `-rgbi` | 0 | off | Resonance-group bond indices. |
| `-rgbi-groups <"i,j,k"> [...]` | 1+ | | Explicit atom groups (comma-separated, per group). |
| `-rgbi_basis nao\|ano` | 1 | `ano` | Orbital basis used for RGBI. |
| `-rgbi_no_sym` | 0 | off | Disable symmetry use in RGBI. |
| `-rgbi_theta` | 0 | off | Report the RGBI mixing angles. |
| `-rgbi_EVs` | 0 | off | Report the RGBI eigenvalues. |

**Every `-nbo_*` and `-nrt*` sub-option must appear *after* the
`-nbo_native`/`-nbo`/`-nbo_parse`/`-convert_to_47` flag it belongs to** — they
are read by a loop starting at that flag's position, so anything earlier on the
line is ignored without a warning.

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
| `-ri_fit <basis> [basis ...]` | 1+ | | RI-fit the density in these auxiliary basis sets, and use RI partitioning. Also `-RI_FIT`. |
| `-ri_fit auto_aux [elements...] [basis...]` | 1+ | | Generate the auxiliary set automatically; element symbols right after `auto_aux` restrict which elements it is generated for, the basis names that follow fill the rest (first match wins). |
| `-multipole_moments <scheme> <lmax>` | 2 | off (`lmax -1`) | Restrain the fit's multipole moments up to `lmax` (0-8), using `Hirshfeld`, `Becke`, `TFVC`, `MBIS` or `EMBIS` moments. Also `-multipole-moments`. |
| `-multipole_strength <w>` | 1 | `1.0` | Weight of the multipole restraints. Must be positive. |
| `-multipole_partition` | 0 | on | Restrain partitioned atomic moments (the default). |
| `-multipole_centre` | 0 | | Restrain moments about a single centre instead. Also `-multipole_center`. |
| `-write_ri_coefs` | 0 | | Fit and write `RI_COEFS.npy`. One-shot. |
| `-ri_cube <...>` | | | Cube of the fitted density. Also `-RI_CUBE`. |
| `-SALTED <model>` | 1 | | Predict the density with a SALTED model — a directory or a `.salted` file. Also `-salted`. |
| `-SALTED_COEFS` | 0 | | Predict and write the coefficients only. Needs `-wfn` (or `-xyz`) first. One-shot. Also `-salted_coefs`. |
| `-SALTED_Training` | 0 | | Generate SALTED training data. Needs `-wfn` and `-ri_fit` first. One-shot. |
| `-salted_charge_constraint` | 0 | off | Constrain the predicted density to the correct electron count. **Use it.** |
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
cost is the two-electron integral ("I") tensor, so the `-XCW_settings` file is
mostly about where that tensor lives.

| Flag | Arguments | Default | Meaning |
| --- | --- | --- | --- |
| `-do_XCW` | 0 | off | Run the XCW lambda scan. |
| `-XCW_settings <file>` | 1 | | Settings file. Keywords include `stream`, `i_tensor_mb`, `save`, `read`, `df_basis`, `i_double`, `i_float`. |
| `-calc_F` | 0 | off | Compute structure factors from the fitted wavefunction. |
| `-xcw_int_precision <eps>` | 1 | `1e-10` | Integral screening threshold. |
| `-xcw_extrapolate` / `-no_xcw_extrapolate` | 0 | on | Extrapolate between lambda steps. |
| `-xcw_incremental` / `-no_xcw_incremental` | 0 | off | Build the Fock matrix incrementally. |
| `-xcw_strong_cutoff <x>` | 1 | `3.0` | Cut-off separating strong from weak reflections. |
| `-xcw_gaussian_halt` | 0 | off | Halt the scan on the Gaussian-statistics criterion. |
| `-convert_XCW <stdout> <step>` | 2 | | Import lambda steps from a Tonto XCW log. One-shot. |

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
| `-cpus <n>` | 1 | all | OpenMP threads. |
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
