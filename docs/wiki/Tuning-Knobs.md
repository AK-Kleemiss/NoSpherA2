# Tuning knobs

Every developer knob NoSpherA2 reads, with its default and what it changes.
They are set on the command line:

```powershell
.\NoSpherA2.exe ... -tune NAME=VALUE -tune OTHER_NAME
```

* `-tune` is repeatable; a bare `NAME` means `NAME=1`.
* Knobs are applied before any other option is parsed, live only in the
  process, and are restored when the run ends, so one in-process run (the DLL,
  the test binary) never leaks a knob into the next. They used to be environment
  variables; setting them in the environment no longer does anything.
* An on/off knob is on whenever it is set, whatever its value, unless the
  table says `0` turns it off. The NAO switches are off for `0`.
* Numeric knobs that fall outside their range are ignored (with a message) and
  the default stays, except the RGBI thresholds, which stop the run.

None of this is a stable interface. Most knobs exist to compare a fast path
against its reference, to time one stage, or to make a rare branch testable.
A knob that changes results says so in the table; a run that sets one is not a
default run.

The source column names the file that reads the knob (under `Src/core/`).
`python scripts/check_tuning_knobs.py` lists any knob in the sources that this
page does not mention, and any name here that no source reads any more.

## XCW (X-ray constrained wavefunctions)

| Knob | Default | Effect | Source |
| --- | --- | --- | --- |
| `NOSPHERA2_XCW_RCHECK=<fchk>` | unset | Scale, extinction and R1 of that converged wavefunction against this model and data, its F_calc in `rcheck_Fcalc.txt`, then stop: no integrals, no SCF. | `XCW_wrapper.cpp` |
| `NOSPHERA2_XCW_ERI_NO_BUDGET` | off | Store the two-electron integrals whatever `-mem` leaves. | `XCW_wrapper.cpp` |
| `NOSPHERA2_ERI_DYNAMIC` | off | The stored-integral Fock build hands rows out dynamically, longest first. J and K then vary at round-off (1e-15 relative) between runs; the default fixed cyclic order is bit for bit reproducible at a given thread count. No gain measured. | `stored_eri.cpp` |
| `NOSPHERA2_XCW_ADJOINT_CHECK` | off | Check the grid transform against its adjoint at a random symmetric D, after the first converged step. | `XCW_wrapper.cpp` |
| `NOSPHERA2_XCW_TSHB_FIT` | off | With `-xcw_partition ri -tshb`: fit the density anew for the `.tshb` instead of reusing the partition's own fit. | `XCW_wrapper.cpp` |
| `NOSPHERA2_XCW_TSHB_TABLE` | off | Write the `.tscb` table beside the `.tshb`, to compare the two. | `XCW_wrapper.cpp` |
| `NOSPHERA2_XCW_TSCB_BY_F` | off | Match the grid route's table columns to the scatterers by the closest f instead of by atom id or label. | `structure_factors.cpp` |
| `NOSPHERA2_DEBUG_ASYMFACT` | off | Print the deleted symmetry operations, the surviving ones and every atom's symmetry factor. | `structure_factors.cpp` |
| `NOSPHERA2_NO_INVERSION_SCREENING` | off | Integrate both h and -h instead of taking -h as the conjugate. | `structure_factors.cpp` |

### XCW grid partition

| Knob | Default | Effect | Source |
| --- | --- | --- | --- |
| `NOSPHERA2_XCW_GRID_DIRECT` | off | Direct-sum transform of the grid instead of Fourier-Bessel. | `xcw_grid.cpp` |
| `NOSPHERA2_XCW_GRID_FB_H` | `0.25` | Node step h r_max of the 4-point Lagrange table of the Fourier-Bessel transform; error ~ (h r)^4. | `xcw_grid.cpp` |
| `NOSPHERA2_XCW_GRID_DTOL` | `0` (off) | Per-block bound on the \|dF\| a kept rho may carry (host only). Off because an SCF moves every element of D each cycle. | `xcw_grid.cpp` |
| `NOSPHERA2_XCW_GRID_TIME` | off | Wall time of each grid stage, printed per call. | `xcw_grid.cpp` |
| `NOSPHERA2_XCW_GRID_AO_TOL` | `-acc` cutoff | AO value below which a block drops an AO; below 1e-20 also turns the Morton reordering off. | `structure_factors.cpp` |
| `NOSPHERA2_XCW_GRID_NO_MORTON` | off | Keep the grid points in shell order instead of Morton order. Sums do not depend on the order; only the speed does. | `structure_factors.cpp` |
| `NOSPHERA2_XCW_GRID_CHUNK` | `1024` | Points per AO block; `0` uses the three radial bands. | `structure_factors.cpp` |
| `NOSPHERA2_XCW_GRID_TSCB_OLD` | off | Build every step's `.tscb` on the tsc route instead of the grid's own per-atom transform. | `structure_factors.cpp` |
| `NOSPHERA2_XCW_GRID_FCHECK` | off | F_calc from the written table against the one the SCF used at the same density. | `XCW_wrapper.cpp` |

### Density-following shares (`-tfvc`, `-mbis`, `-embis` under XCW)

| Knob | Default | Effect | Source |
| --- | --- | --- | --- |
| `NOSPHERA2_XCW_WEIGHT_UPDATE=k` | `-1` | Refit the shares at every k-th evaluation inside the SCF instead of between solves; `0` freezes them per step (the reference). | `XCW_wrapper.cpp` |
| `NOSPHERA2_XCW_WEIGHT_DTOL` | `1e-3` | Density move that refits the shares inside a solve; `0` refits between solves only. | `XCW_wrapper.cpp` |
| `NOSPHERA2_XCW_WEIGHT_FTOL` | `1e-5` | A step is solved again while the refit moves F by more than this fraction of max \|F\|. | `XCW_wrapper.cpp` |
| `NOSPHERA2_XCW_WEIGHT_TOL` | `1e-6` | Convergence of the MBIS/EMBIS fits (never tighter than the grid's own). | `xcw_weights.cpp` |
| `NOSPHERA2_XCW_WEIGHT_MIX` | `6` | Anderson depth of the MBIS/EMBIS fits; `0` the plain fixed point. | `xcw_weights.cpp` |
| `NOSPHERA2_XCW_WEIGHT_CHECK` | off | How far the shares move F and R1 over a step. | `XCW_wrapper.cpp` |
| `NOSPHERA2_XCW_TFVC_OLD` | off | The stand-alone Becke + TFVC weight pass (the reference) instead of the TFVC-only one. | `xcw_weights.cpp` |
| `NOSPHERA2_XCW_CHI_FALLBACK` | off | `make_chi`'s chi = 1 for a bond whose line density has no minimum, instead of d_A/d_B at the line's inflection. Changes the shares. | `xcw_weights.cpp` |

### XCW RI partition (`-xcw_partition ri`)

| Knob | Default | Effect | Source |
| --- | --- | --- | --- |
| `NOSPHERA2_XCW_RI_FULL` | off | Keep (mu nu\|P) as the full square instead of mu <= nu (twice the memory). | `xcw_ri.cpp` |
| `NOSPHERA2_XCW_RI_PACK_CHECK` | off | Compare 2 B^T D of the square and of the packed store at a random symmetric D. | `xcw_ri.cpp` |

## RI fitting and integrals

| Knob | Default | Effect | Source |
| --- | --- | --- | --- |
| `NOS_RI_SCREEN` | `3e-7` | Screening threshold of the three-centre integrals; `0` evaluates every shell triple. | `libCintMain.cpp` |
| `NOS_RI_FAR` | `1e-14` | erfc overlap below which an auxiliary shell is treated as a point multipole (far field); `0` turns the far field off. | `libCintMain.cpp` |
| `NOSPHERA2_RI_ROWS_NOCUT` | off | Evaluate every shell for the RI rows (the reference the equivalence test compares against). | `integrator.cpp` |
| `NOSPHERA2_RI_PENALTY_LOOP` | off | Add the penalty rows one by one (`add_penalty_row`) instead of collected; both give the same H to the last bit. | `integrator.cpp` |

## Grids and scattering factors

| Knob | Default | Effect | Source |
| --- | --- | --- | --- |
| `NOS_SF_DIRECT` | off | Direct sum for every atom instead of Fourier-Bessel on the prototype shells. | `scattering_factors.cpp` |
| `NOS_SF_FB_L=<L>` | D/2 | Fourier-Bessel expansion order instead of half the prototype's largest Lebedev degree; forces every atom onto Fourier-Bessel. | `scattering_factors.cpp` |
| `NOS_SF_FB_EXACT` | off | Fourier-Bessel without the \|k\| interpolation table; forces every atom onto Fourier-Bessel. | `scattering_factors.cpp` |
| `NOS_SF_FB_FORCE` | off | Fourier-Bessel for every atom, also where the direct sum would be faster. | `scattering_factors.cpp` |
| `NOS_GRID_OLD_FLOOR` | off | `lebedev_low` as the innermost angular floor for every element, not only from K on. | `AtomGrid.cpp` |
| `NOS_MBIS_BECKE` | off | MBIS/EMBIS moments summed over every grid with Becke weights instead of each atom's own grid with the atom weights. | `AtomGrid.cpp`, `GridManager.cpp` |
| `NOS_EXP_CRT` | off | The CRT `exp` in MBIS/EMBIS instead of the four-wide `exp_kernel` (within 1 ulp). | `AtomGrid.cpp` |
| `NOSPHERA2_CHI_DEBUG` | off | Print a checksum of the TFVC chi matrix. | `AtomGrid.cpp` |
| `NOS_DENSITY_ACCURACY` | `5e-5` | Accuracy the primitive screening of the density holds, in (0, 1). Looser values move ELI-D basin populations. | `wfn_io.cpp` |

## GPU

| Knob | Default | Effect | Source |
| --- | --- | --- | --- |
| `NOSPHERA2_GPU_BACKEND=cuda\|hip` | probed | Pick the GPU backend instead of probing for one. | `gpu_dispatch.cpp` |
| `NOSPHERA2_GPU_BATCH=<n>` | fits the card | Cap the atoms per structure-factor batch, to exercise batching on a card large enough not to need it. | `sf_gpu.cu` |
| `NOSPHERA2_GRID_LOCAL=1` | off | Force the local (near-centre) grid-weight kernel that only molecules wider than a few hundred centres reach. | `grid_gpu.cu` |
| `NOSPHERA2_BLAS_GPU_MIN_FLOP` | 4e9 x fp64 ratio / 64 | Flop count above which `-gpu_blas` sends a product to the device. `-gflops` re-derives it. | `blas_gpu.cu` |
| `NOS_BASIN_GPU_MIN_US` | `15` | The basin field goes to the device only if one host evaluation costs at least this many microseconds. | `b2c.cpp` |
| `NOS_BASIN_GPU_FIBERS` | `256` | Trajectories each device feed advances together. | `b2c.cpp` |
| `NOS_BASIN_GPU_PROFILE` | off | Time the device basin-field kernels. | `basin_field_gpu.cu` |

## Basin integration (QTAIM, ELI-D)

| Knob | Default | Effect | Source |
| --- | --- | --- | --- |
| `NOS_BASIN_SKIP` | `5e-5` | Electrons per molecule the lightest points may move instead of being climbed; `0` climbs every point. | `b2c.cpp` |
| `NOS_BASIN_STEP_CAP` | `2000` | Steps one trajectory may take before it is given up (at least 100). | `b2c.cpp` |
| `NOS_BASIN_STEP_RELAX` | `1.5` | Factor by which a shortened step grows back toward the base step, in [1, 4]. | `b2c.cpp` |
| `NOS_BASIN_SHRINK` | on | A floor step that fails to rise halves and tries again; `0` turns it off. | `b2c.cpp` |
| `NOS_BASIN_SHRINK_FLOOR` | `0.0625` | Fraction of the base step below which a failing step stops halving, in (0, 1]. 1/64 clears more stalls but splits near-degenerate lone pairs. | `b2c.cpp` |
| `NOS_ADP_CAP` | `8` | Adaptive step: at most this many times the validated floor step. | `b2c.cpp` |
| `NOS_ADP_GROW` | `0.99999` | Adaptive step: midpoint cosine that earns a doubling. | `b2c.cpp` |
| `NOS_ADP_KEEP` | `0.999` | Adaptive step: below this cosine the step is thrown away and retaken at the floor. | `b2c.cpp` |
| `NOS_ADP_REACH` | `0.25` | Adaptive step: fraction of the distance to the nearest maximum. | `b2c.cpp` |
| `NOS_BETA_MARGIN` | `0.9` | Beta spheres stay within this fraction of the smallest sampled safe radius, at most 1. | `b2c.cpp` |
| `NOS_BASIN_MEMO` | on | Streaming runs memoise visited voxels' basins; `0` turns it off. | `b2c.cpp` |
| `NOS_BASIN_MEMO_MB` | `256` | Memory of that memo in MB. | `b2c.cpp` |
| `NOS_BASIN_ROTATE` | on | Rotate the angular grids of the basin-integration grid; `0` turns it off. | `b2c.cpp` |
| `NOS_BASIN_SCAN_REVERSE` | off | Scan seed voxels in the opposite order (diagnostic; results must not change). | `b2c.cpp` |
| `NOS_ELI_SHELL_DIST` | `1.2` | Bohr within which near-degenerate ELI-D shell maxima merge into one basin, in [0, 10); `0` turns the merge off. | `bondwise_analysis.cpp` |
| `NOS_ELI_SHELL_TOL` | `0.05` | Relative value difference up to which those maxima count as degenerate, in [0, 1). | `bondwise_analysis.cpp` |
| `NOS_ELI_SHELL_SADDLE` | `0.99` | Two shell maxima join when the ELI-D between them never drops below this times the lower one, in [0, 1). | `b2c.cpp` |
| `NOS_ELI_SHELL_PERSIST` | `5e-4` | Shell persistence: two shell maxima of one atom join when the arc between them dips less than this; `0` turns it off. Measured safe window 5e-4 to 6e-4. | `b2c.cpp` |
| `NOS_ELI_SHELL_SADDLE_LOG` | off | Log each shell join and persistence decision. | `b2c.cpp` |

## Bonding analysis (RGBI, NAO, NBO)

| Knob | Default | Effect | Source |
| --- | --- | --- | --- |
| `NOS_RGBI_PINV_CUTOFF` | `1e-5` | Singular values below which RGBI pseudo-inverses cut, > 0. Out of range stops the run. | `bondwise_analysis.cpp` |
| `NOS_RGBI_PAIR_TOL` | `1e-4` | Window in which RGBI pairs ionic eigenvalues and spares lone pairs, in (0, 1). | `bondwise_analysis.cpp` |
| `NOS_RGBI_GROUP_LP` | `0.99` | RGBI groups drop ionic pairs as lone pairs above this, in (0, 2). | `bondwise_analysis.cpp` |
| `NOS_RGBI_PAIR_LP` | `0.999` | RGBI atom pairs drop ionic pairs as lone pairs above this, in (0, 2). | `bondwise_analysis.cpp` |
| `NOS_RGBI_NONREL_FREEATOM` | off | Non-relativistic free atoms instead of sf-X2C-1e for all-electron atoms (comparison only). | `bondwise_analysis.cpp` |
| `NOS_RGBI_FREE_ATOM_MAXITER=<n>` | occ's default | Cap the free-atom SCF iterations, so the unconverged warning is testable. | `bondwise_analysis.cpp` |
| `NOS_RGBI_NO_FREEATOM_CACHE` | off | Recompute every free atom, neither reading nor writing either cache. | `bondwise_analysis.cpp` |
| `NOS_FREEATOM_CACHE_DIR=<dir>` | `~/.cache/NoSpherA2/free_atom_cache` | Directory of the free-atom disk cache; empty or `off` disables the disk cache alone. | `bondwise_analysis.cpp` |
| `NOS_RGBI_PARALLEL_FREEATOM` | off | Compute the free-atom SCFs concurrently up front; the serial loop then finds them cached. Opt-in because occ's SCF is not shown re-entrant. | `bondwise_analysis.cpp` |
| `NOS_RGBI_NO_PIN` | off | Do not pin occ to one thread for the free atoms (measures the pin's cost; open shells may then vary). | `bondwise_analysis.cpp` |
| `NOS_RGBI_DEBUG` | off | One line per free atom: cache hits, SCF energies, timings. | `bondwise_analysis.cpp` |
| `NBO_BOND_FLOOR` | `0.15` | Minority-centre share below which a two-centre NBO candidate is no bond. | `nbo.cpp` |
| `NBO_JSON_DIGITS` | `10` | Significant digits in the NBO JSON, 1 to 20. | `nbo_run.cpp` |
| `NAO_PRENAO_NET` | off | Net-density pre-NAO diagnostic. | `nao.cpp` |
| `NAO_LEGACY_CASCADE` | off | Omit the Rydberg re-naturalisation (step 5), for comparison. | `nao.cpp` |
| `NAO_OWSO_OFF` | off | Skip the occupancy-weighted symmetric orthogonalisation. | `nao.cpp` |
| `NAO_CLASS_SPLIT` | off | Valence and Rydberg NAOs in separate blocks. | `nao.cpp` |
| `NAO_CORE_POOLED` | off | Core, valence and Rydberg in one block per (atom, l) (without `NAO_CLASS_SPLIT`). | `nao.cpp` |
| `NAO_DUMP_CPRE` | off | Dump the AO -> pre-NAO transform (NBO 7 unit 32). | `nao.cpp` |
| `NAO_DUMP_STEP3` | off | Dump the occupancies step 4 inherits. | `nao.cpp` |
| `NAO_DUMP_C` | off | Dump the AO -> NAO transform (NBO 7 unit 33). | `nao.cpp` |

## SALTED descriptors

The three overrides give a descriptor that matches no trained model; they are
for timing only and the run says so.

| Knob | Default | Effect | Source |
| --- | --- | --- | --- |
| `NOSPHERA2_SPLINE_ACCURACY` | `1e-6` | Spline accuracy of the SOAP radial basis. | `geometry_aid.cpp` |
| `NOSPHERA2_MAX_RADIAL` | `6` | Radial channels of the descriptor. | `geometry_aid.cpp` |
| `NOSPHERA2_MAX_ANGULAR` | `12` | Angular order of the descriptor. | `geometry_aid.cpp` |
| `NOSPHERA2_TIME_SOAP` | off | Phase timings of the descriptor and SALTED prediction. | `geometry_aid.cpp`, `SALTED_utilities.cpp` |

## Not knobs: the remaining environment variables

Two settings are still read from the environment, because they are needed
before the option parser runs:

| Variable | Effect |
| --- | --- |
| `NOSPHERA2_HIP_RUNTIME` | File of the HIP runtime library to load, ahead of the search through `ROCM_PATH`, `HIP_PATH` and `/opt/rocm`; the kernel registration opens it before `main()` (`hip_runtime_shim.cpp`). |
| `DEBUG_WAIT` | Wait for a debugger to attach at start-up (`app/debug_utils.h`). |
