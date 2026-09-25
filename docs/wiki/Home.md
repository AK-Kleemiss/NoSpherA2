# NoSpherA2 command-line documentation

NoSpherA2 does two things.

1. **It produces aspherical atomic form factors** for Hirshfeld atom refinement
   (HAR). It reads a quantum-chemical wavefunction and a reflection list, cuts
   the molecular density into atomic pieces with a partitioning scheme, Fourier
   transforms each piece, and writes a `.tsc`/`.tscb` table that
   `olex2.refine` uses in place of spherical (IAM) form factors.
2. **It analyses wavefunctions and densities**: cubes (rho, Laplacian, ESP,
   ELF, ELI-D, RDG/NCI, deformation density, MOs), QTAIM and ELI-D basins,
   Hirshfeld surfaces, NBO/E2/NRT bonding analysis, RI density fitting, SALTED
   machine-learned densities, and X-ray constrained wavefunction (XCW) fitting.

There are two ways it gets used:

* **Driven by Olex2.** The NoSpherA2 GUI tab builds the command line for you.
  Everything on these pages is what that tab is doing underneath, which is why
  a few flags (`-group`, `-cif`, `-mtc`, `-acc`) look odd standing alone.
* **Standalone on the command line.** Either as the scattering-factor step of a
  script-driven refinement, or purely as an analysis tool with no
  crystallography involved at all.

## Read this first: where the output goes

**A NoSpherA2 run prints almost nothing to the console.** The very first thing
`run_app` does is redirect `std::cout` into a log file:

```powershell
.\NoSpherA2.exe -wfn mol.gbw -cif mol.cif -hkl mol.hkl -acc 2 -mult 1 -charge 0
Get-Content .\NoSpherA2.log -Tail 40
```

The log is `NoSpherA2.log` in the working directory. `-out <file>` renames it.
(`-out` is read by a dedicated pre-pass before the option parser, so it works
anywhere on the line, and it is not listed in `-h`.)

A handful of commands deliberately restore the console first — `-h`,
`-topology`, `-eli_family`, `-fukui_analysis` — because their whole output is
meant to be read. Everything else, including `-dipole_moments`, writes to the
log.

## Pages

* **[Command reference](Command-Reference.md)** — every flag, by family, with
  arguments and defaults, taken from the option parser.
* **[Use cases](Use-Cases.md)** — numbered, copy-pasteable PowerShell
  invocations for the common jobs.
* **[Troubleshooting](Troubleshooting.md)** — the guards the code actually
  raises, the flag-order rules, and the known traps.

`.\NoSpherA2.exe -h` prints a self-contained help text with the same
information in a single screenful-per-section form. It is generated from the
same source tree as this wiki (`Src/core/convenience.cpp`), so when the two
disagree, trust `-h` on the binary you are running.

## Verbosity

| Flag | Effect |
| --- | --- |
| `-v` | Debug output. |
| `-v2` | Debug output plus per-step progress counts. |
| `-debug` | Same as `-v`. |
| `-h`, `--h`, `-help`, `--help` | Print banner, full help and build date to the console, then exit. |
| `-no_date` | Suppress timings, the build date **and** the "GPU in use" notes. Used by the test suite so goldens do not churn. |
| `-no_date_but_gpu` | As `-no_date`, but keep the GPU notes — use this when you need to know whether the GPU actually ran. |

## Where to cite

NoSpherA2 prints its own literature. Every method that runs emits one or more
lines of the form

```
[QTAIM] Bader, Chem. Rev. 91 (1991) 893, DOI 10.1021/cr00005a013
```

into the log, for exactly the methods that were used in that run — the
partitioning scheme, the grid, the ECP treatment, each property, the NBO/NRT
analysis, SALTED, the dispersion model, and so on. Copy those lines out of your
log rather than guessing from a table here.

The reference list lives in `Src/core/citations.cpp` (40 methods, DOIs
Crossref-checked on 24 Sep 2026) and is the single place it is maintained. If a
method you used did not print a citation, that is a bug worth reporting.

Regardless of the run, NoSpherA2 itself is
`Kleemiss et al., Chem. Sci. 12 (2021) 1675` — printed under the `[NoSpherA2]`
tag.
