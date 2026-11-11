# Correlated wavefunctions

How to give NoSpherA2 a correlated density (MP2, CCSD, CCSD(T), CI, CASSCF), which files it reads from
which program, and what each analysis means once the density is no longer a single determinant.

A correlated density reaches NoSpherA2 as **natural orbitals with fractional occupations**: the eigenvectors and
eigenvalues of the one-particle density matrix (1-RDM). Every reader that keeps per-orbital occupations can carry
them; the density, and everything built from the density alone, is then the correlated one. Nothing in the files
NoSpherA2 reads carries the two-particle density matrix (2-RDM), so quantities that need it are either approximated
from the 1-RDM or refused. The tables below say which.

## Contents

1. [Which files carry a correlated density](#1-which-files-carry-a-correlated-density)
2. [ORCA](#2-orca)
3. [occ inside NoSpherA2 (`-occ`)](#3-occ-inside-nosphera2--occ)
4. [PySCF, Molpro and other molden writers](#4-pyscf-molpro-and-other-molden-writers)
5. [Gaussian](#5-gaussian)
6. [What is exact, what is approximate, what is refused](#6-what-is-exact-what-is-approximate-what-is-refused)
7. [E2(D) and the reference SCF (`-nbo_ref`)](#7-e2d-and-the-reference-scf--nbo_ref)
8. [The natural-occupation header](#8-the-natural-occupation-header)
9. [How good the Müller DI is](#9-how-good-the-müller-di-is)

---

## 1. Which files carry a correlated density

"Verified" means a run in the November 2026 validation (52 `-fba` runs against exact 2-RDM references) or the
October 2026 scattering-factor comparison against ORCA 6.1. Everything else is stated as not verified.

| Source | File given to `-wfn` / `-fba` / `-nbo` | Status |
| --- | --- | --- |
| ORCA relaxed MP2 | `job.mp2nat` | verified (scattering factors, `-fba`) |
| ORCA CCSD / CISD (MDCI) | `job.mdci.nat` | verified (scattering factors, `-fba`) |
| ORCA CASSCF | the CASSCF `.gbw` (with `ActOrbs NatOrbs`) | verified (`-fba`) |
| ORCA full CI (CASSCF over every orbital) | the `.gbw`, with `job.densities` / `job.densitiesinfo` beside it | verified (`-fba`, LiH) |
| ORCA AutoCI CCSD(T), MP2-MP4 | `job.json` from `orca_2json` | verified for scattering factors; `-fba` not run on it |
| occ (`-occ input.toml`) | none: the density stays in memory | verified for scattering factors; not reachable by `-fba`/`-nbo` |
| PySCF natural orbitals | `.molden` from `pyscf.tools.molden.from_mo(..., occ=n)` | verified (`-fba`, FCI references) |
| Molpro, Psi4, other molden writers | `.molden` with `Occup=` per orbital | not verified |
| any program writing `.wfx`/`.wfn` natural orbitals | `.wfx`, `.wfn` | not verified; primitive files, so no NBO/NPA/RGBI/E2(D) |
| Gaussian `Density=Current` | `.fchk` | **not read**: see [Gaussian](#5-gaussian) |
| ORCA DLPNO densities, RI-MP2 natural orbitals, relaxed CCSD | | not verified |

Which tools need what is unchanged by correlation: NBO, NPA, RGBI, E2(D) and `-ibo` need a contracted basis
(`.gbw`, `.molden`, `.fchk`, `.json`); `.wfn`/`.wfx` give the density, cubes, QTAIM and ELI-D only.
`-ibo` is closed-shell SCF only.

## 2. ORCA

The ORCA inputs below are the ones the validation ran (ORCA 6.1, `TightSCF`). Each correlated job also writes
`job.gbw` with the SCF orbitals; leave it beside the natural-orbital file, because NBO uses it as the reference SCF
for E2(D) ([section 7](#7-e2d-and-the-reference-scf--nbo_ref)).

**Relaxed MP2** writes `job.mp2nat`:

```
! MP2 def2-TZVP TightSCF
%mp2
 density relaxed
 natorbs true
end
```

**CCSD** (and **CISD** the same way) writes `job.mdci.nat`:

```
! CCSD def2-TZVP TightSCF
%mdci
 density unrelaxed
 natorbs true
end
```

MDCI's CCSD(T) run still writes the *CCSD* density, without (T). NoSpherA2 reads `job.property.txt` when it lies beside
the file, prints the method and level of the density it found, and warns about this case, about an unrelaxed MP2
density (no orbital response; `%mp2 Density Relaxed` gives it) and about a CC density that is not the unrelaxed one.

**CCSD(T)**, and unrelaxed MP2/MP3/MP4, come from AutoCI through `orca_2json`:

```
! AUTOCI-CCSD(T) NoRI KeepDens
%autoci
 density unrelaxed
end
```

then `orca_2json job.gbw` with `job.json.conf` = `{"MOCoefficients": true, "Densities": ["autocipur"]}` and
`-wfn job.json`. NoSpherA2 takes the `autocipur`, `autocipre` or `mdcip` density from the JSON and diagonalises it in
the orbitals of `job.gbw`; RHF references only. On Windows AutoCI is a separate download
(`Orca6.1.0.Win64_autoci*.zip`, unpacked into the ORCA folder). The CCSD(T)/MPn methods and their accuracy for
scattering factors are in [Command reference: correlated densities](Command-Reference.md#correlated-densities-mp2-ccsd-ccsdt).

**CASSCF** writes its natural orbitals and occupations into `job.gbw` itself:

```
! CASSCF cc-pVDZ TightSCF
%casscf
 nel 6
 norb 6
 mult 1
 ActOrbs NatOrbs
end
```

Because the CASSCF `.gbw` takes the name the RHF reference would have, rename it (e.g. `job.cas.gbw`) and put the RHF
`.gbw` of the same geometry and basis beside it as `job.gbw` if you want E2(D). A `%moinp` / `MORead` chain along a
scan is fine; check the occupations in the [header](#8-the-natural-occupation-header) at every point, since a chain
can keep a σ pair doubly occupied while only the π pairs break.

**Full CI** in ORCA is CASSCF with every orbital active (LiH cc-pVDZ: `nel 4`, `norb 19`). ORCA then skips the orbital
update and leaves the *guess* orbitals in the `.gbw` with diag(γ) as their occupations, which is not the CI density.
NoSpherA2 detects this: when `job.densities` and `job.densitiesinfo` lie beside the `.gbw` it prints a WARNING with
the largest off-diagonal element of γ and uses the natural orbitals of the CI density matrix instead (LiH Li–H DI
0.1437 from the diagonal, 0.1843 corrected, exact 0.1857). Without those files it prints a NOTE and keeps the
diagonal; a densities file that does not match the `.gbw` is refused. Keep `job.densities*` next to the `.gbw` for
this case.

## 3. occ inside NoSpherA2 (`-occ`)

`-occ input.toml` with `method = "mp2"` (relaxed), `"mp3"`, `"mp4(sdq)"`, `"mp4"`, `"ccsd"` or `"ccsd(t)"`
(unrelaxed) runs RHF and the correlated density inside NoSpherA2, and the natural orbitals go straight on to the
scattering factors. Methods, backends and `[scf]` keys: [Command reference](Command-Reference.md#correlated-densities-mp2-ccsd-ccsdt).

These natural orbitals are not written to a file, so the analysis tools (`-fba`, `-nbo`, `-eli_analysis`,
`-topology`) cannot read them. `-fchk` does write the wavefunction, but an fchk has no occupation field and the
occupations come back as integers on reading, so it is not a route to the correlated density. For bonding analysis
of a correlated density use an ORCA or molden file.

## 4. PySCF, Molpro and other molden writers

A molden file with an `Occup=` line per orbital carries natural orbitals. What the validation ran (PySCF FCI, CASCI
and CCSD 1-RDMs):

```python
from pyscf.tools import molden
n, U = numpy.linalg.eigh(dm1_mo)            # 1-RDM in the MO basis, spin-summed
n, U = n[::-1], U[:, ::-1]
molden.from_mo(mol, "fci.molden", mf.mo_coeff @ U, ene=numpy.zeros(len(n)), occ=n)
```

then `NoSpherA2 -fba fci.molden`. Orbital energies of natural orbitals mean nothing; writing zeros is fine, because
NoSpherA2 never uses them on correlated input. A spin-summed (restricted) 1-RDM is split into equal α and β halves.
For E2(D) give the SCF of the same molecule and basis to `-nbo` with `-nbo_ref hf.molden`; that combination has not
been run in the validation, which ran `-fba` without a reference.

Molpro, Psi4, Multiwfn and other molden writers should work the same way if they write fractional `Occup=` values;
none of them has been checked.

## 5. Gaussian

`read_fchk` reads the SCF MO coefficients and does not look at the `Total MP2 Density`, `Total CC Density` or
`Total CI Density` arrays. **A Gaussian fchk of a correlated job is therefore analysed as the SCF it started from**,
with no message. Until that is implemented, run the correlated step in ORCA or PySCF, or write the natural orbitals
with occupations to a molden or wfx file from Gaussian or Multiwfn (not verified).

## 6. What is exact, what is approximate, what is refused

"Exact" means exact for the density in the file: it depends on the 1-RDM only, so the correlated density gives the
correlated value with no further approximation.

| Quantity | On natural-orbital input | Why |
| --- | --- | --- |
| ρ, cubes, ESP, Laplacian, ρ topology (`-topology`) | exact | depend on ρ only |
| scattering factors, `.tsc`/`.tscb` | exact | Fourier transform of ρ |
| QTAIM populations and charges | exact | integrals of ρ |
| NPA, NAO occupancies, NBO occupancies | exact | from the 1-RDM in the NAO basis |
| delocalisation index DI and localisation index λ | **exact** with an ORCA determinant print beside the gbw (below), else **approximate**: Müller pair density, √(n_i n_j) | the exact DI needs the 2-RDM; see [section 9](#9-how-good-the-müller-di-is) |
| ELI-D, ELI-q (and their α/β/triplet members) | **approximate**: the single-determinant formulas evaluated with the natural orbitals and their occupations; the log prints a WARNING when no 2-RDM was found | the correlated ELI needs the 2-RDM |
| RGBI | **approximate**: a 1-RDM bond index on the natural orbitals | no 2-RDM behind it |
| E2 (second-order donor-acceptor energies) | **refused**, with a message | natural orbitals have no orbital energies; replaced by the density transfer and E2(D) ([section 7](#7-e2d-and-the-reference-scf--nbo_ref)) |
| NRT | fed by E2(D) instead of E2 | |
| ELIA | **exact** with an ORCA determinant print beside the gbw, else **not available** | needs the αβ block of the 2-RDM, which no wfn/wfx/molden/fchk/gbw carries; `-eli_family` says so |

**Exact DI/λ and ELIA from an ORCA CASSCF or CI.** No gbw or `.densities` file holds a CI vector or a 2-RDM; the
determinant print in the output does. Run the CAS with

```
%casscf ... PrintWF det
 ci TPrintWF 1e-14 end
end
```

and keep the output as `<stem>.out` (or `.log`) next to `<stem>.gbw`.

**Required input**, all in one folder with one stem:

| File | Needed for | Notes |
| --- | --- | --- |
| `<stem>.gbw` | always | the file given to `-fba` / `-eli_family`; restricted orbitals |
| `<stem>.out` or `<stem>.log` | exact DI/λ, ELIA, ELI-q_s,2 | ORCA's stdout from the same job: `PrintWF det`, `TPrintWF 1e-14`, `nroots 1` (no state average), at most 64 active orbitals |
| `<stem>.densities` + `<stem>.densitiesinfo` | full CI only | the natural orbitals over every orbital; not needed for a CAS with internal orbitals |

Without the `.out`, or with any of the refusals below, the run continues on the Müller DI and `-eli_family` prints
no ELIA column, with the reason in the log.

ORCA has no binary file to use instead. The CASSCF module does write the CI vector
(`<stem>.CIData.tmp.CIVEC_mult-*_root-*.tmp`, in the CSF basis) and the spin-summed active 2-RDM (`<stem>.casg2.tmp`,
n²×n² doubles, trace N(N−1)/2) during the run, but deletes both at the end. No input keyword keeps them, and neither
the gbw, `.densities` nor `orca_2json` carries a 2-RDM. The determinant print is the only persistent source.

`-fba` and `-eli_family` then build the spin 2-RDM from the determinants and log

```
  2-RDM: 10554 determinants from lih_fci.out, 0 internal + 19 active orbitals, |1 - norm| 7.1e-11, 1-RDM vs occupations 1.0e-09
  pair density: exact, from the CI 2-RDM
```

`-eli_family` adds the columns `ELIA` = ρ₂^αβ(r,r)/(ρ_α ρ_β) and `ELI-q_s,2` = 4ρ₂^αβ(r,r)/ρ², the singlet ELI-q
from the 2-RDM; the other members stay on the 1-RDM. The 2-RDM is refused, with the reason in the log and the Müller
path used instead, for a multi-root (state-averaged) print, an unrestricted gbw, more than 64 active orbitals, a CI norm off by more than
1e-6 (TPrintWF too loose), a full CI whose natural orbitals had to be rebuilt while internal orbitals are present, and a CI 1-RDM that misses the
gbw occupations by more than 1e-4 (truncated print, other orbitals). LiH cc-pVDZ full CI: ELIA at six points agrees
with an independent PySCF FCI to 1e-3; DI(Li,H) exact 0.1855 vs Müller 0.1843.

For an SCF wavefunction every one of these is bit-identical to what it was before correlated input was supported:
√(n_i n_j) = n_i n_j for occupations 0 and 1, and the 1-RDM is idempotent.

The Müller DI keeps the sum rule λ_A + ½ Σ_B DI(A,B) = N_A, which the plain product n_i n_j loses on fractional
occupations; the `-fba` log closes it to 1e-4 e on all validation runs. `-fba -v` also writes the atomic overlap
matrices (`qtaim.aom`) for an external 2-RDM DI.

A closed-shell NRT shows a dissociated singlet (H2 at 5 Å) as a 50/50 pair of mirror ionic structures. That is the
right answer for a spin-summed 1-RDM; it is not a diradical structure.

## 7. E2(D) and the reference SCF (`-nbo_ref`)

On a correlated 1-RDM, NBO replaces the E2 table by the electron transfer the density actually shows. For each Lewis
donor k and non-Lewis acceptor j:

- **q(k→j) = D_jk² / (c − D_jj)**, from the 1-RDM D in the NBO basis (c the full occupancy of the spin), with the
  acceptor budget s (how much of each acceptor's occupation the transfers explain) and, per bond, a local 2×2
  natural-occupation pair n₊/n₋ with b_eff = (n₊ − n₋)/2.
- **E2(D) = q (F_jj − F_kk)**, the transfer priced on the Fock matrix of a reference SCF. A bond's own antibond gets
  n₋ (F_aa − F_bb). For a single determinant E2(D) reproduces E(2) to second order; on the validation molecules it
  lies within 0.5–1.5× of the SCF E2.

The reference SCF is

- `-nbo_ref <file>` after `-nbo`, any reader with a contracted basis (`.gbw`, `.molden`, `.fchk`), or
- by default the sibling `.gbw`: NoSpherA2 strips one extension at a time and takes the first `.gbw` it finds
  (`job.mdci.nat` → `job.mdci.gbw` → `job.gbw`; `job.mp2nat` → `job.gbw`). `-fba` has no `-nbo_ref` and always uses
  the sibling.

Its overlap matrix must match the input's to 1e-8 (same basis and geometry) and its spin treatment must be the same;
otherwise E2(D) is skipped and the log says why:

```
E2(D) of the correlated density: Fock of the reference SCF job.gbw
E2(D) of the correlated density: none, no reference SCF (-nbo_ref, or <base>.gbw beside the input)
```

The JSON carries the transfers in its `density_transfer` block, with `e2d_kcal` when a reference was found.

## 8. The natural-occupation header

`-nbo` and `-fba` print a summary of the occupations before the analysis. CCSD water (def2-TZVP):

```
Natural occupation numbers: n from 0.00000 to 2.00000, sum of negative n 0.00000
  N_U (Head-Gordon) 0.25066 e, N_D (Takatsuka) 0.49513 e, HONO/LUNO 1.96246/0.02655, y 0.00053
  idempotent: no, a correlated 1-RDM: DI/lambda/ELIA exact from a printed CI 2-RDM, else Mueller; E2 is refused
```

| Entry | Meaning |
| --- | --- |
| n range, sum of negative n | negative occupations come from unrelaxed or non-variational densities; small ones are harmless, large ones mean the density is not N-representable |
| N_U (Head-Gordon) | Σ min(n, 1 − n) per spin orbital: the number of effectively unpaired electrons |
| N_D (Takatsuka) | Σ 2n(1 − n) per spin orbital: the odd-electron count |
| HONO/LUNO, y | highest occupied / lowest unoccupied natural orbital and Yamaguchi's diradical character y (0 closed shell, 1 pure diradical); restricted input only |
| idempotent | yes when N_D < 1e-6 e (SCF); no switches on the Müller DI, the density transfers and E2(D) |

Weak correlation (single-reference molecules near equilibrium) gives N_U of a few tenths of an electron and y near 0;
a breaking bond or a diradical gives N_U approaching 2 and y approaching 1 (p-benzyne CASSCF: N_U 1.668).

## 9. How good the Müller DI is

Müller DI against the DI from the exact 2-RDM of the same wavefunction, on the same atomic overlap matrices
(PySCF FCI for H2, BH and Li2; frozen-core CASCI for HF; CASSCF(6,6) along the N2 stretch; the CCSD-Lambda 2-RDM for
H2O, NH3, FH, CO, HCN and N2; cc-pVDZ, H2 def2-TZVP):

| Pair type | Müller vs exact |
| --- | --- |
| single bonds near equilibrium (H2, O–H, N–H, F–H, C–O) | within ±3 % |
| triple bonds near equilibrium (N2, HCN C–N), N2 to 1.6 Å | +10 to +25 % |
| stretched bonds (BH 2.5 Å, N2 2.5–3.0 Å) | −7 to −11 %, ≤ 0.04 absolute |
| small non-bonded pairs (H–H in H2O/NH3, Li–Li) | −5 to −63 %, ≤ 0.03 absolute |

The sign of the change from HF is always right. On a triple bond Müller recovers about two thirds of the drop the
exact DI shows (N2 CCSD: HF 3.03, Müller 2.31, exact −32 % from HF). These are limits of the 1-RDM functional, not of
the implementation: closing them needs the 2-RDM.
