# EQC energy decomposition

`-eqc` and `-eqc_wfn` split the energy of a bond-forming reaction,
fragments → parent, into an electronic and an electrostatic part, following
Rahm and Hoffmann's energy decomposition based on the average electron energy
(EQC). The numbers are printed the way Martin Rahm's
[X-analysis](https://github.com/martinrahm/X-analysis) prints them with `-m 2`,
and on the same wavefunctions they agree with it to the printed digits (the
eV factor is cclib's 27.21138505). Flags are listed in the
[command reference](Command-Reference.md#eqc-energy-decomposition).

## The idea

For any wavefunction with n electrons the total energy splits exactly into

    E = n·X̄ + (Vnn − Eee)

* **n·X̄** is the sum of occupation × orbital energy over all occupied spin
  orbitals. X̄ is the average binding energy of an electron in the system; for
  a free atom's valence shell it is the ground-state electronegativity of Rahm,
  Zeng and Hoffmann (2019).
* **Vnn** is the nuclear repulsion.
* **Eee** is the electron–electron repulsion. It is not computed separately:
  the sum of orbital energies counts every electron–electron interaction
  twice, and the decomposition defines Eee as the remainder,
  Eee = Vnn − (E − n·X̄).

For a reaction A + B + … → C every term is taken as products minus reactants:

    ΔE = Δ(n·X̄) + Δ(Vnn − Eee)

Δ(n·X̄) is how much more (negative) or less tightly the electrons are bound on
average after the bond forms. Δ(Vnn − Eee) is the change of the classical
repulsions between the nuclei and between the electrons. Rahm and Hoffmann
(2015, 2016) read a bond dominated by Δ(n·X̄) as covalent in character and one
dominated by Δ(Vnn − Eee) as electrostatic.

## What the output means

Each run writes `<stem>.eqc_log` (and the same sections into
`NoSpherA2.log`). The citations are at the top of it.

**Terms per wavefunction**, in eV, one row each for the parent and every
fragment: n, E, n·X̄, Vnn and Eee, the SCF iterations, the iterations from
occ's own guess (only with `-eqc_cold`) and where the wavefunction came from.

**Fragments → parent**, products minus reactants, in eV:

| Line | Definition |
| --- | --- |
| `Delta E` | ΔE, the reaction energy; negative = the parent is lower. |
| `Delta(nX-bar)` | Δ(n·X̄). |
| `Delta(Vnn-Eee)` | Δ(Vnn − Eee) = ΔE − Δ(n·X̄). |
| `Delta(E/n)`, `Delta(X-bar)`, `Delta(Vnn-Eee)/n` | the three lines above divided by the parent's n. |
| `Delta(Vnn/n)`, `Delta(Eee/n)` | ΔVnn/n and ΔEee/n separately. Both are large and almost cancel; only their difference is meaningful. |
| `Q` | 2 Δ(n·X̄)/ΔE − 1. Q = 1 when all of ΔE is Δ(n·X̄), 0 when the two terms are equal, −1 when Δ(n·X̄) is zero; outside [−1, 1] when the two terms have opposite signs. |
| `Delta\|Eee/E\| (%)` | 100 (\|Eee/E\| of the parent − \|ΣEee/ΣE\| of the fragments). |
| `Covalency (%)` | X-analysis' covalency index, x = Δ(n·X̄), v = Δ(Vnn − Eee): 100\|x\|/(\|v\|+\|x\|) for x < 0, 100 v/(v − x) for x > 0, undefined at x = 0. |
| `Ionicity (%)` | 100 − covalency. |

### Worked example: ethane → two methyl radicals

The `eqc_ethane` integration test (ORCA 6.1 HF/def2-SVP gbw, Mode 1):

    NoSpherA2 -eqc ethane.gbw -eqc_frag 0,2-4 0 2 1,5-7 0 2

    Delta E               -3.694453
    Delta(nX-bar)         -2.197440
    Delta(Vnn-Eee)        -1.497013
    Q                      0.189589
    Covalency (%)         59.479440

The C–C bond is 3.69 eV (HF, no correlation, no relaxation of the methyl
geometry); 2.20 eV of it comes from the electrons being more tightly bound
on average and 1.50 eV from the repulsion terms, so the index puts it at
59 % covalent. X-analysis `-m 2` on the same orbitals gives Δ(n·X̄)
−2.197319 eV.

The heterolytic split of CH3F into CH3⁺ and F⁻ (`eqc_ch3f_pbe0`, PBE0/def2-SVP,
`-eqc_frag 0,2-4 1 1 1 -1 1 -eqc_method pbe0`) is a case where the two terms
have opposite signs: ΔE −14.76 eV is the small difference of Δ(n·X̄)
−79.52 eV and Δ(Vnn − Eee) +64.76 eV, so Q is 9.78 and the covalency 55 %.
The same reaction from ORCA's own wfx files in Mode 2 (`eqc_ch3f_wfx`,
`-eqc_wfn parent.wfx ch3p.wfx fm.wfx`) gives ΔE −14.7564 against −14.7566 eV;
the parent energies of occ and ORCA already differ by 2e-4 hartree (their DFT grids).

## The two modes

**Mode 1** (`-eqc <gbw> -eqc_frag <atoms> <charge> <mult> ...`). Every term
comes from one program, occ, at one level of theory:

1. the ORCA gbw is read and occ re-converges the parent from its orbitals
   (the log prints E(occ) at the gbw orbitals, the iterations needed and, if
   `<stem>.out` sits beside the gbw, E(ORCA) for comparison);
2. each fragment, with the charge and multiplicity you give, is converged
   from the block of the parent's density on its atoms.

The fragments keep the parent's geometry: ΔE is a bond energy at fixed
nuclei, not a reaction energy with relaxed fragments. Starting from the
parent's density puts each fragment in the state that is adiabatically
connected to the bond; `-eqc_cold` additionally converges from occ's own
guess, and a different iteration count there is a hint to check the
fragment state. `-eqc_method` must name the method the gbw was made with
(`hf` or a functional), and `-eqc_basis` is needed when the gbw uses ECPs.

**Mode 2** (`-eqc_wfn <parent> <frag> <frag> ...`). Nothing is converged;
E, the orbital energies and the geometry are read from the files, so the
fragments can be relaxed or computed with any program. `.wfx` and `.fchk`
carry their own total energy; `.wfn`, `.gbw` and `.molden` take it from the
ORCA `.out` of the same name. All files must be at the same level of theory.

## Caveats

* **DFT.** With a functional the code, like X-analysis, uses the plain sum of
  Kohn–Sham orbital energies as n·X̄. Racioppi, Lolur, Hyldgaard and Rahm
  (2023) show that the average electron energy in DFT is not that sum and
  derive the correct expression; it is not implemented here. HF numbers are
  the reference; DFT numbers are comparable only with other DFT numbers
  computed the same way.
* **ECPs.** Vnn uses the full atomic numbers and n counts the core electrons
  the ECP replaces, as cclib does; the core terms cancel in a reaction only
  when the same atoms carry the same ECP on both sides.
* **The index branches.** The covalency formula changes form with the sign of
  Δ(n·X̄) and is undefined at zero (printed as nan); near zero it jumps.

## Background

* M. Rahm, R. Hoffmann, "Toward an Experimental Quantum Chemistry: Exploring
  a New Energy Partitioning", *J. Am. Chem. Soc.* **137** (2015) 10282.
  [10.1021/jacs.5b05600](https://doi.org/10.1021/jacs.5b05600)
* M. Rahm, R. Hoffmann, "Distinguishing Bonds", *J. Am. Chem. Soc.* **138**
  (2016) 3731. [10.1021/jacs.5b12434](https://doi.org/10.1021/jacs.5b12434)
* M. Rahm, T. Zeng, R. Hoffmann, "Electronegativity Seen as the Ground-State
  Average Valence Electron Binding Energy", *J. Am. Chem. Soc.* **141** (2019)
  342. [10.1021/jacs.8b10246](https://doi.org/10.1021/jacs.8b10246)
* S. Racioppi et al. (last author M. Rahm), "In-Situ Electronegativity and the Bridging of
  Chemical Bonding Concepts", *Chem. Eur. J.* **27** (2021) 18156.
  [10.1002/chem.202103477](https://doi.org/10.1002/chem.202103477)
* S. Racioppi, P. Lolur, P. Hyldgaard, M. Rahm, "A Density Functional Theory
  for the Average Electron Energy", *J. Chem. Theory Comput.* **19** (2023)
  799. [10.1021/acs.jctc.2c00899](https://doi.org/10.1021/acs.jctc.2c00899)
* X-analysis, the reference implementation:
  <https://github.com/martinrahm/X-analysis>

The first three are printed in every `.eqc_log`.
