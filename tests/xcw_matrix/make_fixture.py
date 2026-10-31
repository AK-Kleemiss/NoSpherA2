"""Generate the tiny XCW matrix fixture: NH3 in P2_12_12_1 (4 ops, non-centrosymmetric).

The ammonia structure (P2_13, a = 5.1305 A, tests/cytidine_tonto/nh3_XCW) described in its subgroup
P2_12_12_1, so the whole molecule is the asymmetric unit (4 atoms, 8 sto-3g functions) and 4 operations
remain. The grown model adds the neighbour that N1-H1 points at (8 atoms).

Truth: a def2-SVP Hartree-Fock density (XCW lambda 0, ungrown) with Gram-Charlier C and D on N1 and
SHELX extinction x = EXT_TRUE, its sqrt(y)|F_calc| taken from NA2_0000000_Fcalc.txt. I_obs is
(K * that)^2 plus Gaussian noise of sigma(I) = sqrt(G*I + (P*I)^2 + B^2), fixed seed.

    python make_fixture.py --exe D:/git/NoSpherA2/spaces/3f316b01-r1sym/NoSpherA2.exe
"""
import argparse
import itertools
import math
import pathlib
import shutil
import subprocess

import numpy as np

A = 5.1305
WAVELENGTH = 0.71073
STL_MAX = 0.80
EXT_TRUE = 0.2          # SHELX x on NoSpherA2's absolute |F_calc| scale
K, G, P, B = 10.0, 1.0, 0.01, 2.0
SEED = 31

OPS = ["+X,+Y,+Z", "1/2-X,-Y,1/2+Z", "1/2+X,1/2-Y,-Z", "-X,1/2+Y,1/2-Z"]
ATOMS = [("N1", "N", (0.2103, 0.2103, 0.2103)),
         ("H1", "H", (0.3722, 0.2627, 0.1113)),
         ("H2", "H", (0.1113, 0.3722, 0.2627)),
         ("H3", "H", (0.2627, 0.1113, 0.3722))]
U_N = (0.0372, 0.0372, 0.0372, -0.0009, -0.0009, -0.0009)  # U11 U22 U33 U23 U13 U12
U_H = 0.053
# Gram-Charlier on N1 only, CIF order 111 112 113 122 123 133 222 223 233 333 / 1111 ... 3333
C_N = (2.0e-6, 4.0e-7, -3.0e-7, -5.0e-7, 6.0e-7, 3.0e-7, -1.5e-6, -4.0e-7, 5.0e-7, 1.8e-6)
D_N = (5.0e-8, 6.0e-9, -8.0e-9, 1.2e-8, 4.0e-9, 1.0e-8, -7.0e-9, 9.0e-9, -5.0e-9, 6.0e-9,
       4.0e-8, -6.0e-9, 1.1e-8, 8.0e-9, 3.0e-8)


def ops_matrices():
    out = []
    for s in OPS:
        R, t = np.zeros((3, 3)), np.zeros(3)
        for i, comp in enumerate(s.split(",")):
            for term in comp.replace("-", "+-").split("+"):
                if not term:
                    continue
                if term.lstrip("-") in "XYZ":
                    R[i, "XYZ".index(term.lstrip("-"))] = -1 if term.startswith("-") else 1
                else:
                    t[i] = float(eval(term))
        out.append((R, t))
    return out


def write_cif(path, anharm):
    lines = ["data_nh3_p212121",
             "_space_group_crystal_system        'orthorhombic'",
             "_space_group_name_H-M_alt          'P 21 21 21'",
             "_space_group_name_Hall             ' P 2ac 2ab'",
             "loop_", "  _space_group_symop_id", "  _space_group_symop_operation_xyz"]
    lines += [f" {i + 1} {op}" for i, op in enumerate(OPS)]
    lines += ["", f"_cell_length_a                     {A}", f"_cell_length_b                     {A}",
              f"_cell_length_c                     {A}", "_cell_angle_alpha                  90",
              "_cell_angle_beta                   90", "_cell_angle_gamma                  90",
              f"_cell_volume                       {A ** 3:.4f}", "_cell_formula_units_Z              4",
              f"_diffrn_radiation_wavelength       {WAVELENGTH}", "loop_"]
    lines += ["  _atom_site_" + k for k in ("label", "type_symbol", "fract_x", "fract_y", "fract_z",
                                            "U_iso_or_equiv", "adp_type", "occupancy",
                                            "refinement_flags_posn", "site_symmetry_order", "disorder_group")]
    for lab, el, (x, y, z) in ATOMS:
        u, kind = (sum(U_N[:3]) / 3, "Uani") if el == "N" else (U_H, "Uiso")
        lines.append(f" {lab} {el} {x:.4f} {y:.4f} {z:.4f} {u:.5f} {kind} 1 . 1 .")
    lines += ["", "loop_"] + ["  _atom_site_aniso_" + k for k in
                              ("label", "U_11", "U_22", "U_33", "U_23", "U_13", "U_12")]
    lines.append(" N1 " + " ".join(f"{u:.5f}" for u in U_N))
    if anharm in ("C", "CD"):
        lines += ["", "loop_", "  _atom_site_anharm_GC_C_label"]
        lines += [f"  _atom_site_anharm_GC_C_{''.join(map(str, c))}"
                  for c in itertools.combinations_with_replacement((1, 2, 3), 3)]
        lines.append("  N1 " + " ".join(f"{c:.3g}" for c in C_N))
    if anharm == "CD":
        lines += ["", "loop_", "  _atom_site_anharm_GC_D_label"]
        lines += [f"  _atom_site_anharm_GC_D_{''.join(map(str, c))}"
                  for c in itertools.combinations_with_replacement((1, 2, 3), 4)]
        lines.append("  N1 " + " ".join(f"{d:.3g}" for d in D_N))
    lines += ["", "loop_", "  _atom_type_symbol", "  _atom_type_scat_dispersion_real",
              "  _atom_type_scat_dispersion_imag", "  N 0 0", "  H 0 0", ""]
    path.write_text("\n".join(lines))


def grown_xyz(path):
    """Asymmetric unit plus the image molecule whose N is nearest to H1 (the N-H...N contact)."""
    frac = np.array([a[2] for a in ATOMS])
    best = None
    for o, (R, t) in enumerate(ops_matrices()):
        for shift in itertools.product((-1, 0, 1), repeat=3):
            if o == 0 and shift == (0, 0, 0):
                continue
            img = frac @ R.T + t + np.array(shift)
            d = A * np.linalg.norm(img[0] - frac[1])
            if best is None or d < best[0]:
                best = (d, img)
    allf = np.vstack([frac, best[1]])
    els = [a[1] for a in ATOMS] * 2
    body = [f"{e} {x * A:.6f} {y * A:.6f} {z * A:.6f}" for e, (x, y, z) in zip(els, allf)]
    path.write_text(f"{len(els)}\nNH3 + neighbour (H1...N {best[0]:.3f} A), Angstrom\n" + "\n".join(body) + "\n")
    return best[0]


def hkl_list():
    hmax = int(2 * STL_MAX * A) + 1
    out = []
    for h, k, l in itertools.product(range(hmax + 1), repeat=3):
        if (h, k, l) == (0, 0, 0) or math.sqrt(h * h + k * k + l * l) / (2 * A) > STL_MAX:
            continue
        # P2_12_12_1 absences
        if (k == 0 and l == 0 and h % 2) or (h == 0 and l == 0 and k % 2) or (h == 0 and k == 0 and l % 2):
            continue
        out.append((h, k, l))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--work", default="E:/Temp/xcw_matrix_scratch/gen")
    ap.add_argument("--truth-basis", default="def2-svp")
    args = ap.parse_args()
    here = pathlib.Path(__file__).resolve().parent
    work = pathlib.Path(args.work)
    work.mkdir(parents=True, exist_ok=True)

    for tag, name in (("none", "nh3.cif"), ("C", "nh3_C.cif"), ("CD", "nh3_CD.cif")):
        write_cif(here / name, tag)
    print(f"grown neighbour H1...N {grown_xyz(here / 'nh3_grown.xyz'):.3f} A")

    refl = hkl_list()
    shutil.copy(here / "nh3_CD.cif", work / "truth.cif")
    (work / "dummy.hkl").write_text("".join(f"{h:4d}{k:4d}{l:4d}{100.0:8.2f}{1.0:8.2f}\n" for h, k, l in refl))
    (work / "truth_settings.txt").write_text(
        f"normal normal_conv params 22 basis_set {args.truth_basis} F charge 0 mult 1 rhf "
        f"start 0 step_size 0.01 end 0 extinction shelx fixed {EXT_TRUE}\n")
    subprocess.run([args.exe, "-do_XCW", "-no_gpu", "-no_date", "-print_F", "-cpus", "8", "-cif", "truth.cif",
                    "-hkl", "dummy.hkl", "-XCW_settings", "truth_settings.txt"], cwd=work, check=True,
                   stdout=subprocess.DEVNULL)
    fc_lines = (work / "NA2_0000000_Fcalc.txt").read_text().splitlines()
    scale = float(fc_lines[0].split("scale =")[1])
    F_true = {tuple(map(int, l.split()[:3])): float(l.split()[5]) / scale for l in fc_lines[1:] if l.strip()}

    rng = np.random.default_rng(SEED)
    out = []
    for h, k, l in refl:
        I = (K * F_true[(h, k, l)]) ** 2
        s = math.sqrt(G * I + (P * I) ** 2 + B * B)
        out.append(f"{h:4d}{k:4d}{l:4d}{I + rng.normal(0.0, s):8.2f}{s:8.2f}\n")
    (here / "nh3.hkl").write_text("".join(out))
    print(f"{len(out)} reflections, truth scale-free |F| from {args.truth_basis}, K = {K}, x_ext = {EXT_TRUE}")


if __name__ == "__main__":
    main()
