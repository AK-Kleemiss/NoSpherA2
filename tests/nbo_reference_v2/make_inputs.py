"""Rebuild the ORCA inputs for the 22 NBO reference molecules from idealised starts.

    python make_inputs.py <out_dir> [--nprocs N]

The wavefunctions behind tests/nbo_reference/ are gone; only index.json's record of what
the programs produced survived.  This script writes one <name>/<name>.inp per molecule at
exactly the keyword line, charge and multiplicity index.json recorded, starting from a
textbook geometry.  The geometry does not have to match the original start: the acceptance
test in accept.py does, and it is what makes a regenerated wavefunction usable as a stand-in
(basis_functions exactly equal, final energy within tolerance).  A start that relaxes to a
different minimum fails that test loudly instead of quietly becoming a different molecule.

ch3 is the exception and the calibration point: tests/nbo_reference/ch3_orca_reference.47
survived, so ch3 starts from the *optimised* geometry of the original run, read out of that
archive's $COORD block.  Its regenerated energy therefore says how well the rest of the
pipeline (ORCA version, basis, grid defaults) reproduces the original at all.

%pal goes on a line after the keyword line, so ORCA's echo of line 1 - which is what
collect_reference.py's `keyword_line` regex reads - stays byte-identical to the recorded one.
"""
import argparse
import json
import math
import os

KEYWORD_LINE = "! B3LYP def2-TZVP TightSCF Opt"

# charge, multiplicity: from index.json, asserted below rather than trusted here.
CHARGE_MULT = {
    "acetylene": (0, 1), "ammonia": (0, 1), "benzene": (0, 1), "ch3": (0, 2),
    "ethane": (0, 1), "ethene": (0, 1), "formaldehyde": (0, 1), "formate": (-1, 1),
    "hcn": (0, 1), "lif": (0, 1), "n2": (0, 1), "ni_co_4": (0, 1),
    "nitromethane": (0, 1), "no": (0, 2), "o2": (0, 3), "ozone": (0, 1),
    "pf5": (0, 1), "pyridine": (0, 1), "sf6": (0, 1), "so2": (0, 1),
    "ticl4": (0, 1), "water": (0, 1),
}

D2R = math.pi / 180.0


def ring(n, r, start_deg=90.0):
    """n points on a circle in the xy plane, counter-clockwise from start_deg."""
    return [(r * math.cos((start_deg + 360.0 * i / n) * D2R),
             r * math.sin((start_deg + 360.0 * i / n) * D2R), 0.0) for i in range(n)]


def td(r):
    """Four vertices of a tetrahedron at distance r from the origin."""
    c = r / math.sqrt(3.0)
    return [(c, c, c), (c, -c, -c), (-c, c, -c), (-c, -c, c)]


def oh(r):
    return [(r, 0, 0), (-r, 0, 0), (0, r, 0), (0, -r, 0), (0, 0, r), (0, 0, -r)]


def bent(r1, r2, ang, e_c, e_a, e_b):
    """A bent triatomic e_a-e_c-e_b, apex at the origin, bisector along +y."""
    h = ang / 2.0 * D2R
    return [(e_c, (0.0, 0.0, 0.0)),
            (e_a, (r1 * math.sin(h), r1 * math.cos(h), 0.0)),
            (e_b, (-r2 * math.sin(h), r2 * math.cos(h), 0.0))]


def pyramidal(e_c, e_x, r, ang, n=3):
    """XHn with n equivalent bonds at bond angle `ang` between them (NH3-like)."""
    # cos(ang) = cos^2(t) * cos(120) + sin^2(t) ... solve for the polar angle t from the C3 axis
    ca = math.cos(ang * D2R)
    # for n=3: cos(ang) = 1 - 1.5 * sin^2(t)  =>  sin^2 t = (1 - cos ang) * 2/3
    st = math.sqrt((1.0 - ca) * 2.0 / 3.0)
    ct = math.sqrt(max(0.0, 1.0 - st * st))
    at = [(e_c, (0.0, 0.0, 0.0))]
    for i in range(n):
        p = 360.0 * i / n * D2R
        at.append((e_x, (r * st * math.cos(p), r * st * math.sin(p), -r * ct)))
    return at


def ch3_from_archive():
    """The optimised ch3 geometry of the original run, out of ch3_orca_reference.47."""
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, os.pardir, "nbo_reference", "ch3_orca_reference.47")
    z2e = {1: "H", 6: "C"}
    atoms, inside = [], False
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("$COORD"):
                inside = True
                next(f)          # the job-title line
                continue
            if inside:
                if s.startswith("$END"):
                    break
                p = s.split()
                atoms.append((z2e[int(p[0])], (float(p[2]), float(p[3]), float(p[4]))))
    if len(atoms) != 4:
        raise SystemExit("ch3_orca_reference.47: expected 4 atoms, got %d" % len(atoms))
    return atoms


def methyl_on(anchor, r_ch=1.09, ang=109.5):
    """Three H on the atom at `anchor`, pointing down -z, at bond angle `ang` to +z."""
    st = math.sin((180.0 - ang) * D2R)
    ct = math.cos((180.0 - ang) * D2R)
    out = []
    for i in range(3):
        p = 360.0 * i / 3 * D2R
        out.append(("H", (anchor[0] + r_ch * st * math.cos(p),
                          anchor[1] + r_ch * st * math.sin(p),
                          anchor[2] - r_ch * ct)))
    return out


def geometry(name):
    """Heavy atoms first, then hydrogens - the order the surviving ch3 archive uses."""
    if name == "ch3":
        return ch3_from_archive()

    if name == "water":
        return bent(0.958, 0.958, 104.5, "O", "H", "H")
    if name == "ozone":
        return bent(1.272, 1.272, 116.8, "O", "O", "O")
    if name == "so2":
        return bent(1.432, 1.432, 119.0, "S", "O", "O")
    if name == "ammonia":
        return pyramidal("N", "H", 1.012, 107.2)
    if name == "n2":
        return [("N", (0, 0, 0.550)), ("N", (0, 0, -0.550))]
    if name == "o2":
        return [("O", (0, 0, 0.604)), ("O", (0, 0, -0.604))]
    if name == "no":
        return [("N", (0, 0, 0.577)), ("O", (0, 0, -0.577))]
    if name == "lif":
        return [("Li", (0, 0, 0.0)), ("F", (0, 0, 1.580))]
    if name == "hcn":
        return [("C", (0, 0, 0.0)), ("N", (0, 0, 1.160)), ("H", (0, 0, -1.070))]
    if name == "acetylene":
        return [("C", (0, 0, 0.600)), ("C", (0, 0, -0.600)),
                ("H", (0, 0, 1.663)), ("H", (0, 0, -1.663))]
    if name == "ethene":
        at = [("C", (0, 0, 0.665)), ("C", (0, 0, -0.665))]
        s, c = math.sin(60.0 * D2R), math.cos(60.0 * D2R)
        for sz in (1.0, -1.0):
            for sx in (1.0, -1.0):
                at.append(("H", (sx * 1.085 * s, 0.0, sz * (0.665 + 1.085 * c))))
        return at
    if name == "ethane":
        at = [("C", (0, 0, 0.765)), ("C", (0, 0, -0.765))]
        st = math.sin((180.0 - 111.0) * D2R)
        ct = math.cos((180.0 - 111.0) * D2R)
        for i in range(3):          # staggered: the second set offset by 60 degrees
            p = (360.0 * i / 3) * D2R
            at.append(("H", (1.09 * st * math.cos(p), 1.09 * st * math.sin(p),
                             0.765 - 1.09 * ct)))
        for i in range(3):
            p = (60.0 + 360.0 * i / 3) * D2R
            at.append(("H", (1.09 * st * math.cos(p), 1.09 * st * math.sin(p),
                             -0.765 + 1.09 * ct)))
        return at
    if name == "formaldehyde":
        at = [("C", (0, 0, 0.0)), ("O", (0, 1.208, 0.0))]
        s, c = math.sin(61.0 * D2R), math.cos(61.0 * D2R)
        at += [("H", (1.111 * s, -1.111 * c, 0.0)), ("H", (-1.111 * s, -1.111 * c, 0.0))]
        return at
    if name == "formate":
        h = 62.0 * D2R
        at = [("C", (0, 0, 0.0)),
              ("O", (1.255 * math.sin(h), 1.255 * math.cos(h), 0.0)),
              ("O", (-1.255 * math.sin(h), 1.255 * math.cos(h), 0.0)),
              ("H", (0.0, -1.120, 0.0))]
        return at
    if name == "benzene":
        at = [("C", p) for p in ring(6, 1.394)]
        at += [("H", p) for p in ring(6, 2.478)]
        return at
    if name == "pyridine":
        # a regular hexagon is a fine start; the optimiser breaks the symmetry itself
        heavy = ring(6, 1.394)
        at = [("N", heavy[0])] + [("C", p) for p in heavy[1:]]
        outer = ring(6, 2.478)
        at += [("H", p) for p in outer[1:]]
        return at
    if name == "nitromethane":
        # C at the origin, N along +z, the two O on N with ONO 125 deg, the methyl H on -z
        h = 125.0 / 2.0 * D2R
        zn = 1.490
        at = [("C", (0.0, 0.0, 0.0)), ("N", (0.0, 0.0, zn)),
              ("O", (1.223 * math.sin(h), 0.0, zn + 1.223 * math.cos(h))),
              ("O", (-1.223 * math.sin(h), 0.0, zn + 1.223 * math.cos(h)))]
        at += methyl_on((0.0, 0.0, 0.0))
        return at
    if name == "pf5":
        at = [("P", (0, 0, 0.0)), ("F", (0, 0, 1.577)), ("F", (0, 0, -1.577))]
        at += [("F", p) for p in ring(3, 1.534)]
        return at
    if name == "sf6":
        return [("S", (0, 0, 0.0))] + [("F", p) for p in oh(1.564)]
    if name == "ticl4":
        return [("Ti", (0, 0, 0.0))] + [("Cl", p) for p in td(2.181)]
    if name == "ni_co_4":
        at = [("Ni", (0, 0, 0.0))]
        for p in td(1.838):
            at.append(("C", p))
        for p in td(1.838 + 1.141):
            at.append(("O", p))
        # heavy-atom order Ni, 4 C, 4 O - the C/O pairs share a direction, so the
        # tetrahedron helper is called twice at two radii rather than interleaved
        return at
    raise SystemExit("no start geometry for " + name)


def write_inp(out_dir, name, nprocs):
    q, mult = CHARGE_MULT[name]
    at = geometry(name)
    d = os.path.join(out_dir, name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name + ".inp"), "w") as f:
        f.write(KEYWORD_LINE + "\n")
        if nprocs > 1:
            f.write("%%pal nprocs %d end\n" % nprocs)
        f.write("%maxcore 3000\n")
        f.write("* xyz %d %d\n" % (q, mult))
        for e, (x, y, z) in at:
            f.write("  %-3s %14.8f %14.8f %14.8f\n" % (e, x, y, z))
        f.write("*\n")
    with open(os.path.join(d, name + "_start.xyz"), "w") as f:
        f.write("%d\n%s start geometry, charge %d multiplicity %d\n" % (len(at), name, q, mult))
        for e, (x, y, z) in at:
            f.write("%-3s %14.8f %14.8f %14.8f\n" % (e, x, y, z))
    return len(at)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--nprocs", type=int, default=8)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    idx = json.load(open(os.path.join(here, os.pardir, "nbo_reference", "index.json")))
    rec = {m["name"]: m for m in idx["molecules"]}
    if sorted(rec) != sorted(CHARGE_MULT):
        raise SystemExit("molecule list drifted from index.json")
    for name in sorted(rec):
        m = rec[name]
        assert (m["orca"]["charge"], m["orca"]["multiplicity"]) == CHARGE_MULT[name], name
        assert m["orca"]["keyword_line"] == KEYWORD_LINE, name
        n = write_inp(a.out_dir, name, a.nprocs)
        flag = "" if n == m["atoms"] else "  <-- ATOM COUNT MISMATCH, index says %d" % m["atoms"]
        print("%-13s %2d atoms%s" % (name, n, flag))


if __name__ == "__main__":
    main()
