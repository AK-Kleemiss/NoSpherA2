"""ORCA inputs for three open-shell molecules the old reference set never contained.

    python make_radicals.py <out_dir> [--nprocs N]

WHY THESE THREE.  ch3, no and o2 are the only open shells in the 22, and gennbo's alpha and
beta NRT tables agree on all three, so they cannot test the claim that native doubles each
spin's NRT value - a factor of 2 on two equal numbers is invisible.  allyl, HCO and NO2 are
doublets whose spin density is genuinely lopsided: allyl carries it on the two terminal
carbons and not the middle one, HCO on the carbon, NO2 on the nitrogen.  If alpha and beta
differ there, the doubling either shows up or it does not.

WHAT IS DIFFERENT ABOUT THEM.  index.json is the record of the lost original run, and these
three were never in it, so accept.py's identity test - basis_functions equal to a recorded
value, energy within 2.0e-5 Eh of a recorded value - has nothing to compare against and is
not available.  That is a real weakening and it is stated here rather than papered over: for
these three the acceptance is only that the optimisation converged, the SCF converged, the
multiplicity and basis dimension are what this script asked for, and <S**2> is near 0.75.
The comparison they feed is still two-sided - ORCA + gennbo 7 on one side, -nbo_native on the
other, both reading the same .47 - which is the property that matters for arbitration.  What
they cannot do is prove the wavefunction reproduces somebody else's earlier wavefunction.

Level of theory, charge convention and file layout are make_inputs.py's, imported rather than
copied, so a change there cannot silently skip these three.
"""
import argparse
import math
import os
import socket
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_inputs as mi

D2R = mi.D2R

CHARGE_MULT = {"allyl": (0, 2), "hco": (0, 2), "no2": (0, 2)}

# Experimental / textbook starting geometries.  The optimiser moves them; these only have to
# land in the right basin, which is what the converged-Opt check verifies.
ALLYL_CC = 1.388     # C-C, delocalised, between ethane's 1.53 and ethene's 1.33
ALLYL_CCC = 124.0
ALLYL_CH = 1.09


def rot_z(v, deg):
    c, s = math.cos(deg * D2R), math.sin(deg * D2R)
    return (c * v[0] - s * v[1], s * v[0] + c * v[1], 0.0)


def sp2_hydrogens(center, neighbour, r, angles=(120.0, -120.0)):
    """Place H on a planar atom by rotating the direction to its one known neighbour."""
    d = [neighbour[i] - center[i] for i in range(3)]
    n = math.sqrt(sum(x * x for x in d))
    u = tuple(x / n for x in d)
    out = []
    for a in angles:
        w = rot_z(u, a)
        out.append(("H", tuple(center[i] + r * w[i] for i in range(3))))
    return out


def geometry(name):
    """Heavy atoms first, then hydrogens - the order the surviving ch3 archive uses."""
    if name == "no2":
        return mi.bent(1.193, 1.193, 134.1, "N", "O", "O")
    if name == "hco":
        # C apex, O and H arms: C-O 1.175, C-H 1.125, HCO 124.4 deg
        return mi.bent(1.175, 1.125, 124.4, "C", "O", "H")
    if name == "allyl":
        c = mi.bent(ALLYL_CC, ALLYL_CC, ALLYL_CCC, "C", "C", "C")   # C1 centre, C2 and C3 ends
        at = list(c)
        centre, end2, end3 = c[0][1], c[1][1], c[2][1]
        # one H on the central carbon, on the bisector away from the two ends (-y here)
        at.append(("H", (0.0, -ALLYL_CH, 0.0)))
        for end in (end2, end3):
            at += sp2_hydrogens(end, centre, ALLYL_CH)
        return at
    raise SystemExit("no start geometry for " + name)


def demo():
    """The written geometries must have the bond lengths and angles asked for."""
    def dist(a, b):
        return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))

    def angle(a, b, c):
        u = [a[i] - b[i] for i in range(3)]
        v = [c[i] - b[i] for i in range(3)]
        cu = sum(u[i] * v[i] for i in range(3)) / (dist(a, b) * dist(c, b))
        return math.degrees(math.acos(max(-1.0, min(1.0, cu))))

    at = geometry("allyl")
    assert [e for e, _ in at] == ["C", "C", "C", "H", "H", "H", "H", "H"], at
    p = [q for _, q in at]
    assert abs(dist(p[0], p[1]) - ALLYL_CC) < 1e-9, dist(p[0], p[1])
    assert abs(angle(p[1], p[0], p[2]) - ALLYL_CCC) < 1e-6, angle(p[1], p[0], p[2])
    for h in (p[4], p[5]):                       # the two H on the first terminal carbon
        assert abs(dist(p[1], h) - ALLYL_CH) < 1e-9, dist(p[1], h)
    # an sp2 terminal carbon: C1 and its two H spread 120 deg apart, all in the z=0 plane
    assert abs(angle(p[0], p[1], p[4]) - 120.0) < 1e-6, angle(p[0], p[1], p[4])
    assert abs(angle(p[4], p[1], p[5]) - 120.0) < 1e-6, angle(p[4], p[1], p[5])
    assert all(abs(q[2]) < 1e-12 for q in p), "allyl must be planar"
    # no closer contact than a real bond anywhere, which is what a bad rotation would give
    worst = min(dist(p[i], p[j]) for i in range(len(p)) for j in range(i + 1, len(p)))
    assert worst > 1.0, worst
    h = geometry("hco")
    assert [e for e, _ in h] == ["C", "O", "H"], h
    assert abs(angle(h[1][1], h[0][1], h[2][1]) - 124.4) < 1e-6
    n = geometry("no2")
    assert abs(angle(n[1][1], n[0][1], n[2][1]) - 134.1) < 1e-6
    print("demo OK on %s, 1 thread" % socket.gethostname())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--nprocs", type=int, default=4)
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    if a.demo:
        demo()
        return
    demo()                                    # the check runs before it writes anything
    mi.CHARGE_MULT.update(CHARGE_MULT)        # write_inp reads this dict
    mi.geometry = geometry                    # ... and calls the module-global geometry()
    print("make_radicals on %s, 1 thread, out %s, keyword line: %s" %
          (socket.gethostname(), a.out_dir, mi.KEYWORD_LINE))
    for name in sorted(CHARGE_MULT):
        n = mi.write_inp(a.out_dir, name, a.nprocs)
        print("%-8s %d atoms, charge %d multiplicity %d" % ((name, n) + CHARGE_MULT[name]))


if __name__ == "__main__":
    main()
