"""ORCA inputs for the 12 HELD-OUT molecules: material the NAO change was never tuned on.

    python make_holdout.py <out_dir> [--nprocs N] [--demo]

WHY THIS SET EXISTS.  `renat5` was selected on eight molecules - lif, water, ammonia, ethane,
benzene, pf5, so2, sf6 - out of an enumerated space of 120 definitions in which nothing closed
a single molecule.  It was nominated by name before it was measured and the sweep failed to
displace it, but every number attached to it comes from those eight, and eight second-row
closed shells plus two third-row hypervalents is a narrow view of the cascade.  These twelve
are chosen so that the classes the eight never contained are represented, and so that at least
one member is a plausible way for the change to LOSE.  None of them is in the 25 directories of
data/ or data/ radicals; every one is small on purpose, because the point is coverage of the
class structure, not of molecular size.

THE ACCEPTANCE IS WEAKER THAN accept.py's, AND THAT IS STATED, NOT PAPERED OVER.  index.json is
the record of the lost original run; none of these twelve is in it, so the exact
basis-function identity test has nothing to compare against - exactly as for allyl, hco and
no2 in make_radicals.py.  accept_holdout.py records the criterion that is available (converged
optimisation, converged SCF, charge / multiplicity / basis dimension as asked, <S**2> near
S(S+1)) and RECORDS the basis-function count so that the next regeneration of this set does
have the exact test.

Level of theory, charge convention and file layout are make_inputs.py's, imported rather than
copied, so a change there cannot silently skip these twelve.

PER MEMBER, WHAT IT IS FOR (the pre-registration lives in the report; this is the short form):

  sih4       third row, Si: def2-TZVP's d set on Si is a genuinely populated diffuse shell.
  cl2        two equally polarisable third-row centres - the inter-atomic Val<->Ryd channel,
             which is where 39-66 % of gennbo's own transformation was found to sit.
  clf3       hypervalent beyond sf6/pf5: T-shaped, and the central atom keeps two lone pairs.
  bh3        electron-deficient: an empty valence p on boron, so the Val/Ryd boundary has no
             occupancy to guide it.  The one member with almost no Rydberg population.
  nh4_plus   a cation: the density is contracted, so the Rydberg set is squeezed.
  bf4_minus  an anion: the excess charge pushes population INTO the diffuse set.
  nacl       the third-row homologue of lif, which is one of the eight.  A row transfer on an
             otherwise identical bonding situation.
  mg_h2      group 2, which the eight never contained: Mg 3s valence with 3p/3d above it.
  hs         an open shell with the spin on a THIRD-row atom.
  cf3        an open shell whose spin sits on carbon with three fluorines pulling on it.
  ch2        a TRIPLET carbene - the only non-O2 triplet anywhere in this tree.
  zn_cl2     a transition metal beyond ni_co_4 / ticl4: d10, so 3d is full valence and 4d is
             Rydberg, and def2-TZVP gives Zn several d shells plus an f - a case where the
             component-major -> shell-major permutation is certainly not the identity.
  hi         an ECP-bearing heavy atom (def2-ECP, 28 core electrons on iodine).  Its first job
             is to find out whether gennbo and the .47 writer handle an ECP at all.

hi is counted as the thirteenth and is deliberately allowed to fail: whether it can be built
is itself the measurement, and if it cannot, accept_holdout.py records why instead of dropping
it.
"""
import argparse
import math
import os
import socket
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_inputs as mi

D2R = mi.D2R

CHARGE_MULT = {
    "sih4": (0, 1), "cl2": (0, 1), "clf3": (0, 1), "bh3": (0, 1),
    "nh4_plus": (1, 1), "bf4_minus": (-1, 1), "nacl": (0, 1), "mg_h2": (0, 1),
    "hs": (0, 2), "cf3": (0, 2), "ch2": (0, 3), "zn_cl2": (0, 1),
    "hi": (0, 1),
}

# Textbook / experimental starting geometries.  The optimiser moves them; a start only has to
# land in the right basin, which is what the converged-Opt check verifies.
R = {
    "sih4": 1.480, "cl2": 1.988, "nacl": 2.361, "mg_h2": 1.703, "hs": 1.341,
    "hi": 1.609, "zn_cl2": 2.072, "bh3": 1.190, "nh4_plus": 1.021,
    "bf4_minus": 1.390, "cf3": 1.320, "ch2": 1.075,
}
CLF3_AX = 1.598      # the unique short bond
CLF3_EQ = 1.698      # the two long ones
CLF3_ANG = 87.5      # F(long)-Cl-F(short)
CF3_FCF = 111.0      # pyramidal, not planar - CF3 is the classic non-planar carbon radical
CH2_HCH = 133.9      # triplet methylene, not the 102 deg of the singlet


def diatomic(a, b, r):
    return [(a, (0.0, 0.0, 0.0)), (b, (0.0, 0.0, r))]


def linear_xy2(x, y, r):
    return [(x, (0.0, 0.0, 0.0)), (y, (0.0, 0.0, r)), (y, (0.0, 0.0, -r))]


def geometry(name):
    """Heavy atoms first, then hydrogens - the order the surviving ch3 archive uses."""
    if name == "sih4":
        return [("Si", (0.0, 0.0, 0.0))] + [("H", p) for p in mi.td(R["sih4"])]
    if name == "nh4_plus":
        return [("N", (0.0, 0.0, 0.0))] + [("H", p) for p in mi.td(R["nh4_plus"])]
    if name == "bf4_minus":
        return [("B", (0.0, 0.0, 0.0))] + [("F", p) for p in mi.td(R["bf4_minus"])]
    if name == "bh3":
        return [("B", (0.0, 0.0, 0.0))] + [("H", p) for p in mi.ring(3, R["bh3"])]
    if name == "cf3":
        return mi.pyramidal("C", "F", R["cf3"], CF3_FCF)
    if name == "ch2":
        return mi.bent(R["ch2"], R["ch2"], CH2_HCH, "C", "H", "H")
    if name == "cl2":
        return diatomic("Cl", "Cl", R["cl2"])
    if name == "nacl":
        return diatomic("Na", "Cl", R["nacl"])
    if name == "hs":
        return diatomic("S", "H", R["hs"])
    if name == "hi":
        return diatomic("I", "H", R["hi"])
    if name == "mg_h2":
        return linear_xy2("Mg", "H", R["mg_h2"])
    if name == "zn_cl2":
        return linear_xy2("Zn", "Cl", R["zn_cl2"])
    if name == "clf3":
        h = CLF3_ANG * D2R
        return [("Cl", (0.0, 0.0, 0.0)),
                ("F", (0.0, CLF3_AX, 0.0)),
                ("F", (CLF3_EQ * math.sin(h), CLF3_EQ * math.cos(h), 0.0)),
                ("F", (-CLF3_EQ * math.sin(h), CLF3_EQ * math.cos(h), 0.0))]
    raise SystemExit("no start geometry for " + name)


def dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def angle(a, b, c):
    u = [a[i] - b[i] for i in range(3)]
    v = [c[i] - b[i] for i in range(3)]
    cu = sum(u[i] * v[i] for i in range(3)) / (dist(a, b) * dist(c, b))
    return math.degrees(math.acos(max(-1.0, min(1.0, cu))))


def demo():
    """Every written geometry must have the bonds, angles and shape it was asked for."""
    # every molecule: the element list is heavy-first, and nothing is closer than 0.95 A
    for name in sorted(CHARGE_MULT):
        at = geometry(name)
        p = [q for _, q in at]
        worst = min(dist(p[i], p[j]) for i in range(len(p)) for j in range(i + 1, len(p)))
        assert worst > 0.95, (name, worst)
        assert len(at) == len(set(tuple(round(x, 8) for x in q) for _, q in at)), name

    # the four tetrahedral / trigonal centres: all ligand distances equal, and the right angles
    for name, nlig, want in (("sih4", 4, 109.4712), ("nh4_plus", 4, 109.4712),
                             ("bf4_minus", 4, 109.4712), ("bh3", 3, 120.0)):
        at = geometry(name)
        assert len(at) == nlig + 1, (name, len(at))
        c, ligs = at[0][1], [q for _, q in at[1:]]
        for q in ligs:
            assert abs(dist(c, q) - R[name]) < 1e-9, (name, dist(c, q))
        assert abs(angle(ligs[0], c, ligs[1]) - want) < 1e-3, (name, angle(ligs[0], c, ligs[1]))
    assert all(abs(q[2]) < 1e-12 for _, q in geometry("bh3")), "bh3 must be planar"

    # cf3 pyramidal at the angle asked for, and NOT planar: the carbon must sit out of the F3 plane
    at = geometry("cf3")
    c, f = at[0][1], [q for _, q in at[1:]]
    assert len(f) == 3, at
    for q in f:
        assert abs(dist(c, q) - R["cf3"]) < 1e-9, dist(c, q)
    assert abs(angle(f[0], c, f[1]) - CF3_FCF) < 1e-6, angle(f[0], c, f[1])
    assert abs(sum(q[2] for q in f) / 3.0 - c[2]) > 0.2, "cf3 came out planar"

    # ch2 at the triplet angle, and the two diatomic-style triatomics genuinely linear
    at = geometry("ch2")
    assert abs(angle(at[1][1], at[0][1], at[2][1]) - CH2_HCH) < 1e-6
    for name in ("mg_h2", "zn_cl2"):
        at = geometry(name)
        assert abs(angle(at[1][1], at[0][1], at[2][1]) - 180.0) < 1e-6, name
        assert abs(dist(at[0][1], at[1][1]) - R[name]) < 1e-9, name

    # clf3 T-shaped: one short bond, two long ones, all in a plane, and the two long bonds
    # symmetric about the short one - a rotation error would break exactly this
    at = geometry("clf3")
    cl, ax, e1, e2 = (q for _, q in at)
    assert abs(dist(cl, ax) - CLF3_AX) < 1e-9, dist(cl, ax)
    assert abs(dist(cl, e1) - CLF3_EQ) < 1e-9 and abs(dist(cl, e2) - CLF3_EQ) < 1e-9
    assert abs(angle(ax, cl, e1) - CLF3_ANG) < 1e-6, angle(ax, cl, e1)
    assert abs(angle(ax, cl, e2) - CLF3_ANG) < 1e-6, angle(ax, cl, e2)
    assert abs(angle(e1, cl, e2) - 2.0 * CLF3_ANG) < 1e-6, angle(e1, cl, e2)
    assert all(abs(q[2]) < 1e-12 for _, q in at), "clf3 must be planar"

    # the diatomics, including the ECP case
    for name, els in (("cl2", ["Cl", "Cl"]), ("nacl", ["Na", "Cl"]),
                      ("hs", ["S", "H"]), ("hi", ["I", "H"])):
        at = geometry(name)
        assert [e for e, _ in at] == els, (name, at)
        assert abs(dist(at[0][1], at[1][1]) - R[name]) < 1e-9, name

    # none of these may collide with a molecule that already has a directory
    existing = {"acetylene", "allyl", "ammonia", "benzene", "ch3", "ethane", "ethene",
                "formaldehyde", "formate", "hcn", "hco", "lif", "n2", "ni_co_4",
                "nitromethane", "no", "no2", "o2", "ozone", "pf5", "pyridine", "sf6",
                "so2", "ticl4", "water"}
    clash = sorted(set(CHARGE_MULT) & existing)
    assert not clash, "held-out molecule already exists: %s" % clash
    print("demo OK on %s, 1 thread, %d molecules" % (socket.gethostname(), len(CHARGE_MULT)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", nargs="?")
    ap.add_argument("--nprocs", type=int, default=4)
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    if a.demo:
        demo()
        return
    if not a.out_dir:
        ap.error("out_dir is required unless --demo")
    demo()                                    # the check runs before it writes anything
    mi.CHARGE_MULT.update(CHARGE_MULT)        # write_inp reads this dict
    mi.geometry = geometry                    # ... and calls the module-global geometry()
    print("make_holdout on %s, 1 thread, out %s, keyword line: %s" %
          (socket.gethostname(), a.out_dir, mi.KEYWORD_LINE))
    for name in sorted(CHARGE_MULT):
        n = mi.write_inp(a.out_dir, name, a.nprocs)
        print("%-10s %d atoms, charge %2d multiplicity %d" % ((name, n) + CHARGE_MULT[name]))


if __name__ == "__main__":
    main()
