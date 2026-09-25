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

THE SUPPLEMENT: si2h6, thiophene, c2h5, added 25 Sep AFTER the baseline column of the first
thirteen was measured and BEFORE any number from the changed binary existed for any member.
The reason is a measured one and it weakens the first thirteen rather than the change: on the
baseline binary their total Rydberg-population error spans 0.00002 e (sih4) to 0.00607 e (cf3),
mean intra-atomic leak 0.00063 e per atom, while the 22 already in data/ span 0.00119 e (lif)
to 0.31605 e (benzene) at 0.01229 e per atom.  The defect this change acts on scales with
molecular SIZE - benzene, pyridine, ethane, nitromethane are the large ones - and not with row,
charge, diffuseness or spin, which is what the first thirteen were built to vary.  A set whose
baseline defect is 20x smaller than the tuned set's cannot test a claim of a 5x reduction: 63x
of 0.006 e is not a measurement of anything.  The three added members put the instrument back
inside its range while staying new material, and each is one homologue step from a member of
the eight or of the 22, so what is being tested stays a single transfer:

  si2h6      ethane is one of the EIGHT (baseline 0.12985 e) - the same skeleton one row down.
             If the change helps ethane and not disilane, it is second-row-specific.
  thiophene  benzene is the WORST tuned member (baseline 0.31605 e) - an aromatic ring with a
             third-row heteroatom in it, so the aromatic case and the third row at once.
  c2h5       ch3 is in data/ (baseline 0.03421 e) - a hydrocarbon open shell large enough for
             the instrument to have range, where cf3 and hs sit near its floor.
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
    # the supplement: size, because the defect scales with it (see the docstring)
    "si2h6": (0, 1), "thiophene": (0, 1), "c2h5": (0, 2),
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


# the supplement's internals.  si2h6: experimental Si-Si 2.331, Si-H 1.483, H-Si-Si 110.4.
# c2h5: ethane's skeleton with one H removed, C-C shortened to 1.490 as in the radical.
# thiophene: S-C 1.714, C2=C3 1.370, C3-C4 1.423, C-S-C 92.2, C-H 1.080 - the ring closes by
# construction from those, which is what the demo checks rather than trusting typed coordinates.
SI2H6 = (2.331, 1.483, 110.4)
C2H5 = (1.490, 1.090, 111.0)
TH_SC, TH_C2C3, TH_C3C4, TH_CSC, TH_CH = 1.714, 1.370, 1.423, 92.2, 1.080


def diatomic(a, b, r):
    return [(a, (0.0, 0.0, 0.0)), (b, (0.0, 0.0, r))]


def staggered_x2h6(el, r_xx, r_xh, ang, drop=0):
    """X2H6 along z, the second XH3 rotated 60 degrees.  `drop` removes trailing hydrogens,
    which is all a C2H5 radical is before the optimiser relaxes it."""
    h = r_xx / 2.0
    at = [(el, (0.0, 0.0, h)), (el, (0.0, 0.0, -h))]
    st = math.sin((180.0 - ang) * D2R)
    ct = math.cos((180.0 - ang) * D2R)
    for k, off in ((1.0, 0.0), (-1.0, 60.0)):
        for i in range(3):
            p = (off + 120.0 * i) * D2R
            # + ct, so H points AWAY from the other heavy atom: the angle asserted in demo() is
            # H-X-X, and make_inputs' ethane happens to write its hydrogens tucked the other way
            at.append(("H", (r_xh * st * math.cos(p), r_xh * st * math.sin(p),
                             k * (h + r_xh * ct))))
    return at[:len(at) - drop] if drop else at


def ext_h(a, n1, n2, r):
    """One H on ring atom `a`, on the external bisector of its two ring neighbours."""
    def unit(p, q):
        d = dist(p, q)
        return [(p[i] - q[i]) / d for i in range(3)]
    u, v = unit(a, n1), unit(a, n2)
    w = [u[i] + v[i] for i in range(3)]
    n = math.sqrt(sum(c * c for c in w))
    return ("H", tuple(a[i] + r * w[i] / n for i in range(3)))


def thiophene():
    """Planar C4H4S built from its internals: S at the origin, the ring in the xy plane."""
    h = TH_CSC / 2.0 * D2R
    s = (0.0, 0.0, 0.0)
    c2 = (TH_SC * math.sin(h), TH_SC * math.cos(h), 0.0)
    x3 = TH_C3C4 / 2.0
    dy = math.sqrt(TH_C2C3 ** 2 - (c2[0] - x3) ** 2)
    c3 = (x3, c2[1] + dy, 0.0)      # ACROSS the ring from S, not back towards it
    c5 = (-c2[0], c2[1], 0.0)
    c4 = (-c3[0], c3[1], 0.0)
    at = [("S", s), ("C", c2), ("C", c3), ("C", c4), ("C", c5)]
    at.append(ext_h(c2, s, c3, TH_CH))
    at.append(ext_h(c3, c2, c4, TH_CH))
    at.append(ext_h(c4, c3, c5, TH_CH))
    at.append(ext_h(c5, c4, s, TH_CH))
    return at


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
    if name == "si2h6":
        return staggered_x2h6("Si", *SI2H6)
    if name == "c2h5":
        return staggered_x2h6("C", *C2H5, drop=1)
    if name == "thiophene":
        return thiophene()
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

    # the supplement: the bonds and shapes it was asked for, and the right atom counts
    at = geometry("si2h6")
    assert [e for e, _ in at] == ["Si", "Si"] + ["H"] * 6, at
    assert abs(dist(at[0][1], at[1][1]) - SI2H6[0]) < 1e-9, dist(at[0][1], at[1][1])
    for q in [p for _, p in at[2:5]]:
        assert abs(dist(at[0][1], q) - SI2H6[1]) < 1e-9, dist(at[0][1], q)
        assert abs(angle(at[1][1], at[0][1], q) - SI2H6[2]) < 1e-6, angle(at[1][1], at[0][1], q)
    at = geometry("c2h5")
    assert [e for e, _ in at] == ["C", "C"] + ["H"] * 5, at      # one H short of ethane
    assert abs(dist(at[0][1], at[1][1]) - C2H5[0]) < 1e-9
    at = geometry("thiophene")
    assert [e for e, _ in at] == ["S"] + ["C"] * 4 + ["H"] * 4, at
    assert all(abs(q[2]) < 1e-12 for _, q in at), "thiophene must be planar"
    s, c2, c3, c4, c5 = (q for _, q in at[:5])
    for c in (c2, c5):
        assert abs(dist(s, c) - TH_SC) < 1e-9, dist(s, c)
    assert abs(angle(c2, s, c5) - TH_CSC) < 1e-6, angle(c2, s, c5)
    assert abs(dist(c2, c3) - TH_C2C3) < 1e-9, dist(c2, c3)      # the ring must CLOSE
    assert abs(dist(c5, c4) - TH_C2C3) < 1e-9, dist(c5, c4)
    assert abs(dist(c3, c4) - TH_C3C4) < 1e-9, dist(c3, c4)
    for hq, ring_c in zip([q for _, q in at[5:]], (c2, c3, c4, c5)):
        assert abs(dist(hq, ring_c) - TH_CH) < 1e-9, dist(hq, ring_c)
        # external: every H must be FURTHER from the ring centroid than its carbon is
        cen = [sum(q[i] for _, q in at[:5]) / 5.0 for i in range(3)]
        assert dist(hq, cen) > dist(ring_c, cen), "an H points into the ring"

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
