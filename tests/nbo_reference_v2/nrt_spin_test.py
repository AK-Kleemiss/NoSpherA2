"""The failing test for the open-shell NRT factor of two, written before the fix.

    python nrt_spin_test.py <root> [--only ch3,no,o2]

For every open-shell molecule under <root> this asserts that each NRT quantity native prints
for ONE spin equals gennbo's value for THAT SAME spin: valency, covalency, electrovalency and
electron count per atom, and total/covalent/ionic per bond. It fails on today's binary, which
is the point of committing it now - ratio_probe.py measured 16 of 16 matched ch3 pairs at
exactly 2.0000 with residual 0.00000, so the assertion below is expected to report a factor
2 and nothing else.

Why not simply divide by two and be done
----------------------------------------
Because a spin-summed quantity would also be "about twice" a per-spin one and is a different
bug. On ch3 gennbo's alpha electron count is 4 and its beta is 3; their sum is 7. Native
prints 8 and 6. A code that had summed the spins would print 7 in BOTH blocks and would still
satisfy any check on the total number of electrons, while a code that doubles each spin's own
value prints 8 and 6. The two are only distinguishable on a molecule whose alpha and beta
values differ, so `discriminates()` asserts up front that the molecule under test has that
property and skips it otherwise - a test that cannot tell the two candidate causes apart
would pass a wrong fix.

Tolerance is the engine's own NRT tolerance (compare_nbo.TOL["bond_order"], 5.0e-3), so this
test and the comparator cannot disagree about what "equal" means.
"""
import argparse
import json
import os
import socket
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir, "nbo_reference"))
import compare_nbo  # noqa: E402  for TOL only

TOL = compare_nbo.TOL["bond_order"]

ATOM_FIELDS = ("valency", "covalency", "electrovalency", "electron_count")
BOND_FIELDS = ("total", "covalent", "ionic")


def by_spin(rows, key):
    """({spin: {key: row}}, [duplicate labels]) for one NRT table, composite rows dropped.

    Both halves of the return value exist because of the same defect. NBO prints THREE valency
    tables for an open shell - alpha, beta and "(composite alpha+beta)" - and the composite one
    comes after the last "Beta spin orbitals" header, so a parser that inherits the spin from
    that header files it as a second beta table. `spin == "composite"` is therefore dropped here
    on purpose, and any (spin, key) that still arrives twice is REPORTED rather than overwritten:
    a plain dict kept whichever row came last, which on ch3 replaced the real beta C valency
    (1.5000, N=3) with the composite one (3.0000, N=7) and made a doubled native value look
    like agreement. A test that silently drops half its input cannot fail for the right reason.
    """
    out, dupes = {}, []
    for r in rows:
        spin = r.get("spin") or ""
        if spin == "composite":
            continue
        k = key(r)
        if k in out.setdefault(spin, {}):
            dupes.append("%s %s" % (spin, k))
            continue
        out[spin][k] = r
    return out, dupes


def discriminates(g_atoms):
    """True when this molecule's alpha and beta values differ, so doubling and summing give
    different answers. Returns the pair that proves it, for the report."""
    spins = sorted(g_atoms)
    if len(spins) != 2:
        return None
    a, b = g_atoms[spins[0]], g_atoms[spins[1]]
    for k in sorted(set(a) & set(b)):
        for f in ATOM_FIELDS:
            if f in a[k] and f in b[k] and abs(a[k][f] - b[k][f]) > 10 * TOL:
                return (k, f, a[k][f], b[k][f])
    return None


def check(g, n, label, fields, failures):
    """One table. Every (spin, key) gennbo printed must carry gennbo's own value in native."""
    for spin in sorted(g):
        if spin not in n:
            failures.append("%s: native printed no %s block at all" % (label, spin))
            continue
        for k in sorted(g[spin]):
            if k not in n[spin]:
                # Say how big the entry was. gennbo prints its NRT bond-order matrix in full,
                # zeros included, so an absent native counterpart to a 0.00000 row is a
                # printing convention and an absent counterpart to a real number is not - and
                # a bare "did not print it" cannot be told apart from the other.
                big = max([abs(float(g[spin][k][f])) for f in fields if f in g[spin][k]] or [0.0])
                failures.append("%s %s %s: native did not print it (gennbo's largest of %s "
                                "is %.5f)" % (label, spin, k, "/".join(fields), big))
                continue
            for f in fields:
                if f not in g[spin][k]:
                    continue
                ref, got = float(g[spin][k][f]), float(n[spin][k].get(f, float("nan")))
                if not abs(got - ref) <= TOL:
                    ratio = got / ref if ref else float("nan")
                    failures.append("%s %s %s %s: gennbo %.5f native %.5f  ratio %.4f"
                                    % (label, spin, k, f, ref, got, ratio))
        for k in sorted(set(n[spin]) - set(g[spin])):
            failures.append("%s %s %s: native printed it, gennbo did not" % (label, spin, k))


def one(root, mol):
    d = os.path.join(root, mol)
    try:
        with open(os.path.join(d, mol + ".gennbo.nbo.json")) as f:
            g = json.load(f)
        with open(os.path.join(d, mol + ".native.nbo.json")) as f:
            n = json.load(f)
    except IOError as e:
        return {"molecule": mol, "status": "no data: %s" % e}
    if not (g.get("open_shell") and n.get("open_shell")):
        return {"molecule": mol, "status": "closed shell, nothing per-spin to compare"}

    akey = lambda r: "atom %d" % r["atom"]
    bkey = lambda r: "%d-%d" % (r["atom1"], r["atom2"])
    failures = []
    tables = {}
    for side, r in (("gennbo", g), ("native", n)):
        for name, rows, key in (("valency", r.get("nrt", {}).get("valencies", []), akey),
                                ("bond order", r.get("nrt", {}).get("bond_orders", []), bkey)):
            tables[(side, name)], dupes = by_spin(rows, key)
            for d in dupes:
                failures.append("%s %s: %s printed two rows for the same (spin, key) - a "
                                "mis-tagged composite table looks exactly like this"
                                % (name, d, side))
    proof = discriminates(tables[("gennbo", "valency")])
    if not proof:
        return {"molecule": mol,
                "status": "SKIP: alpha and beta agree, so doubling and summing are "
                          "indistinguishable here"}

    for name, fields in (("valency", ATOM_FIELDS), ("bond order", BOND_FIELDS)):
        check(tables[("gennbo", name)], tables[("native", name)], name, fields, failures)
    return {"molecule": mol, "status": "PASS" if not failures else "FAIL",
            "discriminating_pair": proof, "failures": failures}


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv[1:])
    mols = ([m for m in a.only.split(",") if m] or
            sorted(d for d in os.listdir(a.root)
                   if os.path.isdir(os.path.join(a.root, d))))

    print("nrt_spin_test on %s, python %s, tol %g (compare_nbo.TOL['bond_order'])"
          % (socket.gethostname(), sys.version.split()[0], TOL))
    print("each spin's native NRT value against gennbo's value for the SAME spin\n")
    bad = 0
    for mol in mols:
        r = one(a.root, mol)
        print("%-14s %s" % (r["molecule"], r["status"]))
        if r.get("discriminating_pair"):
            k, f, va, vb = r["discriminating_pair"]
            print("    discriminates: gennbo %s %s is %.5f in one spin and %.5f in the other, "
                  "so a spin sum (%.5f) is not a doubling (%.5f/%.5f)"
                  % (k, f, va, vb, va + vb, 2 * va, 2 * vb))
        for line in r.get("failures", []):
            print("    " + line)
        if r["status"] == "FAIL":
            bad += 1
    print("\n%d molecule(s) FAIL" % bad)
    return 1 if bad else 0


def demo():
    """`python nrt_spin_test.py --demo` - the test's own discrimination, on synthetic tables."""
    g_atoms = {"alpha": {"atom 1": {"valency": 1.5, "electron_count": 4.0}},
               "beta": {"atom 1": {"valency": 1.5, "electron_count": 3.0}}}
    assert discriminates(g_atoms), "4 vs 3 must count as discriminating"
    flat = {"alpha": {"atom 1": {"valency": 1.5, "electron_count": 4.0}},
            "beta": {"atom 1": {"valency": 1.5, "electron_count": 4.0}}}
    assert not discriminates(flat), "equal spins cannot tell doubling from summing"

    # today's binary: each spin doubled
    doubled = {s: {"atom 1": {k: 2 * v for k, v in r["atom 1"].items()}}
               for s, r in g_atoms.items()}
    f = []
    check(g_atoms, doubled, "valency", ATOM_FIELDS, f)
    assert len(f) == 4, f
    assert all("ratio 2.0000" in x for x in f), f

    # the wrong fix: spins summed. Must also fail, and not at ratio 2.
    summed = {s: {"atom 1": {"valency": 3.0, "electron_count": 7.0}} for s in g_atoms}
    f = []
    check(g_atoms, summed, "valency", ATOM_FIELDS, f)
    assert len(f) == 4, f
    assert any("1.7500" in x for x in f), f        # 7/4 in the alpha block
    assert any("2.3333" in x for x in f), f        # 7/3 in the beta block

    # the fix
    f = []
    check(g_atoms, json.loads(json.dumps(g_atoms)), "valency", ATOM_FIELDS, f)
    assert not f, f

    # by_spin: the composite table goes, a genuine duplicate is reported, and neither one is
    # allowed to overwrite the row the comparison depends on. ch3's real beta C valency is
    # 1.5/N=3 and its composite is 3.0/N=7, which is the pair that aliased.
    rows = [{"spin": "alpha", "atom": 1, "valency": 1.5, "electron_count": 4.0},
            {"spin": "beta", "atom": 1, "valency": 1.5, "electron_count": 3.0},
            {"spin": "composite", "atom": 1, "valency": 3.0, "electron_count": 7.0}]
    t, dupes = by_spin(rows, lambda r: "atom %d" % r["atom"])
    assert sorted(t) == ["alpha", "beta"], t
    assert t["beta"]["atom 1"]["electron_count"] == 3.0, t["beta"]
    assert not dupes, dupes
    mistagged = [dict(r, spin="beta") if r["spin"] == "composite" else r for r in rows]
    t, dupes = by_spin(mistagged, lambda r: "atom %d" % r["atom"])
    assert dupes == ["beta atom 1"], dupes
    assert t["beta"]["atom 1"]["electron_count"] == 3.0, "the first row must survive"
    print("nrt_spin_test demo OK on %s" % socket.gethostname())
    return 0


if __name__ == "__main__":
    sys.exit(demo() if "--demo" in sys.argv else main(sys.argv))
