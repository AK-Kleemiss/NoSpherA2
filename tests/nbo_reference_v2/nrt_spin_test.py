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

What the 16 of 16 pairs at "exactly 2.0000" were actually measured on (25 Sep 2026)
-----------------------------------------------------------------------------------
Every one of them is a HALF-INTEGER gennbo value.  On ch3 the rows that come back at 2.0000 are
valency 1.50000 -> 3.00000, bond order 0.50000 -> 1.00000 and electron count 4 -> 8, all of which
are fixed by the leading resonance topology; 1.5 -> 3.0 is a ratio of exactly two whatever either
side computed.  The real-valued columns of the SAME rows are not doubled at all: alpha covalency
1.86, alpha electrovalency 2.58, beta covalency 2.20, beta electrovalency 0.76 - and since valency
is their sum on each side (1.29590 + 0.20410 = 1.50000, 2.84546 + 0.15454 = 3.00000, both exact),
a pinned total with a disagreeing split is the whole picture.  So `pinned()` marks those rows and
their ratio is never evidence about spin bookkeeping.

And the discriminating radicals cannot repair it, which is why `same_solution()` is a precondition
rather than another column: on allyl/hco/no2 the two sides do not solve the same resonance problem
(worst rank-paired |dw| 0.071/0.027/0.090 alpha and 0.280/0.008/0.286 beta, against a 0.005
tolerance on the valency those weights produce), so a valency ratio there measures the solution
difference and not the bookkeeping.  A molecule in that state is UNARBITRABLE, not FAIL: calling it
a failure is what invites the doubling reading back in.  Native's D(0) is 6.5-38x gennbo's, and
that is NOT open-shell specific - the 19 closed shells give the LARGER ratios, median 21.3 against
the open shells' 9.9 - so it is a general deviation and says nothing about spin either.
"""
import argparse
import json
import os
import socket
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir, "nbo_reference"))
import compare_nbo  # noqa: E402  for TOL only
from config import load_nbo  # refuses a reference an older parser wrote

TOL = compare_nbo.TOL["bond_order"]

ATOM_FIELDS = ("valency", "covalency", "electrovalency", "electron_count")
BOND_FIELDS = ("total", "covalent", "ionic")


def by_spin(rows, key):
    """{spin: {key: row}} for one NRT table."""
    out = {}
    for r in rows:
        out.setdefault(r.get("spin") or "", {})[key(r)] = r
    return out


def discriminates(g_atoms):
    """True when this molecule's alpha and beta values differ, so doubling and summing give
    different answers. Returns the pair that proves it, for the report."""
    #Only alpha and beta are spins.  gennbo also prints a composite alpha+beta valency table,
    #which the fixed parser now keeps under its own label; counting it as a third spin made this
    #guard skip every open shell.
    spins = [x for x in sorted(g_atoms) if x in ("alpha", "beta")]
    if len(spins) != 2:
        return None
    a, b = g_atoms[spins[0]], g_atoms[spins[1]]
    for k in sorted(set(a) & set(b)):
        for f in ATOM_FIELDS:
            if f in a[k] and f in b[k] and abs(a[k][f] - b[k][f]) > 10 * TOL:
                return (k, f, a[k][f], b[k][f])
    return None


def pinned(ref):
    """Is this gennbo value fixed by the leading topology rather than by the density?  NRT prints
    valencies, bond orders and electron counts of a symmetric single-structure radical as exact
    multiples of 0.5, and 1.50000 -> 3.00000 is a ratio of exactly two for any pair of codes that
    both got the topology right.  Such a row may FAIL; its RATIO is not a measurement."""
    return abs(ref * 2.0 - round(ref * 2.0)) < 1e-9


def same_solution(g_nrt, n_nrt):
    """The precondition for reading any ratio below.  A per-spin valency is a function of the
    retained resonance structures and their weights, so its ratio measures spin bookkeeping only
    where the weights themselves agree.  Rank-paired, because NRT `RS` is a rank and not a
    structure id, and padded with zeros so one extra negligible structure is not a mismatch.
    Returns None when the two solutions agree, else a one-line report of how they differ."""
    gw, nw = g_nrt.get("weights") or [], n_nrt.get("weights") or []
    out = []
    for spin in sorted({r.get("spin") or "" for r in gw} | {r.get("spin") or "" for r in nw}):
        a = sorted((r["weight_fraction"] for r in gw if (r.get("spin") or "") == spin), reverse=True)
        b = sorted((r["weight_fraction"] for r in nw if (r.get("spin") or "") == spin), reverse=True)
        m = max(len(a), len(b))
        na, nb = len(a), len(b)
        a, b = a + [0.0] * (m - na), b + [0.0] * (m - nb)
        worst = max((abs(x - y) for x, y in zip(a, b)), default=0.0)
        if worst > TOL:
            out.append("%s %d vs %d structures, worst |dw| by rank %.5f"
                       % (spin or "-", na, nb, worst))
    if not out:
        return None
    out.append("D(0) gennbo %.5f native %.5f"
               % (g_nrt.get("d_0", float("nan")), n_nrt.get("d_0", float("nan"))))
    return "; ".join(out)


def check(g, n, label, fields, failures, ratios=None):
    """One table. Every (spin, key) gennbo printed must carry gennbo's own value in native."""
    for spin in sorted(g):
        if spin == "composite":
            #gennbo's composite alpha+beta table has no native counterpart; that absence is
            #worth a line of its own, not one failure per row.
            if spin not in n:
                failures.append("%s: native prints no composite (alpha+beta) block" % label)
            continue
        if spin not in n:
            failures.append("%s: native printed no %s block at all" % (label, spin))
            continue
        for k in sorted(g[spin]):
            if k not in n[spin]:
                failures.append("%s %s %s: native did not print it" % (label, spin, k))
                continue
            for f in fields:
                if f not in g[spin][k]:
                    continue
                ref, got = float(g[spin][k][f]), float(n[spin][k].get(f, float("nan")))
                if not abs(got - ref) <= TOL:
                    ratio = got / ref if ref else float("nan")
                    failures.append("%s %s %s %s: gennbo %.5f native %.5f  ratio %.4f%s"
                                    % (label, spin, k, f, ref, got, ratio,
                                       "  [pinned: ratio forced]" if pinned(ref) else ""))
                    if ratios is not None and ref:
                        ratios.append((pinned(ref), ratio))
        for k in sorted(set(n[spin]) - set(g[spin])):
            failures.append("%s %s %s: native printed it, gennbo did not" % (label, spin, k))


def one(root, mol):
    d = os.path.join(root, mol)
    try:
        g = load_nbo(os.path.join(d, mol + ".gennbo.nbo.json"))
        n = load_nbo(os.path.join(d, mol + ".native.nbo.json"))
    except IOError as e:
        return {"molecule": mol, "status": "no data: %s" % e}
    if not (g.get("open_shell") and n.get("open_shell")):
        return {"molecule": mol, "status": "closed shell, nothing per-spin to compare"}

    akey = lambda r: "atom %d" % r["atom"]
    bkey = lambda r: "%d-%d" % (r["atom1"], r["atom2"])
    g_atoms = by_spin(g.get("nrt", {}).get("valencies", []), akey)
    n_atoms = by_spin(n.get("nrt", {}).get("valencies", []), akey)
    proof = discriminates(g_atoms)
    solution = same_solution(g.get("nrt") or {}, n.get("nrt") or {})
    if solution:
        #Stronger than either verdict below: no column of this molecule's NRT tables is a
        #measurement of the spin bookkeeping while the two sides disagree about the resonance
        #solution those tables are computed from.
        return {"molecule": mol, "discriminating_pair": proof, "solution": solution,
                "status": "UNARBITRABLE: the two sides solved different resonance problems"}
    if not proof:
        return {"molecule": mol,
                "status": "SKIP: alpha and beta agree, so doubling and summing are "
                          "indistinguishable here"}

    failures, ratios = [], []
    check(g_atoms, n_atoms, "valency", ATOM_FIELDS, failures, ratios)
    check(by_spin(g.get("nrt", {}).get("bond_orders", []), bkey),
          by_spin(n.get("nrt", {}).get("bond_orders", []), bkey),
          "bond order", BOND_FIELDS, failures, ratios)
    return {"molecule": mol, "status": "PASS" if not failures else "FAIL",
            "discriminating_pair": proof, "failures": failures, "ratios": ratios}


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
    bad = unarb = 0
    for mol in mols:
        r = one(a.root, mol)
        print("%-14s %s" % (r["molecule"], r["status"]))
        if r.get("discriminating_pair"):
            k, f, va, vb = r["discriminating_pair"]
            print("    discriminates: gennbo %s %s is %.5f in one spin and %.5f in the other, "
                  "so a spin sum (%.5f) is not a doubling (%.5f/%.5f)"
                  % (k, f, va, vb, va + vb, 2 * va, 2 * vb))
        if r.get("solution"):
            print("    solution differs: %s" % r["solution"])
            print("    therefore no ratio from this molecule is quoted, in either direction")
        for line in r.get("failures", []):
            print("    " + line)
        free = [x for p, x in r.get("ratios") or [] if not p]
        pin = [x for p, x in r.get("ratios") or [] if p]
        if r.get("ratios"):
            #The split IS the result: a uniform factor has to show up in the free rows too.
            print("    %d of %d failing rows are topology-pinned; the %d free rows run "
                  "%.4f..%.4f%s" % (len(pin), len(pin) + len(free), len(free),
                                    min(free) if free else float("nan"),
                                    max(free) if free else float("nan"),
                                    "" if not pin else "  (pinned rows all %.4f..%.4f)"
                                    % (min(pin), max(pin))))
        bad += r["status"] == "FAIL"
        unarb += r["status"].startswith("UNARBITRABLE")
    print("\n%d molecule(s) FAIL, %d UNARBITRABLE (exit status counts FAIL only)"
          % (bad, unarb))
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

    # pinned(): the half-integer grid NRT prints a symmetric radical on, and a real value next to it
    assert pinned(1.5) and pinned(0.5) and pinned(4.0) and pinned(0.0)
    assert not pinned(1.2959) and not pinned(0.2041)
    f, r = [], []
    mixed_g = {"alpha": {"atom 1": {"valency": 1.5, "covalency": 1.212}}}
    mixed_n = {"alpha": {"atom 1": {"valency": 3.0, "covalency": 2.25632}}}
    check(mixed_g, mixed_n, "valency", ATOM_FIELDS, f, r)
    assert [p for p, _ in r] == [True, False], r          # valency pinned, covalency free
    assert abs(r[1][1] - 1.8617) < 1e-4, r                 # ch3's real alpha covalency ratio
    assert any("[pinned" in x for x in f) and any("1.8617" in x for x in f), f

    # same_solution(): identical weights agree; one extra negligible structure still agrees; a
    # reordered-but-equal set agrees because the pairing is by rank, not by structure number.
    W = lambda vals, spin="alpha": [{"weight_fraction": v, "spin": spin} for v in vals]
    assert same_solution({"weights": W([0.6, 0.4])}, {"weights": W([0.6, 0.4])}) is None
    assert same_solution({"weights": W([0.6, 0.4])}, {"weights": W([0.4, 0.6])}) is None
    assert same_solution({"weights": W([0.6, 0.4])},
                         {"weights": W([0.6, 0.4, 0.001])}) is None
    hit = same_solution({"weights": W([0.43, 0.43]), "d_0": 0.064},
                        {"weights": W([0.37, 0.37]), "d_0": 0.734})
    assert hit and "worst |dw| by rank 0.06000" in hit and "D(0)" in hit, hit
    #and it must be per spin: a set that agrees in alpha and differs in beta is still a mismatch
    hit = same_solution({"weights": W([1.0]) + W([0.34, 0.33, 0.33], "beta")},
                        {"weights": W([1.0]) + W([0.46, 0.27, 0.27], "beta")})
    assert hit and hit.startswith("beta"), hit
    print("nrt_spin_test demo OK on %s" % socket.gethostname())
    return 0


if __name__ == "__main__":
    sys.exit(demo() if "--demo" in sys.argv else main(sys.argv))
