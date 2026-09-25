"""Where does native's NAO population actually differ from NBO 7's - between atoms, or between
classes on the same atom?

The NPA charge comparison put the median disagreement at 0.0117 e, which looks like a small error.
Splitting the same NAO tables by class (core / valence / Rydberg) instead of by atom tells a
different story: on ethane the valence set is 0.12989 e short and the Rydberg set is 0.12986 e
long.  If that swap happens WITHIN each atom, the atomic totals stay nearly right and the charge
comparison cannot see a leak ten times its own size - which is exactly the kind of cancellation
that made the earlier "agreement" meaningless.

This measures both decompositions on the same rows:

    intra-atomic leak   L(atom) = d(Ryd on that atom), against  d(total on that atom)
    ratio               |d(q)| / L  - small ratio means the leak is invisible to NPA charges

    python nao_class_leak.py <root> [molecule ...]
"""
import json
import os
import socket
import sys

CLASSES = ("Cor", "Val", "Ryd")


def by_atom_class(table):
    """{(element, atom): {class: population}} from one printed NAO table."""
    out = {}
    for n in table:
        key = (n["element"], n["atom"])
        cls = n["type"] if n["type"] in CLASSES else "Ryd"
        out.setdefault(key, dict.fromkeys(CLASSES, 0.0))[cls] += n["occupancy"]
    return out


def paired(g, n):
    """Pair row by row, BY RANK inside each (atom, l) block - never by the printed shell label.

    Both sides print an (atom, l) block in descending occupancy, but they number the Rydberg
    shells differently: ethane's C1 s block is 1s 2s 3s 5s 4s in NBO 7 and 1s 2s 3s 4s 5s in
    native, same five orbitals.  Demanding equal shell labels threw away 11 of 22 molecules.
    Element, atom and l must still agree, and a class (Cor/Val/Ryd) disagreement is reported
    rather than silently summed into the wrong bucket.
    """
    if len(g) != len(n):
        return None
    reclassified = 0
    for a, b in zip(g, n):
        if (a["element"], a["atom"], a["lang"]) != (b["element"], b["atom"], b["lang"]):
            return None
        if a["type"] != b["type"]:
            reclassified += 1
    return by_atom_class(g), by_atom_class(n), reclassified


def demo():
    """A pure valence->Rydberg swap on one atom must show as a leak with zero charge change."""
    g = [{"element": "C", "atom": 1, "lang": "s", "shell": "2s", "type": "Val", "occupancy": 1.0},
         {"element": "C", "atom": 1, "lang": "dxy", "shell": "3d", "type": "Ryd", "occupancy": 0.0}]
    n = [dict(g[0], occupancy=0.9), dict(g[1], occupancy=0.1)]
    gc, nc, recl = paired(g, n)
    assert recl == 0, recl
    k = ("C", 1)
    leak = nc[k]["Ryd"] - gc[k]["Ryd"]
    dq = sum(nc[k][c] - gc[k][c] for c in CLASSES)
    assert abs(leak - 0.1) < 1e-12 and abs(dq) < 1e-12, (leak, dq)
    print("demo OK on %s" % socket.gethostname())


if __name__ == "__main__":
    if sys.argv[1:2] == ["--demo"]:
        demo()
        raise SystemExit(0)
    root = sys.argv[1]
    only = sys.argv[2:]
    print("nao_class_leak on %s, 1 thread (pure JSON work), root %s" % (socket.gethostname(), root))
    print("%-12s %9s %9s %9s %9s %9s %7s" %
          ("molecule", "d(Cor)", "d(Val)", "d(Ryd)", "worst L", "worst dq", "dq/L"))
    tot_leak, tot_dq, rows = 0.0, 0.0, 0
    for mol in sorted(os.listdir(root)):
        if only and mol not in only:
            continue
        gp = os.path.join(root, mol, "%s.gennbo.nbo.json" % mol)
        np_ = os.path.join(root, mol, "%s.native.nbo.json" % mol)
        if not (os.path.isfile(gp) and os.path.isfile(np_)):
            continue
        g, n = json.load(open(gp)), json.load(open(np_))
        pair = paired(g.get("nao") or [], n.get("nao") or [])
        if pair is None:
            print("%-12s NAO tables not comparable - skipped" % mol)
            continue
        gc, nc, recl = pair
        if recl:
            print("%-12s %d rows classified into a different class by the two sides" % (mol, recl))
        d_cls = {c: sum(nc[k][c] - gc[k][c] for k in gc) for c in CLASSES}
        leaks = {k: nc[k]["Ryd"] - gc[k]["Ryd"] for k in gc}
        dqs = {k: sum(nc[k][c] - gc[k][c] for c in CLASSES) for k in gc}
        wl = max(leaks, key=lambda k: abs(leaks[k]))
        wq = max(dqs, key=lambda k: abs(dqs[k]))
        print("%-12s %+9.5f %+9.5f %+9.5f %+9.5f %+9.5f %7.2f" %
              (mol, d_cls["Cor"], d_cls["Val"], d_cls["Ryd"], leaks[wl], dqs[wq],
               abs(dqs[wq]) / abs(leaks[wl]) if leaks[wl] else float("nan")))
        for k in gc:
            tot_leak += abs(leaks[k])
            tot_dq += abs(dqs[k])
            rows += 1
    print("\n%d atoms: mean |intra-atomic Rydberg leak| %.5f e, mean |charge deviation| %.5f e" %
          (rows, tot_leak / rows, tot_dq / rows))
    print("The charge comparison sees the second number.  The first is what the NAO construction")
    print("actually gets wrong, and a valence->Rydberg swap on one atom barely moves the first.")
