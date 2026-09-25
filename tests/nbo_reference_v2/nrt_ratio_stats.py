"""Is the open-shell NRT error a factor of two, or only a factor of two on ch3?

ratio_probe.py measured 16 of 16 matched ch3 pairs at exactly 2.0000 with residual 0.00000, and
"native doubles each spin's own value" was written down on the strength of it.  ch3 is the most
symmetric open shell in the set.  This prints the distribution of native/gennbo per spin-resolved
NRT quantity for every open shell under the roots given, so the claim can be checked on a molecule
where the NRT structure set is not trivial.

A real factor of two is a spike at 2.0000 with nothing else in it.  Anything with a spread means
halving the numbers would leave the residue behind and there is a second error underneath.

    python nrt_ratio_stats.py <root> [<root> ...]
"""
import json
import os
import socket
import sys
from config import load_nbo  # refuses a reference an older parser wrote

ATOM_FIELDS = ("valency", "covalency", "electrovalency", "electron_count")
BOND_FIELDS = ("total", "covalent", "ionic")
SPINS = ("alpha", "beta")


def pairs(mol_dir, mol):
    """Yield (label, gennbo, native) for every spin-resolved NRT number both sides printed."""
    g = load_nbo(os.path.join(mol_dir, "%s.gennbo.nbo.json" % mol)).get("nrt") or {}
    n = load_nbo(os.path.join(mol_dir, "%s.native.nbo.json" % mol)).get("nrt") or {}
    for table, fields, key in (("valencies", ATOM_FIELDS, lambda r: r["atom"]),
                               ("bond_orders", BOND_FIELDS,
                                lambda r: (r.get("atom1"), r.get("atom2")))):
        gi, ni = {}, {}
        for src, out in ((g.get(table) or [], gi), (n.get(table) or [], ni)):
            for r in src:
                if r.get("spin") in SPINS:
                    out[(r["spin"], key(r))] = r
        for k in sorted(set(gi) & set(ni), key=repr):
            for f in fields:
                if f in gi[k] and f in ni[k]:
                    yield ("%s %s %s %s" % (table, k[0], k[1], f), float(gi[k][f]), float(ni[k][f]))


def demo():
    """A clean doubling has zero spread; a scattered one does not."""
    assert spread([2.0, 2.0, 2.0]) == 0.0
    assert spread([1.1, 3.0]) > 1.0
    print("demo OK on %s" % socket.gethostname())


def spread(rs):
    return max(rs) - min(rs)


if __name__ == "__main__":
    if sys.argv[1:2] == ["--demo"]:
        demo()
        raise SystemExit(0)
    print("nrt_ratio_stats on %s, 1 thread" % socket.gethostname())
    print("%-12s %5s %8s %8s %8s %8s %7s" %
          ("molecule", "n", "min", "median", "max", "spread", "at 2.000"))
    for root in sys.argv[1:]:
        for mol in sorted(os.listdir(root)):
            d = os.path.join(root, mol)
            if not os.path.isfile(os.path.join(d, "%s.native.nbo.json" % mol)):
                continue
            if not os.path.isfile(os.path.join(d, "%s.gennbo.nbo.json" % mol)):
                continue
            rs = [nv / gv for _, gv, nv in pairs(d, mol) if abs(gv) > 1e-3]
            if not rs:
                continue
            rs.sort()
            exact = sum(1 for r in rs if abs(r - 2.0) < 5e-4)
            print("%-12s %5d %8.4f %8.4f %8.4f %8.4f %4d/%-4d" %
                  (mol, len(rs), rs[0], rs[len(rs) // 2], rs[-1], spread(rs), exact, len(rs)))
    print("\nA spike at 2.0000 with no spread is a factor of two.  A spread is a second error.")
