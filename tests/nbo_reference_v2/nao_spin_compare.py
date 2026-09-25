"""Compare native's per-spin NAO occupancies against gennbo's, orbital by orbital.

The spin-summed NAO table cannot localise the open-shell error: ch3's carbon is +0.0839 e in alpha
and -0.1279 e in beta, so the sum shows 0.044 e and the difference shows 0.212 e. This walks
gennbo's nao_alpha / nao_beta tables against native's and prints, per spin, the worst orbital and
the per-atom totals - which is the first per-orbital view of that error either side has ever had.

    python nao_spin_compare.py <root> [molecule ...]
"""
import json
import os
import socket
import sys
from config import load_nbo  # refuses a reference an older parser wrote


def label(n):
    return "%s%d %s %s %s" % (n["element"], n["atom"], n["lang"], n["type"], n["shell"])


def compare(g, n, spin, name):
    """Pair by printed index, but refuse to pair rows whose labels disagree."""
    gt, nt = g.get("nao_%s" % spin) or [], n.get("nao_%s" % spin) or []
    if not gt or not nt:
        print("  %-5s missing: gennbo %d rows, native %d rows" % (spin, len(gt), len(nt)))
        return None
    if len(gt) != len(nt):
        print("  %-5s LENGTH MISMATCH gennbo %d vs native %d - not comparable" % (spin, len(gt), len(nt)))
        return None
    worst, mismatched, per_atom = (0.0, None), 0, {}
    for a, b in zip(gt, nt):
        if label(a) != label(b):
            mismatched += 1
            continue
        d = b["occupancy"] - a["occupancy"]
        if abs(d) > abs(worst[0]):
            worst = (d, label(a))
        k = (a["element"], a["atom"])
        per_atom[k] = per_atom.get(k, 0.0) + d
    print("  %-5s %d rows, worst %+0.5f e at %s%s" %
          (spin, len(gt), worst[0], worst[1],
           ", %d label mismatches SKIPPED" % mismatched if mismatched else ""))
    for k in sorted(per_atom, key=lambda x: x[1]):
        print("        atom %s %-2d  sum of per-orbital deviations %+0.5f e" % (k[0], k[1], per_atom[k]))
    return per_atom


if __name__ == "__main__":
    root = sys.argv[1]
    only = sys.argv[2:]
    print("nao_spin_compare on %s, 1 thread (pure JSON work), root %s" % (socket.gethostname(), root))
    for mol in sorted(os.listdir(root)):
        if only and mol not in only:
            continue
        d = os.path.join(root, mol)
        gp, np_ = os.path.join(d, "%s.gennbo.nbo.json" % mol), os.path.join(d, "%s.native.nbo.json" % mol)
        if not (os.path.isfile(gp) and os.path.isfile(np_)):
            continue
        g, n = load_nbo(gp), load_nbo(np_)
        if not (g.get("nao_alpha") or n.get("nao_alpha")):
            continue
        print("\n%s" % mol)
        a = compare(g, n, "alpha", mol)
        b = compare(g, n, "beta", mol)
        if a and b:
            print("  per-atom check against the spin-summed numbers:")
            print("        %-6s %12s %12s %12s %12s" % ("atom", "d(alpha)", "d(beta)", "sum=d(q)", "diff=d(spin)"))
            for k in sorted(a, key=lambda x: x[1]):
                print("        %s %-4d %12.5f %12.5f %12.5f %12.5f" %
                      (k[0], k[1], a[k], b[k], a[k] + b[k], a[k] - b[k]))
