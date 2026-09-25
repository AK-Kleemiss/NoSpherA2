"""Where does the Rydberg excess sit: in the l blocks that contain a valence shell, or the others?

The baseline Rydberg population is 6.2x gennbo's on ethane and 3.7x on benzene, but 0.93-0.98x on
pf5, so2 and sf6 - the hypervalent three are essentially right. Both sides agree on which shell is
Rydberg (nao_class_leak.py reports 0 class reclassifications once rows are paired by rank), so the
question is which (atom, l) blocks the excess lands in.

  excess in an l block that also holds a VALENCE shell  ->  the valence/Rydberg orthogonalisation
      (step 3's Schmidt projection and OWSO, or step 4's re-diagonalisation) is leaking population
      from the valence shell into the Rydberg shell next to it.
  excess in a block with NO valence shell (C's d and f) ->  it cannot be a valence->Rydberg leak
      within a block at all, and the population must be arriving from another l, which only the
      Schmidt step can do.

    python ryd_excess_by_l.py <root>
"""
import json
import os
import socket
import sys
from collections import defaultdict

LANG = "spdfghi"


def by_block(rows):
    """{(atom, l): [rows in rank order]} - the l comes from the printed lang letter."""
    out = defaultdict(list)
    for r in rows:
        lang = (r.get("lang") or "").strip().lower()
        out[(r["atom"], LANG.index(lang[0]) if lang and lang[0] in LANG else -1)].append(r)
    return out


def demo():
    """A block holding a Val row must be reported as valence-bearing, one without it must not."""
    b = by_block([{"atom": 1, "lang": "s", "type": "Val", "occupancy": 1.0},
                  {"atom": 1, "lang": "s", "type": "Ryd", "occupancy": 0.01},
                  {"atom": 1, "lang": "d", "type": "Ryd", "occupancy": 0.02}])
    assert any(r["type"] == "Val" for r in b[(1, 0)])
    assert not any(r["type"] == "Val" for r in b[(1, 2)])
    print("demo OK on %s" % socket.gethostname())


if __name__ == "__main__":
    if sys.argv[1:2] == ["--demo"]:
        demo()
        raise SystemExit(0)
    root = sys.argv[1]
    print("ryd_excess_by_l on %s, 1 thread, root %s" % (socket.gethostname(), root))
    print("%-13s %14s %14s   %s" %
          ("molecule", "with valence", "no valence", "worst single block"))
    tw = tn = 0.0
    for mol in sorted(os.listdir(root)):
        d = os.path.join(root, mol)
        fg = os.path.join(d, "%s.gennbo.nbo.json" % mol)
        fn = os.path.join(d, "%s.native.nbo.json" % mol)
        if not (os.path.isfile(fg) and os.path.isfile(fn)):
            continue
        g, n = by_block(json.load(open(fg))["nao"]), by_block(json.load(open(fn))["nao"])
        withv = nov = 0.0
        worst = (0.0, "")
        for k in sorted(set(g) & set(n)):
            if len(g[k]) != len(n[k]):
                continue
            has_val = any(r["type"] == "Val" for r in g[k])
            dd = sum(b["occupancy"] - a["occupancy"]
                     for a, b in zip(g[k], n[k]) if a["type"] == "Ryd")
            if has_val:
                withv += dd
            else:
                nov += dd
            if abs(dd) > abs(worst[0]):
                worst = (dd, "atom %d %s%s" % (k[0], LANG[k[1]] if k[1] >= 0 else "?",
                                               " (has valence)" if has_val else ""))
        tw += withv
        tn += nov
        print("%-13s %+14.5f %+14.5f   %+.5f %s" % (mol, withv, nov, worst[0], worst[1]))
    print("\ntotals over the set: %+.5f e in valence-bearing l blocks, %+.5f e in blocks with no"
          % (tw, tn))
    print("valence shell.  The second number cannot be an intra-block leak.")
