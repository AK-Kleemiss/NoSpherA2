"""Is the Rydberg excess a physical misplacement, or just a constant per Rydberg function?

The coordinator's caution: 0.014-0.016 e missing from every ethane C 2s and 2p, reappearing across
78 Rydberg functions, is also what a per-orbital normalisation or weighting error would look like.
Those two possibilities make opposite predictions and this measures both:

  per-function constant  ->  d(Ryd on an atom) proportional to the NUMBER of Rydberg functions on
                             it, the excess spread thinly and evenly, and d/n_ryd roughly equal
                             across elements and molecules.
  misdrawn partition     ->  the excess concentrated in the FEW Rydberg functions that lie in the
                             same (atom, l) block as an occupied valence shell, and proportional
                             to that atom's valence population, not to its function count.

    python ryd_per_function.py <root>
"""
import json
import os
import socket
import sys

CLASSES = ("Cor", "Val", "Ryd")


def rows(path):
    return json.load(open(path))["nao"]


def demo():
    """A constant-per-function error must give a flat d/n, a concentrated one must not."""
    flat = [0.001] * 20
    conc = [0.02] + [0.0] * 19
    assert abs(max(flat) / sum(flat) - 0.05) < 1e-12
    assert max(conc) / sum(conc) == 1.0
    print("demo OK on %s" % socket.gethostname())


if __name__ == "__main__":
    if sys.argv[1:2] == ["--demo"]:
        demo()
        raise SystemExit(0)
    root = sys.argv[1]
    print("ryd_per_function on %s, 1 thread, root %s" % (socket.gethostname(), root))
    print("%-12s %-5s %6s %10s %10s %10s %7s %7s" %
          ("molecule", "atom", "n_ryd", "d(Ryd)", "d per fn", "val pop", "top1 %", "top3 %"))
    per_fn, by_el = [], {}
    for mol in sorted(os.listdir(root)):
        gp = os.path.join(root, mol, "%s.gennbo.nbo.json" % mol)
        np_ = os.path.join(root, mol, "%s.native.nbo.json" % mol)
        if not (os.path.isfile(gp) and os.path.isfile(np_)):
            continue
        g, n = rows(gp), rows(np_)
        if len(g) != len(n):
            continue
        # per atom: the Rydberg deviations, and how concentrated they are
        atoms = {}
        for a, b in zip(g, n):
            if (a["element"], a["atom"], a["lang"]) != (b["element"], b["atom"], b["lang"]):
                atoms = None
                break
            k = (a["element"], a["atom"])
            d = atoms.setdefault(k, {"ryd": [], "val": 0.0})
            if a["type"] == "Ryd":
                d["ryd"].append(b["occupancy"] - a["occupancy"])
            elif a["type"] == "Val":
                d["val"] += a["occupancy"]
        if not atoms:
            continue
        for k in sorted(atoms, key=lambda x: x[1]):
            d = atoms[k]
            if not d["ryd"]:
                continue
            tot = sum(d["ryd"])
            mag = sorted((abs(x) for x in d["ryd"]), reverse=True)
            s = sum(mag) or 1.0
            per_fn.append(tot / len(d["ryd"]))
            by_el.setdefault(k[0], []).append((tot, len(d["ryd"]), d["val"]))
            if abs(tot) > 0.004:
                print("%-12s %-5s %6d %10.5f %10.6f %10.5f %7.1f %7.1f" %
                      (mol, "%s%d" % k, len(d["ryd"]), tot, tot / len(d["ryd"]), d["val"],
                       100.0 * mag[0] / s, 100.0 * sum(mag[:3]) / s))
    print("\nper-element: is d(Ryd) set by the function count or by the valence population?")
    print("%-4s %5s %10s %10s %10s %10s" % ("el", "n", "mean d", "mean n_fn", "d per fn", "d/val"))
    for el in sorted(by_el):
        v = by_el[el]
        md = sum(x[0] for x in v) / len(v)
        mn = sum(x[1] for x in v) / len(v)
        mv = sum(x[2] for x in v) / len(v)
        print("%-4s %5d %10.5f %10.1f %10.6f %10.5f" %
              (el, len(v), md, mn, md / mn, md / mv if mv else float("nan")))
    lo, hi = min(per_fn), max(per_fn)
    print("\nd per Rydberg function over %d atoms: %.6f to %.6f e, a factor %.0f" %
          (len(per_fn), lo, hi, abs(hi / lo) if lo else float("inf")))
    print("A constant-per-function error would make that range narrow and the concentration")
    print("columns near 100/n; a misdrawn valence boundary puts most of it in the top few.")
