"""Is the committed column reproduced by the new renat5 arm, and what is d(Ryd)'s real floor?

    python3 storedchk.py <holdout_root> <twoarm_out_dir>

d(Ryd) is a SUM over atoms of a quantity gennbo prints to five decimals, so its floor is
N_atoms * 0.5e-5, not 1e-5.  A four-atom molecule cannot resolve 2e-5 and a nine-atom one
cannot resolve 4.5e-5; quoting a flat 1e-5 on a sum is a threshold below the printing floor.

The stored column in ../data_holdout was produced by the SAME binary with no env var, i.e. the
renat5 arm, so stored and new-renat5 must agree to that floor.  If they do not, the difference
is in the binary or the run, not in the arithmetic.
"""
import json
import os
import sys

root, out = sys.argv[1], sys.argv[2]
MOLS = sorted(m for m in os.listdir(root)
              if os.path.isfile(os.path.join(root, m, m + ".gennbo.nbo.json")))


def npa(path):
    return {(r.get("atom"), r.get("element")): r
            for r in (json.load(open(path)).get("npa") or []) if isinstance(r, dict)}


def dryd(gen, nat):
    keys = set(gen) & set(nat)
    return (sum((nat[k].get("rydberg") or 0.0) - (gen[k].get("rydberg") or 0.0) for k in keys),
            len(keys))


print("%-11s %5s %9s %10s %10s %10s %10s" % (
    "molecule", "atoms", "floor", "stored", "new_ren5", "new-stored", "in_floors"))
worst = 0.0
for m in MOLS:
    gen = npa(os.path.join(root, m, m + ".gennbo.nbo.json"))
    st, n = dryd(gen, npa(os.path.join(root, m, m + ".native.nbo.json")))
    nw, _ = dryd(gen, npa(os.path.join(out, "renat5", m, m + ".native.nbo.json")))
    floor = n * 0.5e-5
    print("%-11s %5d %9.5f %+10.5f %+10.5f %+10.5f %10.2f"
          % (m, n, floor, st, nw, nw - st, abs(nw - st) / floor))
    worst = max(worst, abs(nw - st) / floor)
print("\nworst stored-vs-new disagreement: %.2f floors" % worst)
print("<1 floor everywhere means the stored column and the new renat5 arm are the same"
      " measurement and the committed numbers were mislabelled, not miscomputed")
