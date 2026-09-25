"""Verify the committed tree reproduces the README table from its own files.

Reads each molecule's gennbo JSON, the renat5 native JSON and the new baseline native JSON out
of data_holdout/ and prints d(Ryd) for both arms with the molecule's own floor (N_atoms*0.5e-5).
If this does not match the README, the README is quoting numbers the tree cannot produce.
"""
import io
import json
import os

root = r"D:\git\nos_nbo_ref\tests\nbo_reference_v2\data_holdout"


def npa(p):
    d = json.loads(io.open(p, encoding="utf-8", errors="replace").read())
    return {(r.get("atom"), r.get("element")): r
            for r in (d.get("npa") or []) if isinstance(r, dict)}


print("%-11s %5s %9s %10s %10s %8s" % ("molecule", "atoms", "floor", "base", "renat5", "factor"))
bad = 0
for m in sorted(os.listdir(root)):
    d = os.path.join(root, m)
    g = os.path.join(d, m + ".gennbo.nbo.json")
    if not os.path.isfile(g):
        continue
    gen = npa(g)
    ren = npa(os.path.join(d, m + ".native.nbo.json"))
    bas = npa(os.path.join(d, m + ".baseline.native.nbo.json"))
    keys = set(gen) & set(ren) & set(bas)
    sg = sum(gen[k].get("rydberg") or 0.0 for k in keys)
    dr = sum(ren[k].get("rydberg") or 0.0 for k in keys) - sg
    db = sum(bas[k].get("rydberg") or 0.0 for k in keys) - sg
    floor = len(keys) * 0.5e-5
    if abs(db - dr) <= floor:
        bad += 1
        note = "  <- ARMS INDISTINGUISHABLE"
    else:
        note = ""
    print("%-11s %5d %9.5f %+10.5f %+10.5f %8.1f%s"
          % (m, len(keys), floor, db, dr, abs(db) / abs(dr) if dr else float("inf"), note))
print("\n%d molecules whose two arms differ by less than one floor (want 0)" % bad)
assert bad == 0, "the two committed arms are not distinguishable on every molecule"
print("OK: every molecule's baseline and renat5 files are distinguishable above its own floor")
