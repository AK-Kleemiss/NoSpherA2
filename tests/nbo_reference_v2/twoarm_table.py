"""The held-out table, both arms, against gennbo as the arbiter.

    python3 twoarm_table.py <holdout_root> <twoarm_out_dir>

Per molecule and per arm: d(Ryd) = (sum of native's per-atom Rydberg populations) minus
(sum of gennbo's), and the largest per-atom |charge| difference.  gennbo prints five decimals,
so 1e-5 is the floor and nothing below it is read as a number.

Every field is read by its own key out of the parsed JSON; no column offsets, no inherited
context.  Atoms are matched on (atom index, element) and a molecule whose two sides do not
carry the same atom set is reported as such rather than silently summed over a partial set.
"""
import json
import os
import sys

root, out = sys.argv[1], sys.argv[2]
MOLS = sys.argv[3].split() if len(sys.argv) > 3 else sorted(
    m for m in os.listdir(root) if os.path.isfile(os.path.join(root, m, m + ".gennbo.nbo.json")))


def npa(path):
    rows = json.load(open(path)).get("npa") or []
    return {(r.get("atom"), r.get("element")): r
            for r in rows if isinstance(r, dict)}


def stats(gen, nat):
    keys = set(gen) & set(nat)
    miss = (set(gen) ^ set(nat))
    sg = sum(gen[k].get("rydberg") or 0.0 for k in keys)
    sn = sum(nat[k].get("rydberg") or 0.0 for k in keys)
    dq = max((abs((nat[k].get("charge") or 0.0) - (gen[k].get("charge") or 0.0)) for k in keys),
             default=float("nan"))
    return sn - sg, dq, len(keys), len(miss)


print("%-11s %10s %10s   %9s %9s   %9s %9s  %s" % (
    "molecule", "dRyd_base", "dRyd_ren5", "|dq|_base", "|dq|_ren5", "ratio", "atoms", "note"))
rows = []
for m in MOLS:
    g = os.path.join(root, m, m + ".gennbo.nbo.json")
    pb = os.path.join(out, "baseline", m, m + ".native.nbo.json")
    pr = os.path.join(out, "renat5", m, m + ".native.nbo.json")
    if not all(os.path.exists(p) for p in (g, pb, pr)):
        print("%-11s %s" % (m, "MISSING " + " ".join(
            os.path.basename(p) for p in (g, pb, pr) if not os.path.exists(p))))
        continue
    gen = npa(g)
    db, qb, nb, mb = stats(gen, npa(pb))
    dr, qr, nr, mr = stats(gen, npa(pr))
    ratio = (abs(db) / abs(dr)) if abs(dr) > 0 else float("inf")
    note = "" if (mb == 0 and mr == 0) else "atom-set mismatch b=%d r=%d" % (mb, mr)
    print("%-11s %+10.5f %+10.5f   %9.5f %9.5f   %9.2f %9d  %s"
          % (m, db, dr, qb, qr, ratio, nb, note))
    rows.append((m, db, dr, qb, qr, ratio))

if rows:
    wb = max(abs(r[1]) for r in rows)
    wr = max(abs(r[2]) for r in rows)
    rs = sorted(r[5] for r in rows)
    med = rs[len(rs) // 2] if len(rs) % 2 else 0.5 * (rs[len(rs) // 2 - 1] + rs[len(rs) // 2])
    worse = [r[0] for r in rows if abs(r[2]) > 1.5 * abs(r[1])]
    print("\n%d molecules. worst |dRyd| baseline %.5f, renat5 %.5f, factor %.2f"
          % (len(rows), wb, wr, (wb / wr) if wr else float("inf")))
    print("median per-molecule improvement factor %.2f  (H1 wants >=5 on the worst, H2 wants no"
          " member above 1.5x its baseline)" % med)
    print("members worse than 1.5x their baseline: %s" % (", ".join(worse) or "none"))
    print("floor 1e-5 e (gennbo prints five decimals); anything at 0.00001 is at the floor")
