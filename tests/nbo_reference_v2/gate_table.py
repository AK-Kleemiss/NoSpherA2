"""Did the core-block partition move AGREEMENT WITH THE REFERENCE, and in which direction?

    python3 gate_table.py <A.summary.json> <B.summary.json>

The arm JSONs already say how many points failed per quantity per molecule, so this needs no
rerun.  A is the candidate (arm_core), B the baseline.  Reported per quantity: total failed and
total compared on each side, and the number of molecules where A is better / worse / equal, so a
net improvement made of one molecule gaining and another losing cannot hide inside a sum.
"""
import json
import sys
from collections import defaultdict


def load(p, view):
    #rows keyed by molecule -> quantity name -> (failed, compared, max_dev)
    out = {}
    d = json.load(open(p))
    for r in d["rows"]:
        q = {x["name"]: (x["failed"], x["compared"], x.get("max_dev")) for x in r[view]["quantities"]}
        out[r["molecule"]] = q
    return out, d


def main(argv):
    view = "gated" if "--gated" in argv else "ungated"
    A, dA = load(argv[1], view)
    B, dB = load(argv[2], view)
    print("A = %s\nB = %s\nview = %s" % (argv[1], argv[2], view))
    names = []
    for m in A:
        for q in A[m]:
            if q not in names:
                names.append(q)
    print("%-16s %18s %18s   %s" % ("quantity", "A failed/compared", "B failed/compared",
                                    "mols better/worse/equal (A vs B)"))
    tot = [0, 0, 0, 0]
    for q in names:
        fa = ca = fb = cb = 0
        better = worse = equal = 0
        for m in sorted(set(A) & set(B)):
            a, b = A[m].get(q), B[m].get(q)
            if a is None or b is None:
                continue
            fa += a[0]; ca += a[1]; fb += b[0]; cb += b[1]
            if a[0] < b[0]:
                better += 1
            elif a[0] > b[0]:
                worse += 1
            else:
                equal += 1
        tot = [tot[0] + fa, tot[1] + ca, tot[2] + fb, tot[3] + cb]
        print("%-16s %10d/%7d %10d/%7d   %d/%d/%d" % (q, fa, ca, fb, cb, better, worse, equal))
    print("%-16s %10d/%7d %10d/%7d   net failed %+d" % (
        "TOTAL", tot[0], tot[1], tot[2], tot[3], tot[0] - tot[2]))

    #Where a molecule moved, say which way, so "mostly better" is a count and not an impression.
    moved = defaultdict(list)
    for m in sorted(set(A) & set(B)):
        da = sum(v[0] for v in A[m].values()) - sum(v[0] for v in B[m].values())
        if da:
            moved["better" if da < 0 else "worse"].append("%s %+d" % (m, da))
    print("\nper molecule, total failed points A-B:")
    for k in ("better", "worse"):
        print("  %-7s %2d: %s" % (k, len(moved[k]), "  ".join(moved[k]) or "none"))
    print("  equal   %2d" % (len(set(A) & set(B)) - len(moved["better"]) - len(moved["worse"])))
    for tag, d in (("A", dA), ("B", dB)):
        print("%s: %d PASS %d FAIL, %d points, gated %d points (%d dropped), host %s" % (
            tag, d["n_pass"], d["n_fail"], d["n_points"], d["gated"]["n_points"],
            d["gated"]["n_dropped"], d.get("host")))
    return 0


def demo():
    #A row that is better on one molecule and worse on another must not read as "no change".
    import io, os, tempfile
    def w(f, rows):
        json.dump({"host": "h", "n_pass": 0, "n_fail": 2, "n_points": 9,
                   "gated": {"n_pass": 0, "n_points": 4, "n_dropped": 5},
                   "rows": [{"molecule": m, "ungated": {"quantities": [
                       {"name": "npa charge", "failed": f_, "compared": 4, "max_dev": 0.1}]}}
                            for m, f_ in rows]}, open(f, "w"))
    d = tempfile.mkdtemp()
    a, b = os.path.join(d, "a.json"), os.path.join(d, "b.json")
    w(a, [("x", 1), ("y", 3)])
    w(b, [("x", 3), ("y", 1)])
    buf, old = io.StringIO(), sys.stdout
    sys.stdout = buf
    main(["gate_table.py", a, b])
    sys.stdout = old
    out = buf.getvalue()
    assert "net failed +0" in out, out            # the sum alone would say nothing happened
    assert "1/1/0" in out, out                    # one better, one worse, none equal
    assert "better   1: x -2" in out and "worse    1: y +2" in out, out
    print("demo ok")
    return 0


if __name__ == "__main__":
    sys.exit(demo() if "--demo" in sys.argv else main(sys.argv))
