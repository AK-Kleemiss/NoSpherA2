"""Per-molecule legacy -> renat5 table from two compare_all.py summary JSONs.

    py -3.12 table22.py <legacy.summary.json> <renat5.summary.json>

Prints, per molecule, the two quantities that are populations rather than labels:
NPA charge (max |dq| vs gennbo, e) and NAO occupancy (max |dpop| vs gennbo, e), each with its
failed/compared pair, then the totals.  A quantity whose `compared` count differs between the
arms is flagged, because a failure that vanished by being dropped is not a failure that agreed.
"""
import json
import sys

Q = ("npa charge", "nao occupancy")


def rows(p):
    d = json.load(open(p))
    return {r["molecule"]: {q["name"]: q for q in r["quantities"]} for r in d["rows"]}, d


A, dA = rows(sys.argv[1])
B, dB = rows(sys.argv[2])
mols = sorted(set(A) & set(B))

hdr = "%-13s | %-32s | %-32s | %s" % (
    "molecule", "NPA charge max|dq| e  fail/cmp", "NAO occ max|dpop| e   fail/cmp", "tot failed")
print(hdr)
print("-" * len(hdr))
agg = {q: [0.0, 0.0, 0, 0, 0, 0] for q in Q}   # maxA maxB failA cmpA failB cmpB
tf = [0, 0]
cnt = {q: [0, 0, 0] for q in Q}                # improved flat regressed, on max_dev
for m in mols:
    cells = []
    for q in Q:
        a, b = A[m][q], B[m][q]
        agg[q][0] = max(agg[q][0], a["max_dev"]); agg[q][1] = max(agg[q][1], b["max_dev"])
        agg[q][2] += a["failed"]; agg[q][3] += a["compared"]
        agg[q][4] += b["failed"]; agg[q][5] += b["compared"]
        flag = "" if a["compared"] == b["compared"] else "*"
        d = b["max_dev"] - a["max_dev"]
        cnt[q][0 if d < -1e-12 else (2 if d > 1e-12 else 1)] += 1
        cells.append("%9.6f -> %9.6f %3d/%-3d->%3d/%-3d%s" % (
            a["max_dev"], b["max_dev"], a["failed"], a["compared"],
            b["failed"], b["compared"], flag))
    fa = sum(x["failed"] for x in A[m].values()); fb = sum(x["failed"] for x in B[m].values())
    tf[0] += fa; tf[1] += fb
    print("%-13s | %s | %s | %4d ->%4d %s" % (
        m, cells[0], cells[1], fa, fb, "REGRESSED" if fb > fa else ("flat" if fb == fa else "")))
print("-" * len(hdr))
for q in Q:
    g = agg[q]
    print("%-14s worst over 22: %.6f -> %.6f e   failed %d/%d -> %d/%d   "
          "max_dev improved/flat/worse %d/%d/%d" % (
              q, g[0], g[1], g[2], g[3], g[4], g[5], cnt[q][0], cnt[q][1], cnt[q][2]))
print("total failed points %d -> %d" % (tf[0], tf[1]))
print("verdict A: %d PASS %d FAIL   verdict B: %d PASS %d FAIL" % (
    dA["n_pass"], dA["n_fail"], dB["n_pass"], dB["n_fail"]))
print("* = the compared-point count moved between arms, so read failed/compared as a ratio")
