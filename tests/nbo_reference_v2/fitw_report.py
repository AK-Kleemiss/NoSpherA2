"""Merge every fitw_result*.json under the ops directory and print fitw.py's report tables.

    py -3.12 fitw_report.py <ops dir>

The ops directory is the one fitw.py was run against: it holds ops_<mol>/ per molecule and the
fitw_result*.json fitw.py wrote there.  Nothing here fits anything - every number is either read back
out of those json files or recomputed from them with fitw's own instruments, so this script cannot
change a verdict, only tabulate it.
"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fitw import MOLS, compose, prepare                              # noqa: E402
from dmpnao_step4 import sines                                       # noqa: E402

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "ops")
rows = {}
for p in sorted(glob.glob(os.path.join(D, "fitw_result*.json"))):
    for r in json.load(open(p)):
        rows[r["mol"]] = r
order = [m for m in MOLS if m in rows]
# recheck_recovery.py's certified recovery, where control 1 was re-run at the treatment's own budget
RE = json.load(open(os.path.join(D, "fitw_recheck.json"))) if os.path.exists(
    os.path.join(D, "fitw_recheck.json")) else {}


def attained(m, a, v):
    """What the solver demonstrably reaches on THIS arm: the certified recovery if one was run."""
    return RE.get(m, {}).get(a, {}).get("worst", v["recover"]["worst"])
print("molecules present: %s (%d of 8)\n" % (" ".join(order), len(order)))

print("== the fit, per molecule against its own measured floor ==")
print("  fit/attain is the CONSERVATIVE margin: the fitted worst over max(floor, this arm's own")
print("  recovery residual), so an arm whose solver does not itself reach the floor is judged")
print("  against what the solver demonstrably attains, not against what the metric could resolve.")
print("%-8s %-6s %4s %7s %5s %4s | %-9s | %-9s %-9s %-9s | %-9s %-9s %-9s | %-12s %-9s %s" % (
    "mol", "arm", "n", "eq", "unk", "pin", "floor", "base Cor", "base Val", "base Ryd",
    "fit Cor", "fit Val", "fit Ryd", "call", "fit/floor", "fit/attain"))
for m in order:
    r = rows[m]
    for a in ("base", "renat"):
        if a not in r["arms"]:
            continue
        v = r["arms"][a]
        b, f = v["base"], v["vfit"]
        att = max(v["floor"], attained(m, a, v))
        print("%-8s %-6s %4d %7d %5d %4d | %9.2e | %9.3e %9.3e %9.3e | %9.3e %9.3e %9.3e | %-12s %9.1e %.1e" % (
            m, a, r["n"], r["neq"], r["n"] - 3, v["fit"]["pinned"], v["floor"],
            b["Cor"], b["Val"], b["Ryd"], f["Cor"], f["Val"], f["Ryd"], v["call"],
            f["worst"] / v["floor"], f["worst"] / att))

print("\n== class sizes, so that no maximum above is an empty one ==")
for m in order:
    g = prepare(m, os.path.join(D, "ops_" + m))
    c = [int(sum(1 for k in g["cls"] if k == j)) for j in (0, 1, 2)]
    print("%-8s n=%-4d Cor %-4d Val %-4d Ryd %-4d degenerate columns %d" % (
        m, g["n"], c[0], c[1], c[2], g["ndeg"]))

print("\n== controls ==")
print("%-8s %-6s | %-9s %-9s %-9s | %-9s %-9s | %-9s %-9s | %4s %-9s %-9s" % (
    "mol", "arm", "rec res", "rec worst", "rec/floor", "neg res", "neg worst", "WRONG_id",
    "WRONG_low", "ndeg", "moves", "|Sp-1|"))
for m in order:
    r = rows[m]
    for a in ("base", "renat"):
        if a not in r["arms"]:
            continue
        v = r["arms"][a]
        print("%-8s %-6s | %9.2e %9.2e %9.1e | %9.2e %9.2e | %9.3e %9.3e | %4d %9.2e %9.2e" % (
            m, a, v["recover"]["res"], v["recover"]["worst"], v["recover"]["worst"] / v["floor"],
            v["negative"]["res"], v["negative"]["worst"], r["WRONG_identity"], r["WRONG_lowdin"],
            r["ndeg"], v["moves"], r["sp_off"]))

if RE:
    print("\n== control 1 re-run at the treatment's own multistart budget (recheck_recovery.py) ==")
    print("  fitw.py fits the planted target from ONE start while the real fit gets 3, so a recovery")
    print("  that misses is not yet an incompetent solver - the two had different budgets.")
    for m in order:
        for a in ("base", "renat"):
            if m in RE and a in RE[m]:
                e, v = RE[m][a], rows[m]["arms"][a]
                print("  %-8s %-6s recovery worst  1 start %.4e  ->  %d starts %.4e   (fit %.4e, "
                      "conservative margin %.1e)" % (m, a, e["worst_1start"], e["starts"], e["worst"],
                                                     v["vfit"]["worst"],
                                                     v["vfit"]["worst"] / max(v["floor"], e["worst"])))

print("\n== floor draws (the 4 re-quantised verdicts each floor is 3x the spread of) ==")
for m in order:
    for a in ("base", "renat"):
        if a in rows[m]["arms"]:
            print("%-8s %-6s %s" % (m, a, " ".join("%.6e" % x for x in rows[m]["arms"][a]["draws"])))

print("\n== REPRODUCTION control: worst over molecules, base arm at the shipped weights ==")
for a in ("base", "renat"):
    for k in ("Cor", "Val", "Ryd"):
        vals = [(rows[m]["arms"][a]["base"][k], m) for m in order if a in rows[m]["arms"]]
        if vals:
            x, m = max(vals)
            print("  %-6s %-4s %.4e  (%s)" % (a, k, x, m))
print("  V2.5 section 15.6 published, base +split : Cor 4.460e-05 (so2)  Val 1.661e-01 (benzene)"
      "  Ryd 9.984e-01 (ethane)")
print("  V2.5 section 15.6 published, renat+split :                                           "
      "  Ryd 9.143e-01 (lif)")
print("  (V2.5's third arm spec/renat5, Ryd 6.545e-01 on benzene, is NOT one of the two backbones")
print("   fitted here - it is quoted only so the two that are cannot be mistaken for it.)")

print("\n== where the fitted arm's worst column sits, and the fitted weights' form ==")
for m in order:
    g = None
    for a in ("base", "renat"):
        if a not in rows[m]["arms"]:
            continue
        v = rows[m]["arms"][a]
        if g is None:
            g = prepare(m, os.path.join(D, "ops_" + m))
        C = compose(np.array(v["fit"]["w"]), g, a == "renat")
        s = sines(C, g["C33"], g["S"], g["lbl"], g["cl"])
        wd = s["worst"]
        who = " ".join("%s:atom%d l%d occ%.3g n%d sine%.3e" % (
            k, d["atom"], d["l"], d["occ"], d["size"], d["sine"]) for k, d in sorted(wd.items()))
        print("%-8s %-6s %s" % (m, a, who))
        print("%-8s %-6s ryd_by_occ %s | forms %s" % (
            m, a, " ".join("%s=%.2e" % kv for kv in sorted(s["ryd_by_occ"].items())),
            " ".join("%s=%.2f" % kv for kv in sorted(v["forms"].items()))))

calls = {}
for m in order:
    for a, v in rows[m]["arms"].items():
        calls.setdefault(a, []).append(v["call"])
print()
for a, cs in sorted(calls.items()):
    print("%-6s FIT %d/%d  MISS %d/%d  INCONCLUSIVE %d/%d" % (
        a, cs.count("FIT"), len(cs), cs.count("MISS"), len(cs), cs.count("INCONCLUSIVE"), len(cs)))


# ---------------------------------------------------------------------------------------------
# Where the un-fittable mass sits: X = C32^-1 C split into the four (atom, class) blocks.  This
# reuses V2.5 section 15.7's instrument and applies it at the FITTED weights too, so it says which
# block of gennbo's transformation the composition cannot reach at any w.  Diagnostic for the next
# step, NOT a second acceptance test - and note that the TOTAL |X|^2 = tr(Sp^-1) is fixed by
# construction (X = Sp^-1/2 O gives X^t X = O^t Sp^-1 O), so only the FRACTIONS are a measurement.
# ---------------------------------------------------------------------------------------------

def blocks(X, g):
    at = np.array([t[1] for t in g["lbl"]])
    cl = np.array(g["cls"])
    same_a = at[:, None] == at[None, :]
    same_c = cl[:, None] == cl[None, :]
    tot = float((X ** 2).sum())
    out = {}
    for nm, m in (("atom/class", same_a & same_c), ("atom/xclass", same_a & ~same_c),
                  ("offat/class", ~same_a & same_c), ("offat/xclass", ~same_a & ~same_c)):
        out[nm] = float((X[m] ** 2).sum()) / tot
    out["tot"] = tot
    return out


print("Frobenius mass fraction of X = C32^-1 C, by (same atom?, same class?) block")
print("%-8s %-6s %-10s | %-11s %-11s %-11s %-11s | %s" % (
    "mol", "arm", "which", "atom/class", "atom/xclass", "offat/class", "offat/xclass", "|X|^2"))
for m in [x for x in MOLS if x in rows]:
    g = prepare(m, os.path.join(D, "ops_" + m))
    Xt = np.linalg.solve(g["C32"], g["C33"])
    bt = blocks(Xt, g)
    print("%-8s %-6s %-10s | %11.4f %11.4f %11.4f %11.4f | %.3e" % (
        m, "-", "gennbo", bt["atom/class"], bt["atom/xclass"], bt["offat/class"],
        bt["offat/xclass"], bt["tot"]))
    for a in ("base", "renat"):
        if a not in rows[m]["arms"]:
            continue
        for lab, w in (("shipped w", np.array(rows[m]["arms"][a]["fit"]["w"]) * 0 + np.maximum(g["w0"], 1e-06)),
                       ("fitted w", np.array(rows[m]["arms"][a]["fit"]["w"]))):
            X = np.linalg.solve(g["C32"], compose(w, g, a == "renat"))
            b = blocks(X, g)
            print("%-8s %-6s %-10s | %11.4f %11.4f %11.4f %11.4f | %.3e" % (
                m, a, lab, b["atom/class"], b["atom/xclass"], b["offat/class"],
                b["offat/xclass"], b["tot"]))
