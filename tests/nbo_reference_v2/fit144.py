"""Fit the step-4 block definition to gennbo's PRINTED final NAO matrix: the whole discrete space.

    py -3.12 fit144.py --demo
    py -3.12 fit144.py <dir of per-molecule subdirs>            # the grid
    py -3.12 fit144.py <dir> --holdout <dir2>                   # only if something closes 8/8

WHY A SEARCH IS LEGITIMATE HERE, AND WHERE IT WOULD STOP BEING ONE.  spec_steps.py read the
published cascade and measured its four differences one at a time; nothing closed the acceptance
gate.  The honest statement that follows is that the two codes define the (atom, l) block on
DIFFERENT VECTORS, and the measurement that would settle it directly is gennbo's matrix at its own
step-3/step-4 boundary - which it does not print.  But gennbo DOES print the final NAO matrix
(AONAO=W, unit 33), and the space of cascade definitions is small and DISCRETE, so fitting the
definition to that printed matrix is a finite computation rather than a guess.

THE GRID, fixed before the run and not extended afterwards.  Six binary/ternary axes:

    block key       3   (atom, l) pooled | (atom, l) + core split off | (atom, l) + full class split
    class grouping  2   native's three classes | the spec's two (NMB = core+valence)
    Rydberg weights 3   pre-NAO occupancies | post-Schmidt occupancies | step-5 eigenvalues
    step 5          2   absent (native) | present (spec)
    Rydberg ortho   2   one OWSO + Loewdin (native) | heavy/light (spec)
    Schmidt order   2   NMB first, NRB against it (both codes as read) | NRB first

    3 * 2 * 3 * 2 * 2 * 2 = 144 enumerated.  24 of them are INVALID BY CONSTRUCTION - step-5
    weights with no step 5 - so 120 are measured.  That is arithmetic, not a redesign: the invalid
    cell is named here so the count cannot be quietly adjusted later.

AND THE GRID IS NARROWER THAN ITS OWN ARITHMETIC, measured after the run and reported rather than
corrected: step 5 diagonalises the Rydberg block, so its eigenvalues ARE the post-Schmidt
occupancies, and ryd_w="post" == ryd_w="step5" identically whenever step 5 runs (1e-11 on all 8
molecules).  The weights axis is therefore 3-valued only without step 5.  Same lesson as the
retracted sigma instrument in this lane: a quantity fixed by construction measures the
construction, so the DISTINCT-measurement count is printed beside the nominal 120.

THE TRAP THIS IS DESIGNED AGAINST.  A 144-way search scored on 8 molecules WILL produce a
best-fitting member by chance.  So, fixed in advance:

  * ACCEPTANCE IS THE FLOOR RULE ON 8 OF 8 AT ONCE, never best-of-144.  A definition CLOSES only if
    every molecule's core, valence AND Rydberg class sines all sit at that molecule's own instrument
    floor.  If two definitions tie, BOTH are reported.  The grid's spread is printed as context and
    labelled as context; the minimum over the grid is NOT a result and is not quotable as one.
  * ANY HIT IS CONFIRMED ON A HOLDOUT of the reference molecules outside this set, and the holdout
    result is reported whether or not it agrees.
  * THE NEGATIVE IS A RESULT.  If nothing in the 144 reaches the floor, that eliminates the whole
    discrete "native picked the wrong switch" family in one run, and what remains is CONTINUOUS -
    the functional form of the weights, or which operator is diagonalised.  The grid is then NOT
    extended: a 145th combination invented after the fact is a fit, not a measurement.
  * Cor is printed beside Val and Ryd everywhere, and the complementarity equality (Val = Ryd
    whenever the cores agree, because Cor+Val+Ryd is the whole space on both sides) is ASSERTED, so
    the double count retracted in spec_steps.py cannot come back.

CONTROLS (demo(), all at the REAL 1.5e-08 floor - a principal sine is sqrt(1 - sigma^2), so 1e-16
of coefficient is 1e-08 of angle):

  rediag == step4       this file's arbitrary-key re-diagonalisation must reproduce the shipped
                        step4() exactly on both of the keys step4() can express, on real data.
                        Without that the grid would be measuring my replica, not native's step 4.
  identity              the grid member that IS native's shipped definition must reproduce the
                        shipped path to the floor.
  discrimination        a 30 degree valence/Rydberg mix must read sin(30) = 0.5.
  non-degenerate input  the synthetic case must have a non-orthogonal S and pre-NAOs that are NOT
                        already density eigenvectors, or every orthogonalisation is the identity and
                        every arm is a no-op - which is exactly how a discrimination check passed
                        vacuously the first time it was written in this lane.
  the grid moves        the 120 valid members must not collapse to a handful of identical numbers;
                        the count of DISTINCT worst-sine values is printed.
"""

import json
import os
import socket
import sys
from itertools import product

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aonao_compare import CLASS, principal_sines
from spec_steps import FLOOR_SINE, cascade
from step3_span import class_cols, load_sides, orthonormal
from step4_core_block import step4

KEYS = ("pooled", "core", "class")
AXES = (("key", KEYS), ("two_class", (False, True)), ("ryd_w", ("pre", "post", "step5")),
        ("renat", (False, True)), ("heavy_light", (False, True)), ("nrb_first", (False, True)))
NATIVE = dict(key="core", two_class=False, ryd_w="pre", renat=False, heavy_light=False,
              nrb_first=False)
#the candidate spec_steps.py left on the table, named so the grid can be read against it by name
#rather than by rank: native's classes and key plus the spec's step 5 and its step-5 weights.
RENAT5 = dict(NATIVE, ryd_w="step5", renat=True)


def grid():
    """The 144 enumerated definitions, with the 24 invalid ones flagged rather than dropped.

    >>> g = grid()
    >>> len(g), sum(1 for d in g if d["valid"]), sum(1 for d in g if not d["valid"])
    (144, 120, 24)
    >>> sum(1 for d in g if all(d[k] == v for k, v in NATIVE.items()))   # native is in the grid
    1
    """
    out = []
    for vals in product(*[v for _, v in AXES]):
        d = dict(zip([k for k, _ in AXES], vals))
        d["valid"] = not (d["ryd_w"] == "step5" and not d["renat"])
        out.append(d)
    return out


def blocks(lbl, key):
    """The final step's blocks under each of the three keys, shell-major.

    >>> lbl = [(0, 0, 0, 0, 0, 0, 2.0), (1, 0, 0, 0, 1, 1, 1.0), (2, 0, 0, 0, 2, 2, 0.0)]
    >>> blocks(lbl, "pooled"), blocks(lbl, "core"), blocks(lbl, "class")
    ([[0, 1, 2]], [[0], [1, 2]], [[0], [1], [2]])
    """
    assert key in KEYS, key
    keyed = {}
    for i, t in enumerate(lbl):
        keyed.setdefault((t[1], t[2]), []).append(i)
    out = []
    for k in sorted(keyed):
        cols, nm = keyed[k], 2 * k[1] + 1
        if key == "pooled":
            out.append(cols)
            continue
        by = {}
        for sh in range(len(cols) // nm):
            grp = cols[sh * nm:(sh + 1) * nm]
            cls = lbl[grp[0]][5]
            by.setdefault(cls if key == "class" else (0 if cls == 0 else 1), []).extend(grp)
        out.extend(by[c] for c in sorted(by))
    return out


def rediag(C, SPS, lbl, key):
    """step4() with an arbitrary block key.  Verified against step4() itself in demo()."""
    Porb = C.T @ SPS @ C
    out = C.copy()
    for cols in blocks(lbl, key):
        nm = 2 * lbl[cols[0]][2] + 1
        ns = len(cols) // nm
        if ns <= 1:
            continue
        B = np.zeros((ns, ns))
        for m in range(nm):
            idx = [cols[s * nm + m] for s in range(ns)]
            B += Porb[np.ix_(idx, idx)] / nm
        V = np.linalg.eigh(B)[1][:, ::-1]
        for sh in range(ns):
            for m in range(nm):
                out[:, cols[sh * nm + m]] = sum(
                    V[s, sh] * C[:, cols[s * nm + m]] for s in range(ns))
    return out


def sines(C, x, key):
    """Per-class sine against gennbo's printed unit-33 matrix, plus the Rydberg population."""
    S, ncls, gcls = x["S"], x["ncls"], x["gcls"]
    C4 = rediag(C, x["SPS"], x["lbl"], key)
    out = {}
    for c in (0, 1, 2):
        a, b = class_cols(ncls, c), class_cols(gcls, c)
        if a and len(a) == len(b):
            out[CLASS[c]] = float(principal_sines(orthonormal(C4[:, a], S)[0],
                                                  orthonormal(x["Cg"][:, b], S)[0], S).max())
    ryd = class_cols(ncls, 2)
    out["pop"] = float(np.einsum("ij,ij->j", C4[:, ryd], x["SPS"] @ C4[:, ryd]).sum())
    return out


def one(x, d):
    """One definition on one molecule."""
    kw = {k: d[k] for k in ("two_class", "ryd_w", "renat", "heavy_light", "nrb_first")}
    C = cascade(x["Cnpre"], x["pre_occ"], x["ncls"], x["S"], x["SPS"], x["lbl"], **kw)
    r = sines(C, x, d["key"])
    #COMPLEMENTARITY, asserted: three classes with the cores at the floor force sin(Val) = sin(Ryd).
    if not d["two_class"] and r.get("Cor", 1.0) <= 10.0 * x["floor"] and "Val" in r:
        assert abs(r["Val"] - r["Ryd"]) < 1e-03, (d, r)
    return r


def label(d):
    return "%s/%s/%s/%s/%s/%s" % (d["key"], "2cls" if d["two_class"] else "3cls", d["ryd_w"],
                                  "s5" if d["renat"] else "--", "hl" if d["heavy_light"] else "--",
                                  "NRB1" if d["nrb_first"] else "NMB1")


def run(root, mols, tag):
    """Every valid definition on every molecule.  Returns (rows, xs)."""
    xs = {}
    for m in mols:
        try:
            xs[m] = load_sides(m, os.path.join(root, m), core_own_block=True)
            with open(os.path.join(root, m, m + ".aonao.nbo.json")) as f:  # the arbiter's own figure
                xs[m]["gpop"] = float(sum(e["occupancy"] for e in json.load(f)["nao"]
                                          if e["type"] == "Ryd"))
        except (OSError, AssertionError, KeyError) as e:
            print("SKIP %-9s %s" % (m, e))
    rows = []
    for d in grid():
        if not d["valid"]:
            continue
        per = {m: one(xs[m], d) for m in xs}
        worst = max(max(v[c] for c in ("Cor", "Val", "Ryd") if c in v) for v in per.values())
        closed = [m for m in xs if all(per[m].get(c, 1.0) <= 10.0 * xs[m]["floor"]
                                      for c in ("Cor", "Val", "Ryd") if c in per[m])]
        dpop = max(abs(per[m]["pop"] - xs[m]["gpop"]) for m in xs)
        fp = tuple(round(per[m][c], 9) for m in sorted(xs) for c in ("Cor", "Val", "Ryd", "pop"))
        rows.append(dict(d, label=label(d), per=per, worst=worst, dpop=dpop, fp=fp,
                         closed=sorted(closed)))
    print("%s: %d definitions x %d molecules, %d distinct worst-sine values"
          % (tag, len(rows), len(xs), len({round(r["worst"], 9) for r in rows})))
    #ONE AXIS IS FIXED BY CONSTRUCTION, so the grid is narrower than its own arithmetic: step 5
    #diagonalises the Rydberg block, hence its eigenvalues ARE the post-Schmidt occupancies and
    #ryd_w="post" == ryd_w="step5" whenever renat is on (verified to 1e-11 on all 8 molecules,
    #scratchpad degen.py).  Reported, not corrected - the grid was pre-registered as enumerated.
    print("%s: %d of the %d valid definitions are distinct measurements; the rest coincide, chiefly"
          % (tag, len({r["fp"] for r in rows}), len(rows)))
    print("%s: because post-Schmidt occupancies and step-5 eigenvalues are the same numbers by"
          " construction." % tag)
    return rows, xs


def report(rows, xs, tag):
    n = len(xs)
    full = [r for r in rows if len(r["closed"]) == n]
    print("")
    print("=== %s: ACCEPTANCE, the floor rule on %d of %d at once ===" % (tag, n, n))
    if full:
        print("  CLOSES on %d of %d: %d definition(s), ALL reported, no tie broken:" % (n, n, len(full)))
        for r in full:
            print("    %s" % r["label"])
    else:
        print("  CLOSES on %d of %d: NONE of the %d valid definitions." % (n, n, len(rows)))
        best = max(len(r["closed"]) for r in rows)
        print("  the most any definition closes is %d of %d molecule(s)." % (best, n))
    #CONTEXT, NOT A VERDICT: the spread of the grid, so the shape of the space is visible.  The
    #minimum is NOT an identification and is not quotable as one - that is the best-of-144 trap.
    srt = sorted(rows, key=lambda r: r["worst"])
    print("")
    print("  context only, NOT a verdict and NOT quotable as an identification - the grid's spread")
    print("  of the worst class sine, and the worst Rydberg population error, over %d molecules:" % n)
    for r in (srt[0], srt[len(srt) // 2], srt[-1]):
        print("    %-30s worst %.4f  |dpop| %.4f e  closes %d of %d"
              % (r["label"], r["worst"], r["dpop"], len(r["closed"]), n))
    nat = [r for r in rows if all(r[k] == v for k, v in NATIVE.items())][0]
    print("    %-30s worst %.4f  |dpop| %.4f e  <- native as shipped"
          % (nat["label"], nat["worst"], nat["dpop"]))
    r5 = [r for r in rows if all(r[k] == v for k, v in RENAT5.items())]
    if r5:
        print("    %-30s worst %.4f  |dpop| %.4f e  <- renat5, the candidate already on the table"
              % (r5[0]["label"], r5[0]["worst"], r5[0]["dpop"]))
    print("")
    print("  per molecule, Cor / Val / Ryd against each molecule's own floor:")
    print("    %-9s %-9s %-24s %-24s" % ("mol", "floor", "native as shipped", srt[0]["label"][:24]))
    for m in sorted(xs):
        f = xs[m]["floor"]
        a, b = nat["per"][m], srt[0]["per"][m]
        print("    %-9s %-9.2e %-24s %-24s" % (
            m, f, "%.4f/%.4f/%.4f" % (a.get("Cor", -1), a["Val"], a["Ryd"]),
            "%.4f/%.4f/%.4f" % (b.get("Cor", -1), b["Val"], b["Ryd"])))
    return full


def controls(root, mols):
    """Every control, on real data where the shipped code path is the thing being reproduced."""
    #rediag == step4 on both keys step4 can express: without this the grid measures my replica.
    worst = 0.0
    for m in mols[:3]:
        x = load_sides(m, os.path.join(root, m), core_own_block=True)
        for key, own in (("pooled", False), ("core", True)):
            a = rediag(x["C3"], x["SPS"], x["lbl"], key)
            b = step4(x["C3"], x["SPS"], x["lbl"], own)
            worst = max(worst, float(np.abs(a - b).max()))
    assert worst < 1e-12, worst  # same code path, so this one IS a 1e-12 quantity: no sine involved
    print("control rediag == step4 on both expressible keys: %.1e (max |dC|, 3 molecules)" % worst)


def demo():
    """The synthetic controls, on an input where the cascade is not degenerate."""
    from spec_steps import demo as spec_demo
    spec_demo()  # identity at 1.5e-08, discrimination 0.5000, block-key span move, no no-op arms
    g = grid()
    assert len(g) == 144 and sum(1 for d in g if d["valid"]) == 120
    rng = np.random.default_rng(11)
    n = 6
    A = rng.normal(size=(n, n))
    S = np.eye(n) + 0.15 * (A + A.T) / np.abs(A).max()
    S = S / np.sqrt(np.outer(np.diag(S), np.diag(S)))
    B0 = rng.normal(size=(n, n))
    SPS = S @ (B0 @ np.diag([1.9, 0.9, 0.05, 0.02, 0.005, 0.001]) @ B0.T) @ S
    lbl = [(0, 0, 0, 0, 0, 0, 2.0), (1, 0, 0, 0, 1, 1, 1.4)] + [
        (i, 0, 0, 0, i, 2, 0.02 / (i + 1)) for i in range(2, n)]
    #the three keys must be three different block structures on a case that HAS all three classes
    assert blocks(lbl, "pooled") != blocks(lbl, "core") != blocks(lbl, "class"), "keys collapse"
    assert blocks(lbl, "core") != blocks(lbl, "class")
    #and the reversed Schmidt order must be a real change, not a relabelling
    occ = np.array([t[6] for t in lbl])
    cls = [t[5] for t in lbl]
    Cpre = rng.normal(size=(n, n))
    Cpre /= np.sqrt(np.einsum("ij,ij->j", Cpre, S @ Cpre))
    a = cascade(Cpre, occ, cls, S, SPS, lbl)
    b = cascade(Cpre, occ, cls, S, SPS, lbl, nrb_first=True)
    d = float(np.abs(np.abs(np.einsum("ij,ij->j", a, S @ b)) - 1.0).max())
    assert d > 1e-10, d
    print("demo: grid 144 (120 valid), three keys distinct, reversed Schmidt order moves %.1e" % d)


def main(argv):
    if not argv or argv[0] == "--demo":
        demo()
        return 0
    root = argv[0]
    rest = argv[1:]
    hold = None
    if "--holdout" in rest:
        i = rest.index("--holdout")
        hold = rest[i + 1]
        rest = rest[:i] + rest[i + 2:]
    mols = [a for a in rest if not a.startswith("--")] or sorted(
        m for m in os.listdir(root) if os.path.isdir(os.path.join(root, m)))
    print("host=%s numpy=%s threads=1 root=%s molecules=%d"
          % (socket.gethostname(), np.__version__, root, len(mols)))
    demo()
    controls(root, mols)
    rows, xs = run(root, mols, "grid")
    full = report(rows, xs, "grid")
    if not full:
        print("")
        print("NEGATIVE RESULT, and it is a result: no member of the discrete space reaches the")
        print("floor on 8 of 8, so the difference is NOT a wrong switch in the cascade.  What")
        print("remains is continuous - the functional form of the weights, or which operator is")
        print("diagonalised.  The grid is NOT extended to rescue this; a 145th combination invented")
        print("after the fact would be a fit, not a measurement.")
    elif hold:
        hmols = sorted(m for m in os.listdir(hold) if os.path.isdir(os.path.join(hold, m)))
        print("")
        print("HOLDOUT: %d molecules from %s, reported whether or not it agrees" % (len(hmols), hold))
        hrows, hxs = run(hold, hmols, "holdout")
        keep = [r for r in hrows if r["label"] in {f["label"] for f in full}]
        report(keep, hxs, "holdout")
    else:
        print("")
        print("A definition CLOSES: the holdout is now REQUIRED before this is an identification.")
        print("Re-run with --holdout <dir of the reference molecules outside this set>.")
    if "--json" in argv:
        print(json.dumps([{k: v for k, v in r.items() if k not in ("per", "fp")} for r in rows]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
