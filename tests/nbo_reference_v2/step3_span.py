"""Which stage creates the valence/Rydberg subspace disagreement: the cascade, or step 4?

    py -3.12 step3_span.py --demo
    py -3.12 step3_span.py <dir of per-molecule subdirs with .47 .32 .33 .naoc.txt .naocpre.txt
                            .aonao.nbo.json>

WHY THIS, AND WHY IT IS NOT THE THREE THINGS THE BRIEF NAMED.  The open item is the valence and
Rydberg class subspaces differing by a principal sine of 0.086 .. 0.245 (--space, OUTCOME H), which
is 15-40x more than core bleed can explain and is the only thing still holding the 0.5 e.  The three
routes named as open were the orthogonalisation itself (Schmidt vs OWSO), the occupancy weights OWSO
takes, and the class assignment.  Two of the three cannot produce a class-SUBSPACE difference at
all, and that is not an opinion:

  a class's span after step 3 is fixed by the pre-NAO set and the class partition alone.  Step 3
  Schmidt-projects the class out of the span of the higher classes and then right-multiplies by an
  invertible matrix (OWSO, then a Loewdin clean-up).  Right-multiplication by anything invertible
  does not move a span, and the projector depends on the higher classes' SPAN, not on which vectors
  were chosen inside it.  So no choice of weights and no choice between Schmidt and OWSO can change
  the valence span by one part in 1e-15.

  That is already measured, not just argued: NAO_OWSO_OFF (job 586936) drops the weighting entirely
  and every molecule's class total after step 3 is unchanged to all five printed decimals.  The
  weighting moves the vectors inside the class and the DISTRIBUTION over (atom, l) blocks; it cannot
  move the class span.  nao.cpp:363-370 records that run.

  The third route, the class assignment, is instrumented and agrees: assert_same_classes() in
  aonao_compare.py fails the run unless both sides put the same number of each (atom, l) into each
  class, and it passes on 8/8; and the pre-NAOs themselves agree at 1.000000 over 379 shells rank by
  rank (OUTCOME B).  Same pre-NAO vectors, same partition of them -> same step-3 spans on both sides.

  So if the two sides' step-3 spans are the same and their final spans differ by 0.245, the
  difference is made by what comes AFTER step 3 on one side or the other.  Step 4 is the only thing
  there, and it has the mechanism ready: it re-diagonalises the m-averaged density inside each
  (atom, l) block with valence and Rydberg POOLED (nao.cpp:489-505), and writes the eigenvectors
  back into the same column slots in eigenvalue order.  The valence slot of a pooled block gets the
  top eigenvector of a mixed valence/Rydberg block, which is in general not in the step-3 valence
  span.  The core half of exactly this mechanism was measured last week and fixed (core sine
  7.9e-04..6.3e-03 -> 1.1e-05..6.3e-05).

THE MEASUREMENT.  Per side and per class, the principal sine between

    that side's own step-3 class span   and   that side's own final class span

so no cross-code gauge, no pairing and no sign convention enters any number that carries the
verdict.  Native's step-3 span comes from the numpy replica of steps 3 and 4 that already reproduces
native's shipped matrix (step4_core_block.py, gate A); gennbo's comes from its own unit-32 pre-NAOs
and its own class labels, cascaded the same way.  Two cross-code numbers are printed alongside as
context only: the step-3 spans against each other (the premise) and the final spans against each
other (the published 0.086..0.245, so this script can be checked against it).

PRE-REGISTERED, written before the run.

  premise.  sin(native step-3 span, gennbo step-3 span) is at the floor for all three classes on
    8/8.  If it is NOT, then the pre-NAO agreement or the class agreement does not do what the
    paragraph above says it does, this whole argument is wrong, and THAT is the finding - report it
    first and stop.

  H1, step 4 makes it.  Native's Val and Ryd spans move from step 3 to final by about the same
    0.086..0.245, and gennbo's do not move (its step-3 and final Val spans agree to the floor).
    Then the disagreement is native's step-4 valence/Rydberg pooling, gennbo does not pool, and the
    candidate fix is a valence-only block - the half of NAO_CLASS_SPLIT that the refuted arm ran
    together with everything else.  Note what this would NOT settle: NAO_CLASS_SPLIT's arm measured
    the intra-atomic leak 5.6x WORSE, so H1 would put the subspace evidence and the leak evidence in
    direct conflict and both have to be reported, with no fix shipped on the strength of one.

  H2, both pool.  Both sides' Val spans move by a comparable amount and the two finals still differ.
    Then the pooling is common and the difference is in the operator being diagonalised or in the
    m-averaging - which is a different suspect from all three named routes, and the next instrument
    is the m-averaged block spectrum, not the cascade.

  H3, neither moves.  Both sides preserve their step-3 spans to the floor, yet the finals differ by
    0.245.  That is arithmetically impossible together with the premise, so it would mean the floor
    is being read wrong somewhere in this script; NO OUTCOME, and the script is the bug.

CONTROLS, both directions, before any row is read (demo()):
  must call DIFFERENT   one 30-degree Givens rotation mixing a valence column into a Rydberg one
                        has to come out at sin 30 = 0.5, and a pooled step 4 on a block with
                        non-degenerate m-averaged occupancies has to move the valence span;
  must call IDENTICAL   a random orthogonal mixing INSIDE the valence columns, and a random
                        permutation of them, must leave every sine at zero - those are gauges;
  not fixed by construction: the DIFFERENT control is the demonstration.  A class span CAN move
                        across step 4, so a floor reading is a fact about the code and not about
                        the algebra.  The one quantity here that IS fixed by construction is
                        labelled as such in the output: native's step-3 span from the replica
                        versus the same span built from the pre-NAOs alone must agree, because
                        that is the span-invariance theorem this script rests on - it is a check of
                        my algebra, not evidence about NBO.

OUTCOME, 8 molecules, HEAD (core in its own block), run on Flowoffice in numpy - and it refutes two
claims made above, which is why they are left standing rather than edited away.

  the premise HOLDS, and it is the one solid result here: the two sides' step-3 class spans agree to
  0 .. 3.3e-06 against a floor of 5.7e-05 .. 7.5e-05, on all three classes and all 8 molecules.
  Both codes enter step 4 with the same core, valence and Rydberg spaces.

  OUTCOME H2, not H1: BOTH sides' step 4 moves the valence span, native by 0.154 .. 0.305 and
  gennbo by 0.183 .. 0.316.  So gennbo pools valence with Rydberg in that block too, and giving
  valence its own block would move native AWAY from gennbo rather than towards it.  That is
  independent corroboration of NAO_CLASS_SPLIT's refutation from a different instrument, and it
  dissolves the conflict H1 was written to warn about instead of creating one.

  REFUTED, first: "the next instrument is the m-averaged block spectrum".  The m-average is trivial
  for l = 0, so an s block's operator is the density restricted to the block's span and its
  eigenvectors follow from that span alone - yet the two sides' final Val spans disagree at l = 0
  by up to 0.153 (median 0.090), as much as at l = 1.  The premise that made that argument work is
  false: step 4's block is picked out by COLUMN LABELS, and after an OWSO that mixes atoms and l
  values inside a class, the columns labelled (atom a, l) no longer span anything step 3 pinned.
  Step 4 consumes the VECTORS, not the span.

  REFUTED, second, and it is this file's own headline claim: "no choice of weights and no choice
  between Schmidt and OWSO can change the valence span".  True of the step-3 span, and measured so
  by the premise - but false of the FINAL one, because of the sentence above.  Replacing the OWSO
  weights with equal ones (which turns the OWSO into a plain Loewdin exactly) moves the final Val
  span against gennbo by -0.101 to +0.132 per molecule, the same size as the disagreement being
  chased.  So the orthogonalisation and its weights are a live route after all - they just reach
  the class subspace through step 4's block definition rather than through the span.  This is also
  the mechanism of the already-recorded NAO_OWSO_OFF arm, whose step-3 class totals were flat to
  five decimals while every final number moved.

  NOT AN IDENTIFICATION: the unweighted arm is closer to gennbo on 5 of 8 and further on 3,
  including lif (+0.132) and sf6 (+0.045).  A 5/3 split is a lever, not an answer, and the
  unweighted form is refuted as a candidate here as well as by the leak gate it already failed.
  What would identify the form is an arm that closes the valence span on 8/8 at once.

  NOT QUOTABLE from this run: the Rydberg per-(atom, l) figures (0.16 .. 0.96).  They are screened
  only for carrying at least 1e-04 e, and the determinacy criterion is the step-4 block eigenvalue
  GAP, not the population.  A group holding 0.001 e has a near-degenerate spectrum and its span is
  a gauge in exactly the way the 0.37253 -> 0.36744 restatement was.
"""

import json
import os
import socket
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aonao_compare import (CLASS, LANG_L, ORDER, assert_same_classes, principal_sines, read_47,
                           read_lfn32, read_lfn33, read_naoc)
from step4_core_block import step3, step4, sym_power

try:
    from config import load_nbo
except ImportError:  # pragma: no cover - only the stamp check needs it
    load_nbo = None

NUM = {"Cor": 0, "Val": 1, "Ryd": 2}
SEP = 10.0  # a span counts as moved only above this many times the floor


def orthonormal(B, S):
    """S-orthonormal basis of span(B), span-preserving, plus the Gram's smallest eigenvalue.

    >>> S2 = np.eye(2)
    >>> Q, lo = orthonormal(np.array([[1.0, 1.0], [0.0, 1e-8]]), S2)
    >>> bool(np.abs(Q.T @ Q - np.eye(2)).max() < 1e-9), round(lo, 12)
    (True, 0.0)
    """
    G = B.T @ S @ B
    return B @ sym_power(G, -0.5), float(np.linalg.eigvalsh(G).min())


def cascade_spans(Cpre, cls, S):
    """Each class's span after the step-3 Schmidt cascade, from the pre-NAOs alone.

    The OWSO and the weights are deliberately absent: they right-multiply each block by an
    invertible matrix and cannot move its span, which is the claim this script rests on and which
    span_selfcheck() below verifies against the real step 3 rather than assuming.

    >>> S3 = np.eye(3)
    >>> sp, lo = cascade_spans(np.eye(3), [0, 1, 1], S3)
    >>> sorted(sp), float(np.abs(sp[1].T @ sp[0]).max()) < 1e-12
    ([0, 1], True)
    """
    out, low, done = {}, {}, np.zeros((Cpre.shape[0], 0))
    for c in (0, 1, 2):
        cols = [i for i, k in enumerate(cls) if k == c]
        if not cols:
            continue
        B = Cpre[:, cols].copy()
        if done.shape[1]:
            B = B - done @ (done.T @ S @ B)
        out[c], low[c] = orthonormal(B, S)
        done = np.hstack([done, out[c]])
    return out, low


def class_cols(cls, c):
    return [i for i, k in enumerate(cls) if k == c]


def side(Cspan, Cfin, cls, S):
    """sin(step-3 span, final span) per class, for ONE side - no cross-code quantity involved."""
    out = {}
    for c in (0, 1, 2):
        cols = class_cols(cls, c)
        if not cols or c not in Cspan:
            continue
        Q, _ = orthonormal(Cfin[:, cols], S)
        out[CLASS[c]] = float(principal_sines(Cspan[c], Q, S).max())
    return out


def load_sides(mol, d, core_own_block=True):
    """Everything both sides need, read from the same .47."""
    n, S, P = read_47(os.path.join(d, mol + ".47"))
    lpre, Cnpre = read_naoc(os.path.join(d, mol + ".naocpre.txt"), tag="NAOCPRE")
    lnao, Cnat = read_naoc(os.path.join(d, mol + ".naoc.txt"))
    Cg, _, e33 = read_lfn33(os.path.join(d, mol + ".33"), n, S)
    Cgpre = read_lfn32(os.path.join(d, mol + ".32"), n, S)[0]
    p = os.path.join(d, mol + ".aonao.nbo.json")
    naos = sorted((load_nbo(p) if load_nbo is not None else json.load(open(p)))["nao"],
                  key=lambda e: e["index"])
    assert_same_classes(lnao, naos)  # a Val space on one side is the Val space on the other
    assert [t[1:6] for t in lpre] == [t[1:6] for t in lnao], mol
    ncls = [t[5] for t in lnao]
    gcls = [NUM[e["type"]] for e in naos]
    nkey = [(t[1], t[2], t[5]) for t in lnao]
    gkey = [(e["atom"] - 1, LANG_L[e["lang"][0].lower()], NUM[e["type"]]) for e in naos]
    pre_occ = np.array([t[6] for t in lpre])
    SPS = S @ P @ S
    C3 = step3(Cnpre, pre_occ, ncls, S)
    C4 = step4(C3, SPS, lnao, core_own_block)
    e_n = float(np.abs(Cnat.T @ S @ Cnat - np.eye(n)).max())
    floor = float(np.sqrt(2.0 * max(e_n, e33)))
    return dict(mol=mol, n=n, S=S, SPS=SPS, floor=floor, ncls=ncls, gcls=gcls, nkey=nkey,
                gkey=gkey, lbl=lnao, pre_occ=pre_occ, core_own_block=core_own_block, C3=C3, C4=C4,
                Cnat=Cnat, Cg=Cg, Cnpre=Cnpre, Cgpre=Cgpre)


def row(mol, d, core_own_block=True):
    x = load_sides(mol, d, core_own_block)
    S, ncls, gcls = x["S"], x["ncls"], x["gcls"]
    nspan, nlow = cascade_spans(x["Cnpre"], ncls, S)
    gspan, glow = cascade_spans(x["Cgpre"], gcls, S)
    #FIXED BY CONSTRUCTION, printed as a check of the algebra above and not as evidence: the real
    #step 3, weights and OWSO included, must span exactly what the projection alone spans.
    selfchk = {CLASS[c]: float(principal_sines(
        nspan[c], orthonormal(x["C3"][:, class_cols(ncls, c)], S)[0], S).max())
        for c in nspan}
    premise = {}
    for c in nspan:
        if c in gspan and nspan[c].shape[1] == gspan[c].shape[1]:
            premise[CLASS[c]] = float(principal_sines(nspan[c], gspan[c], S).max())
    fin = {}
    for c in (0, 1, 2):
        a, b = class_cols(ncls, c), class_cols(gcls, c)
        if a and len(a) == len(b):
            fin[CLASS[c]] = float(principal_sines(
                orthonormal(x["C4"][:, a], S)[0], orthonormal(x["Cg"][:, b], S)[0], S).max())
    return dict(mol=mol, n=x["n"], floor=x["floor"], selfchk=selfchk, premise=premise,
                native=side(nspan, x["C4"], ncls, S), gennbo=side(gspan, x["Cg"], gcls, S),
                final=fin, gram=dict(native=min(nlow.values()), gennbo=min(glow.values())),
                dims={CLASS[c]: len(class_cols(ncls, c)) for c in (0, 1, 2)},
                al_occ=group_occ(x["C4"], x["SPS"], x["nkey"]),
                al=by_al(x["C4"], x["Cg"], S, x["nkey"], x["gkey"]),
                al_self=by_al(x["C4"], x["C4"], S, x["nkey"], x["nkey"]),
                arms=arms(x))


def arms(x):
    """Two step-3 forms with the SAME class spans, scored on the final valence span vs gennbo's.

    Equal weights turn the OWSO into a plain Loewdin exactly - T = w (w^2 S)^-1/2 = S^-1/2 - so the
    unweighted arm needs no second code path.  Both arms leave every class span untouched (that is
    the invariance the premise measures), so anything they move is moved by step 4 reading the
    VECTORS they chose, which is the mechanism under test.
    """
    S, ncls, lbl = x["S"], x["ncls"], x["lbl"]
    gv = orthonormal(x["Cg"][:, class_cols(x["gcls"], 1)], S)[0]
    out = {}
    for name, w in (("owso", x["pre_occ"]), ("loewdin", np.ones_like(x["pre_occ"]))):
        C4 = step4(step3(x["Cnpre"], w, ncls, S), x["SPS"], lbl, x["core_own_block"])
        nv = orthonormal(C4[:, class_cols(ncls, 1)], S)[0]
        out[name] = float(principal_sines(nv, gv, S).max())
    return out


def by_al(Cn, Cg, S, nkey, gkey):
    """sin between the two sides' final spans of each (atom, l, class) group, matched by that key.

    Legitimate as a cross-code comparison for the same reason the class version is: a span carries no
    column order and no sign.  It is NOT legitimate for the step-3 matrices, and that is the point of
    the whole measurement - after step 3 an (atom, l) subset of a class has no well-defined span,
    because the OWSO mixes atoms and l values inside the class.  Only the class span survives it.

    >>> S2 = np.eye(2)
    >>> C = np.eye(2)
    >>> k = [(0, 0, 0), (0, 1, 1)]
    >>> {a: round(b, 12) for a, b in by_al(C, C, S2, k, k).items()}
    {(0, 0, 0): 0.0, (0, 1, 1): 0.0}
    """
    out = {}
    for key in sorted(set(nkey)):
        ia = [i for i, k in enumerate(nkey) if k == key]
        ib = [i for i, k in enumerate(gkey) if k == key]
        if ia and len(ia) == len(ib):
            out[key] = float(principal_sines(orthonormal(Cn[:, ia], S)[0],
                                             orthonormal(Cg[:, ib], S)[0], S).max())
    return out


def group_occ(C, SPS, key):
    """Population of each (atom, l, class) group, so an empty group's span can be refused.

    An empty shell's step-4 eigenvalue is degenerate with its neighbours' at about 1e-12, so its
    eigenvectors are not determined by anything - the same trap that turned the Rydberg per-shell
    defect 0.37253 into 0.36744.  A span made of such columns is a gauge, not a measurement.

    >>> C = np.eye(2)
    >>> {k: round(v, 6) for k, v in group_occ(C, np.diag([1.5, 0.0]), [(0, 0, 1), (0, 0, 2)]).items()}
    {(0, 0, 1): 1.5, (0, 0, 2): 0.0}
    """
    d = np.einsum("ij,ij->j", C, SPS @ C)
    out = {}
    for i, k in enumerate(key):
        out[k] = out.get(k, 0.0) + float(d[i])
    return out


def verdict(rows, core_own_block):
    f = lambda d, k: d.get(k, float("nan"))
    print("\nstep 4 partition: %s" % ("core in its own block (HEAD)" if core_own_block
                                     else "fully pooled (NAO_CORE_POOLED=1, the old form)"))
    print("\nFIXED BY CONSTRUCTION - the real step 3 against the projection alone, i.e. a check of")
    print("this script's algebra (span invariance under OWSO), not evidence about NBO:")
    print("  " + "  ".join("%s %.0e" % (r["mol"], max(r["selfchk"].values())) for r in rows))
    print("\nPREMISE - the two sides' step-3 spans against each other.  Both sides' pre-NAOs and")
    print("class partitions agree, so these have to be at the floor or the argument is wrong:")
    print("%-9s %9s %9s %9s %9s" % ("mol", "floor", "Cor", "Val", "Ryd"))
    for r in rows:
        print("%-9s %9.1e %9.2e %9.2e %9.2e" % (r["mol"], r["floor"], f(r["premise"], "Cor"),
                                                f(r["premise"], "Val"), f(r["premise"], "Ryd")))
    bad = [r["mol"] for r in rows
           if max([v for v in r["premise"].values()] or [0.0]) > SEP * r["floor"]]
    if bad:
        print("\nPREMISE FAILS on %s.  Then the step-3 spans are NOT common, the deduction in the"
              % ", ".join(bad))
        print("header is void, and this is the finding - nothing below is interpretable.  STOP.")
        return
    print("\nthe measurement, one side at a time: sin(own step-3 span, own final span).")
    print("%-9s %26s %26s %20s" % ("", "native step3 -> final", "gennbo step3 -> final",
                                   "native vs gennbo"))
    print("%-9s %8s %8s %8s %8s %8s %8s %8s %8s %8s" % (
        "mol", "Cor", "Val", "Ryd", "Cor", "Val", "Ryd", "Cor", "Val", "Ryd"))
    for r in rows:
        print("%-9s %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e" % (
            r["mol"], f(r["native"], "Cor"), f(r["native"], "Val"), f(r["native"], "Ryd"),
            f(r["gennbo"], "Cor"), f(r["gennbo"], "Val"), f(r["gennbo"], "Ryd"),
            f(r["final"], "Cor"), f(r["final"], "Val"), f(r["final"], "Ryd")))
    print("\nworst Gram eigenvalue of a cascaded class block (an ill-conditioned span would make")
    print("the numbers above meaningless): native %.2e  gennbo %.2e" % (
        min(r["gram"]["native"] for r in rows), min(r["gram"]["gennbo"] for r in rows)))
    #WHERE, by l.  Step 4 diagonalises the m-AVERAGED density of the block.  For l = 0 there is
    #nothing to average: the block operator is the density restricted to the block's span, its
    #eigenvectors depend on that span alone, and both sides enter with the same span (the premise
    #above) - so an s block MUST come out the same on both sides whatever each side did inside its
    #class.  For l >= 1 the m-average is taken over the columns' component labels, which the
    #within-class orthogonalisation has already mixed, so there the block operator depends on the
    #VECTORS and not only on the span.  That makes the l split a real test with a predicted shape.
    print("\nwhere, by l: sin between the two sides' FINAL spans of each (atom, l, class) group.")
    print("the same instrument against native's own matrix is the IDENTICAL control, in-run:")
    print("  worst self-comparison over every group and molecule: %.1e" % max(
        max(r["al_self"].values()) for r in rows))
    OCC_MIN = 1e-04
    dropped = 0
    for cls_i in (0, 1, 2):
        per_l = {}
        for r in rows:
            for key, v in r["al"].items():
                if key[2] != cls_i:
                    continue
                if r["al_occ"].get(key, 0.0) < OCC_MIN:
                    dropped += 1
                    continue
                per_l.setdefault(key[1], []).append(v)
        if per_l:
            print("  %-3s %s" % (CLASS[cls_i], "  ".join(
                "l=%d n=%-3d worst %.2e median %.2e" % (
                    l, len(v), max(v), float(np.median(v))) for l, v in sorted(per_l.items()))))
    print("  (%d groups refused for carrying under %g e: an empty group's step-4 eigenvalues are"
          % (dropped, OCC_MIN))
    print("   degenerate at ~1e-12, so its vectors are a gauge and its span is not a measurement.)")
    #The arm the finding above points at: two step-3 forms with identical class spans, scored on
    #native's final valence span against gennbo's.  Not a fix and not gated - a build and the
    #22-molecule acceptance gate decide that, and the recorded NAO_OWSO_OFF arm already fails the
    #pre-registered leak gate in step3_stages.py (pf5 0.93 -> 1.12 x gennbo).  This says only
    #whether the within-class vector choice is where the subspace difference is made.
    print("\nsin(native final Val span, gennbo final Val span), two step-3 forms, same spans:")
    print("%-9s %10s %10s %9s" % ("mol", "owso (now)", "loewdin", "change"))
    for r in rows:
        a, b = r["arms"]["owso"], r["arms"]["loewdin"]
        print("%-9s %10.4f %10.4f %+9.4f" % (r["mol"], a, b, b - a))
    better = sum(1 for r in rows if r["arms"]["loewdin"] < r["arms"]["owso"] - 1e-06)
    worse = sum(1 for r in rows if r["arms"]["loewdin"] > r["arms"]["owso"] + 1e-06)
    print("  unweighted is closer to gennbo on %d/%d, further on %d, equal on %d" % (
        better, len(rows), worse, len(rows) - better - worse))
    moved = lambda r, sd, c: f(r[sd], c) > SEP * r["floor"]
    nv = [r["mol"] for r in rows if moved(r, "native", "Val")]
    gv = [r["mol"] for r in rows if moved(r, "gennbo", "Val")]
    fv = [r["mol"] for r in rows if moved(r, "final", "Val")]
    print("\nVal span moved across step 4 at %.0fx floor: native %d/%d, gennbo %d/%d; finals differ"
          " %d/%d" % (SEP, len(nv), len(rows), len(gv), len(rows), len(fv), len(rows)))
    if nv and not gv:
        print("OUTCOME H1: native's step 4 moves the valence span and gennbo's does not.  The")
        print("valence/Rydberg disagreement is made by native's pooled (atom, l) block, downstream")
        print("of the cascade - so Schmidt-vs-OWSO and the OWSO weights are not merely unrefuted,")
        print("they are irrelevant to it.  The candidate is a valence-only block.  This CONFLICTS")
        print("with NAO_CLASS_SPLIT's measured arm (intra-atomic leak 0.0123 -> 0.0690 e); report")
        print("both, ship neither.")
    elif nv and gv:
        print("OUTCOME H2: both sides' step 4 moves the valence span.  The pooling is common, so")
        print("giving valence its own block would move native AWAY from gennbo - corroborating")
        print("NAO_CLASS_SPLIT's refutation from a second instrument.  The m-averaging is NOT the")
        print("suspect: l = 0 has nothing to average and disagrees as much as l = 1 (see the l")
        print("table).  What step 4 consumes is the VECTORS - its blocks are picked out by column")
        print("labels, and after an OWSO that mixes atoms and l inside a class those labels span")
        print("nothing step 3 pinned.  So the within-class orthogonalisation is the live route, and")
        print("the arm above shows it is a lever of the right size without identifying the form.")
    elif not nv and not gv and fv:
        print("OUTCOME H3: neither side moves its span, yet the finals differ.  With the premise")
        print("above that is arithmetically impossible, so this script is reading a floor wrong.")
        print("NO OUTCOME.")
    else:
        print("OUTCOME: gennbo moves its valence span where native does not (%s vs %s)." % (
            ",".join(gv) or "none", ",".join(nv) or "none"))
        print("Then gennbo pools where native does not, and the fix direction is the opposite one.")


def givens(C, i, j, ang):
    out = C.copy()
    out[:, i] = np.cos(ang) * C[:, i] + np.sin(ang) * C[:, j]
    out[:, j] = -np.sin(ang) * C[:, i] + np.cos(ang) * C[:, j]
    return out


def demo():
    """Both controls, and the demonstration that a moved span is not forbidden by the algebra."""
    rng = np.random.default_rng(7)
    n = 8
    S = np.eye(n)
    C3 = np.linalg.qr(rng.normal(size=(n, n)))[0]
    cls = [0, 0, 1, 1, 1, 2, 2, 2]
    span = {c: orthonormal(C3[:, class_cols(cls, c)], S)[0] for c in (0, 1, 2)}

    #IDENTICAL control 1: the same matrix.  The bound is 1e-06 and not 1e-12 for the reason this
    #lane has now hit three times: a principal sine is sqrt(1 - sigma**2), so double precision's own
    #1e-16 in sigma surfaces as 1.5e-08 of angle.  1e-06 is still two decades under anything the
    #table can report and four under the 0.086 at issue.
    s = side(span, C3, cls, S)
    assert max(s.values()) < 1e-06, s

    #IDENTICAL control 2: an arbitrary orthogonal mixing INSIDE valence, plus a permutation of the
    #valence columns.  Both are gauges - a span does not know about them - so the metric must not
    #move.  If it does, every number in the table is about column order.
    vc = class_cols(cls, 1)
    Cg_ = C3.copy()
    Cg_[:, vc] = C3[:, vc] @ np.linalg.qr(rng.normal(size=(len(vc), len(vc))))[0]
    Cg_[:, vc] = Cg_[:, list(rng.permutation(vc))]
    s2 = side(span, Cg_, cls, S)
    assert max(s2.values()) < 1e-06, s2

    #DIFFERENT control: one 30-degree valence/Rydberg Givens rotation.  sin 30 = 0.5 exactly, on
    #both classes - Ryd is Val's complement inside Val+Ryd, so that pair is one fact twice.
    s3 = side(span, givens(C3, vc[0], class_cols(cls, 2)[0], np.pi / 6), cls, S)
    assert abs(s3["Val"] - 0.5) < 1e-12 and abs(s3["Ryd"] - 0.5) < 1e-12, s3
    assert s3["Cor"] < 1e-06, s3  # an untouched class must stay put, or the metric leaks

    #NOT FIXED BY CONSTRUCTION, the part that matters: a pooled step 4 on a real (atom, l) block
    #with non-degenerate m-averaged occupancies DOES move the valence span, and a class-split one
    #does not.  So a floor reading in the table is a fact about the code, not about the algebra.
    lbl = [(i, 0, 0, 0, sh, c, 0.0) for i, (sh, c) in enumerate([(0, 1), (1, 2), (2, 2)])]
    C = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    #P's eigenvectors must NOT be C's columns, or step 4 has nothing to rotate and this control
    #passes for the wrong reason - which is how it failed the first time it was run.
    Q = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    P = Q @ np.diag([1.6, 0.3, 0.02]) @ Q.T
    sp = {1: orthonormal(C[:, [0]], np.eye(3))[0], 2: orthonormal(C[:, [1, 2]], np.eye(3))[0]}
    pooled = side(sp, step4(C, P, lbl, core_own_block=False), [1, 2, 2], np.eye(3))
    assert pooled["Val"] > 1e-03, pooled
    print("demo ok: gauges 0, Val/Ryd Givens 0.5, pooled step 4 moves Val by %.3f" % pooled["Val"])


def main(argv):
    if not argv or argv[0] == "--demo":
        demo()
        return 0
    root = argv[0]
    pooled = "--pooled" in argv[1:]
    mols = [a for a in argv[1:] if not a.startswith("--")] or sorted(
        m for m in os.listdir(root) if os.path.isdir(os.path.join(root, m)))
    print("host=%s  numpy=%s  cwd-independent root=%s  molecules=%d" % (
        socket.gethostname(), np.__version__, root, len(mols)))
    rows = []
    for m in mols:
        try:
            rows.append(row(m, os.path.join(root, m), core_own_block=not pooled))
        except (OSError, AssertionError) as e:
            print("SKIP %-9s %s" % (m, e))
    if rows:
        verdict(rows, core_own_block=not pooled)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
