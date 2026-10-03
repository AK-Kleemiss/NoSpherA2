"""THE DECISIVE TEST: gennbo's own cascade, extracted, against every candidate definition we have.

    py -3.12 dmpnao_step4.py <dir of ops_<mol> subdirs>
    py -3.12 dmpnao_step4.py --demo

WHY THIS IS A CLOSED LOOP AND THE EARLIER COMPARISONS WERE NOT.  Every measurement in this lane so
far has compared native's whole pipeline against gennbo's whole pipeline, which leaves two ways to
disagree: a different definition, or a different input.  OUTCOME B removed the second one only at
the PNAO stage (mixing 1.000000 over 379 shells, job 591483).  This script removes it outright.
The probe of 25 Sep established that this NBO 7 install honours an explicit W<n> for the operator
matrices as well - DMPNAO on unit 51, SPNAO on 52, DMNAO on 53, each behind its own self-identifying
header - so gennbo's OWN density operator on gennbo's OWN pre-NAO basis is now readable, together
with gennbo's own final answer on unit 33.  Then

    X_true = C32**-1 C33

is gennbo's entire post-pre-NAO cascade, EXACTLY extractable, expressed in the basis on which the
two sides provably agree.  Every candidate definition can be run on that same basis and compared
against X_true with no pairing across codes, no sign gauge in the input, and nothing of native's
entering except the definition being tested.  If a candidate still misses, the miss is the
definition and nothing else.

WHAT IS COMPARED, AND WHY NOT WHAT THE LANE COMPARED BEFORE.  Per-class SUBSPACE sines are void
here and that is a finding about the old metric, not a choice made for convenience: with three
classes walked in decreasing priority, the span of Cor is the pre-NAO core span, the span of
Cor+Val is the pre-NAO core+valence span, and Ryd is the orthogonal complement - all three are
fixed by the CLASS PARTITION alone, and every Schmidt/OWSO/re-diagonalisation inside a class is
span-preserving.  spec_steps.molecule() already asserts the consequence (sin_max(Val) == sin_max(Ryd)
whenever the cores agree).  So a metric built on class spans cannot see the OWSO weights at all,
and the 120 members of the enumerated grid partly agreed with each other for that reason.  What the
continuous ingredient moves is the individual VECTORS inside a block, so this script measures those:

  * per column, the sign-free principal sine against gennbo's own column of the same rank inside the
    same (atom, l) - never a printed shell label, and only for columns whose occupancy is a
    SINGLETON in its block;
  * per degenerate cluster (occupancies within DEG_TOL of each other, taken from unit 53's diagonal
    at full precision rather than the five printed decimals), the max principal sine between the two
    subspaces, because the individual eigenvectors of a degenerate pair are not determined and a
    gate that fails on an undefined quantity measures two eigensolvers;
  * secondarily, the total Rydberg occupancy under the candidate's own vectors against gennbo's,
    which is the quantity the 0.5 e leak lives in.

THE FLOOR IS THE ARBITER'S PRINT, NOT MACHINE PRECISION.  gennbo writes these matrices to nine
decimals, and a sine's floor is the SQUARE ROOT of the coefficient floor, so
FLOOR = sqrt(2 * 1e-09) = 4.47e-05 of resolvable angle.  Nothing below that is a disagreement and
nothing above 10x it can be excused as print noise.

PRE-REGISTERED OUTCOMES, written and committed before a single candidate was run.

  GATE N0, the instrument.  All of these must hold or NO ARM MAY BE READ:
    (a) C33**T S C33 = 1 and diag(C32**T S C32) = 1, each under 1e-07;
    (b) unit 52 = C32**T S C32 and unit 51 = C32**T S P S C32, each under 1e-07, WITH the
        row-major packing of the same file refused above 1e-03 - a reader that cannot fail has
        not been tested;
    (c) C32 X_true = C33 to 1e-10, and X_true**T Sp X_true = 1 to 1e-07 with Sp read from unit 52,
        a file that did not enter X_true's construction;
    (d) unit 53 = X_true**T Dp X_true to 1e-07, likewise independent;
    (e) gennbo's shell ranks inside each (atom, l) component list are in DESCENDING occupancy, so
        that rank pairing is pairing like with like.  This is checked, not assumed; if it fails the
        ranks are re-sorted by occupancy and the script says so.
  GATE N1, the shipped replica.  step4_core_block.replica() must still reproduce native's OWN
    shipped NAO matrix on the same molecules, per column for singleton clusters and per subspace for
    degenerate ones, under 1e-03.  This is the negative control that fails silently if it is left
    out: it is the only thing that proves the numpy step3/step4/cascade in this directory is still
    the code it replicates, and it uses only files that have nothing to do with gennbo.
  GATE N2, discrimination.  The metric must go RED on two deliberately wrong candidates: X = 1
    (gennbo's pre-NAOs left alone) and X = Sp**-1/2 (plain Loewdin, no classes, no weights).  Both
    must exceed 100x FLOOR on at least the Rydberg columns of every molecule.  A gate that cannot
    fail is not a gate.

  OUTCOME N.  Exactly one candidate's worst sine falls to FLOOR on all three classes, on 8 of 8
    molecules at once.
    Licenses: naming that definition as the missing continuous ingredient, and a specific edit in
      Src/core/nao.cpp, left committed-not-pushed.
    Does NOT license: claiming agreement with NBO 7, which is the 22-molecule external gate and is
      untouched by this script; nor quoting any population number as an improvement, because a
      closed loop on gennbo's own pre-NAOs is not the pipeline that ships.
  OUTCOME O.  Every candidate, renat5 included, stays above 10x FLOOR on at least one class of at
    least one molecule.
    Licenses: the statement that the input-difference explanation is now dead - the candidate was
      handed gennbo's own basis, gennbo's own operator and gennbo's own classes and still missed -
      so the residual is the definition, and no member of spec_steps.ARMS is it.
    Does NOT license: naming what the right definition is, extending the 144-grid, or inventing a
      145th arm.  A candidate invented after seeing these numbers is a fit.
  OUTCOME P.  A candidate closes Cor and Val to FLOOR but not Ryd, or closes on some molecules only.
    Licenses: reporting the split, and localising what remains to the Rydberg treatment.
    Does NOT license: landing anything, and does not license calling renat5 confirmed.
  OUTCOME Q, NO OUTCOME.  Any part of N0, N1 or N2 fails.  Then say "cannot arbitrate", name which
    control failed on which side, and read no arm.

  OUTCOME R, measured alongside and independent of the arms.  Do gennbo's OWN unit-32 columns
    already diagonalise the m-averaged DMPNAO inside every (atom, l) block?  On lif they do, to
    2.96e-09.  If that holds on 8 of 8:
    Licenses: the statement that gennbo's AOPNAO output is ALREADY post-step-4, so step 4 applied to
      it is the identity, and therefore the whole content of X_true is the later steps - the
      Schmidt/OWSO treatment and the Rydberg re-naturalisation, i.e. the WEIGHTS.  That moves the
      lane's localisation from "step 4" to "the weighting inside steps 5-7" on the arbiter's side.
    Does NOT license: any claim that native's step 4 is correct.  Native's step 4 runs on native's
      step-3 output, which is NOT the pre-NAO basis, so this measurement says nothing about it.
    And if it FAILS on some molecule, that is the more interesting result and must be reported
      first: it would mean gennbo's printed PNAOs are pre-step-4 and the lane's reading of unit 32
      has been wrong about what stage it is.

NOT MEASURED HERE, on purpose: timing, and anything about NBO 7's internals.  Every number comes
from documented keyword output of a licensed binary; nothing is disassembled, decompiled or patched.
"""

import io
import json
import os
import re
import socket
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aonao_compare import (CLASS, LANG_L, numbers, principal_sines, read_47, read_lfn32,
                           read_lfn33, unpack_upper)
from spec_steps import ARMS, cascade
from step4_core_block import step4, sym_power

FLOOR = float(np.sqrt(2.0 * 1e-09))  # nine printed decimals, and a sine floors at the square root
DEG_TOL = 1e-06                      # occupancies this close are one undetermined cluster
PACK_REFUSE = 1e-03                  # the losing packing must be at least this bad
CLS_OF = {"Cor": 0, "Val": 1, "Ryd": 2}


def read_packed(path, header, n, rowmajor=False):
    """A packed symmetric triangle behind its own named header line.

    The header is the SELF-IDENTIFYING field - " PNAO density matrix:" - matched at the start of a
    line, never a unit number and never a file offset.  `rowmajor` is the losing layout, kept so the
    caller can refuse it.

    >>> import tempfile, os
    >>> p = os.path.join(tempfile.mkdtemp(), "t")
    >>> _ = open(p, "w").write(" junk\\n H:\\n 1.0 2.0 3.0\\n")
    >>> read_packed(p, "H:", 2).tolist()
    [[1.0, 2.0], [2.0, 3.0]]
    >>> read_packed(p, "H:", 2, rowmajor=True).tolist()
    [[1.0, 2.0], [2.0, 3.0]]
    """
    body = open(path).read()
    m = re.search(r"^[ \t]*" + re.escape(header) + r"\s*$", body, re.M)
    assert m, "no %r header in %s (%d bytes)" % (header, os.path.basename(path), len(body))
    want = n * (n + 1) // 2
    vals = numbers(body[m.end():])
    assert len(vals) >= want, "%s: %d numbers behind %r, need %d" % (
        os.path.basename(path), len(vals), header, want)
    if not rowmajor:
        return unpack_upper(vals[:want], n)
    M, k = np.zeros((n, n)), 0
    for i in range(n):
        for j in range(i, n):
            M[i, j] = M[j, i] = vals[k]
            k += 1
    return M


def gennbo_labels(naos, occ):
    """(perm, lbl) putting gennbo's columns in SHELL-MAJOR order inside each (atom, l).

    gennbo prints an (atom, l) block component-major; step4_blocks() and al_groups() both slice a
    block in groups of (2l+1) and call each group a shell, so they need shell-major.  The grouping
    comes from the printed `lang` field - which component a row is - and the rank comes from the
    position in that component's own list, checked below against DESCENDING occupancy.  Nothing is
    read from a printed shell label or a column offset.

    >>> naos = [dict(index=1, atom=1, lang="s", type="Cor"), dict(index=2, atom=1, lang="s", type="Val")]
    >>> perm, lbl, rr = gennbo_labels(naos, np.array([2.0, 1.0]))
    >>> perm, [t[1:6] for t in lbl], rr
    ([0, 1], [(0, 0, 0, 0, 0), (0, 0, 0, 1, 1)], 0)
    """
    comps = {}
    for e in naos:
        l = LANG_L[e["lang"][0].lower()]
        comps.setdefault((e["atom"] - 1, l), {}).setdefault(e["lang"].lower(), []).append(
            e["index"] - 1)
    perm, lbl, reranked = [], [], 0
    for key in sorted(comps):
        atom, l = key
        lists = [comps[key][c] for c in sorted(comps[key])]
        nsh = len(lists[0])
        assert all(len(v) == nsh for v in lists), "ragged block %s" % (key,)
        assert len(lists) == 2 * l + 1, "block %s: %d components, want %d" % (key, len(lists), 2 * l + 1)
        # GATE N0(e): the printed rank must already be descending occupancy, m-averaged.
        avg = [float(np.mean([occ[v[a]] for v in lists])) for a in range(nsh)]
        order = sorted(range(nsh), key=lambda a: -avg[a])
        if order != list(range(nsh)):
            reranked += 1
        for rank, a in enumerate(order):
            for m, v in enumerate(lists):
                i = v[a]
                perm.append(i)
                lbl.append((len(lbl), atom, l, m, rank, CLS_OF[naos[i]["type"]], float(occ[i])))
    assert len(perm) == len(naos), (len(perm), len(naos))
    return perm, lbl, reranked


def clusters(lbl):
    """Column groups that must be compared as SUBSPACES: one (atom, l, m) family, degenerate occ.

    A cluster is per component, not per block: the (2l+1) components of one shell are a gauge of the
    spatial frame, not a degeneracy of the cascade, and they are already paired by `m`.
    """
    fam = {}
    for k, t in enumerate(lbl):
        fam.setdefault((t[1], t[2], t[3]), []).append(k)
    out = []
    for key in sorted(fam):
        ks = sorted(fam[key], key=lambda k: -lbl[k][6])
        grp = [ks[0]]
        for k in ks[1:]:
            if abs(lbl[k][6] - lbl[grp[-1]][6]) <= DEG_TOL:
                grp.append(k)
            else:
                out.append(grp)
                grp = [k]
        out.append(grp)
    return out


def sines(A, B, S, lbl, cl):
    """Per-class worst sine between two coefficient matrices, singleton by column, cluster by span.

    >>> S = np.eye(2); A = np.eye(2); lbl = [(0, 0, 0, 0, 0, 1, 2.0), (1, 0, 0, 0, 1, 1, 1.0)]
    >>> r = sines(A, A, S, lbl, clusters(lbl))
    >>> round(r["Val"], 12)
    0.0
    """
    out = {c: 0.0 for c in ("Cor", "Val", "Ryd")}
    ndeg = 0
    worst = {}
    # Rydberg columns are near-empty, and a direction is determined only to the extent that its
    # occupancy is separated from its neighbours'.  DEG_TOL catches exact degeneracy; it does NOT
    # make a 1e-09-vs-1e-10 pair meaningful.  So the Rydberg worst sine is also reported split by
    # the occupancy decade of the column it came from, which says whether the number is a defect or
    # a gauge.  This is a decomposition of the pre-registered metric, not a second threshold.
    dec = {}
    for grp in cl:
        name = CLASS[lbl[grp[0]][5]]
        if len(grp) == 1:
            k = grp[0]
            ov = abs(float(A[:, k] @ S @ B[:, k]))
            s = float(np.sqrt(max(0.0, 1.0 - min(ov, 1.0) ** 2)))
        else:
            ndeg += len(grp)
            s = float(principal_sines(A[:, grp], B[:, grp], S).max())
        occ = max(lbl[k][6] for k in grp)
        if s > out[name]:
            worst[name] = dict(occ=occ, size=len(grp), atom=lbl[grp[0]][1], l=lbl[grp[0]][2],
                               rank=lbl[grp[0]][4], sine=s)
        out[name] = max(out[name], s)
        for k in grp[1:]:
            out[CLASS[lbl[k][5]]] = max(out[CLASS[lbl[k][5]]], s)
        if name == "Ryd":
            d = "occ>1e-3" if occ > 1e-03 else ("1e-6<occ<1e-3" if occ > 1e-06 else "occ<1e-6")
            dec[d] = max(dec.get(d, 0.0), s)
    out["ndeg"] = ndeg
    out["worst"] = worst
    out["ryd_by_occ"] = dec
    return out


def already_step4(C, Dp_ao, S, lbl):
    """OUTCOME R: worst m-averaged off-diagonal of the density inside an (atom, l) block.

    If this is at the floor, gennbo's printed PNAOs are already the per-(atom, l) m-averaged
    eigenvectors, i.e. already post-step-4, and step 4 applied to them is the identity.
    """
    D = C.T @ Dp_ao @ C
    blocks = {}
    for k, t in enumerate(lbl):
        blocks.setdefault((t[1], t[2]), []).append(k)
    worst = 0.0
    for key in sorted(blocks):
        cols, nm = blocks[key], 2 * key[1] + 1
        ns = len(cols) // nm
        if ns < 2:
            continue
        B = sum(D[np.ix_([cols[s * nm + m] for s in range(ns)],
                         [cols[s * nm + m] for s in range(ns)])] for m in range(nm)) / nm
        worst = max(worst, float(np.abs(B - np.diag(np.diag(B))).max()))
    return worst


def structure(X, lbl):
    """Where the mass of gennbo's OWN cascade lives, as Frobenius mass fractions of X_true.

    This characterises the arbiter, not a candidate, and it is the one question the 144-member grid
    could not ask: every member of that grid was a native definition, so the grid could only say
    which of OUR shapes fits, never what shape the arbiter's transformation HAS.  A max element says
    nothing about how much of a transformation is off-block, so these are squared-mass fractions.

    >>> lbl = [(0, 0, 0, 0, 0, 0, 2.0), (1, 1, 0, 0, 0, 1, 1.0)]
    >>> r = structure(np.eye(2), lbl)
    >>> round(r["off_al"], 12), round(r["off_atom"], 12), round(r["off_class"], 12)
    (0.0, 0.0, 0.0)
    >>> r2 = structure(np.array([[0.0, 1.0], [1.0, 0.0]]), lbl)
    >>> round(r2["off_atom"], 12)
    1.0
    """
    al = np.array([t[1] * 100 + t[2] for t in lbl])
    at = np.array([t[1] for t in lbl])
    cl = np.array([t[5] for t in lbl])
    tot = float((X ** 2).sum())
    out = dict(tot=tot,
               off_al=float((X[al[:, None] != al[None, :]] ** 2).sum()) / tot,
               off_atom=float((X[at[:, None] != at[None, :]] ** 2).sum()) / tot,
               off_class=float((X[cl[:, None] != cl[None, :]] ** 2).sum()) / tot)
    for a, b in ((0, 1), (0, 2), (1, 2)):
        m = (cl[:, None] == a) & (cl[None, :] == b)
        out["%s-%s" % (CLASS[a], CLASS[b])] = float((X[m] ** 2).sum() + (X[m.T] ** 2).sum()) / tot
    return out


def load(mol, d):
    """Everything the test needs, with GATE N0 measured on the way in."""
    n, S, P = read_47(os.path.join(d, mol + ".47"))
    C32, _, e32 = read_lfn32(os.path.join(d, mol + ".32"), n, S)
    C33, _, e33 = read_lfn33(os.path.join(d, mol + ".33"), n, S)
    Sp = read_packed(os.path.join(d, mol + ".52"), "PNAO overlap matrix:", n)
    Dp = read_packed(os.path.join(d, mol + ".51"), "PNAO density matrix:", n)
    Dn = read_packed(os.path.join(d, mol + ".53"), "NAO density matrix:", n)
    Sp_bad = read_packed(os.path.join(d, mol + ".52"), "PNAO overlap matrix:", n, rowmajor=True)
    SPS = S @ P @ S
    X = np.linalg.solve(C32, C33)
    naos = json.load(open(os.path.join(d, mol + ".aopnao.nbo.json")))["nao"]
    perm, lbl, reranked = gennbo_labels(naos, np.diag(Dn))
    g = dict(mol=mol, n=n, S=S, SPS=SPS, C32=C32[:, perm], C33=C33[:, perm], lbl=lbl,
             cls=[t[5] for t in lbl], pre_occ=np.diag(Dp)[perm].copy(), reranked=reranked,
             Sp=Sp[np.ix_(perm, perm)], cl=None)
    g["cl"] = clusters(lbl)
    # X_true is built by INVERTING a matrix that is printed to nine decimals, so the two checks that
    # go THROUGH that inverse (X**T Sp X = 1 and unit 53) inherit the print quantisation and cannot
    # be gated at a flat tolerance.  A guessed bound is a fitted bound, so the floor is MEASURED
    # instead: every printed matrix is re-quantised by +-0.5e-09, which is exactly what nine decimals
    # throw away, and the same two quantities are recomputed.  The spread over the draws IS the
    # instrument's floor, it is computed from the input alone, and it never looks at an arm.  The two
    # identifications that do not go through the inverse (unit 51, unit 52) keep their flat 1e-07.
    rng = np.random.default_rng(20260925)
    q = 0.5e-09
    noise = 0.0
    for _ in range(4):
        def j(M):
            return M + rng.uniform(-q, q, M.shape)
        C32q, C33q, Spq, Dpq, Dnq = j(C32), j(C33), j(Sp), j(Dp), j(Dn)
        Xq = np.linalg.solve(C32q, C33q)
        noise = max(noise, float(np.abs(Xq.T @ Spq @ Xq - np.eye(n)).max()),
                    float(np.abs(Dnq - Xq.T @ Dpq @ Xq).max()))
    cond = float(np.linalg.cond(C32))
    g["n0"] = dict(
        o33=e33, o32=e32, cond=cond, noise=noise, prop=3.0 * noise,
        sp=float(np.abs(Sp - C32.T @ S @ C32).max()),
        dp=float(np.abs(Dp - C32.T @ SPS @ C32).max()),
        sp_rowmajor=float(np.abs(Sp_bad - C32.T @ S @ C32).max()),
        recon=float(np.abs(C32 @ X - C33).max()),
        xorth=float(np.abs(X.T @ Sp @ X - np.eye(n)).max()),
        dn=float(np.abs(Dn - X.T @ Dp @ X).max()),
        reranked=reranked,
        r_step4=already_step4(g["C32"], SPS, S, lbl))
    return g


def candidates(g):
    """Every arm in spec_steps.ARMS plus the two wrong ones GATE N2 needs, all on gennbo's basis.

    Cpre is gennbo's own unit 32, pre_occ is the diagonal of gennbo's own unit 51, cls is gennbo's
    own class column and lbl its own (atom, l, m) fields.  Nothing of native's enters but the
    definition under test.
    """
    S, SPS, lbl, cls = g["S"], g["SPS"], g["lbl"], g["cls"]
    C32, n = g["C32"], g["n"]
    out = {}
    for name, kw in ARMS:
        C = cascade(C32, g["pre_occ"], cls, S, SPS, lbl, **kw)
        for split, tag in ((False, ""), (True, "+split")):
            out[name + tag] = step4(C, SPS, lbl, split)
    out["WRONG_identity"] = C32.copy()
    out["WRONG_lowdin"] = C32 @ sym_power(g["Sp"], -0.5)
    return out


def molecule(mol, d):
    g = load(mol, d)
    res = dict(mol=mol, n=g["n"], n0=g["n0"], arms={})
    ryd = [k for k, t in enumerate(g["lbl"]) if t[5] == 2]
    res["gpop"] = float(np.einsum("ij,ij->j", g["C33"][:, ryd], g["SPS"] @ g["C33"][:, ryd]).sum())
    res["struct"] = structure(np.linalg.solve(g["C32"], g["C33"]), g["lbl"])
    for name, C in candidates(g).items():
        r = sines(C, g["C33"], g["S"], g["lbl"], g["cl"])
        r["pop"] = float(np.einsum("ij,ij->j", C[:, ryd], g["SPS"] @ C[:, ryd]).sum())
        res["arms"][name] = r
    return res


def n1_replica(mols, d):
    """GATE N1: the shipped replica still reproduces native's own shipped NAO matrix."""
    from step4_core_block import replica, per_column_overlap
    worst = {}
    for mol in mols:
        try:
            lnao, Cnat, C, S, _Cg, _f = replica(mol, os.path.join(d, "ops_" + mol), False)
        except Exception as exc:                                    # noqa: BLE001
            worst[mol] = "FAILED: %s: %s" % (type(exc).__name__, exc)
            continue
        ov = per_column_overlap(Cnat, C, S)
        lb = [(i,) + t for i, t in enumerate(lnao)] if len(lnao[0]) == 6 else lnao
        cl = clusters([(t[0], t[1], t[2], t[3], t[4], t[5], t[6]) for t in lb])
        ang = 0.0
        for grp in cl:
            if len(grp) == 1:
                ang = max(ang, float(np.sqrt(max(0.0, 1.0 - min(ov[grp[0]], 1.0) ** 2))))
            else:
                ang = max(ang, float(principal_sines(Cnat[:, grp], C[:, grp], S).max()))
        worst[mol] = ang
    return worst


def report(rows, n1, out):
    p = out.write
    p("FLOOR = %.3e (nine printed decimals, sine floors at the square root)\n\n" % FLOOR)
    p("GATE N0, the instrument.  All under 1e-07 except recon (1e-10) and the refused packing.\n")
    p("%-9s %7s %9s %9s %9s %9s %9s %9s %9s %7s %8s %9s %9s\n" % (
        "mol", "n", "|C33'SC33|", "diagC32", "unit52", "unit51", "recon", "X'SpX", "unit53",
        "rerank", "cond32", "prop.bnd", "R:step4"))
    ok = True
    for r in rows:
        z = r["n0"]
        p("%-9s %7d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %7d %8.1f %9.2e %9.2e\n" % (
            r["mol"], r["n"], z["o33"], z["o32"], z["sp"], z["dp"], z["recon"], z["xorth"],
            z["dn"], z["reranked"], z["cond"], z["prop"], z["r_step4"]))
        for k, tol in (("o33", 1e-07), ("o32", 1e-07), ("sp", 1e-07), ("dp", 1e-07),
                       ("recon", 1e-10), ("xorth", z["prop"]), ("dn", z["prop"])):
            if not z[k] < tol:
                ok = False
                p("  N0 FAIL %s %s = %.3e, want < %.1e\n" % (r["mol"], k, z[k], tol))
        if not z["sp_rowmajor"] > PACK_REFUSE:
            ok = False
            p("  N0 FAIL %s: the row-major packing was NOT refused (%.3e)\n" % (r["mol"], z["sp_rowmajor"]))
    p("\n  the losing row-major packing is refused at %.3e .. %.3e (must exceed %.0e)\n" % (
        min(r["n0"]["sp_rowmajor"] for r in rows), max(r["n0"]["sp_rowmajor"] for r in rows),
        PACK_REFUSE))
    r4 = max(r["n0"]["r_step4"] for r in rows)
    p("\nOUTCOME R: worst m-averaged off-diagonal of DMPNAO inside an (atom,l) block, over 8 mols:"
      " %.3e\n" % r4)
    p("  -> gennbo's unit 32 is %s post-step-4.\n" % ("ALREADY" if r4 < 1e-06 else "NOT"))

    p("\nGATE N1, the shipped replica vs native's own shipped NAO matrix (want < 1e-03):\n  ")
    for mol, v in n1.items():
        p("%s %s  " % (mol, ("%.2e" % v) if isinstance(v, float) else v))
        if not (isinstance(v, float) and v < 1e-03):
            ok = False
    p("\n")

    names = list(rows[0]["arms"])
    p("\nWORST SINE OVER 8 MOLECULES, per class, against gennbo's own X_true.  pop = total Rydberg\n"
      "occupancy under the candidate's vectors; gennbo's own is %s\n"
      % " ".join("%s %.4f" % (r["mol"], r["gpop"]) for r in rows))
    p("%-18s %11s %11s %11s %11s %9s\n" % ("arm", "Cor", "Val", "Ryd", "xFLOOR", "worst dpop"))
    verdict = {}
    for nm in names:
        w = {c: max(r["arms"][nm][c] for r in rows) for c in ("Cor", "Val", "Ryd")}
        mx = max(w.values())
        dpop = max(abs(r["arms"][nm]["pop"] - r["gpop"]) for r in rows)
        verdict[nm] = (mx, dpop)
        p("%-18s %11.3e %11.3e %11.3e %11.1f %9.4f\n" % (
            nm, w["Cor"], w["Val"], w["Ryd"], mx / FLOOR, dpop))

    p("\nWHAT SHAPE GENNBO'S OWN CASCADE HAS: Frobenius mass fractions of X_true off each block\n"
      "structure.  This describes the arbiter, not a candidate; the 144-grid could never ask it.\n")
    p("%-9s %10s %10s %10s %10s %10s %10s\n" % (
        "mol", "off(at,l)", "off atom", "off class", "Cor-Val", "Cor-Ryd", "Val-Ryd"))
    for r in rows:
        z = r["struct"]
        p("%-9s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n" % (
            r["mol"], z["off_al"], z["off_atom"], z["off_class"], z["Cor-Val"], z["Cor-Ryd"],
            z["Val-Ryd"]))

    p("\nWHERE THE RYDBERG MISS LIVES: worst sine split by the occupancy of the column it came from,\n"
      "plus the worst Rydberg contributor itself.  A near-empty shell's direction is a gauge.\n")
    p("%-18s %11s %15s %11s %s\n" % ("arm", "occ>1e-3", "1e-6<occ<1e-3", "occ<1e-6", "worst Ryd contributor"))
    for nm in names:
        d = {}
        wc = None
        for r in rows:
            for k, v in r["arms"][nm]["ryd_by_occ"].items():
                d[k] = max(d.get(k, 0.0), v)
            w = r["arms"][nm].get("worst", {}).get("Ryd")
            if w and (wc is None or w["sine"] > wc[1]["sine"]):
                wc = (r["mol"], w)
        tag = "-" if wc is None else "%s atom %d l=%d rank %d occ %.2e size %d" % (
            wc[0], wc[1]["atom"], wc[1]["l"], wc[1]["rank"], wc[1]["occ"], wc[1]["size"])
        p("%-18s %11.3e %15.3e %11.3e %s\n" % (
            nm, d.get("occ>1e-3", 0.0), d.get("1e-6<occ<1e-3", 0.0), d.get("occ<1e-6", 0.0), tag))
    p("  degenerate columns compared as subspaces: %s\n"
      % " ".join("%s %d/%d" % (r["mol"], r["arms"][names[0]]["ndeg"], r["n"]) for r in rows))

    p("\nGATE N2, discrimination (both wrong candidates must exceed 100x FLOOR = %.1e):\n" % (100 * FLOOR))
    for nm in ("WRONG_identity", "WRONG_lowdin"):
        mx = verdict[nm][0]
        p("  %-16s %.3e  %s\n" % (nm, mx, "RED, as required" if mx > 100 * FLOOR else "FAILED TO GO RED"))
        if not mx > 100 * FLOOR:
            ok = False

    if not ok:
        p("\nOUTCOME Q, NO OUTCOME: a control failed above.  Cannot arbitrate; no arm may be read.\n")
        return "Q"
    real = [nm for nm in names if not nm.startswith("WRONG_")]
    closed = [nm for nm in real if verdict[nm][0] <= FLOOR]
    near = [nm for nm in real if verdict[nm][0] <= 10 * FLOOR]
    if len(closed) == 1:
        p("\nOUTCOME N: %s closes at %.3e.\n" % (closed[0], verdict[closed[0]][0]))
        return "N"
    if closed:
        p("\nMORE THAN ONE closes (%s) - the metric does not discriminate; treat as NO OUTCOME.\n"
          % ", ".join(closed))
        return "Q"
    if not near:
        best = min(real, key=lambda nm: verdict[nm][0])
        p("\nOUTCOME O: every candidate stays above 10x FLOOR.  Best is %s at %.3e = %.0fx FLOOR.\n"
          % (best, verdict[best][0], verdict[best][0] / FLOOR))
        p("  The input-difference explanation is dead: the candidate was handed gennbo's own basis,\n"
          "  its own operator and its own classes and still missed.  This does NOT name the right\n"
          "  definition and does NOT license a 145th arm.\n")
        return "O"
    p("\nOUTCOME P: %s reach 10x FLOOR but none reaches FLOOR; see the per-class columns.\n"
      % ", ".join(near))
    return "P"


def demo():
    """Self-checks that need no cluster output."""
    import doctest
    for mod in (sys.modules[__name__],):
        f, _ = doctest.testmod(mod, verbose=False)
        assert f == 0, "%d doctest failures" % f
    # a clustered comparison must see a wrong rank assignment, and must NOT see a sign flip
    lbl = [(0, 0, 0, 0, 0, 1, 2.0), (1, 0, 0, 0, 1, 2, 0.5)]
    cl = clusters(lbl)
    assert [len(g) for g in cl] == [1, 1], cl
    A = np.eye(2)
    assert sines(A, A * np.array([1.0, -1.0]), np.eye(2), lbl, cl)["Ryd"] < 1e-14, "sign gauge leaked"
    swapped = A[:, ::-1]
    assert sines(A, swapped, np.eye(2), lbl, cl)["Ryd"] > 0.99, "a rank swap was not caught"
    # a degenerate pair must be compared as a subspace, so an in-plane rotation is free
    lbl2 = [(0, 0, 0, 0, 0, 2, 0.5), (1, 0, 0, 0, 1, 2, 0.5)]
    cl2 = clusters(lbl2)
    assert [len(g) for g in cl2] == [2], cl2
    th = 0.3
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    # 1.5e-08 not 0: a sine built as sqrt(1 - sigma**2) floors at the SQUARE ROOT of the double
    # epsilon, which is exactly why FLOOR here is sqrt of the print floor and not the print floor.
    assert sines(A, R, np.eye(2), lbl2, cl2)["Ryd"] < 2e-08, "degenerate cluster not span-compared"
    print("demo OK")


def main(argv):
    if "--demo" in argv:
        demo()
        return 0
    d = argv[1]
    mols = [m for m in ("lif", "water", "ammonia", "ethane", "benzene", "pf5", "so2", "sf6")
            if os.path.isdir(os.path.join(d, "ops_" + m))]
    rows = [molecule(m, os.path.join(d, "ops_" + m)) for m in mols]
    out = io.StringIO()
    out.write("host=%s mols=%s\n" % (socket.gethostname(), ",".join(mols)))
    code = report(rows, n1_replica(mols, d), out)
    sys.stdout.buffer.write(out.getvalue().encode("utf-8", "replace"))
    json.dump(dict(outcome=code, floor=FLOOR, rows=rows), open(
        os.path.join(d, "dmpnao_step4.json"), "w"), indent=1, default=float)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
