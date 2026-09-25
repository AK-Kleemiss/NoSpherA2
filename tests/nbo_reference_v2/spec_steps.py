"""The published NAO cascade has SEVEN steps; native has four.  Do the missing three close the gap?

    py -3.12 spec_steps.py --demo
    py -3.12 spec_steps.py <dir of per-molecule subdirs with .47 .32 .33 .naoc.txt .naocpre.txt
                            .aonao.nbo.json>

WHY THIS EXISTS, AND WHY IT IS READING BEFORE MEASURING.  The localisation is settled: both codes
enter native's step 4 with the SAME class spans (0 .. 3.3e-06 against a 5.7e-05 floor, 8/8), and
their FINAL valence/Rydberg class subspaces differ by a principal sine of 0.086 .. 0.245.  So the
last step owns it, and what it eats is the VECTORS its (atom, l) blocks are built from.  The open
question was which block/vector definition gennbo uses - and that question is published, so it was
read rather than bisected.

WHAT THE SPEC SAYS.  Reed/Weinstock/Weinhold's NAO construction is a SEVEN-step cascade.  The NBO 7
manual itself does not restate the algorithm (it is a program manual: NMB/NRB are defined in A.1 and
the matrices are named in B.2.5, nothing more), so the arbiter for the algorithm is an independent,
published, peer-reviewed implementation of the same paper - JANPA (Nikolaienko, Bulavin, Hovorun,
Comput. Theor. Chem. 1050, 15 (2014)), read here in its Python port's npa.py, whose own step
banners name the seven steps:

    1 pre-NAOs: per (atom, l), m-averaged, (S P S) c = w S c            native step 1  IDENTICAL
    2 split into NMB and NRB                                            native step 2  DIFFERENT
    3 OWSO within the NMB, weights = pre-NAO occupancies                native step 3  DIFFERENT
    4 Schmidt the NRB against the new NMB                               native step 3  same
    5 INTRACENTER NATURALISATION OF THE NRB, per (atom, l)              native          ABSENT
    6 OWSO within the NRB, weights = STEP 5's eigenvalues               native step 3  DIFFERENT
    7 final intracenter natural transformation, full (atom, l) block    native step 4  same key

Four differences, each a definite implementation change rather than a search direction:

  (i)   THE MISSING STEP.  The spec re-diagonalises the m-averaged density inside each (atom, l) of
        the RYDBERG set after the Schmidt and BEFORE the Rydberg OWSO, as a generalised eigenproblem
        in the non-orthogonal Schmidt-projected vectors.  Native goes straight from the Schmidt to
        the OWSO.  This is exactly a "which vectors reach the final block" difference, which is the
        localisation's own category.
  (ii)  THE RYDBERG WEIGHTS.  The spec weights the Rydberg OWSO with the step-5 eigenvalues - the
        occupancies of the Schmidt-projected, re-naturalised Rydberg orbitals.  Native weights it
        with the PRE-NAO occupancies, which are the occupancies those orbitals had before anything
        was projected out of them.  nao.cpp:349-352 argues for pre-NAO weights, but the thing it
        argues against is re-running the whole cascade with updated weights, which is not this.
  (iii) TWO CLASSES, NOT THREE.  The spec's OWSO is over the whole natural minimal basis, core and
        valence together.  Native runs a three-class cascade and Schmidt-projects the valence out of
        the core.  nao.cpp:334-336 records that pooling them moved epoxide/nh3bh3/benzene by
        0.002 e of NPA charge - but an NPA charge is now known to be blind to this defect (0.695 e
        of benzene's Rydberg set moved with all 111 charges agreeing to 1e-10), so that test did not
        measure what it was read as measuring and this arm is open again.
  (iv)  HEAVY/LIGHT RYDBERG.  The spec's default Rydberg orthogonalisation is not one OWSO: the
        orbitals above 1e-4 are OWSO'd, the rest are Schmidt-projected against them and then
        Loewdin'd.  Native does one OWSO plus a Loewdin clean-up.

WHAT THE SPEC ALSO SAYS, AND IT CONTRADICTS SOMETHING ALREADY LANDED-BUT-GATED.  Step 7's block is
keyed on (atom, l) ALONE, with no class in the key: the spec pools core, valence and Rydberg in one
final re-diagonalisation.  The core-block fix (7a466b0, gated) splits the core off.  That fix was not
argued from the spec, it was measured - native's core subspace moved onto the instrument floor on
8/8, by a factor of 138 on lif - so the two are in real tension and this script reports the spec key
and the split key side by side rather than choosing.  Both are run on every arm.

ENUMERATION, WRITTEN BEFORE ANY ARM WAS MEASURED.  The five questions the coordinator named, with
the spec's answer and which arm tests it:

  does the block key carry the class?          spec: NO, (atom, l) only          arm key=pooled/split
  pre-NAO columns or already-OWSO'd columns?   spec: already-OWSO'd (step 6 out) all arms
  OWSO within class then across, or across?    spec: WITHIN, never across; and the
                                               NMB is ONE class, not two          arm two_class
  what are the occupancy weights?              spec: pre-NAO occ for the NMB,
                                               STEP-5 eigenvalues for the NRB     arm ryd_w
  what is the Schmidt ordering?                spec: NMB first, then NRB against
                                               it; inside the NRB, heavy before
                                               light                              arm heavy_light

The arms, and the order they are tried in, is the order the spec implies - the full spec first, then
its parts, so that a hit on a part is read as a part of a definition that is written down and not as
a knob that happened to help:

  base      native as shipped                                        (reference, not an arm)
  spec      all four differences at once, spec key (pooled)
  spec_hl   spec with the heavy/light Rydberg orthogonalisation (the spec's own default)
  renat     native + the missing step 5 only (Rydberg weights still pre-NAO)
  weights   native + Rydberg OWSO weights from the post-Schmidt density (no step 5)
  twoclass  native + core and valence pooled in one OWSO

ACCEPTANCE, FIXED BEFORE THE RUN.  An arm CLOSES only if, on 8 OF 8 MOLECULES AT ONCE, the valence
AND Rydberg class-subspace sines against gennbo fall to that molecule's own instrument floor.  No
mean is reported as a verdict and no failed-point count is: a mean over 8 would have called the
equal-weights lever an improvement.  Every arm prints, per molecule, the signed change in both
classes and the count of molecules closer and further WITH THE MAGNITUDES.

CONTROLS, both directions, on every arm:

  IDENTITY (must read the floor).  The arm with every switch set to native's own choice must
    reproduce the base numbers to 1.5e-08 - the REAL floor, the square root of a 1e-16 coefficient
    floor, because a principal sine is sqrt(1 - sigma^2).  A 1e-12 assertion here is unphysical and
    has already fired once in this lane.
  DISCRIMINATION (must read different).  A 30-degree Givens rotation mixing one valence column into
    one Rydberg column of the same atom must move the class sine by about sin(30 deg) = 0.5.
  NOT FIXED BY CONSTRUCTION.  A within-class rotation does not move a class span, so the metric
    would be void if the class membership were fixed - it is not: the final step re-diagonalises a
    block that POOLS valence and Rydberg and re-assigns the classes by eigenvalue rank, so which
    vectors end up in the valence span depends on the block.  demo() shows a case where changing
    only the block definition moves the class span by 0.5 while the total span is untouched.

NOT MEASURED HERE, on purpose: timings, NPA charges (blind to this by measurement), and anything
about NRT.

--------------------------------------------------------------------------------------------------
SECOND ROUND.  Three additions, and the first of them corrects this script's own first-round report.

(A) THE CLASS SINES ARE NOT TWO INDEPENDENT NUMBERS, and I reported them as if they were.  In round 1
    the Val and Ryd sines were EQUAL to every printed digit on all 8 molecules in base, renat and
    weights (0.1025/0.1025, 0.2445/0.2445, ...) and differed only in the arms that pool the core.
    That is complementarity, not coincidence: Cor + Val + Ryd is the whole space on both sides, so
    when the core subspaces agree to the floor the largest principal angle between the two Val spaces
    EQUALS the largest between the two Ryd spaces.  So "closer on 8 of 8 for BOTH classes" was one
    quantity counted twice, and round 2 prints the Cor sine next to them so that the reading is
    visible rather than inferred, and asserts the equality where the theorem says it must hold.  The
    arms are unaffected; the strength of the claim about them is.

(B) ONE POST-HOC ARM, declared as post-hoc, with its prediction fixed before it runs.  Round 1
    refuses exactly two of the spec's four elements - (iii) two classes (Val closer 1 / further 7)
    and step 7's class-free block key (Val closer 2 / further 6) - while (i) and (ii) are closer on
    8/8.  The arm that holds the accepted parts together with native's choices on the refused ones is
    not in the round-1 grid:

      renat5    native's THREE classes, native's core-split key, plus spec step 5 AND the step-5
                Rydberg weights                                         (= spec minus (iii))
      renat5_hl renat5 with the spec's heavy/light Rydberg scheme

    PREDICTION, written before the run: if the four spec elements are separable, renat5 beats spec on
    Val and renat on Ryd, on 8 of 8.  If it does not - if it lands between them or splits the
    molecules - the four interact, no single element is "the" missing piece, and that is the finding.
    This is a combination of items already enumerated above, chosen by which ones the measurement
    refused; it is not a new knob and it is not scored as a discovery if it merely helps.

(C) THE 0.5 e ITSELF, as a secondary both codes print.  Every arm's total Rydberg population
    (sum of the diagonal of C^T S P S C over native's Rydberg columns after the final step) is
    printed against gennbo's own Rydberg occupancy sum from its .aonao.nbo.json.  The subspace sine
    is the acceptance metric; this says whether an arm moves the QUANTITY the lane exists for.  An
    arm could close the sine and leave the population, or the reverse, and both would be findings.
"""

import json
import os
import socket
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aonao_compare import CLASS, principal_sines
from step3_span import class_cols, load_sides, orthonormal
from step4_core_block import owso, step4, sym_power

HEAVY = 1.0e-04  # JANPA's own -heavyNRBthreshold default
FLOOR_SINE = 1.5e-08  # sqrt of a 1e-16 coefficient floor: the real identity-control floor


def al_groups(lbl, cols):
    """(atom, l) -> that group's columns, in the order given, shell-major with m contiguous.

    >>> lbl = [(0, 0, 0, 0, 0, 0, 2.0), (1, 0, 0, 0, 1, 2, 0.0), (2, 0, 1, 0, 0, 2, 0.0)]
    >>> sorted(al_groups(lbl, [0, 1, 2]).items())
    [((0, 0), [0, 1]), ((0, 1), [2])]
    """
    out = {}
    for c in cols:
        out.setdefault((lbl[c][1], lbl[c][2]), []).append(c)
    return out


def m_average(M, cols, nm):
    """The (ns x ns) m-average of M over the shells of one (atom, l), cols shell-major."""
    ns = len(cols) // nm
    B = np.zeros((ns, ns))
    for m in range(nm):
        idx = [cols[s * nm + m] for s in range(ns)]
        B += M[np.ix_(idx, idx)] / nm
    return B


def naturalize(C, SPS, S, lbl, cols):
    """Spec step 5: per (atom, l), the m-averaged GENERALISED eigenproblem in non-orthogonal vectors.

    Returns (C_new, weights) with weights[col] the eigenvalue of the shell that column belongs to -
    m-averaged, so the (2l+1) components of a shell share one weight and stay degenerate.  This is
    native's step 1 run again on the Schmidt-projected Rydberg set, which is what the spec's
    "intracenter naturalization of new NRBs" is: same code shape, different input vectors.

    >>> S2 = np.eye(2); C = np.eye(2)
    >>> lbl = [(0, 0, 0, 0, 0, 2, 0.0), (0, 0, 0, 0, 1, 2, 0.0)]
    >>> Cn, w = naturalize(C, np.diag([0.5, 0.1]), S2, lbl, [0, 1])
    >>> [round(x, 10) for x in w]            # descending, and the vectors are already natural
    [0.5, 0.1]
    """
    Cn = C.copy()
    w = np.zeros(C.shape[1])
    Sc, Pc = C.T @ S @ C, C.T @ SPS @ C
    for key, grp in al_groups(lbl, cols).items():
        nm = 2 * key[1] + 1
        ns = len(grp) // nm
        Sb, Pb = m_average(Sc, grp, nm), m_average(Pc, grp, nm)
        X = sym_power(Sb, -0.5)
        ev, V = np.linalg.eigh(X @ Pb @ X)
        for sh in range(ns):
            k = ns - 1 - sh  # descending occupancy
            c = X @ V[:, k]
            for m in range(nm):
                Cn[:, grp[sh * nm + m]] = sum(c[s] * C[:, grp[s * nm + m]] for s in range(ns))
                w[grp[sh * nm + m]] = ev[k]
    return Cn, w


def post_weights(C, SPS, lbl, cols):
    """The m-averaged diagonal occupancy of each column in the CURRENT vectors (arm `weights`)."""
    P = C.T @ SPS @ C
    w = np.zeros(C.shape[1])
    for key, grp in al_groups(lbl, cols).items():
        nm = 2 * key[1] + 1
        ns = len(grp) // nm
        B = m_average(P, grp, nm)
        for sh in range(ns):
            for m in range(nm):
                w[grp[sh * nm + m]] = B[sh, sh]
    return w


def orthogonalize(B, S, w, heavy_light):
    """One class's within-class orthogonalisation: one OWSO, or the spec's heavy/light scheme.

    >>> S2 = np.eye(2); B = np.array([[1.0, 0.4], [0.0, 1.0]])
    >>> Q = B @ orthogonalize(B, S2, np.array([2.0, 1e-9]), True)
    >>> bool(np.abs(Q.T @ Q - np.eye(2)).max() < 1e-10)
    True
    >>> bool(abs(np.abs(Q[:, 0] @ B[:, 0]) / np.linalg.norm(B[:, 0]) - 1.0) < 1e-10)  # heavy kept
    True
    """
    G = B.T @ S @ B
    if not heavy_light:
        T = owso(G, w)
        return T @ sym_power(T.T @ G @ T, -0.5)
    hv = [i for i in range(len(w)) if w[i] > HEAVY]
    lt = [i for i in range(len(w)) if w[i] <= HEAVY]
    T = np.eye(len(w))
    if hv:
        T[np.ix_(hv, hv)] = owso(G[np.ix_(hv, hv)], w[hv])
    G2 = T.T @ G @ T
    if hv and lt:  # Schmidt the light ones against the heavy ones
        Sc = np.eye(len(w))
        Sc[np.ix_(hv, lt)] = -G2[np.ix_(hv, lt)]
        T = T @ Sc
        G2 = T.T @ G @ T
    if lt:  # plain Loewdin among the light ones
        L = np.eye(len(w))
        L[np.ix_(lt, lt)] = sym_power(G2[np.ix_(lt, lt)], -0.5)
        T = T @ L
        G2 = T.T @ G @ T
    return T @ sym_power(G2, -0.5)


def cascade(Cpre, pre_occ, cls, S, SPS, lbl, two_class=False, renat=False, ryd_w="pre",
            heavy_light=False):
    """Steps 3-6 as one parametrised cascade; native's own choices are the defaults.

    two_class    pool core and valence into one OWSO (spec step 3) instead of native's three classes
    renat        run spec step 5 on the Rydberg set after the Schmidt
    ryd_w        "pre" (native), "post" (current m-averaged occupancies), "step5" (spec)
    heavy_light  the spec's heavy/light Rydberg scheme instead of one OWSO
    """
    assert ryd_w in ("pre", "post", "step5"), ryd_w
    assert ryd_w != "step5" or renat, "step5 weights need step 5"
    C = Cpre.copy()
    groups = ([[0, 1], [2]] if two_class else [[0], [1], [2]])
    done = np.zeros((C.shape[0], 0))
    for g in groups:
        cols = [i for i, k in enumerate(cls) if k in g]
        if not cols:
            continue
        B = C[:, cols].copy()
        if done.shape[1]:
            B -= done @ (done.T @ S @ B)
            B /= np.sqrt(np.maximum(np.einsum("ij,ij->j", B, S @ B), 1e-300))
        w = np.maximum(pre_occ[cols], 0.0)
        if g == [2]:
            C[:, cols] = B
            if renat:
                C, wn = naturalize(C, SPS, S, lbl, cols)
                B = C[:, cols].copy()
                if ryd_w == "step5":
                    w = np.maximum(wn[cols], 0.0)
            if ryd_w == "post":
                w = np.maximum(post_weights(C, SPS, lbl, cols)[cols], 0.0)
        B = B @ orthogonalize(B, S, w, heavy_light)
        C[:, cols] = B
        done = np.hstack([done, B])
    return C


ARMS = (
    ("base", {}),
    ("spec", dict(two_class=True, renat=True, ryd_w="step5")),
    ("spec_hl", dict(two_class=True, renat=True, ryd_w="step5", heavy_light=True)),
    ("renat", dict(renat=True)),
    ("weights", dict(ryd_w="post")),
    ("twoclass", dict(two_class=True)),
    ("renat5", dict(renat=True, ryd_w="step5")),  # post-hoc, see (B): spec minus two_class
    ("renat5_hl", dict(renat=True, ryd_w="step5", heavy_light=True)),
)


def class_sines(C, x, key_split):
    """Per-class final-subspace sine against gennbo, plus the Rydberg population (secondary C)."""
    S, ncls, gcls = x["S"], x["ncls"], x["gcls"]
    C4 = step4(C, x["SPS"], x["lbl"], key_split)
    out = {}
    for c in (0, 1, 2):
        a, b = class_cols(ncls, c), class_cols(gcls, c)
        if a and len(a) == len(b):
            out[CLASS[c]] = float(principal_sines(orthonormal(C4[:, a], S)[0],
                                                  orthonormal(x["Cg"][:, b], S)[0], S).max())
    ryd = class_cols(ncls, 2)
    out["pop"] = float(np.einsum("ij,ij->j", C4[:, ryd], x["SPS"] @ C4[:, ryd]).sum())
    return out


def molecule(mol, d):
    """Every arm x both block keys for one molecule, plus the identity control."""
    x = load_sides(mol, d, core_own_block=True)
    res = {"mol": mol, "n": x["n"], "floor": x["floor"], "arms": {}}
    for name, kw in ARMS:
        C = cascade(x["Cnpre"], x["pre_occ"], x["ncls"], x["S"], x["SPS"], x["lbl"], **kw)
        res["arms"][name] = {"split": class_sines(C, x, True),
                             "pooled": class_sines(C, x, False)}
    #IDENTITY CONTROL on real data: the shipped code path, reached through this script's own
    #cascade with native's switches, must agree with step3_span's step 3 to the real floor.
    from step3_span import step3
    C0 = step3(x["Cnpre"], x["pre_occ"], x["ncls"], x["S"])
    Cb = cascade(x["Cnpre"], x["pre_occ"], x["ncls"], x["S"], x["SPS"], x["lbl"])
    res["identity"] = float(principal_sines(orthonormal(C0, x["S"])[0],
                                            orthonormal(Cb, x["S"])[0], x["S"]).max())
    res["identity_cols"] = float(np.abs(np.abs(np.einsum(
        "ij,ij->j", C0, x["S"] @ Cb)) - 1.0).max())
    #(C) gennbo's own Rydberg population, from the arbiter's printed occupancies
    res["gpop"] = float(sum(e["occupancy"] for e in json.load(
        open(os.path.join(d, mol + ".aonao.nbo.json")))["nao"] if e["type"] == "Ryd"))
    #(A) COMPLEMENTARITY: with three classes and the cores agreeing at the floor, sin_max(Val) must
    #EQUAL sin_max(Ryd) - the two are not independent measurements.  Asserted where it must hold so
    #that the round-1 double count cannot be repeated silently.
    for name, kw in ARMS:
        a = res["arms"][name]["split"]
        if not kw.get("two_class") and a.get("Cor", 1.0) <= 10.0 * x["floor"] and "Val" in a:
            assert abs(a["Val"] - a["Ryd"]) < 1e-03, (mol, name, a["Val"], a["Ryd"])
    return res


def verdict(rows):
    print("")
    print("floor (instrument) per molecule:")
    for r in rows:
        print("  %-9s n=%-4d floor=%.2e  identity(span)=%.2e identity(col)=%.2e"
              % (r["mol"], r["n"], r["floor"], r["identity"], r["identity_cols"]))
    bad = [r["mol"] for r in rows if r["identity_cols"] > FLOOR_SINE]
    print("IDENTITY CONTROL at the real floor %.1e: %s" % (
        FLOOR_SINE, "PASS" if not bad else "FAIL on " + ",".join(bad)))
    for key in ("split", "pooled"):
        print("")
        print("=== final block key: %s ===" % ("(atom, l, core-split)" if key == "split"
                                               else "(atom, l) only, the spec's key"))
        base = {r["mol"]: r["arms"]["base"][key] for r in rows}
        for cl in ("Cor", "Val", "Ryd"):
            print("  %s sine vs gennbo" % cl)
            head = "    %-9s %-9s" % ("mol", "floor")
            print(head + "".join("%-10s" % n for n, _ in ARMS))
            for r in rows:
                line = "    %-9s %-9.2e" % (r["mol"], r["floor"])
                for n, _ in ARMS:
                    v = r["arms"][n][key].get(cl)
                    line += "%-10s" % ("-" if v is None else "%.4f" % v)
                print(line)
            for n, _ in ARMS:
                if n == "base":
                    continue
                d = [(r["mol"], r["arms"][n][key].get(cl, float("nan")) - base[r["mol"]].get(cl, 0.0))
                     for r in rows if cl in base[r["mol"]]]
                closer = [m for m, v in d if v < -1e-06]
                further = [m for m, v in d if v > 1e-06]
                closed = [r["mol"] for r in rows
                          if cl in r["arms"][n][key]
                          and r["arms"][n][key][cl] <= 10.0 * r["floor"]]
                print("    %-9s %s closer %d / further %d / flat %d, worst %+0.4f best %+0.4f"
                      " | AT FLOOR %d of %d" % (
                          n, cl, len(closer), len(further), len(d) - len(closer) - len(further),
                          max((v for _, v in d), default=0.0), min((v for _, v in d), default=0.0),
                          len(closed), len(d)))
    print("")
    print("(A) COUPLING: Val and Ryd are ONE measurement wherever Cor sits at the floor - the count")
    print("    'closer on 8/8 for both classes' is a single coupled quantity, not two.")
    print("")
    print("(C) total Rydberg population, e (secondary; gennbo's own printed occupancy sum):")
    print("    %-9s %-9s" % ("mol", "gennbo") + "".join("%-10s" % n for n, _ in ARMS))
    for r in rows:
        line = "    %-9s %-9.4f" % (r["mol"], r["gpop"])
        for n, _ in ARMS:
            line += "%-10.4f" % r["arms"][n]["split"]["pop"]
        print(line)
    for n, _ in ARMS:
        d = [abs(r["arms"][n]["split"]["pop"] - r["gpop"]) for r in rows]
        b = [abs(r["arms"]["base"]["split"]["pop"] - r["gpop"]) for r in rows]
        print("    %-9s |dpop| worst %.4f e, closer than base on %d of %d"
              % (n, max(d), sum(1 for i in range(len(d)) if d[i] < b[i] - 1e-06), len(d)))
    print("")
    print("ACCEPTANCE (8/8 at once, both Val and Ryd at the floor):")
    for key in ("split", "pooled"):
        for n, _ in ARMS:
            ok = [r["mol"] for r in rows
                  if all(r["arms"][n][key].get(cl, 1.0) <= 10.0 * r["floor"]
                         for cl in ("Val", "Ryd"))]
            print("  %-8s %-7s CLOSED %d of %d %s" % (
                n, key, len(ok), len(rows), "<<< CLOSES" if len(ok) == len(rows) else ""))


def demo():
    """Controls before data: the metric must read the floor on an identity and 0.5 on a mix."""
    rng = np.random.default_rng(7)
    n = 6
    #a NON-orthogonal metric and pre-NAOs that are NOT already natural: with S = 1 and eigenvectors
    #of the density as inputs every orthogonalisation in the cascade is the identity and no arm can
    #be told from another, so the discrimination check below would pass vacuously.
    A = rng.normal(size=(n, n))
    S = np.eye(n) + 0.15 * (A + A.T) / np.abs(A).max()
    S = S / np.sqrt(np.outer(np.diag(S), np.diag(S)))
    B0 = rng.normal(size=(n, n))
    P = B0 @ np.diag([1.9, 0.9, 0.05, 0.02, 0.005, 0.001]) @ B0.T
    SPS = S @ P @ S
    #one atom, one l, six shells: a core, a valence and four Rydberg
    lbl = [(0, 0, 0, 0, 0, 0, 2.0), (1, 0, 0, 0, 1, 1, 1.4)] + [
        (i, 0, 0, 0, i, 2, 0.02 / (i + 1)) for i in range(2, n)]
    cls = [t[5] for t in lbl]
    Cpre = rng.normal(size=(n, n))
    Cpre /= np.sqrt(np.einsum("ij,ij->j", Cpre, S @ Cpre))
    occ = np.array([t[6] for t in lbl])

    #IDENTITY: the arm with native's switches is the shipped cascade, to the REAL floor 1.5e-08
    #and not to 1e-12 - a principal sine is sqrt(1 - sigma^2), so 1e-16 of coefficient is 1e-08
    #of angle.  Asserting 1e-12 here has already failed once in this lane for that reason.
    from step3_span import step3
    a = cascade(Cpre, occ, cls, S, SPS, lbl)
    b = step3(Cpre, occ, cls, S)
    assert np.abs(np.abs(np.einsum("ij,ij->j", a, S @ b)) - 1.0).max() < FLOOR_SINE, "identity control"

    #DISCRIMINATION: a 30 degree valence/Rydberg mix must move the valence span by sin(30).
    val, ryd = [i for i, k in enumerate(cls) if k == 1], [i for i, k in enumerate(cls) if k == 2]
    M = a.copy()
    c, s = np.cos(np.pi / 6), np.sin(np.pi / 6)
    M[:, val[0]] = c * a[:, val[0]] + s * a[:, ryd[0]]
    M[:, ryd[0]] = -s * a[:, val[0]] + c * a[:, ryd[0]]
    moved = float(principal_sines(orthonormal(a[:, val], S)[0],
                                  orthonormal(M[:, val], S)[0], S).max())
    assert abs(moved - 0.5) < 1e-06, moved

    #NOT FIXED BY CONSTRUCTION: the class span is not invariant under the final step's block
    #definition.  Pooling atom 0's core 1s with its valence 2s lets the valence slot take the second
    #eigenvector of a mixed block; splitting the core off keeps the step-3 valence vector.  The
    #TOTAL span is identical either way (the final step is a rotation inside the block), which is
    #exactly why the CLASS span is the thing measured and the total span would be void.
    lbl2 = [(0, 0, 0, 0, 0, 0, 2.0), (1, 0, 0, 0, 1, 1, 1.0)] + [
        (i, 1, 0, 0, i - 2, 2, 0.0) for i in range(2, n)]
    cls2 = [t[5] for t in lbl2]
    base = np.linalg.qr(rng.normal(size=(n, n)))[0]
    mix = np.eye(n)
    mix[:2, :2] = [[c, -s], [s, c]]
    SPS2 = (base @ mix) @ np.diag([1.9, 0.6, 0.01, 0.005, 0.002, 0.001]) @ (base @ mix).T
    v1 = [i for i, k in enumerate(cls2) if k == 1]
    p = step4(base, SPS2, lbl2, False)
    q = step4(base, SPS2, lbl2, True)
    span_moved = float(principal_sines(orthonormal(p[:, v1], S)[0],
                                       orthonormal(q[:, v1], S)[0], S).max())
    assert span_moved > 0.1, span_moved  # the block definition alone moves the class span
    tot = float(principal_sines(orthonormal(p, S)[0], orthonormal(q, S)[0], S).max())
    #while the total span cannot move: that is what makes the class version a real test.  The
    #threshold is 1e-06 and not 1e-08 for the reason this lane has now learned three times - a
    #principal sine is sqrt(1 - sigma^2), so double precision in the vectors is ~1e-08 in the sine
    #and this control reads 4.9e-08 on an exact rotation.
    assert tot < 1e-06, tot

    #the parametrised arms must not be no-ops: each one must change the cascade somewhere
    for name, kw in ARMS:
        if not kw:
            continue
        d = np.abs(np.abs(np.einsum(
            "ij,ij->j", cascade(Cpre, occ, cls, S, SPS, lbl, **kw), S @ a)) - 1.0).max()
        assert d > 1e-10, (name, d)  # every arm must be a real change to the cascade
    print("demo: identity %.1e, discrimination %.4f, block-definition span move %.4f - all controls"
          " pass" % (np.abs(np.abs(np.einsum("ij,ij->j", a, S @ b)) - 1.0).max(), moved, span_moved))


def main(argv):
    if not argv or argv[0] == "--demo":
        demo()
        return 0
    root = argv[0]
    mols = [a for a in argv[1:] if not a.startswith("--")] or sorted(
        m for m in os.listdir(root) if os.path.isdir(os.path.join(root, m)))
    print("host=%s numpy=%s threads=1 root=%s molecules=%d"
          % (socket.gethostname(), np.__version__, root, len(mols)))
    demo()
    rows = []
    for m in mols:
        try:
            rows.append(molecule(m, os.path.join(root, m)))
        except (OSError, AssertionError, KeyError) as e:
            print("SKIP %-9s %s" % (m, e))
    if rows:
        verdict(rows)
        if "--json" in argv[1:]:
            print(json.dumps(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
