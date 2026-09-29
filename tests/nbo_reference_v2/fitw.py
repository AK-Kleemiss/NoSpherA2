"""SOLVE for gennbo's OWSO weights instead of enumerating switches: fit w to X_true.

    py -3.12 -m doctest fitw.py
    py -3.12 fitw.py --demo                      # the three controls alone, synthetic
    py -3.12 fitw.py <dir of ops_<mol> subdirs>  # the fit, 8 molecules x 2 backbones

================================================================================================
THIS FILE IS THE PRE-REGISTRATION.  It is committed BEFORE a single fitted number exists, the way
da1594e was, and nothing below is edited afterwards except to append the run's own output path.
================================================================================================

WHY A FIT AND NOT A 145th ARM.  The discrete space is exhausted: 144 cascade definitions were
enumerated (block key x class grouping x Rydberg weights x step 5 x Rydberg orthogonalisation x
Schmidt order), 24 invalid by construction, 120 measured, and NOT ONE closes a single molecule
(V2.4 section 14).  V2.5 then explained that negative structurally rather than recording it: gennbo's
whole post-pre-NAO cascade is exactly extractable as X_true = C32^-1 C33, and its Frobenius mass is
39-66 % OFF-ATOM, essentially all valence<->Rydberg.  No block-local redefinition can reach mass
that is half off-atom, and every one of the 120 varied a block-local definition.  What CAN reach it
is the inter-atomic OWSO across a class - and the only continuous ingredient in an OWSO is its
WEIGHT VECTOR, which native takes from the pre-NAO occupancies at Src/core/nao.cpp:349-352.

So the weights stop being three named vectors to choose between and become n unknowns to solve for.
Extending the grid to a 145th combination would be a fit dressed as a measurement; solving for a
continuous parameter with the fit's own residual as the verdict is a measurement, and it has a
decisive answer in BOTH directions (see ACCEPTANCE).

THE OBJECT THAT IS FITTED, and why it is not fixed by construction.  With S the AO overlap, the
composition is code that is already committed, evaluated at a weight vector w:

    C(w) = step4( cascade(C32, w, cls, S, SPS, lbl, renat=<backbone>, ryd_w="pre"),
                  SPS, lbl, core_own_block=True )

and the comparison against the arbiter is the per-column overlap matrix in the AO metric

    M(w) = C33^t S C(w).

M is orthogonal for every w, because both column sets are S-orthonormal - so the gauge warning of
this lane applies in full and is quotiented out rather than measured:

  * NOT MEASURED: T T^t (identically G^-1 for any orthonormalising T), G^-1 itself, or the SINGULAR
    VALUES of M (all 1.000000 for any orthogonal M, on any input - a random permutation reproduced
    that interval exactly, which is why it was retracted).  Equivalently, writing T = Sp^-1/2 O with
    Sp = C32^t S C32, everything that is common to every orthonormalising T is divided out and only
    the ORTHOGONAL FACTOR O is compared; M = O_true^t O(w) is exactly that comparison.
  * MEASURED: how far M is from DIAGONAL.  That is not fixed by construction, and the proof is that
    the same instrument reads 0.994 on WRONG_identity and 0.994 on WRONG_lowdin (V2.5 section 15.10)
    and reproduces V2.5's whole arm table to every printed digit (see the REPRODUCTION control).

RESIDUAL.  r(w) = [ M_ij for i != j ; 1 - |M_jj| ], i.e. n^2 equations for n unknowns, heavily
overdetermined, with the per-molecule count printed in the output.  Minimised by
scipy.optimize.least_squares over theta with w = exp(theta).

GAUGES, all quotiented out explicitly, NONE of them fitted:

  1. The documented component-major -> shell-major permutation between gennbo's column order and
     native's is CONSTRUCTED by gennbo_labels() from the printed `lang` field and ASSERTED
     (`reranked == 0`, descending m-averaged occupancy), never fitted.  Handing already_step4() an
     unpermuted C32 with permuted labels was the seventh read-by-position defect in this lane and it
     failed on exactly the molecules where the permutation is not the identity.
  2. Column SIGN.  A NAO with a flipped sign is the same NAO and eigh's sign is arbitrary, so the
     diagonal enters as 1 - |M_jj| and never as 1 - M_jj.  Sign fixing cannot manufacture agreement:
     it caps a wrong column's contribution at sqrt(2) instead of 2.
  3. DEGENERATE occupancy clusters, where the individual vectors are undetermined and only the
     subspace is defined.  Those off-diagonal entries are MASKED OUT of the residual (their count is
     printed per molecule) and the verdict metric compares them as subspaces, which is what
     dmpnao_step4.sines() already does.
  4. OWSO is invariant under an overall positive rescaling of w WITHIN one class -
     T = W (W S W)^-1/2 is unchanged by W -> aW - so the reported solution is normalised to
     max(w) = 1 per class and the unknown count is stated as n - (number of classes present).
     Reporting n unknowns would overstate the fit's freedom; reporting n^2 equations against
     n - 3 unknowns is the honest count.

THE BOX, and why it is exactly the shipped code's own dynamic range.  owso() clamps w to
>= 1e-6 * max(w) within each class, so a weight below that is not representable by the shipped path
at all and the objective is flat in it.  The fit therefore searches w in [1e-6, 1] after the per-class
normalisation of gauge 4, which is precisely the range owso() can express.  Entries that converge
pinned to the lower bound are counted and reported as unidentifiable, because a parameter the
objective cannot see is not a fitted parameter.

TWO BACKBONES, both already named members of the enumerated grid - no new switch is invented:

    base   = native as shipped              (renat=False, core-split step-4 key: the C++ default)
    renat  = native + the published step 5  (renat=True,  same key) - the arm renat5 differs from
             this only in taking its Rydberg weights FROM step 5, which is exactly the choice this
             fit replaces by an unknown, so `renat` is that arm with its weights set free.

FLOOR - MEASURED, NEVER GUESSED.  A flat pre-registered 1e-07 fired on 4 of 8 molecules earlier in
this lane and a cond(C32)-based replacement was itself refuted (cond is 9.6 .. 10.0 for every
molecule while the error plainly grows with n).  So: re-quantise C32 and C33 by +-0.5e-09 (they are
printed to nine decimals) and the .47's S and P by +-0.5e-12 (twelve-digit E-format, three orders
finer, so they are not expected to set the floor - measured anyway rather than assumed), four draws,
recompute the verdict metric at the shipped weights, and take 3x the spread.  A sine's floor is the
SQUARE ROOT of the coefficient floor - that has bitten this lane three times - and the re-quantisation
recipe captures that automatically because the metric it re-evaluates IS a sine.

ACCEPTANCE, per molecule, fixed here:

    FIT           worst per-column sine at the converged w  <=  that molecule's measured floor
    MISS          worst per-column sine at the converged w  >=  10 x that floor
    INCONCLUSIVE  strictly between - reported as INCONCLUSIVE, never rounded into either

    The weights are an IDENTIFICATION only on FIT for 8 of 8 MOLECULES AT ONCE.  Best-of-N is not
    acceptance.  7 of 8 is a MISS overall and is reported as one.  A printed 0.00e+00 over 0/N
    columns is an EMPTY MAXIMUM, not an agreement, so every maximum prints the count it is over.

    The verdict is ONE number per molecule - the worst per-column sine over all classes.  Cor, Val
    and Ryd are printed only as its decomposition, because Val and Ryd are complementary subspaces
    and their principal angles are EQUAL whenever the cores agree: that identity was the twelfth
    retired metric-reading in this lane and it is not going to be reported as two confirmations
    again.

WHAT EACH OUTCOME LICENSES, written before the run:

  OUTCOME W-FIT (8 of 8 FIT).  The fitted w IS gennbo's weight vector up to the per-class scale
    gauge.  Licenses: reading the weights off and reporting the functional form they correspond to -
    tested against the pre-NAO occupancies native uses, the post-Schmidt occupancies, the step-5
    eigenvalues, net rather than gross populations, and integer powers of each, by correlation and by
    ratio; and a known one-line edit at nao.cpp:349-352.  Does NOT license quoting -nbo_native
    against NBO 7 until the 22-molecule external gate is re-run, nor skipping the holdout: any hit
    is confirmed on reference molecules outside these 8 before it is called an identification.

  OUTCOME W-MISS (any molecule MISSes).  OWSO with a positive diagonal weight vector is the wrong
    FORM: no weight choice can close the gap, on the arbiter's own basis with the arbiter's own
    operators and the arbiter's own classes.  Licenses: retiring the entire "functional form of the
    weights" branch and forcing the remaining branch - which operator is diagonalised - which is
    what the next month would otherwise be spent on.  Does NOT license any claim about WHICH
    operator, nor any statement that native's OWSO is wrong in some other respect.

  A split outcome (some FIT, some MISS, or anything INCONCLUSIVE) licenses NEITHER and is reported
    as the inconclusive result it is.

CONTROLS, all three mandatory, reported with their numbers whatever they say:

  1. RECOVERY / DISCRIMINATION.  Plant w* drawn log-uniform over [1e-4, 1] - every entry strictly
     inside the box, so every entry is identifiable - generate C(w*) through THIS file's own
     composition, and fit it starting from a different point.  PASS iff the residual returns to the
     floor AND max|log10(w_fit/w*)| after the per-class gauge is removed is small.  If this fails,
     the solver is broken and nothing else in the run means anything, and the run says so first.

  2. NEGATIVE.  Feed the fit a C that is provably NOT of the composed form: C_neg = C32 Sp^-1/2 Q,
     a plain Loewdin times an INDEPENDENT random orthogonal Q from a fixed-seed Gaussian QR.  Q is
     built with no reference to Sp, to Dp, to the labels or to any eigenvector of the problem - a
     pooled-block control in this lane first passed for the wrong reason because the planted mixing
     was built from the basis's OWN orthogonal matrix and therefore never left the span it was meant
     to leave.  The fit MUST NOT reach the floor, and its converged number is printed.

  3. NON-DEGENERATE INPUT, four readings, because a degenerate synthetic case once made every arm a
     no-op and an S = 1 case with already-natural pre-NAOs made every orthogonalisation the identity
     so the discrimination passed VACUOUSLY:
       (a) |Sp - 1| per molecule, so the orthogonalisation is not the identity;
       (b) the composition MOVES when w moves - worst sine between C at the shipped weights and at
           the shipped weights times a random per-column factor in [0.5, 2];
       (c) the instrument goes RED on two wrong answers, WRONG_identity and WRONG_lowdin;
       (d) the degenerate-cluster count per molecule, printed beside the masked-entry count.

  REPRODUCTION (a fourth, free): the verdict metric at the shipped weights must reproduce V2.5
  section 15.6's table - base+split worst Cor 4.460e-05, Val 1.661e-01, Ryd 9.984e-01 - from an
  instrument in this file.  A number that cannot reproduce the published one is not measuring the
  published thing.

NOT MEASURED HERE, on purpose: timings, NPA charges (measured blind to this defect), NRT, and
anything about the 22-molecule external gate, which this file does not touch.

AMENDMENTS AFTER THE PRE-REGISTRATION COMMIT (1409901).  Everything above is the frozen text; these
are the changes made to the code afterwards, each one forced by a CONTROL FAILURE and not by any
real-data result, and NONE of them touches the acceptance rule, the floor recipe or the backbones:

  a. resid() sign-canonicalises each column of M by the sign of its own diagonal element, and drops
     the diagonal residual for members of a degenerate cluster.  Forced by control 1: the recovery
     fit returned 1.2e+00 on lif and moved zero steps, because the raw residual JUMPED by 0.75 for a
     1e-10 step along a direction the cascade is exactly invariant under (eigh's column sign).  Gauge
     2 of the pre-registration already declared the sign to be a gauge; taking |M_jj| quotients it on
     the diagonal only, and the off-diagonal needed the same treatment.
  b. x0 for both fits is log(max(w0, 1e-6)) rather than log(w0): a pre-NAO occupancy of exactly zero
     gave -inf, i.e. outside the box of gauge 4.
  c. Two reporting-layer bugs: a KeyError in the per-row print line and a json writer that could not
     serialise an array.
  d. The solver's finite-difference step is DIFF_STEP = 1e-4, not scipy's default, and the real fit
     runs from STARTS = 3 starting points.  Also forced by control 1 and tuned ONLY on the recovery
     target, whose answer is known by construction: the default step left the optimiser differentiating
     eigen-decomposition noise (lif 2.076e-04 against a 1.26e-07 floor, terminating after 13 steps),
     and 16 (jac, diff_step, x_scale) combinations were measured there before one was fixed.  Both
     changes can only lower the achieved residual, i.e. only ever make a MISS harder to claim.
"""

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aonao_compare import read_47, read_lfn32, read_lfn33                        # noqa: E402
from dmpnao_step4 import clusters, gennbo_labels, load, read_packed, sines       # noqa: E402
from spec_steps import cascade                                                   # noqa: E402
from step4_core_block import step4, sym_power                                    # noqa: E402

MOLS = ("lif", "water", "ammonia", "ethane", "benzene", "pf5", "so2", "sf6")
BACKBONES = (("base", False), ("renat", True))
WLO, WHI = 1e-06, 1.0     # exactly owso()'s own representable dynamic range, per class
Q32, Q47 = 0.5e-09, 0.5e-12   # the two print quantisations, nine decimals and twelve digits
DRAWS = 4
SPREAD = 3.0
FIT_K, MISS_K = 1.0, 10.0     # acceptance multiples of the measured floor
SEED = 20260925
DIFF_STEP = 1e-04            # measured on the recovery target, not scipy's default
STARTS = 3                   # multistart for the real fit only


def norm_w(w, cls):
    """Divide out OWSO's per-class scale gauge: max(w) = 1 inside every class.

    >>> [round(x, 6) for x in norm_w(np.array([2.0, 1.0, 4.0]), [0, 0, 1])]
    [1.0, 0.5, 1.0]
    """
    out = np.array(w, dtype=float)
    for c in sorted(set(cls)):
        idx = [i for i, k in enumerate(cls) if k == c]
        m = max(out[idx].max(), 1e-300)
        out[idx] = out[idx] / m
    return out


def compose(w, g, renat):
    """C(w): the committed cascade + step 4, with w as the OWSO weight vector and nothing else free."""
    C = cascade(g["C32"], w, g["cls"], g["S"], g["SPS"], g["lbl"], renat=renat, ryd_w="pre")
    return step4(C, g["SPS"], g["lbl"], True)


def mask_of(g):
    """Off-diagonal entries that are an undetermined gauge, not a disagreement: intra-cluster pairs.

    >>> g = dict(n=2, cl=[[0, 1]])
    >>> mask_of(g).tolist()
    [[False, False], [False, False]]
    """
    keep = ~np.eye(g["n"], dtype=bool)
    for grp in g["cl"]:
        if len(grp) > 1:
            keep[np.ix_(grp, grp)] = False
    return keep


def maskd_of(g):
    """Diagonal entries that carry information: those outside every degenerate cluster.

    >>> maskd_of(dict(n=3, cl=[[0], [1, 2]])).tolist()
    [True, False, False]
    """
    keepd = np.ones(g["n"], dtype=bool)
    for grp in g["cl"]:
        if len(grp) > 1:
            keepd[list(grp)] = False
    return keepd


def resid(theta, g, renat, keep, Cref):
    """[off-diagonal M ; 1 - diag M] with M = Cref^t S C(w) sign-canonicalised, w = exp(theta).

    The column sign of C(w) is eigh's arbitrary gauge, not an answer: flipping it flips a whole column
    of M, and the raw off-diagonal residual therefore JUMPED by 0.75 for a 1e-10 step along a direction
    the cascade is exactly invariant under, which left least_squares differentiating noise (it moved
    zero steps, and the recovery control failed at 1.2e+00 on lif).  Each column is canonicalised
    against the reference by the sign of its own diagonal element - near any candidate solution that
    element is +-1, nowhere near the kink at zero.  Members of a degenerate occupancy cluster carry no
    diagonal residual at all: their intra-cluster rotation is undetermined, and with the cluster's
    off-diagonal block already masked out, M's orthogonality makes that block redundant rather than
    informative.  Neither change touches the acceptance rule; both are gauge quotients.
    """
    M = Cref.T @ g["S"] @ compose(np.exp(theta), g, renat)
    d = np.diag(M).copy()
    M = M * np.sign(np.where(d == 0.0, 1.0, d))[None, :]
    return np.concatenate([M[keep], 1.0 - np.diag(M)[g["keepd"]]])


def verdict(C, Cref, g):
    """The one verdict number and its per-class decomposition, from the lane's own sine instrument."""
    r = sines(C, Cref, g["S"], g["lbl"], g["cl"])
    r["worst"] = max(r["Cor"], r["Val"], r["Ryd"])
    return r


def jitter_floor(mol, d, g, renat):
    """The verdict metric's own floor: re-quantise every printed input, DRAWS draws, SPREAD x spread."""
    rng = np.random.default_rng(SEED)
    vals = []
    for _ in range(DRAWS):
        n, S, P = read_47(os.path.join(d, mol + ".47"))
        S = S + rng.uniform(-Q47, Q47, S.shape)
        P = P + rng.uniform(-Q47, Q47, P.shape)
        S, P = (S + S.T) / 2, (P + P.T) / 2
        C32 = read_lfn32(os.path.join(d, mol + ".32"), n, S)[0] + rng.uniform(-Q32, Q32, (n, n))
        C33 = read_lfn33(os.path.join(d, mol + ".33"), n, S)[0] + rng.uniform(-Q32, Q32, (n, n))
        perm = g["perm"]
        h = dict(g)
        h["S"], h["SPS"] = S, S @ P @ S
        h["C32"], h["C33"] = C32[:, perm], C33[:, perm]
        w = norm_w(np.maximum(g["pre_occ"] + rng.uniform(-Q32, Q32, n), 0.0), g["cls"])
        vals.append(verdict(compose(w, h, renat), h["C33"], h)["worst"])
    return SPREAD * (max(vals) - min(vals)), vals


def prepare(mol, d):
    """load() plus what the fit needs on top: the permutation CONSTRUCTED and ASSERTED, w0, the mask.

    The component-major -> shell-major permutation is rebuilt here from the same self-identifying
    printed fields load() reads (the `lang` field of the .aopnao.nbo.json, ranked by the unit-53
    diagonal) and then ASSERTED against load()'s own permuted C32 - never inferred from an offset and
    never fitted.  That assertion is the guard against the seventh read-by-position defect in this
    lane, where already_step4() was handed an unpermuted C32 with permuted labels.
    """
    g = load(mol, d)
    n, S = g["n"], g["S"]
    naos = json.load(open(os.path.join(d, mol + ".aopnao.nbo.json")))["nao"]
    Dn = read_packed(os.path.join(d, mol + ".53"), "NAO density matrix:", n)
    perm, lbl, reranked = gennbo_labels(naos, np.diag(Dn))
    assert reranked == 0, "%s: gennbo's printed ranks are not descending occupancy" % mol
    assert [t[1:6] for t in lbl] == [t[1:6] for t in g["lbl"]], mol
    C32raw = read_lfn32(os.path.join(d, mol + ".32"), n, S)[0]
    assert np.array_equal(C32raw[:, perm], g["C32"]), "%s: permutation does not reproduce load()'s C32" % mol
    g["perm"] = perm
    g["w0"] = norm_w(np.maximum(g["pre_occ"], 0.0), g["cls"])
    g["keep"] = mask_of(g)
    g["keepd"] = maskd_of(g)
    g["ndeg"] = sum(len(x) for x in g["cl"] if len(x) > 1)
    return g


def fit_one(g, renat, Cref, x0, starts=1, max_nfev=200):
    """least_squares on theta, boxed to owso()'s own representable range; best of `starts` starts.

    DIFF_STEP is not scipy's default.  On the recovery target - which is exactly reachable, so every
    failure there is the solver's - the default 2-point step left lif at 2.076e-04 while diff_step=1e-4
    reached 5.960e-08, the sine floor, because the cascade's eigen-decompositions are not smooth at the
    1.5e-08 scale the default probes.  Sixteen (jac, diff_step, x_scale) combinations were measured on
    that target before this one was fixed, and it is then used unchanged for every molecule and arm.

    `starts` > 1 adds the flat vectors and a log-uniform draw as extra starting points and keeps the
    lowest residual.  That is multistart for a plateaued optimiser, not best-of-N acceptance: it can
    only ever make a MISS harder to claim, and the acceptance rule still runs on 8 of 8 molecules.
    """
    from scipy.optimize import least_squares
    lo, hi = np.log(WLO), np.log(WHI)
    rng = np.random.default_rng(SEED + 2)
    xs = [np.clip(x0, lo + 1e-09, hi - 1e-09)]
    for k in range(1, starts):
        xs.append(np.clip(np.log(np.full(g["n"], 0.1)) if k == 1
                          else np.log(norm_w(10.0 ** rng.uniform(-4, 0, g["n"]), g["cls"])),
                          lo + 1e-09, hi - 1e-09))
    t0, best = time.time(), None
    nfev = 0
    for x in xs:
        r = least_squares(resid, x, bounds=(lo, hi), args=(g, renat, g["keep"], Cref), method="trf",
                          jac="2-point", diff_step=DIFF_STEP, x_scale=1.0,
                          xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=max_nfev)
        nfev += int(r.nfev)
        if best is None or r.cost < best.cost:
            best = r
    w = norm_w(np.exp(best.x), g["cls"])
    pinned = int((best.x <= lo + 1e-06).sum())
    return dict(w=w, cost=float(best.cost), nfev=nfev, starts=len(xs), secs=time.time() - t0,
                pinned=pinned, status=int(best.status),
                res=float(np.linalg.norm(best.fun)), neq=int(best.fun.size))


def forms(w, g):
    """Which function of which occupancy does the fitted w look like?  Correlations, per class."""
    out = {}
    pre = np.maximum(g["pre_occ"], 0.0)
    cands = {"pre": pre, "pre^2": pre ** 2, "sqrt(pre)": np.sqrt(pre), "1/pre": 1.0 / np.maximum(pre, 1e-30)}
    for c in sorted(set(g["cls"])):
        idx = [i for i, k in enumerate(g["cls"]) if k == c]
        if len(idx) < 3:
            continue
        lw = np.log(np.maximum(w[idx], 1e-300))
        for name, v in cands.items():
            lv = np.log(np.maximum(v[idx], 1e-300))
            if lv.std() > 0 and lw.std() > 0:
                out["%d/%s" % (c, name)] = float(np.corrcoef(lw, lv)[0, 1])
    return out


def demo():
    """The three controls on synthetic data, where the right answer is known.

    Non-degenerate by construction and checked: Sp is not the identity, the composition moves when w
    moves, and the planted weights come back.

    >>> ok = demo()
    >>> ok["sp_off"] > 0.1                      # (3a) the orthogonalisation is not the identity
    True
    >>> ok["moves"] > 1e-03                     # (3b) the composition is not flat in w
    True
    >>> ok["recover_res"] < 1e-06               # (1) the planted w comes back
    True
    >>> ok["negative_res"] > 1e-02              # (2) an off-form target is refused
    True
    """
    rng = np.random.default_rng(7)
    na, nl = 3, 2
    lbl, cls = [], []
    for a in range(na):
        for l in range(nl):
            for sh in range(2):
                for m in range(2 * l + 1):
                    lbl.append((len(lbl), a, l, m, sh, 0 if (l == 0 and sh == 0) else (1 if sh == 0 else 2), 0.0))
                    cls.append(lbl[-1][5])
    n = len(lbl)
    A = rng.normal(size=(n, n))
    S = A @ A.T / n + np.eye(n)
    d = 1.0 / np.sqrt(np.diag(S))
    S = S * d[:, None] * d[None, :]
    B = rng.normal(size=(n, n))
    P = B @ B.T / n
    SPS = S @ P @ S
    # C32 = 1, so Sp = S and the pre-NAO overlap is NOT the identity.  Taking C32 = S^-1/2 here would
    # make every orthogonalisation in the cascade the identity and the controls would pass VACUOUSLY -
    # that is exactly the degenerate synthetic case this lane has already been caught by once.
    C32 = np.eye(n)
    lbl = [(t[0], t[1], t[2], t[3], t[4], t[5], float(np.abs(rng.normal()) + 0.1)) for t in lbl]
    g = dict(n=n, S=S, SPS=SPS, lbl=lbl, cls=cls, C32=C32, pre_occ=np.array([t[6] for t in lbl]))
    g["cl"] = clusters(lbl)
    g["keep"] = mask_of(g)
    g["keepd"] = maskd_of(g)
    g["w0"] = norm_w(g["pre_occ"], cls)
    out = {"n": n, "sp_off": float(np.abs(C32.T @ S @ C32 - np.eye(n)).max())}
    C0 = compose(g["w0"], g, False)
    out["moves"] = verdict(compose(norm_w(g["w0"] * rng.uniform(0.5, 2.0, n), cls), g, False), C0, g)["worst"]
    wp = norm_w(10.0 ** rng.uniform(-4, 0, n), cls)
    Cp = compose(wp, g, False)
    r = fit_one(g, False, Cp, np.log(np.full(n, 0.1)))
    out["recover_res"] = r["res"]
    out["recover_dlog"] = float(np.abs(np.log10(r["w"] / wp)).max())
    Q = np.linalg.qr(rng.normal(size=(n, n)))[0]
    Cneg = C32 @ sym_power(C32.T @ S @ C32, -0.5) @ Q
    rn = fit_one(g, False, Cneg, np.log(np.full(n, 0.1)), max_nfev=400)
    out["negative_res"] = rn["res"]
    return out


def main(argv):
    if not argv or argv[0] == "--demo":
        d = demo()
        print("DEMO / synthetic controls, n = %d" % d["n"])
        print("  (3a) |Sp-1|                    %.4e   (must be >> 0: not the identity)" % d["sp_off"])
        print("  (3b) composition moves with w  %.4e   (must be >> floor)" % d["moves"])
        print("  (1)  recovery residual         %.4e   max|dlog10 w| %.3e" % (d["recover_res"], d["recover_dlog"]))
        print("  (2)  negative residual         %.4e   (must NOT reach the floor)" % d["negative_res"])
        return 0
    d = argv[0]
    mols = argv[1].split(",") if len(argv) > 1 else list(MOLS)
    rows = []
    for mol in mols:
        sub = os.path.join(d, "ops_" + mol)
        g = prepare(mol, sub)
        row = dict(mol=mol, n=g["n"], ndeg=g["ndeg"], neq=int(g["keep"].sum() + g["keepd"].sum()),
                   sp_off=float(np.abs(g["Sp"] - np.eye(g["n"])).max()), arms={})
        rng = np.random.default_rng(SEED)
        C33 = g["C33"]
        # (3c) the instrument must go red on two wrong answers.
        row["WRONG_identity"] = verdict(g["C32"], C33, g)["worst"]
        row["WRONG_lowdin"] = verdict(g["C32"] @ sym_power(g["Sp"], -0.5), C33, g)["worst"]
        for name, renat in BACKBONES:
            C0 = compose(g["w0"], g, renat)
            base = verdict(C0, C33, g)
            fl, draws = jitter_floor(mol, sub, g, renat)
            # (3b) on real data: the composition moves when w moves.
            moves = verdict(compose(norm_w(g["w0"] * rng.uniform(0.5, 2.0, g["n"]), g["cls"]), g, renat),
                            C0, g)["worst"]
            x0 = np.log(np.maximum(g["w0"], WLO))
            f = fit_one(g, renat, C33, x0, starts=STARTS)
            vf = verdict(compose(f["w"], g, renat), C33, g)
            # (1) recovery, on this molecule's real S/SPS/labels.
            wp = norm_w(10.0 ** rng.uniform(-4, 0, g["n"]), g["cls"])
            rc = fit_one(g, renat, compose(wp, g, renat), np.log(np.full(g["n"], 0.1)))
            # (2) negative, an independent random orthogonal on top of a plain Loewdin.
            Q = np.linalg.qr(np.random.default_rng(SEED + 1).normal(size=(g["n"], g["n"])))[0]
            Cneg = g["C32"] @ sym_power(g["Sp"], -0.5) @ Q
            rn = fit_one(g, renat, Cneg, x0, starts=STARTS)
            neg = dict(res=rn["res"], worst=verdict(compose(rn["w"], g, renat), Cneg, g)["worst"])
            call = ("FIT" if vf["worst"] <= FIT_K * fl else
                    ("MISS" if vf["worst"] >= MISS_K * fl else "INCONCLUSIVE"))
            row["arms"][name] = dict(
                floor=fl, draws=draws, base=base, fit=f, vfit=vf, call=call, moves=moves,
                recover=dict(res=rc["res"], dlog=float(np.abs(np.log10(rc["w"] / wp)).max()),
                             worst=verdict(compose(rc["w"], g, renat), compose(wp, g, renat), g)["worst"]),
                negative=neg,
                forms=forms(f["w"], g))
            print("%-9s %-6s n=%-4d eq=%-6d unk=%-4d floor=%.2e base=%.4e fit=%.4e %-12s "
                  "pin=%-3d rec=%.1e neg=%.4e" % (
                      mol, name, g["n"], row["neq"], g["n"] - len(set(g["cls"])),
                      fl, base["worst"], vf["worst"], call, f["pinned"],
                      rc["res"], neg["worst"]))
            sys.stdout.flush()
        rows.append(row)
    calls = {}
    for r in rows:
        for a, v in r["arms"].items():
            calls.setdefault(a, []).append(v["call"])
    print()
    for a, cs in calls.items():
        print("%-6s FIT %d/%d  MISS %d/%d  INCONCLUSIVE %d/%d" % (
            a, cs.count("FIT"), len(cs), cs.count("MISS"), len(cs), cs.count("INCONCLUSIVE"), len(cs)))
    tag = "" if list(mols) == list(MOLS) else "_" + "-".join(mols)
    out = os.path.join(d, "fitw_result%s.json" % tag)
    json.dump(rows, open(out, "w"), indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
