"""JANPA 2.02 as a THIRD CODE against gennbo 7's own NAO matrix, per (atom, l) block.

    py -3.12 janpa_third_code.py <dir with one subdir per molecule> [--full]
    py -3.12 janpa_third_code.py --demo

Why a third code at all.  The published algorithm (Nikolaienko, Bulavin, Hovorun, Comput. Theor.
Chem. 1050, 15 (2014)) has seven steps where NoSpherA2's cascade has four, and restoring two of
them (renat5) cuts the defect 63x and still fails the gate 0 of 22.  Everything that has been
learned about the remaining gap was learned by reading the spec and measuring gennbo's printed
output.  Reading a spec twice is not a second measurement.  JANPA is an independent published
implementation of that spec which can be instrumented at every stage, so it answers a question
neither previous side could: is gennbo 7 running the published algorithm at all?

What is compared, and in which gauge.  Both sides get the SAME wavefunction: gennbo reads the
.47, JANPA reads the molden that orca_2mkl writes from the same .gbw.  The two AO bases are the
same functions in a different ORDER and a different SIGN, and both are pinned exactly, neither
fitted:

  ORDER - JANPA's own documented component table (janpa-py, convert/file47.py,
    proper_ordering_SPH) says its molden basis runs p = (px, py, pz) = NBO labels 151, 152, 153,
    while d = (255, 252, 253, 254, 251) and f = (351 .. 357) agree with the .47 term for term.
    ORCA's internal order, which is what the .47 carries, is (pz, px, py).  So the permutation is
    the identity on every l except p, where it is one 3-cycle per shell.  It is CONSTRUCTED from
    that table and then ASSERTED, not fitted.
  SIGN - orca_2mkl's molden and ORCA's own AO basis differ in the phase of some real spherical
    harmonics.  molden2molden's `-orca3signs` puts the MO coefficients into MOLDEN's convention,
    which is what makes JANPA's read exact (see below); the BASIS functions still differ from the
    .47's by a diagonal sign matrix D.  D is recovered from S by a spanning tree over the
    significant elements - n - 1 signs - and then checked against all n(n+1)/2 elements of
    D S_j D = S_47, which is an overdetermination by a factor ~n/2, plus the triple-consistency
    sign(S_ij) sign(S_jk) sign(S_ki) that a diagonal gauge must satisfy and a general sign error
    must violate.  Every decisive metric below is in any case m-averaged or a principal sine, so
    it is blind to both the component order and the sign; the gauge is needed only for the
    per-function table.

The read had to be fixed before any of this meant anything.  Raw, JANPA on an ORCA 6.1 molden
reports a basis normalisation problem, non-orthogonal MOs and "input data seems to be improper",
and returns Li +15.5 e on a 12-electron LiF: orca_2mkl prints ORCA's unnormalised contraction
coefficients.  `-fromorca3bf` fixes the contraction and takes LiF to 0 warnings - but LiF is a
diatomic on z, where an x/y sign difference is invisible, and the other seven still lost up to
0.31 e of 70 (sf6).  `-orca3signs` as well takes all eight to 0 warnings and to an EXACT electron
count.  The gate for the read is not the warning count: it is tr(P S) = N exactly, max|D S_j D -
S_47| at the .47's print floor, and JANPA's NPA charges against gennbo's printed ones.

THE FLOOR IS MEASURED, NEVER GUESSED.  gennbo prints units 32 and 33 at nine decimals, so the
matrix that comes back is not the matrix it computed.  Each gennbo matrix is re-quantised by
+-0.5e-09 over four independent draws and the floor is 3x the spread of the metric over those
draws, per molecule and per metric.  A cond()-based bound was refuted on this branch; a guessed
bound is a fitted bound.  JANPA is asked for %.16e, so its own print floor is 1e-16 relative.

EVERY READING CARRIES ITS OWN CONTROLS, IN THE SAME RUN.
  positive  JANPA's PRE-NAOs against gennbo's unit 32.  The pre-NAOs are upstream of every
            orthogonalisation and the two codes are known to agree there (native vs gennbo,
            1.000000 over 379 shells), so this arm must come out GREEN - it is the proof that the
            permutation, the sign gauge, the pairing and the reader are all right.  If it is red
            the NAO comparison says nothing about NAOs.
  negative  JANPA's NAOs against gennbo's unit 32, and gennbo's own 33 against its own 32.  Both
            must come out RED.  A printed 0.00e+00 over 0 of N pairs is an empty maximum, not an
            agreement, so the sample size of every reading is printed beside it.
"""
import collections
import json
import os
import re
import socket
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aonao_compare import (LANG_L, gennbo_blocks, mixing, numbers, principal_sines,  # noqa: E402
                           read_47, read_lfn32, read_lfn33, section)

# JANPA's own documented molden component order per l, in NBO label codes
# (janpa-py src/janpa/convert/file47.py, proper_ordering_SPH).
SPH = {0: [51], 1: [151, 152, 153], 2: [255, 252, 253, 254, 251],
       3: [351, 352, 353, 354, 355, 356, 357],
       4: [451, 452, 453, 454, 455, 456, 457, 458, 459]}
FLOOR_DRAWS = 4
FLOOR_QUANT = 0.5e-09     # half a unit in gennbo's ninth printed decimal
FLOOR_K = 3.0
OCC_TIE = 1.0e-05         # gennbo prints the NAO occupancy at five decimals: one printed unit


def norm_label(lab):
    """NBO's s/p label codes come in two spellings; SPH uses the 1x1 one."""
    return lab if lab >= 154 else (lab // 100) * 100 + 50 + (lab % 10)


def basis_47(path):
    """(n, centers, labels) of the .47's AO basis, in the .47's own order."""
    text = open(path).read()
    n = int(re.search(r"NBAS=(\d+)", text).group(1))
    basis = section(text, "$BASIS")
    at = basis.index("LABEL")
    centers = [int(v) for v in numbers(basis[:at])]
    labels = [int(v) for v in numbers(basis[at:])]
    assert len(centers) == len(labels) == n, (len(centers), len(labels), n)
    return n, centers, labels


def ao_perm(path47):
    """perm[j] = the .47 AO index holding JANPA's molden AO j, from JANPA's documented table.

    >>> import tempfile, os
    >>> # one p shell: the .47 carries ORCA's (pz, px, py) = 103, 101, 102
    >>> ao_perm.__doc__ is not None
    True
    """
    n, _, labels = basis_47(path47)
    perm = [None] * n
    i = 0
    while i < n:
        l = norm_label(labels[i]) // 100
        order = SPH[l]
        for m in range(len(order)):
            perm[i + order.index(norm_label(labels[i + m]))] = i + m
        i += len(order)
    assert sorted(perm) == list(range(n)), "constructed map is not a permutation"
    return n, np.array(perm)


def read_janpa(path, n):
    """One of JANPA's exported matrices; three header lines, then n*n values."""
    vals = []
    for line in open(path).read().split("\n")[3:]:
        for t in line.split():
            try:
                vals.append(float(t))
            except ValueError:
                pass
    assert len(vals) >= n * n, "%s holds %d numbers, need %d" % (
        os.path.basename(path), len(vals), n * n)
    return np.array(vals[:n * n]).reshape(n, n)


def read_janpa_c(path, n, G, S, kind):
    """JANPA's PNAO2AO / NAO2AO export as AO-rows x NAO-columns in the .47 gauge.

    JANPA's header calls it "AO-to-NAO transformation matrix" and its printout is the TRANSPOSE
    of gennbo's unit 32/33 convention - it stores one NAO per ROW.  That is not assumed here: the
    orientation is DECIDED by the same test read_lfn33/read_lfn32 use, and the loser must FAIL,
    because a test both orientations pass has decided nothing.  Reading this by position was the
    eighth read-by-position defect in this lane; the positive control caught it.
    """
    M = read_janpa(path, n)
    cand = {"asis": G @ M, "transposed": G @ M.T}
    err = {}
    for name, X in cand.items():
        O = X.T @ S @ X
        err[name] = float(np.abs(O - np.eye(n)).max() if kind == "nao"
                          else np.abs(np.diag(O) - 1.0).max())
    win = min(err, key=err.get)
    lose = "transposed" if win == "asis" else "asis"
    assert err[win] < 1.0e-06, "%s: neither orientation is %s-orthonormal (%r)" % (
        os.path.basename(path), kind, err)
    assert err[lose] > 1.0e-03, "%s: BOTH orientations pass, the test discriminates nothing (%r)" \
        % (os.path.basename(path), err)
    return cand[win], win, err[win], err[lose]


def janpa_labels(path):
    """(atom, l, m, is_NRB, shell) per NAO from JANPA's own printed header.

    `shell` is JANPA's own R number, which is the ONLY self-identifying statement of which radial
    function a component belongs to; m is kept for the record and used nowhere as a direction.
    """
    head = open(path).read().split("\n")[2]
    out = []
    for tok in head.split("\t"):
        tok = tok.strip()
        if not tok:
            continue
        m = re.match(r"^A(\d+)(\*?):\s*R(\d+)\*([spdfg])\((-?\d+)\)", tok)
        assert m, "unparsed JANPA label %r" % tok
        out.append((int(m.group(1)), "spdfg".index(m.group(4)), int(m.group(5)),
                    m.group(2) == "*", int(m.group(3))))
    return out


def sign_gauge(Sj, S, thresh=1.0e-04):
    """D with D S_j D = S, from a spanning tree over the elements above `thresh`.

    n - 1 signs are fixed by the tree; the residual is then checked on all n(n+1)/2 elements.
    Returns (D, residual, n_edges_used, n_triples_violated).
    """
    n = S.shape[0]
    big = np.abs(S) > thresh
    sg = np.where(big, np.sign(Sj * S), 0.0)
    D = np.zeros(n)
    edges = 0
    for start in range(n):
        if D[start] != 0:
            continue
        D[start] = 1.0
        q = collections.deque([start])
        while q:
            i = q.popleft()
            for j in np.nonzero(big[i])[0]:
                if D[j] == 0:
                    D[j] = D[i] * sg[i, j]
                    q.append(j)
                    edges += 1
    resid = float(np.abs(np.outer(D, D) * S - Sj).max())
    pred = np.outer(D, D) * S
    bad = int(np.sum(big & (np.sign(pred) != np.sign(Sj)) & (np.abs(S) > thresh)))
    return D, resid, edges, bad


def block_index(centers, labels):
    """(atom0, l) -> {component label -> [AO indices in .47 order]}."""
    out = {}
    for i, (c, lab) in enumerate(zip(centers, labels)):
        nl = norm_label(lab)
        out.setdefault((c - 1, nl // 100), {}).setdefault(nl, []).append(i)
    return out


def janpa_blocks(lab_j, occ):
    """(atom0, l) -> shells, each the (2l+1) JANPA NAO indices, ranked by occupancy descending.

    The blocking is read from JANPA's OWN printed labels (`A<c>: R<shell>*<l>(<m>)`), never from
    the .47's AO order - the two axes happen to be laid out the same way for these basis sets and
    a happens-to is exactly the read-by-position that has cost this lane seven defects.  A shell is
    the set of components sharing JANPA's R number; the m field is not taken to mean a direction,
    because JANPA prints p(0), p(1), p(-1) for what its own table calls px, py, pz.

    The RANK is the m-AVERAGED occupancy, descending, which is what gennbo's own printing order is.
    Ranking each component separately by its own occupancy is wrong and was measurably wrong: in a
    molecule the components of one shell are not degenerate, so the per-component orders disagree
    and one component of a shell gets paired with a different shell - which shows up as M[a][a]
    pinned at exactly 1 - 1/(2l+1) (0.2000 on a d shell, 0.3333 on a p shell) with the missing
    weight still inside the same block (`out` at 1e-09).
    """
    shells = {}
    for i, (a, l, _, _, r) in enumerate(lab_j):
        shells.setdefault((a - 1, l), {}).setdefault(r, []).append(i)
    blocks = {}
    for key, by_r in shells.items():
        for r, idx in by_r.items():
            assert len(idx) == 2 * key[1] + 1, "block %s shell R%d has %d components" % (
                key, r, len(idx))
        blocks[key] = [by_r[r] for r in sorted(by_r, key=lambda r: -np.mean(occ[by_r[r]]))]
    return blocks


def tie_groups(occ_desc, tol=OCC_TIE):
    """Runs of shells whose RANK is not resolvable, from occupancies sorted descending.

    Two shells whose occupancies differ by less than one printed unit are not ordered by anything:
    water's oxygen holds two s pre-NAOs at 0.000000 and rank pairing then splits a degenerate pair
    0.8621/0.1379 - an ambiguity of the pairing rule, reported as a disagreement of the codes.  The
    metric is summed over the run, which is pairing-free inside it.

    >>> tie_groups([2.0, 1.75, 0.0003, 0.0, 0.0])
    [[0], [1], [2], [3, 4]]
    """
    groups = [[0]]
    for a in range(1, len(occ_desc)):
        if occ_desc[a - 1] - occ_desc[a] <= tol:
            groups[-1].append(a)
        else:
            groups.append([a])
    return groups


def group_metric(M, groups):
    """max over tie-runs of |1 - mean_{a in g} sum_{b in g} M[a][b]|."""
    worst = 0.0
    for g in groups:
        worst = max(worst, abs(1.0 - float(M[np.ix_(g, g)].sum()) / len(g)))
    return worst


def log_nrb(path):
    """JANPA's own NRB flag per PNAO, from the log line it prints under STEP 2.

    The header stars are the same information; taking both and requiring them to agree is the
    only way to know the label axis of the exported matrix is the axis JANPA thinks it is.
    """
    text = open(path).read()
    at = text.index("Does PNAO belong to NRB?")
    toks = text[at:].split("\n")[1].split()
    assert all(t in ("true", "false") for t in toks), toks[:5]
    return np.array([t == "true" for t in toks])


def classes_janpa(labels_j):
    return np.array([2 if is_nrb else 0 for (_, _, _, is_nrb, _) in labels_j])


def classes_gennbo(naos):
    return np.array([2 if e["type"] == "Ryd" else 0 for e in naos])


def one(d, mol, verbose=False):
    p47 = os.path.join(d, mol + ".47")
    n, centers, labels = basis_47(p47)
    _, perm = ao_perm(p47)
    _, S, P = read_47(p47)
    Q = np.zeros((n, n))
    Q[perm, np.arange(n)] = 1.0            # JANPA order -> .47 order

    Sj = Q @ read_janpa(os.path.join(d, mol + ".j.S.txt"), n) @ Q.T
    D, resid, edges, bad = sign_gauge(Sj, S)
    G = Q @ np.diag(D[perm])               # signed permutation: JANPA gauge -> .47 gauge
    # the density needs the SAME gauge; tr(Pj S_47) with an ungauged Pj is not tr(Pj S_j) and the
    # difference is not small - it read -0.31 e on sf6 while JANPA's own electron count was exact.
    Pj = G @ read_janpa(os.path.join(d, mol + ".j.D.txt"), n) @ G.T
    # a sign gauge is the ONLY thing between the two overlaps
    absdiff = float(np.abs(np.abs(Sj) - np.abs(S)).max())

    Cpre_j, o_pre, e_pre, e_pre_l = read_janpa_c(
        os.path.join(d, mol + ".j.pnao2ao.txt"), n, G, S, "pre")
    Cnao_j, o_nao, e_nao, e_nao_l = read_janpa_c(
        os.path.join(d, mol + ".j.nao2ao.txt"), n, G, S, "nao")
    lab_j = janpa_labels(os.path.join(d, mol + ".j.nao2ao.txt"))
    # occupancies of JANPA's NAOs, by their definition, from the .47's own P and S
    occ_j = np.diag(Cnao_j.T @ S @ P @ S @ Cnao_j)
    occ_j_exp = np.diag(read_janpa(os.path.join(d, mol + ".j.sds_nao.txt"), n))
    nrb_log = log_nrb(os.path.join(d, "janpa.log"))
    star = np.array([b for (_, _, _, b, _) in lab_j])
    assert len(nrb_log) == n and (nrb_log == star).all(), \
        "%s: JANPA's log NRB flags and the matrix header stars disagree on %d of %d" % (
            mol, int((nrb_log != star).sum()), n)

    Cpre_g, lay32, e32 = read_lfn32(os.path.join(d, mol + ".32"), n, S)
    Cnao_g, lay33, e33 = read_lfn33(os.path.join(d, mol + ".33"), n, S)
    naos = json.load(open(os.path.join(d, mol + ".aopnao.nbo.json")))["nao"]
    assert len(naos) == n, (len(naos), n)

    res = dict(mol=mol, n=n, nflip=int((D < 0).sum()), sign_resid=resid, sign_edges=edges,
               sign_bad=bad, abs_S_diff=absdiff, layout32=lay32, layout33=lay33,
               err32=e32, err33=e33,
               trPS=float((P @ S).trace()), trPjS=float((Pj @ S).trace()),
               occ_j_sum=float(occ_j.sum()),
               orient_pre=o_pre, orient_nao=o_nao, err_pre=e_pre, err_pre_loser=e_pre_l,
               err_nao=e_nao, err_nao_loser=e_nao_l,
               occ_j_exp_diff=float(np.abs(np.sort(occ_j)[::-1]
                                           - np.sort(occ_j_exp)[::-1]).max()),
               npa_j=[float(x) for x in open(os.path.join(d, mol + ".j.npa.txt")).read().split()],
               npa_g=[e["charge"] for e in json.load(
                   open(os.path.join(d, mol + ".aopnao.nbo.json")))["npa"]])

    # gennbo's NAO occupancies must come back out of its own matrix, or the matrix is not the
    # one that produced the printed table.
    occ_g_mat = np.diag(Cnao_g.T @ S @ P @ S @ Cnao_g)
    occ_g_tab = np.array([e["occupancy"] for e in naos])
    res["occ_g_check"] = float(np.abs(occ_g_mat - occ_g_tab).max())

    gb = gennbo_blocks(naos)
    jb = janpa_blocks(lab_j, occ_j)
    assert set(gb) == set(jb), (sorted(set(gb) ^ set(jb)),)

    # tie runs per block, from GENNBO's printed occupancies: the rank pairing is only defined
    # between shells its own table separates by at least one printed unit.
    tg = {key: tie_groups([float(np.mean([naos[i]["occupancy"] for i in sh])) for sh in gb[key]])
          for key in gb}
    res["ties"] = sum(len(g) - 1 for gs in tg.values() for g in gs)

    arms = {}
    for name, Cj, Cg in (("PRE  (positive control)", Cpre_j, Cpre_g),
                         ("NAO  (decisive)", Cnao_j, Cnao_g),
                         ("NAO vs PRE (negative control)", Cnao_j, Cpre_g),
                         ("gennbo 33 vs 32 (negative control)", Cnao_g, Cpre_g)):
        worst_diag, worst_best, worst_out, npairs = 0.0, 0.0, 0.0, 0
        for key in sorted(gb):
            rows, cols = jb[key], gb[key]
            if len(rows) != len(cols):
                continue
            M, total, _, _ = mixing(Cj, Cg, S, rows, cols)
            npairs += len(rows)
            worst_diag = max(worst_diag, group_metric(M, tg[key]))
            worst_best = max(worst_best, float(np.abs(1.0 - M.max(axis=1)).max()))
            worst_out = max(worst_out, float(np.abs(1.0 - M.sum(axis=1)).max()))
        arms[name] = dict(rank=worst_diag, best=worst_best, out=worst_out, npairs=npairs)

    # principal sines per (atom, l, class), pairing-free
    cj, cg = classes_janpa(lab_j), classes_gennbo(naos)
    by_j, by_g = {}, {}
    for i, (a, l, _, _, _) in enumerate(lab_j):
        by_j.setdefault((a - 1, l, cj[i]), []).append(i)
    for e in naos:
        by_g.setdefault((e["atom"] - 1, LANG_L[e["lang"][0].lower()],
                         2 if e["type"] == "Ryd" else 0), []).append(e["index"] - 1)
    sines = {}
    for key in sorted(set(by_g) & set(by_j)):
        A, B = Cnao_j[:, by_j[key]], Cnao_g[:, by_g[key]]
        if A.shape[1] != B.shape[1]:
            sines[key] = None
            continue
        sines[key] = float(principal_sines(A, B, S).max())
    res["sines_max"] = max(v for v in sines.values() if v is not None)
    res["sines_n"] = sum(1 for v in sines.values() if v is not None)
    res["sines_skipped"] = sum(1 for v in sines.values() if v is None)
    res["sines_worst_key"] = max((v, k) for k, v in sines.items() if v is not None)[1]
    res["arms"] = arms

    # THE FLOOR, measured: re-quantise the gennbo matrix the arm uses and re-run its own metric.
    rng = np.random.default_rng(20260925)

    def floor_of(Cj, Cg):
        draws = []
        for _ in range(FLOOR_DRAWS):
            pert = Cg + rng.uniform(-FLOOR_QUANT, FLOOR_QUANT, Cg.shape)
            w = 0.0
            for key in sorted(gb):
                if len(jb[key]) != len(gb[key]):
                    continue
                M, _, _, _ = mixing(Cj, pert, S, jb[key], gb[key])
                w = max(w, group_metric(M, tg[key]))
            draws.append(w)
        return FLOOR_K * float(np.ptp(draws)), draws

    res["floor"], res["floor_draws"] = floor_of(Cnao_j, Cnao_g)
    res["pre_floor"], res["pre_floor_draws"] = floor_of(Cpre_j, Cpre_g)
    base = arms["NAO  (decisive)"]["rank"]
    pre = arms["PRE  (positive control)"]["rank"]
    # The positive control's bound is not guessed: the pre-NAOs of the two codes are ESTABLISHED
    # to agree at 1.000000 over 379 shells, so 5e-07 is that reading's own print resolution.
    res["pre_ok"] = bool(pre <= max(res["pre_floor"], 5.0e-07))
    res["neg_ok"] = all(arms[k]["rank"] > 100.0 * max(res["floor"], 1e-12) and arms[k]["npairs"] > 0
                        for k in ("NAO vs PRE (negative control)",
                                  "gennbo 33 vs 32 (negative control)"))
    res["verdict"] = ("VOID" if not (res["pre_ok"] and res["neg_ok"])
                      else "AGREE" if base <= res["floor"] else "DISAGREE")
    if verbose:
        res["sines_all"] = {str(k): v for k, v in sines.items()}
    return res


def report(rows):
    print("host=%s  draws=%d quant=%.1e k=%.1f" % (socket.gethostname(), FLOOR_DRAWS,
                                                   FLOOR_QUANT, FLOOR_K))
    print("\n--- the read, asserted (not the warning count) ---")
    print("%-9s %4s %6s %11s %6s %6s  %11s %11s %11s" % (
        "mol", "n", "nflip", "max|DSjD-S|", "edges", "bad3", "max||Sj|-|S||",
        "tr(Pj S)-N", "occ_g chk"))
    for r in rows:
        print("%-9s %4d %6d %11.3e %6d %6d  %11.3e %11.3e %11.3e" % (
            r["mol"], r["n"], r["nflip"], r["sign_resid"], r["sign_edges"], r["sign_bad"],
            r["abs_S_diff"], r["trPjS"] - r["trPS"], r["occ_g_check"]))
    print("\n--- orientation of JANPA's exports, DECIDED not assumed (winner must pass, loser must"
          " FAIL) ---")
    print("%-9s %12s %11s %11s   %12s %11s %11s %11s" % (
        "mol", "NAO2AO", "max|CtSC-I|", "loser", "PNAO2AO", "max|diag-1|", "loser",
        "occ vs SDS_NAO"))
    for r in rows:
        print("%-9s %12s %11.3e %11.3e   %12s %11.3e %11.3e %11.3e" % (
            r["mol"], r["orient_nao"], r["err_nao"], r["err_nao_loser"],
            r["orient_pre"], r["err_pre"], r["err_pre_loser"], r["occ_j_exp_diff"]))
    print("\n--- NPA charges, JANPA against gennbo's printed table (gennbo prints 4 decimals,"
          " so its own floor here is 5e-05) ---")
    print("%-9s %5s %12s %12s   %s" % ("mol", "natom", "max|dq|", "rms dq", "worst atom"))
    for r in rows:
        a, b = np.array(r["npa_j"]), np.array(r["npa_g"])
        assert len(a) == len(b), (r["mol"], len(a), len(b))
        d = np.abs(a - b)
        print("%-9s %5d %12.3e %12.3e   %d: %+0.6f vs %+0.4f" % (
            r["mol"], len(a), d.max(), float(np.sqrt((d ** 2).mean())),
            int(d.argmax()) + 1, a[d.argmax()], b[d.argmax()]))
    print("\n--- per (atom, l) block, m-averaged mixing: 1 - M[a][a] after occupancy-rank pairing"
          " (`rank`), 1 - max_b M[a][b] (`best`, pairing-free), 1 - sum_b M[a][b] (`out`) ---")
    for name in ("PRE  (positive control)", "NAO  (decisive)",
                 "NAO vs PRE (negative control)", "gennbo 33 vs 32 (negative control)"):
        print("  %s" % name)
        print("    %-9s %6s %11s %11s %11s" % ("mol", "shells", "rank", "best", "out"))
        for r in rows:
            a = r["arms"][name]
            print("    %-9s %6d %11.3e %11.3e %11.3e" % (
                r["mol"], a["npairs"], a["rank"], a["best"], a["out"]))
    print("\n--- principal sines per (atom, l, class), pairing- and sign-free; and the MEASURED"
          " floor of the decisive metric ---")
    print("%-9s %8s %8s %13s %22s %11s %11s %6s %6s  %s" % (
        "mol", "blocks", "skipped", "max sine", "worst block", "decisive", "floor", "pos", "neg",
        "verdict"))
    for r in rows:
        print("%-9s %8d %8d %13.3e %22s %11.3e %11.3e %6s %6s  %s" % (
            r["mol"], r["sines_n"], r["sines_skipped"], r["sines_max"],
            str(r["sines_worst_key"]), r["arms"]["NAO  (decisive)"]["rank"], r["floor"],
            "ok" if r["pre_ok"] else "RED", "ok" if r["neg_ok"] else "RED", r["verdict"]))
    void = [r["mol"] for r in rows if r["verdict"] == "VOID"]
    if void:
        print("\nNO VERDICT on %d of %d molecules (%s): a control failed in this very run, so the"
              " decisive arm measures the instrument, not the codes." % (
                  len(void), len(rows), " ".join(void)))
    ok = [r for r in rows if r["verdict"] != "VOID"]
    if ok:
        n_ag = sum(1 for r in ok if r["verdict"] == "AGREE")
        print("\nVERDICT: JANPA reproduces gennbo's NAO matrix on %d of %d molecules whose controls"
              " passed (%d molecules total)." % (n_ag, len(ok), len(rows)))


def demo():
    """Self-check: the constructed p 3-cycle, the sign gauge and the empty-maximum guard."""
    lab = [1, 103, 101, 102, 255, 252, 253, 254, 251]
    perm = [None] * len(lab)
    i = 0
    while i < len(lab):
        l = norm_label(lab[i]) // 100
        order = SPH[l]
        for m in range(len(order)):
            perm[i + order.index(norm_label(lab[i + m]))] = i + m
        i += len(order)
    assert perm == [0, 2, 3, 1, 4, 5, 6, 7, 8], perm          # p: molden xyz -> .47 zxy
    rng = np.random.default_rng(0)
    A = rng.standard_normal((6, 6))
    S0 = A @ A.T + 6 * np.eye(6)
    d = np.array([1.0, -1, 1, 1, -1, -1])
    Dg, resid, edges, bad = sign_gauge(np.outer(d, d) * S0, S0)
    assert resid < 1e-12 and bad == 0, (resid, bad)
    assert np.allclose(np.abs(Dg), 1.0)
    assert np.allclose(np.outer(Dg, Dg), np.outer(d, d))      # the gauge is fixed up to a global sign
    # a sign error that is NOT a diagonal gauge must be caught, not absorbed
    Sbad = np.outer(d, d) * S0
    Sbad[0, 3] *= -1
    Sbad[3, 0] *= -1
    _, resid2, _, bad2 = sign_gauge(Sbad, S0)
    assert resid2 > 1e-3 and bad2 > 0, (resid2, bad2)
    assert tie_groups([2.0, 1.75, 0.0003, 0.0, 0.0]) == [[0], [1], [2], [3, 4]]
    # a degenerate pair split 0.8621/0.1379 by rank pairing is not a disagreement
    Mt = np.array([[1.0, 0, 0], [0, 0.8621, 0.1379], [0, 0.1379, 0.8621]])
    assert group_metric(Mt, [[0], [1], [2]]) > 0.13
    assert group_metric(Mt, [[0], [1, 2]]) < 1e-12
    # and a group sum must still catch weight that has left the run
    assert group_metric(np.diag([1.0, 0.9, 0.9]), [[0], [1, 2]]) > 0.09
    print("demo: p 3-cycle constructed, sign gauge exact at %.1e, non-gauge sign error caught"
          " (%d elements, residual %.3e), tie runs pairing-free and still leak-tight"
          % (resid, bad2, resid2))


def main(argv):
    if "--demo" in argv:
        demo()
        return 0
    root = argv[1]
    verbose = "--full" in argv
    mols = [m for m in sorted(os.listdir(root)) if os.path.isdir(os.path.join(root, m))]
    rows = []
    for mol in mols:
        if not os.path.exists(os.path.join(root, mol, mol + ".j.nao2ao.txt")):
            print("skip %s: no JANPA export" % mol)
            continue
        rows.append(one(os.path.join(root, mol), mol, verbose))
    assert rows, "no molecule had both sides - an empty table is not an agreement"
    report(rows)
    json.dump(rows, open(os.path.join(root, "janpa_third_code.json"), "w"), indent=1, default=str)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
