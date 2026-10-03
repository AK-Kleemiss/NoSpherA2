"""Does step 4 re-mix the core, and is that why native's cores are right only to 1e-05?

    py -3.12 step4_core_block.py --demo
    py -3.12 step4_core_block.py <dir of per-molecule subdirs with .47 .33 .naoc.txt .naocpre.txt>

WHY THIS EXISTS, and why it needs no rebuild.  Steps 3 and 4 of build_naos are a few lines of
linear algebra on quantities that are already dumped: Cpre and the pre-NAO occupancies from
NAO_DUMP_CPRE, the class and (atom, l, m, shell) labels with them, and S and P from the SAME .47
both codes read.  So the arm can be run here, in numpy, against native's own output as the control.
A replica is worthless unless it IS the thing it replicates, so gate A below reproduces native's
shipped NAO matrix before any arm is read; if it does not, this script says NO OUTCOME and stops.

WHAT THE CODE READING SAID (Src/core/nao.cpp, steps 3 and 4).  The coordinator's narrowed question
was whether gennbo's core priority is a restriction applied BEFORE step 4 rather than a class term
inside it.  Native ALREADY applies it as a restriction before step 4: step 3 walks the three classes
in decreasing priority and Schmidt-projects each out of everything above it, so the core keeps its
shape exactly at that point (nao.cpp:331-345).  The question therefore answers itself in the
opposite direction from the one it was asked in - and in doing so it MOVES the defect, because what
comes next undoes the restriction.  Step 4 keys its blocks on

    l_blocks[{ orbitals[i].atom, orbitals[i].l }]        nao.cpp:342-345

with no class in the key, and then re-diagonalises the m-averaged density inside each block.  An
atom's core 1s and its valence 2s share (atom, l) = (a, 0), so the very last step mixes the core
back into the valence it had just been protected from.  NAO_CLASS_SPLIT=1 adds the class to that
key for ALL THREE classes at once and is refuted - it separates valence from Rydberg, which is the
destructive half.  The core-only variant has never been run, and it is what the code reading points
at.

WHY THE MAGNITUDES FIT, stated before the arm so it cannot be fitted afterwards.  Three numbers
measured earlier in this lane look inconsistent and are not:

  the AONAO comparison found a mean m-averaged mixing defect of 0.00000 over 34 core shells, at
      five printed decimals, and that was reported as the cores being exactly right;
  --cascade found native's core partition deviating by 7.4e-05 .. 8.7e-04 where gennbo deviates by
      at most 1.3e-10, six decades apart;
  --space finds the core SUBSPACES differing by a principal sine of 7.9e-04 .. 6.3e-03.

A principal sine is sqrt(1 - sigma**2), so a coefficient-level deficit of d surfaces as about
sqrt(2 d): 6.3e-03 of sine is 2e-05 of coefficient, which prints as 0.00000 at five decimals.  The
three agree on one statement - native's cores are right to about 1e-05 and gennbo's to 1e-10 - and
that statement is exactly what a small unitary core/valence remix in step 4 would produce.  It also
says what this defect is NOT: 1e-05 in coefficients cannot be the 0.5 e leak, so this is a
correctness defect being localised, not the leak being found.

WHAT THE ARM MUST DO TO ITS OWN CREDIT.  Cauchy interlacing decides the direction in advance.  With
the core in the block the valence shell is assigned the SECOND eigenvalue of the (ns x ns)
m-averaged block; with the core removed it is assigned the FIRST eigenvalue of the submatrix, and
lambda_1(B') >= lambda_2(B) always.  So this arm can only move valence population UP and Rydberg
population DOWN - the correcting direction for the known leak - which is the opposite of
NAO_CLASS_SPLIT, whose own comment records that splitting can only increase a Rydberg excess.  That
makes the leak side of this arm a real test: it is allowed to fail by overshooting.

PRE-REGISTERED OUTCOMES, written before the run.

  GATE A, replica fidelity.  The baseline replica must reproduce native's shipped NAO matrix,
    measured sign-free as an ANGLE per column and per degenerate cluster: the replica's own worst
    angle must stay under 1e-03, and the undetermined (degenerate) directions must not disagree by
    more than ten times the determined ones.  Below that the replica is not native and nothing
    downstream means anything: NO OUTCOME, stop.

    GATE A AS FIRST WRITTEN WAS WRONG, and its failure is a finding rather than a bug in the
    replica.  Column by column it passed 6/8 at 1e-10 and failed water at 2.1e-03 and ammonia at
    3.5e-07 - in both cases on EXACTLY TWO columns, both empty Rydberg s shells of the heavy atom,
    and in both cases the m-averaged block those shells came from has a degenerate eigenvalue pair:
    water's gap is 2.1e-13 and ammonia's 6.1e-10.  Eigenvectors of a degenerate pair are not
    determined, so those columns have no value to agree about; the two SUBSPACES agree to 8.3e-06
    and 1.4e-05, which is native's own coefficient floor.  A gate that fails on an undefined
    quantity is measuring the arbitrary choice of two eigensolvers.
      The repair is not to loosen the tolerance.  It is to compare what is defined: cluster the
    shells of each block by their m-averaged eigenvalue, compare a singleton cluster column by
    column, and compare a degenerate cluster as a subspace.  Note what this deliberately does NOT
    do - it does not compare the whole BLOCK's span, which would be void: step 4 is a rotation
    inside the block, so the block span is invariant under it and the test could not fail.  A
    cluster is a proper subspace of the block, so the clustered test still fails if a shell is
    assigned to the wrong rank.
      This also puts a caveat on a number already reported in this lane.  The AONAO per-shell
    mixing defect (Core 0.00000 / Valence 0.00783 / Rydberg 0.37253, rank-paired) compares
    Rydberg shells COLUMN BY COLUMN, and the degenerate-cluster count below says how many of those
    shells had no defined coefficients on either side.  That part of the Rydberg figure is a gauge,
    not a defect.
  GATE B, discrimination.  The arm must not be a no-op: at least one core column must move by more
    than 1e-09 in the same overlap measure.  A gate that cannot fail is not a gate.

  OUTCOME K.  The arm collapses the core principal sine to the floor (i.e. to gennbo's own 1e-10
    territory) and leaves the valence and Rydberg sines within the core disagreement's own size.
    Then the mechanism is CONFIRMED and the core defect is localised to one line - the l_blocks key
    at nao.cpp:342-345 pooling core columns - with the fix being a core-only block.  Shipping still
    waits on the 22-molecule acceptance gate with pf5/so2/sf6 held flat.
  OUTCOME L.  The arm leaves the core sine where it was.  Then step 4 is not what perturbs the core,
    the mechanism is REFUTED, and the core defect stays located-and-open for Florian.
  OUTCOME M.  The arm collapses the core sine but moves valence or Rydberg by far more than the core
    disagreement can account for.  Then the two are coupled through the shared block; report both
    sides and let the acceptance gate decide, do not call it a fix.

NOT MEASURED HERE, on purpose: anything about timing.  The arm changes a block partition in a step
that is microseconds of a run.
"""

import io
import json
import os
import socket
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aonao_compare import (CLASS, ORDER, class_mass, principal_sines, read_47, read_lfn33,
                           read_naoc)

CORE_SEP = 10.0  # a class counts as moved only above this many times the floor


def sym_power(M, p, rel_floor=1e-10):
    """M**p for symmetric psd M, with the same relative eigenvalue floor nao.cpp uses.

    >>> M = np.diag([4.0, 1e-30])
    >>> float(sym_power(M, -0.5)[0, 0])
    0.5
    >>> float(sym_power(M, -0.5)[1, 1])   # dropped, not 1e+15
    0.0
    """
    w, V = np.linalg.eigh(M)
    cut = rel_floor * max(w.max(), 1e-300)
    f = np.where(w > cut, np.power(np.abs(w), p), 0.0 if p < 0.0 else np.power(np.clip(w, 0, None), p))
    return (V * f) @ V.T


def owso(S, weights):
    """T = W (W S W)^-1/2, with nao.cpp's weight clamp and weight-scaled eigenvalue floor.

    >>> float(np.abs(owso(np.eye(2), np.array([2.0, 1.0])) - np.eye(2)).max()) < 1e-12
    True
    """
    w = np.maximum(weights, 0.0)
    wmax = max(w.max(), 1e-300)
    w = np.maximum(w, 1e-6 * wmax)
    s = w.min() / wmax
    return w[:, None] * sym_power((w[:, None] * S) * w[None, :], -0.5, 1e-10 * s * s)


def pre_occupations(Cpre, SPS, labels):
    diag = np.diag(Cpre.T @ SPS @ Cpre)
    out = diag.copy()
    blocks = {}
    for i, t in enumerate(labels):
        blocks.setdefault((t[1], t[2], t[4]), []).append(i)
    for cols in blocks.values():
        out[cols] = np.mean(diag[cols])
    return out


def step3(Cpre, pre_occ, cls, S):
    """The three-class Schmidt + OWSO cascade, as nao.cpp:337-393."""
    C = Cpre.copy()
    done = np.zeros((C.shape[0], 0))
    for c in (0, 1, 2):
        cols = [i for i, k in enumerate(cls) if k == c]
        if not cols:
            continue
        B = C[:, cols].copy()
        w = np.maximum(pre_occ[cols], 0.0)
        if done.shape[1]:
            B -= done @ (done.T @ S @ B)
            B /= np.sqrt(np.maximum(np.einsum("ij,ij->j", B, S @ B), 1e-300))
        B = B @ owso(B.T @ S @ B, w)
        B = B @ sym_power(B.T @ S @ B, -0.5)
        C[:, cols] = B
        done = np.hstack([done, B])
    return C


def step4_blocks(lbl, core_own_block):
    """(atom, l) -> column lists, shell-major; the core split off into its own block or not.

    >>> lbl = [(0, 0, 0, 0, 0, 0, 2.0), (1, 0, 0, 0, 1, 1, 1.0), (2, 0, 0, 0, 2, 2, 0.0)]
    >>> step4_blocks(lbl, False)
    [[0, 1, 2]]
    >>> step4_blocks(lbl, True)
    [[0], [1, 2]]
    """
    keyed = {}
    for i, t in enumerate(lbl):
        keyed.setdefault((t[1], t[2]), []).append(i)
    out = []
    for key in sorted(keyed):
        cols, nm = keyed[key], 2 * key[1] + 1
        if not core_own_block:
            out.append(cols)
            continue
        core, rest = [], []
        for sh in range(len(cols) // nm):
            grp = cols[sh * nm:(sh + 1) * nm]
            (core if lbl[grp[0]][5] == 0 else rest).extend(grp)
        out.extend(b for b in (core, rest) if b)
    return out


def step4(C, SPS, lbl, core_own_block=False):
    """The m-averaged within-block re-diagonalisation, as nao.cpp:417-481."""
    Porb = C.T @ SPS @ C
    out = C.copy()
    for cols in step4_blocks(lbl, core_own_block):
        nm = 2 * lbl[cols[0]][2] + 1
        ns = len(cols) // nm
        if ns <= 1:
            continue
        B = np.zeros((ns, ns))
        for m in range(nm):
            idx = [cols[s * nm + m] for s in range(ns)]
            B += Porb[np.ix_(idx, idx)] / nm
        V = np.linalg.eigh(B)[1][:, ::-1]  # descending eigenvalue -> shell 0, 1, ...
        for sh in range(ns):
            for m in range(nm):
                out[:, cols[sh * nm + m]] = sum(
                    V[s, sh] * C[:, cols[s * nm + m]] for s in range(ns))
    return out


def replica(mol, d, core_own_block):
    """(labels, C_native, C_replica, S, C_gennbo, floor) for one molecule."""
    n, S, P = read_47(os.path.join(d, mol + ".47"))
    lpre, Cpre = read_naoc(os.path.join(d, mol + ".naocpre.txt"), tag="NAOCPRE")
    lnao, Cnat = read_naoc(os.path.join(d, mol + ".naoc.txt"))
    Cg, _, e33 = read_lfn33(os.path.join(d, mol + ".33"), n, S)
    assert [t[1:6] for t in lpre] == [t[1:6] for t in lnao], mol
    pre_occ = pre_occupations(Cpre, S @ P @ S, lpre)
    cls = [t[5] for t in lnao]
    C = step4(step3(Cpre, pre_occ, cls, S), S @ P @ S, lnao, core_own_block)
    e_n = float(np.abs(Cnat.T @ S @ Cnat - np.eye(Cnat.shape[1])).max())
    return lnao, Cnat, C, S, Cg, float(np.sqrt(2.0 * max(e_n, e33)))


def per_column_overlap(A, B, S):
    """|a_i^T S b_i| per column: sign-free, so no gauge has to be chosen."""
    return np.abs(np.einsum("ij,ij->j", A, S @ B))


def block_eigs(C, SPS, cols, nm):
    """The m-averaged (ns x ns) density of one step-4 block, descending."""
    ns = len(cols) // nm
    B = np.zeros((ns, ns))
    for m in range(nm):
        idx = [cols[s * nm + m] for s in range(ns)]
        B += (C.T @ SPS @ C)[np.ix_(idx, idx)] / nm
    return np.linalg.eigvalsh(B)[::-1]


def clusters_of(ev, tol=1e-08):
    """Shell ranks grouped by degenerate m-averaged eigenvalue.

    >>> clusters_of(np.array([2.0, 1.5, 0.3, 1e-12, 0.0]))
    [[0], [1], [2], [3, 4]]
    >>> clusters_of(np.array([1.0, 1.0]))
    [[0, 1]]
    """
    out = [[0]]
    for i in range(1, len(ev)):
        if abs(ev[i] - ev[out[-1][-1]]) <= tol * max(1.0, abs(ev[i])):
            out[-1].append(i)
        else:
            out.append([i])
    return out


def fidelity(Cnat, Crep, S, SPS, lbl):
    """Worst 1 - agreement between native and the replica, degeneracy-aware.

    A singleton eigenvalue cluster is compared column by column; a degenerate cluster, whose
    eigenvectors are not determined, is compared as a subspace.  Returns
    (worst_column, worst_subspace, n_degenerate_shells, n_shells), BOTH AS ANGLES.

    The first version returned the column half as a cosine deficit and the cluster half as a sine,
    and then thresholded them against two different numbers.  That was incommensurable, and the
    square root is why: 1 - |cos| = 1.05e-10 IS an angle of 1.45e-05.  Read in the same unit the two
    halves swapped places - the degenerate clusters came out at 2.8e-06 and 4.5e-06, an order BETTER
    than the columns they were being failed against.  So both are returned as angles and there is
    one gate, self-calibrated: the cluster disagreement may not exceed ten times the column
    disagreement.  That threshold is a measured quantity rather than a number chosen to pass - it
    says the undetermined directions are no worse than the determined ones, which is exactly what
    gate A has to establish, and it can fail.
    """
    worst, wsub, ndeg, nsh = 0.0, 0.0, 0, 0
    for cols in step4_blocks(lbl, False):
        nm = 2 * lbl[cols[0]][2] + 1
        ns = len(cols) // nm
        nsh += ns
        groups = clusters_of(block_eigs(Cnat, SPS, cols, nm)) if ns > 1 else [[0]]
        for g in groups:
            idx = [cols[sh * nm + m] for sh in g for m in range(nm)]
            if len(g) == 1:
                c = per_column_overlap(Crep[:, idx], Cnat[:, idx], S).min()
                worst = max(worst, float(np.sqrt(max(0.0, 1.0 - min(c, 1.0) ** 2))))
                continue
            ndeg += len(g)
            sv = np.linalg.svd(Cnat[:, idx].T @ S @ Crep[:, idx], compute_uv=False)
            wsub = max(wsub, float(np.sqrt(max(0.0, 1.0 - sv.min() ** 2))))
    return worst, wsub, ndeg, nsh


def run(root):
    mols = sorted(m for m in os.listdir(root) if os.path.isdir(os.path.join(root, m)))
    print("step4_core_block on %s, 1 thread, root %s" % (socket.gethostname(), root))
    rows = []
    for mol in mols:
        d = os.path.join(root, mol)
        if not all(os.path.isfile(os.path.join(d, mol + e))
                   for e in (".47", ".33", ".naoc.txt", ".naocpre.txt")):
            print("%-9s SKIP, missing an input" % mol)
            continue
        lbl, Cnat, Cbase, S, Cg, floor = replica(mol, d, False)
        _, _, Carm, _, _, _ = replica(mol, d, True)
        n47, S47, P47 = read_47(os.path.join(d, mol + ".47"))
        cls = [t[5] for t in lbl]
        core = [i for i, k in enumerate(cls) if k == 0]
        fid, fsub, ndeg, nsh = fidelity(Cnat, Cbase, S, S47 @ P47 @ S47, lbl)
        moved = float(1.0 - per_column_overlap(Cbase[:, core], Carm[:, core], S).min()) if core else 0.0
        sines = {}
        for tag, Cx in (("base", Cbase), ("arm", Carm)):
            idx_n = {c: [i for i, k in enumerate(cls) if CLASS[k] == c] for c in ORDER}
            gcls = [x["type"] for x in sorted(json.load(
                io.open(os.path.join(d, mol + ".aonao.nbo.json")))["nao"], key=lambda x: x["index"])]
            idx_g = {c: [i for i, k in enumerate(gcls) if k == c] for c in ORDER}
            sines[tag] = {c: (float(principal_sines(Cx[:, idx_n[c]], Cg[:, idx_g[c]], S).max())
                              if idx_n[c] and idx_g[c] else 0.0) for c in ORDER}
            sines[tag]["dVR"] = class_mass(Cx, Cg, S, [CLASS[k] for k in cls], gcls)[("Val", "Ryd")]
        rows.append(dict(mol=mol, floor=floor, fid=fid, fsub=fsub, moved=moved, ndeg=ndeg,
                         nsh=nsh, **sines))
    return rows


def verdict(rows):
    print("\n%s" % ("=" * 96))
    print("GATE A, replica fidelity: worst 1 - |c_replica^T S c_native| over all columns.")
    print("GATE B, discrimination:   worst 1 - |c_base^T S c_arm| over the CORE columns only.")
    print("Gate A is degeneracy-aware: a singleton eigenvalue cluster is compared column by column,")
    print("a degenerate one as a subspace, because its eigenvectors are not determined.  The last")
    print("column counts the shells that fell in a degenerate cluster - for those, a COLUMN-wise")
    print("comparison on either side of this lane is a gauge and not a measurement.")
    print("Both halves are ANGLES, so they are commensurable: 1 - |cos| = 1e-10 is an angle of")
    print("1.4e-05.  The gate is self-calibrated - the undetermined directions may not disagree by")
    print("more than 10x the determined ones - plus an absolute ceiling of 1e-03 on the replica.")
    print("%-9s %11s %11s %11s  %8s   %s" % ("mol", "A columns", "A clusters", "gate B", "degen/sh",
                                             "verdict"))
    ok = True
    for r in rows:
        bad = []
        if not r["fid"] < 1e-3:
            bad.append("A FAILS, replica is not native")
        if not r["fsub"] <= 10.0 * max(r["fid"], 1e-12):
            bad.append("A FAILS, degenerate clusters worse than 10x the columns")
        if not r["moved"] > 1e-9:
            bad.append("B FAILS, arm is a no-op")
        ok = ok and not bad
        print("%-9s %11.2e %11.2e %11.2e  %4d/%-4d   %s" % (
            r["mol"], r["fid"], r["fsub"], r["moved"], r["ndeg"], r["nsh"],
            ", ".join(bad) if bad else "both pass"))
    print("degenerate shells over all molecules: %d of %d" % (sum(r["ndeg"] for r in rows),
                                                              sum(r["nsh"] for r in rows)))
    if not ok:
        print("\nNO OUTCOME: a gate failed, so the arm is not read.  A replica that is not native")
        print("measures nothing, and an arm that changes nothing cannot be evidence either way.")
        return
    print("\nthe measurable part, principal sine of each class subspace against gennbo's own;")
    print("floor = sqrt(2 max|C^T S C - 1|), and a class counts as moved only above %.0f x floor."
          % CORE_SEP)
    print("%-9s %9s %-22s %-22s %-22s" % ("mol", "floor", "Cor base -> arm", "Val base -> arm",
                                          "Ryd base -> arm"))
    for r in rows:
        print("%-9s %9.1e %-22s %-22s %-22s" % (
            r["mol"], r["floor"],
            "%.2e -> %.2e" % (r["base"]["Cor"], r["arm"]["Cor"]),
            "%.2e -> %.2e" % (r["base"]["Val"], r["arm"]["Val"]),
            "%.2e -> %.2e" % (r["base"]["Ryd"], r["arm"]["Ryd"])))
    print("\nd_VR, the dimensions of valence space the two sides disagree about (NOT electrons):")
    for r in rows:
        print("  %-9s %.6f -> %.6f   (%+.6f)" % (r["mol"], r["base"]["dVR"], r["arm"]["dVR"],
                                                 r["arm"]["dVR"] - r["base"]["dVR"]))
    core_fixed = [r for r in rows if r["arm"]["Cor"] < CORE_SEP * r["floor"]]
    core_moved = [r for r in rows if r["arm"]["Cor"] < 0.5 * r["base"]["Cor"]]
    vr_moved = [r for r in rows
                if max(abs(r["arm"][c] - r["base"][c]) for c in ("Val", "Ryd")) > r["base"]["Cor"]]
    print("\ncore sine below %.0f x floor after the arm: %d/%d;  at least halved: %d/%d;"
          % (CORE_SEP, len(core_fixed), len(rows), len(core_moved), len(rows)))
    print("valence or Rydberg moved by more than the core disagreement's own size: %d/%d"
          % (len(vr_moved), len(rows)))
    if not core_moved:
        print("\nOUTCOME L: the arm leaves the core subspace where it was, so step 4's shared")
        print("(atom, l) block is NOT what perturbs native's cores.  The mechanism is REFUTED.")
        print("The core defect stays LOCATED AND OPEN, and the next place to look is step 3's own")
        print("Schmidt+OWSO on the core class, not step 4.")
    elif len(vr_moved) > len(rows) // 2:
        print("\nOUTCOME M: the arm does move the core, but valence and Rydberg move with it by more")
        print("than the core disagreement can account for - the classes are coupled through the")
        print("shared block.  Reported as a mechanism, NOT as a fix; the 22-molecule acceptance")
        print("gate with pf5/so2/sf6 held flat is what decides, and it has not been run here.")
    else:
        print("\nOUTCOME K: the core defect is localised to ONE LINE - the l_blocks key at")
        print("nao.cpp:342-345, which pools core columns into step 4's (atom, l) block and so")
        print("re-mixes the core that step 3 had just Schmidt-protected.  Fix = a core-only block.")
        print("Shipping waits on the 22-molecule acceptance gate with pf5/so2/sf6 held flat.")
    print("\n(outcome letter: %s)" % ("L" if not core_moved
                                      else "M" if len(vr_moved) > len(rows) // 2 else "K"))


def demo():
    import doctest
    for f in (sym_power, owso, step4_blocks, clusters_of):
        r = doctest.run_docstring_examples(f, globals(), verbose=False, name=f.__name__)
        del r
    assert doctest.testmod(sys.modules[__name__], verbose=False).failed == 0

    #Step 4 on a two-shell s block whose density is already diagonal must be the identity up to
    #column signs, and on a NON-diagonal one it must order the shells by descending occupancy.
    lbl = [(0, 0, 0, 0, 0, 0, 2.0), (1, 0, 0, 0, 1, 1, 1.0)]
    S2 = np.eye(2)
    C0 = np.eye(2)
    P2 = np.array([[0.3, 0.4], [0.4, 1.7]])  # the SECOND function is the occupied one
    C1 = step4(C0, P2, lbl)
    #0.4 of off-diagonal means the top eigenvector is a MIX, not the second basis function; the
    #check is that shell 0 is dominated by it, which is what "ordered by occupancy" actually claims.
    assert abs(C1[1, 0]) > 0.9, C1
    assert abs(float(np.linalg.det(C1.T @ S2 @ C1)) - 1.0) < 1e-12

    #The gate that matters: with a core in the block the arm MUST change the core column, and with
    #the core alone in its (atom, l) it must not.  Otherwise the arm is measuring nothing.
    lbl3 = [(0, 0, 0, 0, 0, 0, 2.0), (1, 0, 0, 0, 1, 1, 1.0), (2, 0, 0, 0, 2, 2, 0.0)]
    S3, P3 = np.eye(3), np.array([[1.99, 0.05, 0.01], [0.05, 1.2, 0.03], [0.01, 0.03, 0.004]])
    base = step4(np.eye(3), P3, lbl3, False)
    arm = step4(np.eye(3), P3, lbl3, True)
    assert per_column_overlap(base[:, [0]], arm[:, [0]], S3).min() < 1.0 - 1e-06
    #and the core really is untouched by the arm, which is the whole point of the partition
    assert abs(abs(arm[0, 0]) - 1.0) < 1e-12, arm

    #Cauchy interlacing, the direction claim the header makes: the valence shell can only gain.
    Bfull = P3
    lam_full = np.linalg.eigvalsh(Bfull)[::-1]
    lam_sub = np.linalg.eigvalsh(Bfull[1:, 1:])[::-1]
    assert lam_sub[0] >= lam_full[1] - 1e-12, (lam_sub, lam_full)

    #fidelity must SEPARATE a mis-assigned shell from a degenerate gauge, or it is the same void
    #gate again.  Two non-degenerate shells swapped must fail; a rotation inside a degenerate pair
    #must pass.  Both on the same block, so only the degeneracy differs.
    lblf = [(i, 0, 0, 0, i, 1 if i else 0, 0.0) for i in range(3)]
    Sf = np.eye(3)
    Pnd = np.diag([2.0, 1.0, 0.5])            # all three distinct
    Pdg = np.diag([2.0, 1.0, 1.0])            # shells 1 and 2 degenerate
    swap = np.eye(3)[:, [0, 2, 1]]
    th = 0.3
    rot = np.eye(3)
    rot[1:, 1:] = [[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]]
    assert fidelity(np.eye(3), swap, Sf, Pnd, lblf)[0] > 0.5      # mis-assignment caught
    assert fidelity(np.eye(3), rot, Sf, Pnd, lblf)[0] > 0.29      # same rotation, not degenerate;
    #  0.29 and not 0.04 because both halves are angles now: sin(0.3) = 0.2955, where the cosine
    #  deficit was 0.0447.  The unit is the whole point of the repair.
    #the degenerate gauge is forgiven, down to the SINE floor - 1e-16 in the cosine is 1.5e-08 here
    assert fidelity(np.eye(3), rot, Sf, Pdg, lblf)[1] < 1e-7
    assert fidelity(np.eye(3), rot, Sf, Pdg, lblf)[0] < 1e-7      # its column half stays clean
    assert fidelity(np.eye(3), np.eye(3), Sf, Pdg, lblf)[2] == 2  # the pair is counted

    print("demo ok")


if __name__ == "__main__":
    if "--demo" in sys.argv[1:]:
        demo()
    else:
        verdict(run(sys.argv[1]))
