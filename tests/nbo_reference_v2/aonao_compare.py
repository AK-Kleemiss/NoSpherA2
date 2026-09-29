"""Compare NoSpherA2's AO -> NAO transformation with gennbo 7's own, orbital by orbital.

    py -3.12 aonao_compare.py <dir with one subdirectory per molecule> [--full]
    py -3.12 aonao_compare.py --pre <same, with a .32 and a .naocpre.txt per molecule>
    py -3.12 aonao_compare.py --cascade <same, needing both ends: .32/.naocpre.txt and .33/.naoc.txt>
    py -3.12 aonao_compare.py --demo

--pre is the other end of a bisection: `$NBO AOPNAO=W $END` (unit 32 on this install) and
NAO_DUMP_CPRE=1 give the AO -> PRE-NAO matrices, which are upstream of every orthogonalisation, so
the two modes bracket the whole cascade.  It needs its own admissibility gate because neither side
prints a pre-NAO occupancy table - see pre_check - and its verdict is pre-registered in
aopnao_stage.sh rather than chosen once the numbers are on screen.

Every comparison before this one was of a FINAL table - NPA charges, NAO occupancies, class
totals.  A final table can say that an (atom, l) block came out with the wrong population; it
cannot say which orbital has the wrong shape, and a fix needs the second thing.  NBO 7 will
emit the intermediate after all: `$NBO AONAO=W $END` writes its AO -> NAO matrix to lfn 33 at
nine decimals (verified on this install, probe job 588007; AONAO=W48 is refused, 48 is
reserved).  `NAO_DUMP_C=1` writes ours.  Both sides read the SAME .47, so the AO order, the
overlap S and the density P are the same arrays and nothing is converted.

What is measured, per (atom, l) block: the m-averaged mixing matrix

    M[a][b] = (1 / (2l+1)) * sum_{m,m'} ( c_native(a,m)^T S c_gennbo(b,m') )**2

between the block's shells on the two sides, paired BY RANK (position within the block) and
never by shell label - ethane's C s block is printed 1s 2s 3s 5s 4s by NBO 7, and pairing on
those labels silently discarded 11 of 22 molecules once already.  M is m-averaged, so it needs
no mapping between NBO's component order and libcint's and is blind to a sign flip; M[a][a] = 1
means native's shell a IS gennbo's shell a.  Two numbers come out of it:

    in-block leak    1 - M[a][a] summed over b != a inside the block: native's shell a is partly
                     a different shell of the same (atom, l) - a valence/Rydberg mix, which is
                     exactly the leak the final tables see, now resolved per orbital
    out-of-block     1 - sum_b M[a][b] over the block: native's shell a is partly on another
                     atom or another l, which no intra-atomic step can produce

VALIDATION, before any of that is believed (a metric that has not been shown commensurable is
worth nothing - three have been retired on this branch for exactly that):

  1. the reshape of lfn 33 is not documented in the file, so both layouts are tried and exactly
     one must give C^T S C = 1;
  2. the NAO occupancies printed in the same run's .nbo must come back out of the matrix, as
     diag(C_g^T S P S C_g) - that is the NAO-basis density because C^-1 = C^T S, and both it
     and the wrong form C_g^T P C_g are computed so the report names which one matched;
  3. the AONAO run used a reduced keylist (AONAO=W alone, no NRT - minutes on benzene), which
     assumes the NAO stage is upstream of NBO/E2/NRT.  Its NAO table must therefore equal the
     stamped reference JSON's to 1e-5, or the shortcut is wrong and the comparison is void.
"""
import collections
import json
import os
import re
import socket
import sys

import contextlib
import io

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from config import load_nbo
except Exception:  # running from a copy without config.py next to it
    load_nbo = None

LANG_L = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}
CLASS = {0: "Cor", 1: "Val", 2: "Ryd"}


def numbers(text):
    return [float(t) for t in re.findall(r"[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?|[-+]?\d+(?:[EeDd][-+]?\d+)?",
                                         text.replace("D", "E"))]


def section(text, name):
    """The body of a $SECTION ... $END block of a .47 archive."""
    i = text.index(name)
    j = text.index("$END", i)
    return text[i + len(name):j]


def unpack_upper(vals, n):
    """A packed upper triangle in the .47's UPPER order - COLUMN by column - as a full matrix.

    Column by column is (1,1) (1,2) (2,2) (1,3) (2,3) (3,3), not row by row.  The first version of
    this function filled it row by row and its doctest used n = 2, where the two orders coincide, so
    the doctest could not fail; S and P were both scrambled and the layout gate then refused every
    molecule.  n = 3 is the smallest case that tells them apart, which is why it is the example.

    >>> unpack_upper([1., 2., 3., 4., 5., 6.], 3).tolist()
    [[1.0, 2.0, 4.0], [2.0, 3.0, 5.0], [4.0, 5.0, 6.0]]
    """
    m = np.zeros((n, n))
    k = 0
    for j in range(n):
        for i in range(j + 1):
            m[i, j] = m[j, i] = vals[k]
            k += 1
    return m


def read_47(path):
    text = open(path).read()
    n = int(re.search(r"NBAS=(\d+)", text).group(1))
    basis = section(text, "$BASIS")
    label_at = basis.index("LABEL")
    centers = [int(v) for v in numbers(basis[:label_at])]
    packed = n * (n + 1) // 2
    S = unpack_upper(numbers(section(text, "$OVERLAP")), n)
    dens = numbers(section(text, "$DENSITY"))
    # An OPEN archive puts alpha THEN beta in $DENSITY - measured, ratio exactly 2.0000 on
    # ch3 (1225 -> 2450) and o2 (1953 -> 3906), 1.0000 on water and benzene.  This used to
    # read [:packed], which on those two silently returned the ALPHA density wearing a
    # total-density label.  Refuse instead: every caller here compares spin-summed NAO
    # tables, so there is no correct spin for this reader to pick on its own.
    assert len(dens) == packed, (
        "%s: $DENSITY holds %d values, not the %d of one packed triangle (ratio %.4f). An "
        "open-shell archive carries alpha then beta here and this reader is spin-blind - the "
        "caller must read the spin it wants, by name." % (path, len(dens), packed, len(dens) / packed))
    P = unpack_upper(dens, n)
    assert len(centers) == n, (len(centers), n)
    return n, S, P


NAO_HEADER = "NAOs in the AO basis:"
PRE_HEADER = "PNAOs in the AO basis:"


def read_lfn(path, n, block, layout_test, what):
    """One of NBO 7's written transformation matrices, with the layout decided rather than assumed.

    `block` is matched at the start of a line, because "NAOs in the AO basis:" is a substring of
    "PNAOs in the AO basis:" - a plain `in` test would read a PNAO matrix as an NAO one and nothing
    downstream would notice.  `layout_test` returns the deviation from a property the correct layout
    has and the transpose does not; the loser is asserted to FAIL, because a test both layouts pass
    has not decided anything.
    """
    body = open(path).read()
    #A VOID, not a traceback: sf6's lfn 33 came back 0 bytes with rc = 0 from an earlier job.
    m = re.search(r"^[ \t]*" + re.escape(block), body, re.M)
    assert m, "no %r block in %s (%d bytes)" % (block, os.path.basename(path), len(body))
    vals = numbers(body[m.start():])[:n * n]
    assert len(vals) == n * n, "%s holds %d numbers, need %d" % (
        os.path.basename(path), len(vals), n * n)
    a = np.array(vals)
    cands = {"column-major": a.reshape(n, n).T, "row-major": a.reshape(n, n)}
    errs = {k: layout_test(C) for k, C in cands.items()}
    best = min(errs, key=errs.get)
    assert errs[best] < 1e-6, "neither layout of %s is %s: %s" % (
        os.path.basename(path), what, errs)
    other = [k for k in errs if k != best][0]
    assert errs[other] > 1e-6, "both layouts pass, the test cannot decide: %s" % errs
    return cands[best], best, errs[best]


def read_lfn33(path, n, S):
    """gennbo's AO -> NAO matrix, with the layout decided by C^T S C = 1 rather than assumed."""
    return read_lfn(path, n, NAO_HEADER,
                    lambda C: float(np.abs(C.T @ S @ C - np.eye(n)).max()), "S-orthonormal")


def read_lfn32(path, n, S):
    """gennbo's AO -> pre-NAO matrix (`$NBO AOPNAO=W $END`, unit 32 on this install).

    The layout CANNOT be decided by C^T S C = 1 here: pre-NAOs on different atoms are not
    orthogonal, which is the whole point of them, so the NAO test would fail on both candidates.
    What is still true of every pre-NAO is that it is normalised, so diag(C^T S C) = 1 decides it -
    a weaker property, and the assertion that the transpose fails it is what keeps that honest.
    """
    return read_lfn(path, n, PRE_HEADER,
                    lambda C: float(np.abs(np.diag(C.T @ S @ C) - 1.0).max()), "column-normalised")


def read_naoc(path, tag="NAOC"):
    """The native dump: (labels, C) with labels[i] = (atom, l, m, shell, class, occ)."""
    labels, cols = [], []
    for line in open(path):
        if not line.startswith(tag + " "):
            continue
        t = line.split()
        if t[1] == "index":
            continue
        labels.append((int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5]), int(t[6]), float(t[7])))
        cols.append([float(x) for x in t[8:]])
    C = np.array(cols).T  # column i = NAO i
    return labels, C


def gennbo_blocks(naos):
    """(atom0, l) -> list over rank of the (2l+1) gennbo NAO indices that make up that shell.

    gennbo prints an (atom, l) block COMPONENT-major: for water's oxygen p set the order is
    2p_x 3p_x 4p_x 5p_x, then 2p_y 3p_y ..., not 2p_x 2p_y 2p_z.  Slicing the block in groups of
    (2l+1) therefore grouped four shells of ONE component together and called it a shell, which is
    what held M[a][a] at ~1/3 for every p shell - one of three m pairs could match.  The grouping is
    read from the printed `lang` field, which says which component each row is, and the shell RANK is
    the position within that component's own list - never the printed shell label, which is each
    side's private bookkeeping.

    >>> naos = [dict(index=1, atom=1, lang="px"), dict(index=2, atom=1, lang="px"),
    ...         dict(index=3, atom=1, lang="py"), dict(index=4, atom=1, lang="py"),
    ...         dict(index=5, atom=1, lang="pz"), dict(index=6, atom=1, lang="pz")]
    >>> gennbo_blocks(naos)[(0, 1)]
    [[0, 2, 4], [1, 3, 5]]
    """
    comps = {}
    for e in naos:
        l = LANG_L[e["lang"][0].lower()]
        comps.setdefault((e["atom"] - 1, l), {}).setdefault(e["lang"].lower(), []).append(
            e["index"] - 1)
    blocks = {}
    for key, by_comp in comps.items():
        lists = [by_comp[c] for c in sorted(by_comp)]
        nsh = len(lists[0])
        assert all(len(v) == nsh for v in lists), "block %s has ragged components: %s" % (
            key, {c: len(v) for c, v in by_comp.items()})
        assert len(lists) == 2 * key[1] + 1, "block %s has %d components, expected %d" % (
            key, len(lists), 2 * key[1] + 1)
        blocks[key] = [[v[a] for v in lists] for a in range(nsh)]
    return blocks


def mixing(Cn, Cg, S, rows, cols):
    """The m-averaged mixing matrix between two shell lists of one (atom, l) block.

    `rows` and `cols` are lists of shells, each a list of the (2l+1) indices of its components, so
    neither side's component order has to match the other's - the m sum is over all pairs.

    Returns (M, total, w, O3), with w[a][j] the m-averaged weight of native shell a on gennbo NAO j,
    so any other grouping of the columns - same atom other l, other atom - is a sum over w and needs
    no second pass over the matrices, and O3[a][m][j] the SIGNED overlap itself.  w discards the sign
    and the coherence between two gennbo NAOs; O3 keeps both, which is what the signed, occupancy-
    weighted prediction needs - and what makes that prediction checkable exactly rather than to
    leading order.
    """
    flat = [i for sh in rows for i in sh]
    nm = len(rows[0])
    O = Cn[:, flat].T @ S @ Cg  # native block rows against ALL gennbo NAOs
    O3 = O.reshape(len(rows), nm, -1)
    w = (O3 ** 2).sum(axis=1) / nm
    M = np.zeros((len(rows), len(cols)))
    for a in range(len(rows)):
        for b, cb in enumerate(cols):
            M[a, b] = w[a, cb].sum()
    total = w.sum(axis=1)  # 1.0 by completeness
    return M, total, w, O3


LEVELS = [("per NAO, full block", False, False), ("per NAO, class sub-block", False, True),
          ("m-averaged, full block", True, False), ("m-averaged, class sub-block", True, True)]


def block_spectra(Dn, Dg, nshells, gshells, ncls, gcls):
    """Per (atom, l) block: |d eigen| and each side's self-consistency, at four candidate levels.

    "Spectrum or vectors" can only be asked at the level each side actually diagonalises, and that
    level is not obvious - assuming one and reading the answer is how three metrics were retired on
    this branch.  So all four candidates are computed and the one that is quoted is the one where
    BOTH sides reproduce their own reported occupancies as eigenvalues (`self` ~ 0):

      per NAO vs m-averaged   the NAO recipe diagonalises the density averaged over the (2l+1)
                              components of a shell, so that one radial function serves all m.  The
                              m-averaged matrix is nsh x nsh with element sum_m D[(a,m)][(b,m)], its
                              diagonal is the shell population, and a side that m-averages does NOT
                              reproduce its occupancies as eigenvalues of the per-NAO block.
      full block vs class     the recipe's last diagonalisation runs inside the natural minimal and
                              Rydberg parts separately, which is what leaves a valence/Rydberg
                              partition to be wrong about in the first place.

    Returns {level name: (|d eigen|, self native, self gennbo, n class-size mismatches)}.

    >>> D = np.diag([2.0, 2.0, 2.0, 0.1, 0.1, 0.1])            # two p shells, one atom
    >>> D[0, 3] = D[3, 0] = 0.3                                # cancels in the m sum ...
    >>> D[1, 4] = D[4, 1] = -0.3                               # ... so the m-averaged block is
    >>> sh = [[0, 1, 2], [3, 4, 5]]                            #     diagonal and the per-NAO is not
    >>> r = block_spectra(D, D, sh, sh, ["Val", "Ryd"], ["Val", "Ryd"])
    >>> round(r["m-averaged, full block"][1], 12), round(r["per NAO, full block"][1], 3)
    (0.0, 0.185)
    """
    def arm(D, shells, avg):
        if avg:
            A = np.array([[sum(D[sa[i], sb[i]] for i in range(len(sa))) for sb in shells]
                          for sa in shells])
        else:
            ix = [i for sh in shells for i in sh]
            A = D[np.ix_(ix, ix)]
        e = np.sort(np.linalg.eigvalsh(A))[::-1]
        return e, float(np.abs(e - np.sort(np.diag(A))[::-1]).sum())

    out = {}
    for name, avg, per_class in LEVELS:
        if per_class:
            keys = sorted(set(ncls) | set(gcls))
            groups = [([s for s, c in zip(nshells, ncls) if c == k],
                       [s for s, c in zip(gshells, gcls) if c == k]) for k in keys]
        else:
            groups = [(nshells, gshells)]
        deig = sn = sg = 0.0
        mism = 0
        for ns, gs in groups:
            if not ns and not gs:
                continue
            if len(ns) != len(gs):
                mism += 1
                continue
            en, s1 = arm(Dn, ns, avg)
            eg, s2 = arm(Dg, gs, avg)
            deig += float(np.abs(en - eg).sum())
            sn += s1
            sg += s2
        out[name] = (deig, sn, sg, mism)
    return out


def pre_check(C, shells, S, D):
    """A pre-NAO set's own defining property, on ONE side, inside one m-averaged (atom, l) block.

    Neither side prints a pre-NAO occupancy table, so the admissibility gate that made the AO -> NAO
    comparison readable - both sides reproducing their own reported numbers - has no counterpart
    here and something else has to carry it.  What both sides have is the DEFINITION of a pre-NAO:
    inside one (atom, l) block its columns are S-orthonormal, and m-averaged they diagonalise
    S P S.  Native by construction, gennbo by NBO's own definition.  Either failing is outcome C.

    The m-averaged test also checks the component grouping, which is otherwise an assumption: a
    wrong grouping mixes components into one "shell" and the off-diagonal does not vanish.

    Returns (|G - 1| max, |off-diagonal| max, eigenvalues descending).

    >>> S = np.eye(4)
    >>> D = np.diag([2.0, 2.0, 0.1, 0.1])
    >>> D[0, 2] = D[2, 0] = 0.3          # cancels in the m sum ...
    >>> D[1, 3] = D[3, 1] = -0.3         # ... so the m-averaged block is diagonal
    >>> o, off, e = pre_check(np.eye(4), [[0, 1], [2, 3]], S, D)
    >>> round(o, 12), round(off, 12), e.round(6).tolist()
    (0.0, 0.0, [4.0, 0.2])
    """
    flat = [i for sh in shells for i in sh]
    G = C[:, flat].T @ S @ C[:, flat]
    ortho = float(np.abs(G - np.eye(len(flat))).max())
    A = np.array([[sum(D[sa[i], sb[i]] for i in range(len(sa))) for sb in shells] for sa in shells])
    off = float(np.abs(A - np.diag(np.diag(A))).max())
    return ortho, off, np.sort(np.diag(A))[::-1]


def compare_pre(mol, d):
    """The AO -> pre-NAO matrices from both sides, on the same .47.

    Two things are different from the NAO comparison and both make this arm SHARPER, not weaker:

      - there is no out-of-block question.  A pre-NAO is strictly intra-atomic, so both sides' block
        spans the same subspace (that atom's AOs of that l) and the in-block m-averaged mixing sums
        to exactly 1.  That is an identity, so it is a gate: if it misses, one side's pre-NAOs are
        not confined to the block and the comparison is void rather than small.
      - `mixing`'s completeness total is NOT 1 here, because pre-NAOs on different atoms overlap.
        The in-block row sum is what is 1, and it is the one that is checked.

    The saved NAOCPRE occupation column duplicates the final NAOC occupation column on all eight
    cases and is not a pre-NAO measurement.  It is not used: occupations here come from C^T S P S C.
    """
    n, S, P = read_47(os.path.join(d, mol + ".47"))
    Cg, layout, lerr = read_lfn32(os.path.join(d, mol + ".32"), n, S)
    labels, Cn = read_naoc(os.path.join(d, mol + ".naocpre.txt"), tag="NAOCPRE")
    assert Cn.shape == (n, n), "native dump is %s, .47 says %d" % (Cn.shape, n)
    run_path = os.path.join(d, mol + ".aopnao.nbo.json")
    run = load_nbo(run_path) if load_nbo is not None else json.load(open(run_path))
    naos = run["nao"]
    Dg = Cg.T @ S @ P @ S @ Cg
    Dn = Cn.T @ S @ P @ S @ Cn
    rows = []
    for (atom, l), gshells in sorted(gennbo_blocks(naos).items()):
        nm = 2 * l + 1
        nrows = [i for i, lab in enumerate(labels) if lab[1] == atom and lab[2] == l]
        nrows.sort(key=lambda i: (labels[i][4], labels[i][3]))  # shell rank, then m
        assert len(nrows) == nm * len(gshells), "block (%d,%d): %d native vs %d gennbo" % (
            atom, l, len(nrows), nm * len(gshells))
        nshells = [nrows[a * nm:(a + 1) * nm] for a in range(len(gshells))]
        on, offn, en = pre_check(Cn, nshells, S, Dn)
        og, offg, eg = pre_check(Cg, gshells, S, Dg)
        M, _, _, _ = mixing(Cn, Cg, S, nshells, gshells)
        inblock = M.sum(axis=1)
        for a in range(M.shape[0]):
            order = np.argsort(M[a])[::-1]
            #A degenerate eigenvalue fixes its EIGENSPACE, not its eigenvectors: two shells with the
            #same pre-occupancy can be any orthogonal mix of each other on either side, and M[a][a]
            #would then read < 1 with nothing wrong.  So the invariant is the projector - M summed
            #over the columns degenerate with a - and the rank-paired diagonal is reported next to
            #it.  Where a shell is non-degenerate the two are the SAME number, which is why both are
            #carried instead of the looser one replacing the tighter.
            deg = [b for b in range(len(gshells)) if abs(en[b] - en[a]) < DEG_TOL]
            rows.append(dict(mol=mol, atom=atom, l=l, rank=a, nsh=len(gshells),
                             cls=CLASS[labels[nshells[a][0]][5]], gcls=naos[gshells[a][0]]["type"],
                             diag=float(M[a, a]), best=int(order[0]),
                             bestval=float(M[a, order[0]]), inblock=float(inblock[a]),
                             degsum=float(M[a, deg].sum()), ndeg=len(deg),
                             bestdeg=int(order[0]) in deg,
                             ortho_n=on, off_n=offn, ortho_g=og, off_g=offg,
                             deig=float(np.abs(en - eg).sum()),
                             occ_n=float(en[a]), occ_g=float(eg[a])))
    return dict(mol=mol, n=n, layout=layout, lerr=lerr, rows=rows)


#The gate thresholds, argued rather than tuned.  gennbo writes nine decimals, so a coefficient
#carries up to 5e-10 and an n-term dot product of them ~1e-8; a 222-column block accumulates a
#little more.  1e-5 is three decades above that floor and three below the NAO-level leaks this is
#bisecting.
PRE_GATE = 1e-5
#The in-block sum is an exact identity in exact arithmetic, so it is tempting to gate it at 1e-9 -
#and that is a threshold set BELOW the measurement floor: the first run refused all eight molecules
#at 5e-10 to 2.2e-09, which is the nine-decimal roundoff and not a missing dimension.  A span that
#really differs loses O(0.1) of the sum, so 1e-7 separates the two by six decades either way.
PRE_SPAN = 1e-7
PRE_FLOOR = 1e-4
#Two pre-occupancies this close are degenerate for the purpose of fixing eigenvectors: gennbo's own
#tables print five decimals, and the matrices here carry ~1e-8.
DEG_TOL = 1e-6


def pre_verdict(rows_by_mol):
    """The pre-registered A/B/C verdict of the bisection, and nothing beyond what it licenses."""
    rows = [r for rs in rows_by_mol.values() for r in rs]
    if not rows:
        print("\nno PNAO rows: nothing to say")
        return
    print("\n=== PNAO bisection: the AO -> pre-NAO matrices, both sides, same .47")
    bad_self, bad_span = [], []
    for mol, rs in sorted(rows_by_mol.items()):
        sn = max(max(r["ortho_n"], r["off_n"]) for r in rs)
        sg = max(max(r["ortho_g"], r["off_g"]) for r in rs)
        span = max(abs(r["inblock"] - 1.0) for r in rs)
        leak = [1.0 - r["diag"] for r in rs]
        dleak = [1.0 - r["degsum"] for r in rs]
        mism = [r for r in rs if r["best"] != r["rank"] and not r["bestdeg"]]
        print("  %-9s shells %3d  self-check native %.2e gennbo %.2e  in-block sum-1 %.2e"
              "  leak max %.5f mean %.5f  degeneracy-invariant max %.5f (%d degenerate shells)"
              "  rank mismatches %d  |d eigen| max %.2e" % (
                  mol, len(rs), sn, sg, span, max(leak), sum(leak) / len(leak), max(dleak),
                  sum(1 for r in rs if r["ndeg"] > 1), len(mism), max(r["deig"] for r in rs)))
        if sn > PRE_GATE or sg > PRE_GATE:
            bad_self.append((mol, sn, sg))
        if span > PRE_SPAN:
            bad_span.append((mol, span))
    #C first, because a failed self-check makes the rest unquotable however it came out.
    if bad_self or bad_span:
        print("\n  OUTCOME C: the PNAO matrices on disk cannot arbitrate.")
        for mol, sn, sg in bad_self:
            side = "native" if sn > sg else "gennbo"
            print("    %-9s self-check fails on %s (native %.2e, gennbo %.2e) - a pre-NAO set that"
                  " is not S-orthonormal and SPS-diagonal inside its own block is not a pre-NAO"
                  " set, and nothing is quoted against it" % (mol, side, sn, sg))
        for mol, span in bad_span:
            print("    %-9s in-block mixing sums to 1 + %.2e: the two sides' blocks do not span the"
                  " same subspace, so the rotation is not the only difference" % (mol, span))
        return
    worst = min(r["diag"] for r in rows)
    #The verdict runs on the degeneracy invariant, because that is the quantity the two sides are
    #both entitled to; the rank-paired diagonal is printed beside it and they coincide wherever a
    #shell is alone in its eigenvalue.
    wdeg = min(r["degsum"] for r in rows)
    #A "wrong" best partner inside one eigenspace is not an ordering error either: which member of a
    #degenerate pair a side prints first is arbitrary on both sides.
    mism = [r for r in rows if r["best"] != r["rank"] and not r["bestdeg"]]
    print("\n  gate passed on both sides everywhere (< %.0e), and the in-block sums are 1 to %.0e:"
          " the comparison is readable and it is a pure rotation question." % (PRE_GATE, PRE_SPAN))
    if 1.0 - wdeg <= PRE_FLOOR and not mism:
        print("  OUTCOME B: pre-NAOs AGREE (worst eigenspace-projected %.6f, worst rank-paired"
              " M[a][a] %.6f, no rank mismatches)." % (wdeg, worst))
        print("    Licensed: the disagreement is created strictly between the pre-NAO and the"
              " m-averaged NAO block - the orthogonalisation cascade and the class input it takes,"
              " steps 2 to 4, and nothing earlier.")
        print("    NOT licensed: which of Schmidt and OWSO. Two arms over a bracketed interval now,"
              " not a guess.")
    else:
        print("  OUTCOME A: pre-NAOs DISAGREE (worst eigenspace-projected %.6f, worst rank-paired"
              " M[a][a] %.6f, %d rank mismatches of %d)." % (wdeg, worst, len(mism), len(rows)))
        print("    Licensed: the two sides differ already at step 1's atomic (atom, l) eigenproblem,"
              " upstream of every orthogonalisation; no downstream arbitration is interpretable"
              " until that is fixed.")
        print("    NOT licensed: WHICH ingredient of step 1 (m-averaging, gross vs net density, the"
              " S^-1/2 metric, the shell grouping). That is a further bisection with its own arms.")
    #By the invariant first, then by the rank-paired diagonal: with the invariant at 1.0 everywhere
    #the second key is what puts the interesting rows on screen instead of eight arbitrary ones.
    for r in sorted(rows, key=lambda r: (round(r["degsum"], 6), r["diag"]))[:8]:
        print("    %-9s atom %2d l=%d rank %d/%d %s/%-4s  M[a][a]=%.6f  eigenspace(%d)=%.6f"
              "  best rank %d (%.6f)  occ %.5f/%.5f" % (
                  r["mol"], r["atom"], r["l"], r["rank"], r["nsh"], r["cls"], r["gcls"], r["diag"],
                  r["ndeg"], r["degsum"], r["best"], r["bestval"], r["occ_n"], r["occ_g"]))


def main_pre(root, verbose=False):
    print("aonao_compare --pre on %s, 1 thread, root %s" % (socket.gethostname(), root))
    rows_by_mol, void, missing = {}, [], []
    for mol in sorted(os.listdir(root)):
        d = os.path.join(root, mol)
        if not os.path.isdir(d):
            continue
        if not os.path.isfile(os.path.join(d, mol + ".32")):
            missing.append(mol)
            continue
        try:
            res = compare_pre(mol, d)
        except AssertionError as e:
            void.append((mol, str(e)))
            print("%-10s VOID: %s" % (mol, e))
            continue
        print("%-10s n=%-4d lfn32 %s  |diag(C^T S C)-1|=%.1e  blocks/shells %d" % (
            res["mol"], res["n"], res["layout"], res["lerr"], len(res["rows"])))
        rows_by_mol[mol] = res["rows"]
    print("\ndenominator: %d molecules compared, %d VOID, %d without an lfn 32" % (
        len(rows_by_mol), len(void), len(missing)))
    if missing:
        print("  no lfn 32: %s" % " ".join(missing))
    pre_verdict(rows_by_mol)
    return 0


#------------------------------------------------------------------- the cascade's net transform
#Pre-registered before any of these numbers exist; the outcome letters continue the A/B/C of
#aopnao_stage.sh.  The bisection put the disagreement strictly between the pre-NAO level, where the
#two sides agree to 1.000000 on every shell of all eight molecules, and the m-averaged NAO block,
#where they do not.  Between those two levels sits the orthogonalisation cascade - and BOTH of its
#ends are already on disk for all eight: AO -> PNAO (gennbo unit 32, native NAOCPRE) and AO -> NAO
#(unit 33, NAOC), read off the same .47.  So the cascade's ENTIRE net effect is available per side as
#the PNAO -> NAO transform T with Cnao = Cpre T, with no new job, no cluster time and no assumption
#about what the cascade contains.  T is not a guess; the only question is its SHAPE.
#
#  D.  A gate fails.  Quote NOTHING for that molecule: the matrices on disk do not determine its net
#      transform.  Both gates are ONE-SIDED - each side's T must reproduce THAT side's own NAOs from
#      THAT side's own PNAOs to the printing floor, and T^T G T must be the identity with
#      G = Cpre^T S Cpre.  If no molecule survives, that is the whole result and nothing follows.
#  E.  The two transforms have DIFFERENT shape.  Licenses: the cascade's STRUCTURE differs between
#      the codes, so the fault is not a coefficient inside a correct cascade, and the interval closes
#      with no substitution experiment at all.  Does NOT license naming the step - unless the block
#      that differs is one step's own block, and then it is named BY that block and not by the story.
#  F.  SAME shape, differing magnitude.  Licenses: a coefficient error inside a correctly structured
#      cascade, which is what makes the OWSO-weight and class-partition substitution arms worth
#      running.  Does NOT license which coefficient, and does not rank the two candidates.
#
#The two candidate structures have DIFFERENT measurable shapes, so neither has to be assumed:
#
#  A Schmidt step against a priority order is TRIANGULAR in that partition.  Native's step 3 says so
#  in its own comment - core, then valence, then Rydberg, each projected out of everything above it -
#  and that forces T[i, j] = 0 whenever class(i) > class(j): a core NAO is built out of core PNAOs
#  alone, a valence NAO may reach into core, a Rydberg NAO into both.  Measured as the Frobenius mass
#  of every ordered class block of T, as a share of ||T||.
#
#  An occupancy-weighted symmetric orthogonalisation is T = W^-1 M with M symmetric and W the
#  positive diagonal of weights, so SOME positive diagonal makes it symmetric.  The weights are not
#  known and do not need to be: |T_ij| / |T_ji| must factor as d_j / d_i, i.e.
#  log|T_ij| - log|T_ji| = a_j - a_i, a least-squares fit whose residual is the deviation.  It is
#  read on the WITHIN-class blocks, where the OWSO acts and the Schmidt step is exactly the identity
#  (projecting a valence PNAO against core changes its core rows, never its valence rows), and only
#  on pairs from DIFFERENT (atom, l): step 4's per-(atom, l) re-diagonalisation is a further
#  orthogonal factor inside those, so same-(atom, l) pairs are reported apart instead of mixed in.
#
#Both measures are GAUGE-INVARIANT, which is why they are the ones quoted.  Every PNAO and every NAO
#column carries an arbitrary sign on each side, so T -> diag(s) T diag(t) with s, t = +-1 is the SAME
#transform.  max|T - T^T| is not invariant under that and is never quoted; |T_ij| / |T_ji| is, every
#block norm is, and the sign structure is only quoted in the one form that survives it - whether the
#products sign(T_ij) sign(T_ji) FACTOR as u_i u_j, which is the same cocycle condition in the sign
#group (the raw products do not survive: the gauge multiplies them by (s_i t_i)(s_j t_j)).
#
#Thresholds: native dumps 10 fixed decimals and gennbo 9, so an entry carries ~5e-10 and the solve
#inflates that by cond(Cpre), which is printed beside every gate so a failure is diagnosable rather
#than mysterious.  1e-6 is three decades above that floor and three below the structure being looked
#for - and a threshold below the measurement floor is not a gate, which is exactly how the in-block
#identity refused all eight at 1e-9 one level up.
CASC_GATE = 1e-6
#A pair enters the symmetry fit only if BOTH its entries clear this; below it a ratio is roundoff
#over roundoff.  The class-block masses are not floored - a norm needs no floor.
CASC_FLOOR = 1e-3
#exp(resid) - 1 is the worst relative error left in the ratios after the best diagonal: 1e-3 admits
#"some positive diagonal makes it symmetric" at the same decade as the entry floor.
CASC_SYM = 1e-3
#A class block counts as EMPTY below this share of ||T||: six decades under a filled one.
CASC_ZERO = 1e-6
ORDER = ("Cor", "Val", "Ryd")


def cascade_transform(Cpre, Cnao, S):
    """The net PNAO -> NAO transform T with Cnao = Cpre T, plus its two one-sided gates.

    Returns (T, reproduce, ortho, cond).  `reproduce` = max|Cpre T - Cnao| is the printing-floor
    test.  `ortho` = max|T^T G T - 1| with G = Cpre^T S Cpre is an EXACT identity, because the NAOs
    are S-orthonormal - and it is free.  It also says where the magnitude lives: G is built from the
    pre-NAOs, which the two sides already agree on to 1.000000, so whatever T carries in magnitude is
    constrained by a matrix both sides share, and what is left to compare is shape.

    >>> rng = np.random.default_rng(0)
    >>> S = np.eye(5); Cpre = rng.normal(size=(5, 5)); T0 = rng.normal(size=(5, 5))
    >>> T, rep, _, _ = cascade_transform(Cpre, Cpre @ T0, S)
    >>> bool(np.abs(T - T0).max() < 1e-10), bool(rep < 1e-10)
    (True, True)
    >>> _, rep, _, _ = cascade_transform(Cpre[:, :2], Cpre[:, 2:4], S)   # not in the span
    >>> bool(rep > 1e-3)
    True
    >>> Q = np.linalg.qr(rng.normal(size=(5, 5)))[0]      # S-orthonormal NAOs, S = 1
    >>> _, _, orth, _ = cascade_transform(Cpre, Q, S)
    >>> bool(orth < 1e-12)
    True
    """
    T, _, _, _ = np.linalg.lstsq(Cpre, Cnao, rcond=None)
    rep = float(np.abs(Cpre @ T - Cnao).max())
    G = Cpre.T @ S @ Cpre
    ortho = float(np.abs(T.T @ G @ T - np.eye(T.shape[1])).max())
    return T, rep, ortho, float(np.linalg.cond(Cpre))


def sym_shape(T, pairs, floor=CASC_FLOOR):
    """How far T is from being made symmetric by SOME positive diagonal, without knowing which one.

    log|T_ij| - log|T_ji| = a_j - a_i, solved for a by least squares over the pairs whose two entries
    both clear `floor`; `resid` is the worst leftover in log units, so exp(resid) - 1 is a relative
    error in the ratios.  `plain` is max|log|T_ij| - log|T_ji|| on the same pairs - what those ratios
    would give if no diagonal were allowed - carried so that a small residual can never be read as
    "T is symmetric" when it is nothing of the kind.

    >>> M = np.array([[1.0, 0.4, -0.2], [0.4, 1.0, 0.3], [-0.2, 0.3, 1.0]])
    >>> T = np.diag([1.0, 2.0, 0.5]) @ M            # W^-1 M, so diagonally symmetrizable
    >>> p = [(0, 1), (0, 2), (1, 2)]
    >>> r = sym_shape(T, p, floor=0.01)
    >>> r["n"], round(r["resid"], 12), round(r["plain"], 4)
    (3, 0.0, 1.3863)
    >>> g = np.diag([1.0, -1.0, 1.0])               # a sign gauge changes nothing measured
    >>> r2 = sym_shape(g @ T @ g, p, floor=0.01)
    >>> (r2["n"], round(r2["resid"], 12)) == (r["n"], round(r["resid"], 12))
    True
    >>> T[0, 1] = 0.9        # break the ratio structure: the triangle no longer closes, and its
    >>> round(sym_shape(T, p, floor=0.01)["resid"], 4)   # 0.8109 of discrepancy spreads over 3 edges
    0.2703
    >>> sym_shape(T, p, floor=10.0)["n"]            # nothing clears the floor
    0
    """
    idx = [(i, j) for i, j in pairs if min(abs(T[i, j]), abs(T[j, i])) > floor]
    if not idx:
        return dict(n=0, resid=None, plain=None, pairs=[])
    A = np.zeros((len(idx), T.shape[0]))
    L = np.empty(len(idx))
    for k, (i, j) in enumerate(idx):
        A[k, j], A[k, i] = 1.0, -1.0
        L[k] = np.log(abs(T[i, j])) - np.log(abs(T[j, i]))
    a, _, _, _ = np.linalg.lstsq(A, L, rcond=None)
    return dict(n=len(idx), resid=float(np.abs(L - A @ a).max()),
                plain=float(np.abs(L).max()), pairs=idx)


def sign_frustration(pairs, neg):
    """Do the products sign(T_ij) sign(T_ji) FACTOR as u_i u_j?  Returns (edges, violations).

    A W^-1 M form with M symmetric and W positive has sign(T_ij) = sign(T_ji) on every pair, but the
    sign of a PNAO column and of an NAO column are each an arbitrary gauge, so that raw statement is
    not measurable: T -> diag(s) T diag(t) multiplies the product by (s_i t_i)(s_j t_j).  What
    survives the gauge is whether the products factor, which is a two-colouring of the pair graph -
    the same cocycle condition as the magnitudes, in the sign group.  `neg` holds the (i, j), i < j,
    whose product is negative.  The colouring is taken from a spanning tree, so the count is an UPPER
    bound on the minimum frustration; zero is exact, and zero is the only reading quoted.

    >>> sign_frustration([(0, 1), (1, 2), (0, 2)], {(0, 1), (0, 2)})
    (3, 0)
    >>> sign_frustration([(0, 1), (1, 2), (0, 2)], {(0, 1)})
    (3, 1)
    >>> sign_frustration([], set())
    (0, 0)
    """
    adj = {}
    for i, j in pairs:
        adj.setdefault(i, []).append(j)
        adj.setdefault(j, []).append(i)
    val = {}
    for start in adj:
        if start in val:
            continue
        val[start], stack = 0, [start]
        while stack:
            i = stack.pop()
            for j in adj[i]:
                if j not in val:
                    val[j] = val[i] ^ (1 if (min(i, j), max(i, j)) in neg else 0)
                    stack.append(j)
    bad = sum(1 for i, j in pairs
              if val[i] ^ val[j] != (1 if (min(i, j), max(i, j)) in neg else 0))
    return len(pairs), bad


def tri_shape(T, rcls, ccls, order=ORDER):
    """Frobenius mass of every ordered class block of T as a share of ||T||.

    A Schmidt step against the priority order leaves the blocks with class(row) > class(col) empty -
    a core NAO is built from core PNAOs alone.  Which triangle is the empty one is not assumed
    either: both are returned and the report prints them side by side.

    >>> c = ["Cor", "Val", "Ryd"]
    >>> b = tri_shape(np.array([[1.0, 0.5, 0.5], [0.0, 1.0, 0.5], [0.0, 0.0, 1.0]]), c, c)
    >>> b[("Val", "Cor")], b[("Ryd", "Val")]
    (0.0, 0.0)
    >>> round(b[("Cor", "Val")], 6)
    0.258199
    """
    tot = float(np.linalg.norm(T))
    out = {}
    for a in order:
        rows = [i for i, c in enumerate(rcls) if c == a]
        if not rows:
            continue
        for b in order:
            cols = [j for j, c in enumerate(ccls) if c == b]
            if not cols:
                continue
            out[(a, b)] = float(np.linalg.norm(T[np.ix_(rows, cols)])) / tot
    return out


def class_locality(T, cls, al, order=ORDER):
    """Per cross-class block, how its mass splits into same-(atom, l) and cross-(atom, l).

    This is the measure that actually separates the two steps, and it came out of the triangularity
    numbers rather than being designed with them: the priority Schmidt step is a MOLECULAR projection,
    so what it puts into a block is spread over every atom, while the only other thing that can move
    amplitude between classes is the per-(atom, l) re-diagonalisation, which is strictly intra-atomic
    and intra-l by construction on both sides.  A cross-class block whose mass is entirely
    same-(atom, l) was therefore produced by the second, and one carrying cross-(atom, l) mass was
    not.  Both halves are block norms, so both survive the sign gauge.

    >>> c = ["Cor", "Val", "Val"]
    >>> al = [(0, 0), (0, 0), (1, 0)]
    >>> T = np.array([[1.0, 0.0, 0.0], [0.3, 1.0, 0.0], [0.0, 0.0, 1.0]])
    >>> s, x = class_locality(T, c, al)[("Val", "Cor")]
    >>> round(s, 4), round(x, 4)                 # the 0.3 sits on atom 0's own l = 0
    (0.1697, 0.0)
    """
    tot = float(np.linalg.norm(T))
    keys = sorted(set(al))
    code = np.array([keys.index(x) for x in al])
    same = code[:, None] == code[None, :]
    out = {}
    for a in order:
        rows = [i for i, c in enumerate(cls) if c == a]
        for b in order:
            cols = [j for j, c in enumerate(cls) if c == b]
            if a == b or not rows or not cols:
                continue
            blk, m = T[np.ix_(rows, cols)], same[np.ix_(rows, cols)]
            out[(a, b)] = (float(np.linalg.norm(blk[m])) / tot,
                           float(np.linalg.norm(blk[~m])) / tot)
    return out


def cascade_sides(mol, d):
    """Both sides' net PNAO -> NAO transform, measured the same way off the same .47."""
    n, S, P = read_47(os.path.join(d, mol + ".47"))
    Cg_pre, lay32, e32 = read_lfn32(os.path.join(d, mol + ".32"), n, S)
    Cg_nao, lay33, e33 = read_lfn33(os.path.join(d, mol + ".33"), n, S)
    lpre, Cn_pre = read_naoc(os.path.join(d, mol + ".naocpre.txt"), tag="NAOCPRE")
    lnao, Cn_nao = read_naoc(os.path.join(d, mol + ".naoc.txt"))
    assert Cn_pre.shape == Cn_nao.shape == (n, n), "native dumps %s / %s, .47 says %d" % (
        Cn_pre.shape, Cn_nao.shape, n)
    #A row of T and a column of T are the same orbital only if native's two dumps are in one order.
    #Asserted on the self-identifying fields rather than assumed - a permutation would turn every
    #shape number below into some other matrix's.  The CLASS field is compared separately, because a
    #re-labelling between the two dumps would be a finding and not a pairing failure.
    assert [l[:5] for l in lpre] == [l[:5] for l in lnao], (
        "native's NAOCPRE and NAOC are not in the same orbital order")
    recls = sum(1 for a, b in zip(lpre, lnao) if a[5] != b[5])

    def load(suffix):
        p = os.path.join(d, mol + suffix)
        return load_nbo(p) if load_nbo is not None else json.load(open(p))

    naos = sorted(load(".aonao.nbo.json")["nao"], key=lambda e: e["index"])
    pre_naos = sorted(load(".aopnao.nbo.json")["nao"], key=lambda e: e["index"])
    assert len(naos) == n, "gennbo's NAO table has %d rows, .47 says %d" % (len(naos), n)
    #gennbo's rows are the unit-32 columns and its columns the unit-33 ones, which is a POSITIONAL
    #pairing.  It is not free: the --pre arm already found each (atom, l) block of unit 32 at the NAO
    #table's own indices, S-orthonormal inside the block and diagonalising the m-averaged S P S
    #there, which is what says unit 32 is indexed like the NAO table.  The two runs are also required
    #to agree about the classes, which catches a mismatched pair of jobs.
    gcls = [e["type"] for e in naos]
    assert [e["type"] for e in pre_naos] == gcls, (
        "the AONAO and AOPNAO runs disagree about the NAO classes")
    #A class block is only the SAME block on both sides if both sides put the same orbitals into it,
    #so the verdict below is about the transforms rather than about the labels only once this holds.
    #Measured IDENTICAL on all eight, which is what the core-block verdict rests on.
    assert_same_classes(lnao, naos)
    sides = {}
    for side, Cpre, Cnao, cls, al in (
            ("native", Cn_pre, Cn_nao, [CLASS[l[5]] for l in lnao], [(l[1], l[2]) for l in lnao]),
            ("gennbo", Cg_pre, Cg_nao, gcls,
             [(e["atom"] - 1, LANG_L[e["lang"][0].lower()]) for e in naos])):
        T, rep, ortho, cond = cascade_transform(Cpre, Cnao, S)
        cross, same = [], []
        for a in range(n):
            for b in range(a + 1, n):
                if cls[a] == cls[b]:
                    (same if al[a] == al[b] else cross).append((a, b))
        sc, ss = sym_shape(T, cross), sym_shape(T, same)
        neg = {(i, j) for i, j in sc["pairs"] if T[i, j] * T[j, i] < 0}
        edges, frus = sign_frustration(sc["pairs"], neg)
        sides[side] = dict(rep=rep, ortho=ortho, cond=cond, blocks=tri_shape(T, cls, cls),
                           loc=class_locality(T, cls, al), sym=sc, sym_same=ss, edges=edges,
                           frus=frus, cls=cls,
                           dom=sum(1 for j in range(n) if int(np.abs(T[:, j]).argmax()) == j),
                           norm=float(np.linalg.norm(T)))
    #The two pre-NAO sets carry the same metric or they do not, and that is checkable WITHOUT pairing
    #a single column: the spectrum of G is invariant under any permutation of the columns and under
    #any per-column sign, which is exactly the gauge freedom between the two sides.  It matters here
    #because T^T G T = 1 holds on both sides, so if the spectra agree then both transforms are
    #isometries of one metric and no magnitude difference between them can be hiding in G.
    #(The entrywise |G_n| - |G_g| is NOT the test and was void when tried: native orders its columns
    #(atom, l, shell, m) while gennbo orders an (atom, l) block component-major, so entrywise it
    #compares two different orderings and reports a permutation as a disagreement.)
    Gn, Gg = Cn_pre.T @ S @ Cn_pre, Cg_pre.T @ S @ Cg_pre
    deig = float(np.abs(np.sort(np.linalg.eigvalsh(Gn)) - np.sort(np.linalg.eigvalsh(Gg))).max())
    return dict(mol=mol, n=n, lay32=lay32, lay33=lay33, e32=e32, e33=e33, recls=recls, sides=sides,
                deig=deig, dtrace=float(abs(Gn.trace() - Gg.trace())))


def cascade_shape(s):
    """(the set of empty class blocks, is it diagonally symmetrizable) - the two shape statements."""
    empty = frozenset(k for k, v in s["blocks"].items() if k[0] != k[1] and v < CASC_ZERO)
    return empty, bool(s["sym"]["n"] > 0 and s["sym"]["resid"] < CASC_SYM and s["frus"] == 0)


def cascade_verdict(results):
    """The pre-registered D/E/F verdict, and not one word past what it licenses."""
    ok, dropped = [], []
    for r in results:
        bad = ["%s rep=%.1e ortho=%.1e" % (side, s["rep"], s["ortho"])
               for side, s in sorted(r["sides"].items())
               if not (s["rep"] < CASC_GATE and s["ortho"] < CASC_GATE)]
        if bad:
            dropped.append((r, bad))
        else:
            ok.append(r)
    print("\ngate: T reproduces its OWN side's NAOs and satisfies T^T G T = 1, both at %g" %
          CASC_GATE)
    for r, bad in dropped:
        print("  OUTCOME D %-10s quoted nowhere below: %s" % (r["mol"], "; ".join(bad)))
    if not ok:
        print("\nOUTCOME D on every molecule: the matrices on disk do not determine the net")
        print("transform, so nothing about the cascade is quoted, in either direction.")
        return
    print("  %d of %d molecules pass on both sides" % (len(ok), len(results)))
    print("\nclass-block mass as a share of ||T||.  A Schmidt step against Cor < Val < Ryd empties")
    print("every block with class(row) > class(col); the other triangle is what it fills.")
    print("%-9s %-7s %9s %9s %11s %9s" % ("mol", "side", "wrong-tri", "right-tri", "worst-wrong",
                                          "cond"))
    for r in ok:
        for side, s in sorted(r["sides"].items()):
            lo = {k: v for k, v in s["blocks"].items()
                  if ORDER.index(k[0]) > ORDER.index(k[1])}
            hi = {k: v for k, v in s["blocks"].items()
                  if ORDER.index(k[0]) < ORDER.index(k[1])}
            worst = max(lo.items(), key=lambda kv: kv[1], default=(("-", "-"), 0.0))
            print("%-9s %-7s %9.2e %9.2e %11s %9.1e" % (
                r["mol"], side, sum(v * v for v in lo.values()) ** 0.5,
                sum(v * v for v in hi.values()) ** 0.5,
                ("%s>%s %.0e" % (worst[0][0], worst[0][1], worst[1])) if worst[1] > CASC_ZERO
                else "-", s["cond"]))
    #The separation, printed before the cut is used: last time a threshold was set below the
    #measurement floor and refused all eight, so a cut is only a gate if the two populations it
    #divides sit either side of it with room to spare.
    filled = [(v, r["mol"], side, k) for r in ok for side, s in r["sides"].items()
              for k, v in s["blocks"].items() if k[0] != k[1] and v >= CASC_ZERO]
    empty = [(v, r["mol"], side, k) for r in ok for side, s in r["sides"].items()
             for k, v in s["blocks"].items() if k[0] != k[1] and v < CASC_ZERO]
    if empty and filled:
        hi, lo = max(empty), min(filled)
        print("\nseparation: the largest mass called empty is %.1e (%s %s %s>%s), the smallest called"
              % (hi[0], hi[1], hi[2], hi[3][0], hi[3][1]))
        print("filled is %.1e (%s %s %s>%s).  The cut at %g lies between two populations %.0f decades"
              % (lo[0], lo[1], lo[2], lo[3][0], lo[3][1], CASC_ZERO, np.log10(lo[0] / hi[0])))
        print("apart, which is what makes it a gate and not a floor.")
    print("\nthe same metric on both sides, without pairing a column: max|dEig(G)| over the pre-NAO")
    print("Gram matrices, whose spectrum is invariant under the permutation and the signs that are")
    print("the gauge here.  T^T G T = 1 on both sides, so agreement here leaves no magnitude")
    print("difference between the two transforms hiding in G.")
    for r in ok:
        print("  %-9s max|dEig| %.1e   |dTrace| %.1e   (n = %d)" % (
            r["mol"], r["deig"], r["dtrace"], r["n"]))
    print("\nlocality of the cross-class mass: same-(atom, l) / cross-(atom, l), share of ||T||.")
    print("This was pre-registered as a discriminator - a per-(atom, l) re-diagonalisation fills only")
    print("the first column, a molecular Schmidt projection also the second - and the measurement")
    print("REFUTED that reading of it: step 4 mixes step-3 OUTPUTS, which are already molecular, so")
    print("one intra-atomic rotation of a delocalised valence NAO puts other atoms' valence pre-NAOs")
    print("into a core NAO.  Both columns are therefore consistent with step 4 alone and the split is")
    print("reported for the record, not as evidence.")
    print("%-9s %-7s %-19s %-19s %-19s" % ("mol", "side", "Val>Cor same/cross",
                                           "Ryd>Cor same/cross", "Ryd>Val same/cross"))
    for r in ok:
        for side, s in sorted(r["sides"].items()):
            cells = []
            for k in (("Val", "Cor"), ("Ryd", "Cor"), ("Ryd", "Val")):
                sm, cr = s["loc"].get(k, (float("nan"),) * 2)
                cells.append("%8.1e /%8.1e" % (sm, cr))
            print("%-9s %-7s %s" % (r["mol"], side, " ".join(cells)))
    print("\nsymmetrizability on the within-class, cross-(atom, l) pairs above %g.  resid is in log" %
          CASC_FLOOR)
    print("units, so exp(resid) is the worst ratio error left after the best diagonal; `plain` is the")
    print("same pairs with no diagonal allowed, and a small resid next to a LARGE plain is the")
    print("finding - both small would only mean T is nearly symmetric outright.")
    print("%-9s %-7s %6s %9s %9s %7s %10s" % (
        "mol", "side", "pairs", "resid", "plain", "signbad", "same-(a,l)"))
    for r in ok:
        for side, s in sorted(r["sides"].items()):
            q = s["sym"]
            print("%-9s %-7s %6d %9s %9s %7d %10s" % (
                r["mol"], side, q["n"],
                "-" if q["resid"] is None else "%.2e" % q["resid"],
                "-" if q["plain"] is None else "%.2e" % q["plain"], s["frus"],
                "-" if s["sym_same"]["resid"] is None else "%.2e" % s["sym_same"]["resid"]))
    rows, differ = [], []
    for r in ok:
        en, sn = cascade_shape(r["sides"]["native"])
        eg, sg = cascade_shape(r["sides"]["gennbo"])
        rows.append((r["mol"], sorted(en ^ eg), sn, sg))
        if en != eg or sn != sg:
            differ.append(r["mol"])
    print("\nshape per molecule: the empty-block set and diagonal symmetrizability, both")
    print("gauge-invariant.  `blocks differing` is the symmetric difference of the two empty sets.")
    for mol, diffblocks, sn, sg in sorted(rows):
        print("  %-9s symmetrizable native %-5s gennbo %-5s  blocks differing: %s" % (
            mol, sn, sg, ", ".join("%s>%s" % b for b in diffblocks) or "none"))
    #TWO THINGS NOT TO QUOTE, both of the shape the lane has been burnt by.
    #  1. the gates above are near-tautological.  Cpre is square and well conditioned (cond ~ 10), so
    #     T = Cpre^-1 Cnao reproduces Cnao by construction, and T^T G T = 1 reduces to
    #     Cnao^T S Cnao = 1, which is read_lfn33's own layout test.  Passing them says the matrices
    #     are full rank and S-orthonormal, nothing about the cascade.  The one-sided check that is
    #     NOT vacuous was already paid for at the level below: each side's pre-NAOs satisfy the
    #     pre-NAO definition (native 4.91e-10, gennbo 5.78e-09), and it is what says the dumped Cpre
    #     is a pre-NAO set at all.  It also bounds the ambiguity that is left: any other proper
    #     pre-NAO set differs by a rotation inside one (atom, l) block, which cannot create or remove
    #     content on a DIFFERENT atom.
    #  2. ||T|| agrees between the sides to every printed digit on all eight, and that is worth
    #     exactly nothing: T^T G T = 1 with square T gives T T^T = G^-1 identically, so
    #     ||T||^2 = tr(G^-1) is fixed by the metric the two sides already share.  It is the same
    #     pattern as the NRT doubling - a quantity fixed by construction measures the construction -
    #     and it is recorded here so nobody quotes it later.  The identity was CHECKED and not just
    #     argued: max|T T^T - G^-1| is 3.2e-10 to 1.5e-6 over the sixteen sides (at each side's own
    #     printing floor, gennbo's nine decimals being the looser one) and ||T||^2 reproduces tr(G^-1)
    #     in every printed digit, e.g. lif 59.652594 and benzene 145994.7196.  The same identity says what the object
    #     of interest actually is: T = G^-1/2 O with O orthogonal, so the entire cascade is one
    #     orthogonal matrix per side and every difference between them lives in O.
    core = [("Val", "Cor"), ("Ryd", "Cor")]
    core_n = max(max(r["sides"]["native"]["blocks"][k] for k in core) for r in ok)
    core_g = max(max(r["sides"]["gennbo"]["blocks"][k] for k in core) for r in ok)
    print("\nREAD PER PARTITION, because the two partitions do not give the same answer.")
    print("\ncore partition - OUTCOME E on %d of %d, and it is named by the block." % (
        len(differ), len(ok)))
    print("gennbo's core NAOs are built from core pre-NAOs ALONE: Val>Cor and Ryd>Cor at most %.1e" %
          core_g)
    print("over all eight, which is its own numerical zero.  Native's carry %.1e at most and %.1e at" %
          (core_n, min(min(r["sides"]["native"]["blocks"][k] for k in core) for r in ok)))
    print("least, five to six decades above the same floor, on 8 of 8.  So the priority partition")
    print("itself differs: native lets non-core pre-NAOs into a core NAO and gennbo does not.  The")
    print("only step in nao.cpp that can cross classes is step 4, whose blocks are keyed on")
    print("(atom, l) with no class in the key - that is a reading of the code, not a measurement, and")
    print("it is already falsifiable: NAO_CLASS_SPLIT=1 is exactly the arm that empties these two")
    print("blocks, and it is ALREADY measured and rejected (pf5/so2/sf6 go 0.93/0.94/0.98 ->")
    print("1.12/1.05/1.19 x gennbo).  So this is a real, located, structural difference that the")
    print("acceptance gate has already refused as a fix - a second defect, not the one being hunted.")
    print("It cannot be the leak by magnitude either: 1e-4 of ||T|| against 0.5 e of moved charge.")
    print("\nvalence/Rydberg partition - NO OUTCOME.  This is the partition the leak lives in, and")
    print("the instrument does not resolve it.  Ryd>Val is filled on BOTH sides (%.1e to %.1e), so" % (
        min(s["blocks"][("Ryd", "Val")] for r in ok for s in r["sides"].values()),
        max(s["blocks"][("Ryd", "Val")] for r in ok for s in r["sides"].values())))
    print("neither side is triangular there and the difference is magnitude, not shape; and the")
    print("symmetry measure refuses both sides at resid 1.9 to 6.6 log units, which it MUST, because")
    print("both codes end with a per-(atom, l) re-diagonalisation and a right rotation destroys the")
    print("W^-1 M signature it was built to detect.  A measure that refuses both sides by construction")
    print("discriminates nothing, so nothing is quoted from it - the resolution does not pass, and")
    print("that is the result.")
    print("\nthe one thing in this arm whose RANKING tracks the failure, quoted as a correlation and")
    print("nothing more - it is not gated, and pf5/so2/sf6 are not flat in it, which is the")
    print("discriminator the acceptance test uses:")
    ratios = sorted(((r["sides"]["native"]["blocks"][("Ryd", "Val")] /
                      r["sides"]["gennbo"]["blocks"][("Ryd", "Val")], r["mol"]) for r in ok),
                    reverse=True)
    print("  Ryd>Val native/gennbo:  " + "  ".join("%s %.2f" % (m, v) for v, m in ratios))
    print("  the top two are the two molecules whose final Rydberg excess is worst; below them the")
    print("  ordering does not track it, so this licenses a next measurement and no claim.")


def main_cascade(root):
    print("aonao_compare --cascade on %s, 1 thread, root %s" % (socket.gethostname(), root))
    results, void, missing = [], [], []
    for mol in sorted(os.listdir(root)):
        d = os.path.join(root, mol)
        if not os.path.isdir(d):
            continue
        if not all(os.path.isfile(os.path.join(d, mol + e)) for e in (".32", ".33")):
            missing.append(mol)
            continue
        try:
            r = cascade_sides(mol, d)
        except AssertionError as e:
            void.append((mol, str(e)))
            print("%-10s VOID: %s" % (mol, e))
            continue
        print("%-10s n=%-4d lfn32 %s lfn33 %s  ||T|| native %.1f gennbo %.1f  diag-dominant cols "
              "%d/%d, %d/%d  class relabels %d" % (
                  mol, r["n"], r["lay32"], r["lay33"], r["sides"]["native"]["norm"],
                  r["sides"]["gennbo"]["norm"], r["sides"]["native"]["dom"], r["n"],
                  r["sides"]["gennbo"]["dom"], r["n"], r["recls"]))
        results.append(r)
    print("\ndenominator: %d molecules measured, %d VOID, %d without both matrices" % (
        len(results), len(void), len(missing)))
    for mol, why in void:
        print("  VOID %-10s %s" % (mol, why))
    if missing:
        print("  incomplete: %s" % " ".join(missing))
    if results:
        cascade_verdict(results)
    return 0


#------------------------------------------------------------------------------------------------
#--space: the cascade as ONE orthogonal matrix per side, compared where a comparison is possible.
#
#--cascade ended with T = G^-1/2 O, so the whole orthogonalisation cascade is one orthogonal matrix
#per side and every difference between the sides lives in O.  This arm compares the two O's.
#
#WHAT WAS PROPOSED FOR THE FIRST GATE, AND WHY IT CANNOT BE ONE.  The singular values of O_n^T O_g P,
#"all 1.000000 iff they differ by nothing", are all 1.000000 FULL STOP: O_n, O_g and P are each
#orthogonal, a product of orthogonal matrices is orthogonal, and an orthogonal matrix has every
#singular value equal to 1 identically.  The principal angles between the two FULL spaces vanish for
#the same reason - both sides span the whole AO space, so completeness fixes them.  This is the
#||T||^2 = tr(G^-1) trap one level up, and it is DEMONSTRATED rather than argued: the very matrix
#whose singular values come out 1.000000 is the matrix whose class blocks disagree, the value does
#not move when the permutation is replaced by a deliberately wrong one, and demo() rotates valence
#into Rydberg by 30 degrees and still gets 1.000000 out of it.  So it is computed, printed as VOID,
#and never quoted.
#
#WHAT IS MEASURABLE.  A subspace has no column order and no column signs:
#
#    U = Cnao_n^T S Cnao_g           both sides' NAOs in the shared AO basis, S from the shared .47
#    mass[C, C'] = ||U[C, C']||_F^2  dimensions of native's class-C space lying in gennbo's C' space
#
#plus the principal angles between span(Cnao_n[:, C]) and span(Cnao_g[:, C]) per class, and per
#(atom, class) for the "where".  None of this needs the permutation between native's
#(atom, l, shell, m) order and gennbo's component-major order, and none of it needs a sign rule: the
#permutation is not resolved here, it is ELIMINATED.  That is strictly better than pinning it,
#because a quantity that never enters the arithmetic cannot be fitted to make the answer nicer.  The
#class membership per (atom, l) is asserted identical first, or "native's Val space" and "gennbo's
#Val space" are not the same question.
#
#WHAT IS FIXED BY CONSTRUCTION IN THE MASS TABLE, stated before it is read so nobody quotes it as a
#finding: U is orthogonal, so every row block AND every column block of mass sums to that class's
#dimension.  With the cores in near-agreement that FORCES mass[Val, Ryd] ~ mass[Ryd, Val] - the two
#directions of the valence/Rydberg disagreement are near-equal by the sum rule, not by measurement.
#One number carries the content: d_VR = mass[Val, Ryd], the dimensions of valence space the two sides
#disagree about.  It is a subspace mass and NOT electrons; converting it needs occupancies.
#
#THE FLOOR IS MEASURED, AND IT IS THE SQUARE ROOT OF THE COEFFICIENT FLOOR.  A principal sine is
#sqrt(1 - sigma**2) ~ sqrt(2 * (1 - sigma)), so a deficit of delta in a cosine shows up as sqrt(2
#delta) in the angle: the sqrt AMPLIFIES roundoff, and a decade of printed precision buys only half a
#decade of resolvable angle.  This is not theory - the self-comparison in demo() came out at 4e-08 on
#coefficients good to 1e-15 and failed a 1e-12 assertion, which is how the factor was found.  So the
#floor is sqrt(2 * max|C^T S C - 1|) over both sides and a class counts as different only above
#SPACE_SEP times it.  Anyone reading gennbo's nine decimals as a 1e-09 floor here would be a factor
#of ~10000 too optimistic and would call noise a leak.
#A self-comparison is NOT the way to measure the floor (the sines of a space against itself are zero
#by construction); it is a check on this code, and it is labelled as one in demo().
#
#PRE-REGISTERED, before any number was read:
#  G  all class subspaces agree to the floor.  Then the cascade produces the same core, valence and
#     Rydberg SPACES, the leak is not a subspace difference at all, and the working hypothesis that
#     the orthogonalisation cascade is at fault is REFUTED - it would have to be the occupancy and
#     m-averaging step downstream of it.
#  H  the valence and/or Rydberg subspaces differ measurably while the cores agree.  The cascade's
#     valence/Rydberg split is then located as a subspace difference and d_VR quantifies it; pf5,
#     so2 and sf6 must then be checked for flatness, since those three are what the acceptance test
#     discriminates on.
#  I  the core subspaces differ too.  Corroborates the core defect --cascade already located and adds
#     nothing to it - reported as corroboration, never as a second finding.
#  J  the gate fails: the instrument moves under something it must ignore, or fails to move under an
#     injected leak.  NO OUTCOME, quote nothing, exactly as --cascade's valence/Rydberg half did.
SPACE_SEP = 10.0


def cascade_orthogonal(Cpre, Cnao, S):
    """O = G^1/2 T, the whole cascade as one orthogonal matrix, with max|O^T O - 1| beside it.

    >>> rng = np.random.default_rng(3)
    >>> A = rng.normal(size=(6, 6)); S = A @ A.T + 6 * np.eye(6)
    >>> Cp = np.linalg.inv(np.linalg.cholesky(S)).T        # any S-orthonormal set, so G = 1
    >>> Q, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    >>> O, res = cascade_orthogonal(Cp, Cp @ Q, S)
    >>> bool(res < 1e-12), bool(np.abs(O - Q).max() < 1e-10)
    (True, True)
    """
    w, V = np.linalg.eigh(Cpre.T @ S @ Cpre)
    T, _, _, _ = np.linalg.lstsq(Cpre, Cnao, rcond=None)
    O = (V * np.sqrt(w)) @ V.T @ T
    return O, float(np.abs(O.T @ O - np.eye(O.shape[1])).max())


def principal_sines(A, B, S):
    """sin of the principal angles between span(A) and span(B), both S-orthonormal, descending.

    Zero iff the two spans coincide.  Neither a column order nor a column sign enters, which is the
    whole reason this comparison needs no permutation between the sides and no sign convention.

    >>> S6 = np.eye(6)
    >>> A = np.eye(6)[:, :2]
    >>> float(principal_sines(A, A[:, ::-1], S6).max())     # a reordering is the SAME space
    0.0
    >>> B = np.column_stack([A[:, 0] * np.cos(0.3) + np.eye(6)[:, 3] * np.sin(0.3), A[:, 1]])
    >>> bool(abs(principal_sines(A, B, S6).max() - np.sin(0.3)) < 1e-12)
    True
    """
    sv = np.linalg.svd(A.T @ S @ B, compute_uv=False)
    return np.sqrt(np.clip(1.0 - sv ** 2, 0.0, None))[::-1]


def class_mass(Cn, Cg, S, ncls, gcls, order=ORDER):
    """mass[(C, C')] = dimensions of native's class-C NAO space lying in gennbo's class-C' space.

    Row AND column blocks sum to the class dimension because U is orthogonal, so a row sum is not a
    measurement, and the near-equality of mass[Val, Ryd] with mass[Ryd, Val] is forced rather than
    found - the second doctest is one Val/Ryd Givens rotation and both come out 0.64.

    >>> S4, C = np.eye(4), np.eye(4)
    >>> cls = ["Val", "Val", "Ryd", "Ryd"]
    >>> m = class_mass(C, C, S4, cls, cls)
    >>> round(m[("Val", "Val")], 12), m[("Val", "Ryd")]
    (2.0, 0.0)
    >>> g = np.array([[0.6, 0, -0.8, 0], [0, 1, 0, 0], [0.8, 0, 0.6, 0], [0, 0, 0, 1.0]])
    >>> m = class_mass(C, g, S4, cls, cls)
    >>> round(m[("Val", "Ryd")], 6), round(m[("Ryd", "Val")], 6)
    (0.64, 0.64)
    """
    U = Cn.T @ S @ Cg
    idx_n = {c: [i for i, k in enumerate(ncls) if k == c] for c in order}
    idx_g = {c: [i for i, k in enumerate(gcls) if k == c] for c in order}
    return {(a, b): float(np.linalg.norm(U[np.ix_(idx_n[a], idx_g[b])]) ** 2)
            if idx_n[a] and idx_g[b] else 0.0 for a in order for b in order}


def assert_same_classes(lnao, naos):
    """Both sides must put the same number of orbitals of each (atom, l) into each class.

    Counted per (atom, l) because that is invariant under the two sides' differing column orders, so
    no pairing is needed.  Without it, a block or a space called Val on one side is not the same
    question as the one called Val on the other and every number below would be about the labels.
    """
    cnt_n = collections.Counter((l[1], l[2], CLASS[l[5]]) for l in lnao)
    cnt_g = collections.Counter((e["atom"] - 1, LANG_L[e["lang"][0].lower()], e["type"])
                                for e in naos)
    assert cnt_n == cnt_g, ("native and gennbo classify differently per (atom, l): %s - a space "
                            "called Val on one side is then not the same space on the other"
                            % ((cnt_n - cnt_g) + (cnt_g - cnt_n)))


def space_sides(mol, d):
    """Both sides' class subspaces compared in the shared AO basis, plus the VOID metric."""
    n, S, P = read_47(os.path.join(d, mol + ".47"))
    Cg, lay33, e33 = read_lfn33(os.path.join(d, mol + ".33"), n, S)
    lnao, Cn = read_naoc(os.path.join(d, mol + ".naoc.txt"))
    p = os.path.join(d, mol + ".aonao.nbo.json")
    naos = sorted((load_nbo(p) if load_nbo is not None else json.load(open(p)))["nao"],
                  key=lambda e: e["index"])
    assert_same_classes(lnao, naos)
    ncls = [CLASS[l[5]] for l in lnao]
    gcls = [e["type"] for e in naos]
    natom = [l[1] for l in lnao]
    gatom = [e["atom"] - 1 for e in naos]
    e_n = float(np.abs(Cn.T @ S @ Cn - np.eye(n)).max())
    #sqrt, not the residual itself: a principal sine is sqrt(1 - sigma**2), so a cosine deficit of
    #delta surfaces as sqrt(2 delta) in the angle.  See the header - the factor was found by a failed
    #assertion, not assumed.
    floor = float(np.sqrt(2.0 * max(e_n, e33)))
    cls_sines, worst_atom = {}, {}
    for c in ORDER:
        a = Cn[:, [i for i, k in enumerate(ncls) if k == c]]
        b = Cg[:, [i for i, k in enumerate(gcls) if k == c]]
        cls_sines[c] = float(principal_sines(a, b, S).max()) if a.shape[1] else 0.0
        #"where", resolved by atom and still pairing-free: one atom's class space per side.
        per = []
        for at in sorted(set(natom)):
            ia = [i for i, k in enumerate(ncls) if k == c and natom[i] == at]
            ib = [i for i, k in enumerate(gcls) if k == c and gatom[i] == at]
            if ia and len(ia) == len(ib):
                per.append((float(principal_sines(Cn[:, ia], Cg[:, ib], S).max()), at))
        worst_atom[c] = max(per) if per else (0.0, -1)
    mass = class_mass(Cn, Cg, S, ncls, gcls)
    dims = {c: sum(1 for k in ncls if k == c) for c in ORDER}
    #The VOID metric, computed so it can be shown to be void rather than asserted to be.  A
    #deliberately WRONG permutation is used for the second number: if the singular values do not move
    #when the pairing is scrambled, the pairing was never being tested.
    Cn_pre = read_naoc(os.path.join(d, mol + ".naocpre.txt"), tag="NAOCPRE")[1]
    Cg_pre = read_lfn32(os.path.join(d, mol + ".32"), n, S)[0]
    On, res_n = cascade_orthogonal(Cn_pre, Cn, S)
    Og, res_g = cascade_orthogonal(Cg_pre, Cg, S)
    wrong = np.random.default_rng(0).permutation(n)
    sv_id = np.linalg.svd(On.T @ Og, compute_uv=False)
    sv_wrong = np.linalg.svd(On.T @ Og[:, wrong], compute_uv=False)
    return dict(mol=mol, n=n, floor=floor, e_n=e_n, e_g=e33, dims=dims, sines=cls_sines,
                worst_atom=worst_atom, mass=mass, res_n=res_n, res_g=res_g,
                void_sv=(float(sv_id.min()), float(sv_id.max())),
                void_sv_wrong=(float(sv_wrong.min()), float(sv_wrong.max())))


def space_verdict(rows):
    print("\n%s" % ("=" * 96))
    print("the VOID metric first, so it cannot be mistaken for the result.  sigma(O_n^T O_g) is")
    print("bounded below and above by:")
    lo = min(r["void_sv"][0] for r in rows)
    hi = max(r["void_sv"][1] for r in rows)
    wlo = min(r["void_sv_wrong"][0] for r in rows)
    whi = max(r["void_sv_wrong"][1] for r in rows)
    print("  correct pairing  %.12f .. %.12f" % (lo, hi))
    print("  WRONG pairing    %.12f .. %.12f   <- scrambling the pairing changes nothing," % (
        wlo, whi))
    print("     which is the proof that the pairing was never being tested.  A product of orthogonal")
    print("     matrices is orthogonal and every singular value of an orthogonal matrix is 1, so this")
    print("     metric cannot fail and is quoted nowhere.  demo() rotates valence into Rydberg by 30")
    print("     degrees and still gets 1.000000 from it.")
    print("\nthe measurable part.  floor = each side's own max|C^T S C - 1|; a class counts as")
    print("different only above %.0f x floor." % SPACE_SEP)
    print("\n%-9s %9s %9s %9s %9s %7s %7s %7s" % (
        "mol", "floor", "sin Cor", "sin Val", "sin Ryd", "dVR", "dim Val", "dVR/dim"))
    for r in rows:
        print("%-9s %9.1e %9.2e %9.2e %9.2e %7.3f %7d %7.3f" % (
            r["mol"], r["floor"], r["sines"]["Cor"], r["sines"]["Val"], r["sines"]["Ryd"],
            r["mass"][("Val", "Ryd")], r["dims"]["Val"],
            r["mass"][("Val", "Ryd")] / max(r["dims"]["Val"], 1)))
    #A THIRD quantity fixed by construction, and the one most likely to be mis-read off the table
    #above: sin Val and sin Ryd are EQUAL, and have to be.  The complementarity theorem says the
    #nonzero principal angles between two equal-dimensional subspaces equal those between their
    #orthogonal complements; inside the Val+Ryd space, Ryd IS Val's complement on each side, and
    #the class assertion already guarantees both sides give Val the same dimension.  So the two
    #columns are ONE measurement, not two agreeing ones, and "Val 8/8 Ryd 8/8" below is one fact
    #reported twice.  Printed as a measured difference rather than asserted, because the cores are
    #only NEARLY common - the residual is what the core defect leaks into it.
    print("\nforced by complementarity and therefore ONE measurement, not two: sin Val vs sin Ryd")
    print("  " + "  ".join("%s %.1e" % (r["mol"], abs(r["sines"]["Val"] - r["sines"]["Ryd"]))
                           for r in rows))
    print("  (Ryd is Val's complement inside Val+Ryd, so these are equal by theorem.  The residual")
    print("   is bounded by the core disagreement, which is why it is not exactly zero.)")
    print("\nforced by the sum rule and therefore not a finding: mass[Val, Ryd] vs mass[Ryd, Val]")
    for r in rows:
        print("  %-9s %.6f vs %.6f   (difference %.2e)" % (
            r["mol"], r["mass"][("Val", "Ryd")], r["mass"][("Ryd", "Val")],
            abs(r["mass"][("Val", "Ryd")] - r["mass"][("Ryd", "Val")])))
    print("\nwhere, resolved by atom - the worst atom of each class subspace, pairing-free:")
    for r in rows:
        print("  %-9s %s" % (r["mol"], "  ".join(
            "%s worst atom %d at %.2e" % (c, r["worst_atom"][c][1], r["worst_atom"][c][0])
            for c in ORDER)))
    #The pre-registered read, applied in the order it was written.
    diff = {c: [r for r in rows if r["sines"][c] > SPACE_SEP * r["floor"]] for c in ORDER}
    print("\nclasses called different at %.0f x floor: %s" % (
        SPACE_SEP, "  ".join("%s %d/%d" % (c, len(diff[c]), len(rows)) for c in ORDER)))
    if not any(diff[c] for c in ORDER):
        print("OUTCOME G: every class subspace agrees to the floor on every molecule.  The cascade")
        print("produces the SAME core, valence and Rydberg spaces, so the leak is not a subspace")
        print("difference and the hypothesis that the orthogonalisation cascade is at fault is")
        print("REFUTED - it has to be the occupancy/m-averaging step downstream of it.")
    else:
        which = "H" if not diff["Cor"] else "H and I"
        print("OUTCOME %s: the class subspaces themselves differ, so the two cascades do not even"
              % which)
        print("agree about WHICH space is valence and which is Rydberg.")
        print("d_VR = the dimensions of valence space the two sides disagree about: %.3f to %.3f" % (
            min(r["mass"][("Val", "Ryd")] for r in rows),
            max(r["mass"][("Val", "Ryd")] for r in rows)))
        flat = [(r["mol"], r["mass"][("Val", "Ryd")] / max(r["dims"]["Val"], 1)) for r in rows]
        flat.sort(key=lambda t: -t[1])
        print("d_VR per valence dimension, the form that can be compared across molecules:")
        print("  " + "  ".join("%s %.4f" % t for t in flat))
        print("the acceptance test discriminates on pf5/so2/sf6 staying flat, so those three are")
        print("what a candidate fix has to leave alone; they are %s here." % ", ".join(
            "%s %.4f" % t for t in flat if t[0] in ("pf5", "so2", "sf6")))
        if diff["Cor"]:
            print("the core half is OUTCOME I - corroboration of the core defect --cascade already")
            print("located, and NOT a second finding.")
        print("(outcome letter: %s)" % which)


def main_space(root):
    print("aonao_compare --space on %s, 1 thread, root %s" % (socket.gethostname(), root))
    rows, void, missing = [], [], []
    for mol in sorted(os.listdir(root)):
        d = os.path.join(root, mol)
        if not os.path.isdir(d):
            continue
        if not all(os.path.isfile(os.path.join(d, mol + s))
                   for s in (".33", ".naoc.txt", ".32", ".naocpre.txt")):
            missing.append(mol)
            continue
        try:
            r = space_sides(mol, d)
        except AssertionError as e:
            void.append((mol, str(e)))
            print("%-10s VOID: %s" % (mol, e))
            continue
        print("%-10s n=%-4d floor %.1e (native %.1e / gennbo %.1e)  O^T O - 1: %.1e / %.1e" % (
            mol, r["n"], r["floor"], r["e_n"], r["e_g"], r["res_n"], r["res_g"]))
        rows.append(r)
    print("\ndenominator: %d molecules measured, %d VOID, %d incomplete" % (
        len(rows), len(void), len(missing)))
    if missing:
        print("  incomplete: %s" % " ".join(missing))
    if rows:
        space_verdict(rows)
    return 0


def compare(mol, d, verbose=False):
    n, S, P = read_47(os.path.join(d, mol + ".47"))
    Cg, layout, ortho = read_lfn33(os.path.join(d, mol + ".33"), n, S)
    labels, Cn = read_naoc(os.path.join(d, mol + ".naoc.txt"))
    assert Cn.shape == (n, n), "native dump is %s, .47 says %d" % (Cn.shape, n)
    #Through the refusing loader, not json.load: this file is a reference like any other, and the
    #rule on this branch is that a stale one is refused rather than reinterpreted.  The cluster tree
    #was at 99efd0e, which predates the parser_version stamp, so nbo_run.{cpp,h} were synced with
    #nao.cpp before the staging job ran - without that, an unstamped v1 file would have been read
    #here and only the closed-shell NAO occupancies happening to be unaffected would have saved it.
    run_path = os.path.join(d, mol + ".aonao.nbo.json")
    run = load_nbo(run_path) if load_nbo is not None else json.load(open(run_path))
    naos = run["nao"]

    #The NAO-basis density of an S-orthonormal set is C^T S P S C, because C^-1 = C^T S - NOT
    #C^T P C, which is what the first version of this check used.  The two are only equal for an
    #orthogonal AO basis, so rather than assert the algebra both are computed and the one that
    #reproduces the printed table is reported: if the wrong one had been used the validation would
    #have failed loudly instead of passing on a coincidence.
    printed = np.array([e["occupancy"] for e in naos])
    #The full matrix, not just its diagonal: the off-diagonal elements are the coherence between two
    #gennbo NAOs, and they are exactly what separates the signed prediction below from an estimate.
    Dg = Cg.T @ S @ P @ S @ Cg
    forms = {"C^T S P S C": np.diag(Dg), "C^T P C": np.diag(Cg.T @ P @ Cg)}
    occ_form = min(forms, key=lambda k: np.abs(forms[k] - printed).max())
    occ_err = float(np.abs(forms[occ_form] - printed).max())

    ref_err = None
    if load_nbo is not None:
        ref_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", mol,
                                mol + ".gennbo.nbo.json")
        if os.path.exists(ref_path):
            ref = load_nbo(ref_path)
            ro = np.array([e["occupancy"] for e in ref["nao"]])
            ref_err = float(np.abs(ro - printed).max()) if ro.shape == printed.shape else float("inf")

    #The native side's NAO-basis density, by the same identity: needed whole because the block
    #spectra are taken from sub-blocks of it at several groupings.
    Dn = Cn.T @ S @ P @ S @ Cn
    gb = gennbo_blocks(naos)
    #Which atom and which l each GENNBO column belongs to, so the out-of-block weight can be split
    #into the two things it can be.  Same atom, different l is an intra-atomic cross-l mix, which is
    #what the final table's 92 %-inter-l excess would look like one level down; different atom cannot
    #be that and would move charge between atoms instead.
    g_atom = np.array([e["atom"] - 1 for e in naos])
    g_l = np.array([LANG_L[e["lang"][0].lower()] for e in naos])
    rows = []
    for (atom, l), gshells in sorted(gb.items()):
        nm = 2 * l + 1
        nrows = [i for i, lab in enumerate(labels) if lab[1] == atom and lab[2] == l]
        nrows.sort(key=lambda i: (labels[i][4], labels[i][3]))  # shell rank, then m
        assert len(nrows) == nm * len(gshells), "block (%d,%d): %d native vs %d gennbo" % (
            atom, l, len(nrows), nm * len(gshells))
        nshells = [nrows[a * nm:(a + 1) * nm] for a in range(len(gshells))]
        M, total, w, O3 = mixing(Cn, Cg, S, nshells, gshells)
        #"Spectrum" has to mean eigenvalues, not each side's reported diagonal.  The density
        #restricted to a block's own subspace is C_block^T S P S C_block, and its eigenvalues are
        #the populations a perfect diagonalisation of THAT subspace would report - so:
        #  blockeig  the two sides' eigenvalue sets, order-free.  Differing means the two subspaces
        #            do not hold the same density, and no choice of vectors inside the block can
        #            repair it.
        #  self_n    each side's reported diagonal against its OWN eigenvalues.  Non-zero means that
        #            side did not diagonalise its block, which is a different defect and is
        #            localised without reference to the other side.
        spectra = block_spectra(
            Dn, Dg, nshells, gshells,
            [CLASS[labels[sh[0]][5]] for sh in nshells],
            [naos[sh[0]]["type"] if naos[sh[0]]["type"] in CLASS.values() else "Ryd"
             for sh in gshells])
        same_other_l = (g_atom == atom) & (g_l != l)
        other_atom = g_atom != atom
        #Per-shell populations on each side, m-averaged, so a degenerate PAIR inside this block can
        #be recognised.  Both sides' own reported numbers, read from the fields each side printed -
        #no new matrix and nothing inferred from an offset.
        nocc_sh = [sum(labels[i][6] for i in sh) / nm for sh in nshells]
        gocc_sh = [float(printed[gs].sum()) / nm for gs in gshells]
        for a in range(M.shape[0]):
            lab = labels[nshells[a][0]]
            g = naos[gshells[a][0]]
            order = np.argsort(M[a])[::-1]
            worst_b = int(order[0])
            #Same rule as compare_pre above, applied one level later: a shared m-averaged block
            #eigenvalue fixes an EIGENSPACE and nothing inside it, so for two shells of equal
            #population `M[a][a]` charges one eigensolver's arbitrary choice as an error.  The
            #invariant is the projector onto that eigenspace.  BOTH sides must be degenerate - if
            #only one is, the two are not describing the same indeterminacy and the disagreement is
            #real.  Where a shell is alone, degsum IS diag, which is why both are carried.
            deg = [b for b in range(M.shape[1])
                   if abs(nocc_sh[b] - nocc_sh[a]) < DEG_TOL
                   and abs(gocc_sh[b] - gocc_sh[a]) < DEG_TOL]
            #Gap to the runner-up: a pairing that is merely close is worth knowing about, and a
            #systematic off-by-one has a LARGE gap on the wrong column, not a small one.
            gap = float(M[a, order[0]] - M[a, order[1]]) if M.shape[1] > 1 else float("nan")
            #The three column sets must be disjoint and exhaustive, or the split is bookkeeping
            #fiction: own (atom,l) block + same atom other l + other atom = everything = total.
            part = float(M[a].sum() + w[a, same_other_l].sum() + w[a, other_atom].sum())
            assert abs(part - float(total[a])) < 1e-9, "block (%d,%d) rank %d: partition %.12f vs total %.12f" % (
                atom, l, a, part, total[a])
            #The signed arm.  A native shell's population is exactly
            #    sum_m sum_{b,b'} O[m][b] O[m][b'] Dg[b][b']
            #because C_g^-1 = C_g^T S, so `predfull` must reproduce the native run's own printed
            #occupancy - a gate on the two matrices and the density together, not an estimate.
            #`pred` keeps only the b = b' terms, which is the same information the unsigned metric
            #has (w times an occupancy) but with the direction left in: it is what native's
            #occupancy WOULD be if it were a weighted average of gennbo's, and its deviation from
            #gennbo's own shell population therefore has a sign.
            pred = nm * float(w[a] @ printed)
            predfull = float(np.einsum("mb,mc,bc->", O3[a], O3[a], Dg))
            rows.append(dict(spectra=spectra, pred=pred, predfull=predfull,
                             occshell=float(sum(labels[i][6] for i in nshells[a])),
                             goccshell=float(printed[gshells[a]].sum()),
                             mol=mol, atom=atom, l=l, rank=a, cls=CLASS[lab[5]],
                             gcls=g["type"], gshell=g["shell"], occ=lab[6], gocc=g["occupancy"],
                             diag=float(M[a, a]), inblock=float(M[a].sum()),
                             degsum=float(M[a, deg].sum()), ndeg=len(deg),
                             crossl=float(w[a, same_other_l].sum()),
                             otheratom=float(w[a, other_atom].sum()),
                             best=worst_b, bestval=float(M[a, worst_b]), gap=gap,
                             total=float(total[a])))
    return dict(mol=mol, n=n, layout=layout, ortho=ortho, occ_err=occ_err, occ_form=occ_form,
                ref_err=ref_err, rows=rows)


def report(res, verbose=False):
    rows = res["rows"]
    leak = [1.0 - r["diag"] for r in rows]
    out = [1.0 - r["inblock"] for r in rows]
    mism = [r for r in rows if r["best"] != r["rank"]]
    print("%-10s n=%-4d lfn33 %s  |C^T S C-1|=%.1e  occ (%s) vs printed %.1e  vs stamped ref %s" % (
        res["mol"], res["n"], res["layout"], res["ortho"], res["occ_form"], res["occ_err"],
        "n/a" if res["ref_err"] is None else "%.1e" % res["ref_err"]))
    print("   shells %d   in-block leak max %.5f mean %.5f   out-of-block max %.5f   rank mismatches %d" % (
        len(rows), max(leak), sum(leak) / len(leak), max(out), len(mism)))
    worst = sorted(rows, key=lambda r: r["diag"])[:5 if not verbose else len(rows)]
    for r in worst:
        print("     atom %2d l=%d rank %d %s/%s %-4s occ %8.5f/%8.5f  M[a][a]=%.5f  best rank %d (%.5f)  in-block %.5f" % (
            r["atom"], r["l"], r["rank"], r["cls"], r["gcls"], r["gshell"], r["occ"], r["gocc"],
            r["diag"], r["best"], r["bestval"], r["inblock"]))


def by_group(all_rows):
    """Two arms, because a rank pairing is not part of the physics.

    `1 - M[a][a]` is the rank-paired leak and it counts an ordering difference as an error; the
    ordering-insensitive arm `1 - max_b M[a][b]` cannot, because it takes the best partner in the
    block whatever its rank.  Where the two differ the shells are merely ordered differently (LiF's
    2p Rydberg pair: 0.00007 on the rank partner, 0.98861 one rank over - the shape is right); where
    both are large the shape itself is wrong and no relabelling recovers it.
    """
    #A FOURTH view, and the one the step-4 degeneracy by-product forces.  `1 - M[a][a]` is not a
    #measurement on a shell whose population is shared with another shell of the same block: the
    #eigenvalue they share fixes their eigenspace and nothing inside it, so which one native calls
    #rank a is its eigensolver's choice.  `1 - degsum` projects onto that eigenspace instead.  It is
    #the rule compare_pre already uses on the pre-NAOs, at the same DEG_TOL, not a new one - and it
    #is NOT a looser version of the rank-paired arm: on a shell that is alone the two are the same
    #number, and the count of shells where they can differ is printed below so the correction's
    #reach is visible rather than asserted.
    for what, val in (("mean 1-M[a][a], rank-paired", lambda r: 1.0 - r["diag"]),
                      ("mean 1-degsum, degeneracy-projected", lambda r: 1.0 - r["degsum"]),
                      ("mean 1-max_b M[a][b], ordering-insensitive", lambda r: 1.0 - r["bestval"]),
                      ("mean weight outside the shell's own (atom,l) block",
                       lambda r: 1.0 - r["inblock"])):
        print("\nwhere the shape error sits (%s, n shells):" % what)
        for key, name in ((lambda r: r["l"], "l"), (lambda r: r["cls"], "class")):
            groups = {}
            for r in all_rows:
                groups.setdefault(key(r), []).append(val(r))
            print("  by %-5s %s" % (name, "  ".join("%s: %.5f (%d)" % (k, sum(v) / len(v), len(v))
                                                    for k, v in sorted(groups.items(), key=str))))
    ind = [r for r in all_rows if r["ndeg"] > 1]
    cls_all = sorted(set(r["cls"] for r in all_rows))
    print("\nshells whose rank pairing is a GAUGE (degenerate on both sides, so 1-M[a][a] is not a"
          " measurement there): %d of %d   %s" % (
              len(ind), len(all_rows),
              "  ".join("%s %d/%d" % (k, sum(1 for r in ind if r["cls"] == k),
                                      sum(1 for r in all_rows if r["cls"] == k)) for k in cls_all)))

    #The metric above normalises every shell to 1 whatever its occupancy, so an empty Rydberg shell
    #with a badly wrong shape counts the same as a doubly-occupied valence shell that is nearly right.
    #Populations do not work that way, and the failure being chased is a population failure, so the
    #same three arms are also summed with (2l+1)*occ as the weight.  It is an estimate: M is
    #m-averaged, so this assumes the m components of a shell share its occupancy.
    print("\non the charge scale (sum (2l+1)*occ*leak, electrons):")
    for what, val in (("rank-paired", lambda r: 1.0 - r["diag"]),
                      ("degeneracy-projected", lambda r: 1.0 - r["degsum"]),
                      ("ordering-insensitive", lambda r: 1.0 - r["bestval"]),
                      ("out-of-block", lambda r: 1.0 - r["inblock"]),
                      ("  of that, same atom other l", lambda r: r["crossl"]),
                      ("  of that, another atom", lambda r: r["otheratom"]),
                      ("in-block, other shell", lambda r: r["inblock"] - r["diag"])):
        groups = {}
        for r in all_rows:
            groups.setdefault(r["cls"], []).append((2 * r["l"] + 1) * r["occ"] * val(r))
        print("  %-29s %s   total %.5f" % (
            what, "  ".join("%s: %.5f" % (k, sum(v)) for k, v in sorted(groups.items())),
            sum(sum(v) for v in groups.values())))


def final_table(mol):
    """(d(Val), worst |dq| per atom) from the stamped reference pair - the FINAL table this is
    supposed to explain.  Read here rather than copied from nao_class_leak's printout, so the
    bridge cannot quietly compare against a number that has since changed."""
    if load_nbo is None:
        return None
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", mol)
    paths = [os.path.join(d, "%s.%s.nbo.json" % (mol, side)) for side in ("gennbo", "native")]
    if not all(os.path.exists(p) for p in paths):
        return None
    g, n = (load_nbo(p)["nao"] for p in paths)
    if len(g) != len(n):
        return None

    def by_atom(table):
        out = {}
        for e in table:
            cls = e["type"] if e["type"] in ("Cor", "Val", "Ryd") else "Ryd"
            out.setdefault(e["atom"], {"Cor": 0.0, "Val": 0.0, "Ryd": 0.0})[cls] += e["occupancy"]
        return out

    gc, nc = by_atom(g), by_atom(n)
    dval = sum(nc[a]["Val"] - gc[a]["Val"] for a in gc)
    dq = max(abs(sum(nc[a][c] - gc[a][c] for c in gc[a])) for a in gc)
    return dval, dq


def bridge(rows_by_mol, final=None):
    """Does the intermediate account for the final table, molecule by molecule?

    The final tables say native's valence set is short by d(Val) and its Rydberg set long by the
    same amount, INTRA-ATOMICALLY - benzene 0.31611 e, ethane 0.12989 e, and pf5/so2/sf6 the other
    way by ~0.02 e.  The intermediate says native's valence NAOs are the wrong shape.  Those are
    only the same finding if the mis-shape is large enough to produce the mis-population, so it is
    put on the same scale and the same molecules:

        intra   sum over VALENCE shells of (2l+1)*occ*(weight on the same atom but a different
                shell or a different l) - the part that can move population between classes of one
                atom without moving charge off it
        inter   the same for weight on ANOTHER atom - which cannot do that, and would show up as a
                charge error instead

    intra is an ESTIMATE of an upper bound, not an upper bound: the population a mis-shaped shell
    actually moves is w*(occ_a - occ_b) to leading order, and occ_b >= 0, so w*occ_a is on the high
    side - but M is m-averaged, so the weight assumes a shell's m components share its occupancy,
    and both classes' rows contribute to d(Val), which is why the Rydberg rows' intra weight is
    printed alongside and the ratio is taken on the sum.  A ratio near 1 on a molecule whose d(Val)
    is tenths of an electron is therefore a real quantitative bridge; a ratio of 0.3-0.8 on one
    whose d(Val) is 0.001-0.02 e is at the estimate's own resolution and says only that the two
    numbers are the same size.  Nothing here can check the SIGN: a mixing weight is positive, and
    pf5/so2/sf6 have native's valence set too LARGE while the other five have it too small.
    """
    print("\nbridge to the final table (electrons):")
    print("  %-10s %10s %10s %10s %7s %10s %10s" % (
        "molecule", "d(Val)", "intra Val", "intra Ryd", "ratio", "inter Val", "worst dq"))
    short, inter_tot, dq_tot = [], 0.0, 0.0
    for mol, rows in sorted(rows_by_mol.items()):
        wt = lambda r: (2 * r["l"] + 1) * r["occ"]
        def intra_of(cls):
            return sum(wt(r) * (r["inblock"] - r["diag"] + r["crossl"])
                       for r in rows if r["cls"] == cls)
        iv, ir = intra_of("Val"), intra_of("Ryd")
        inter = sum(wt(r) * r["otheratom"] for r in rows if r["cls"] == "Val")
        inter_tot += inter
        ft = (final or final_table)(mol)
        if ft is None:
            print("  %-10s %10s %10.5f %10.5f %7s %10.5f %10s" % (
                mol, "n/a", iv, ir, "n/a", inter, "n/a"))
            continue
        dval, dq = ft
        dq_tot = max(dq_tot, dq)
        ratio = (iv + ir) / abs(dval) if dval else float("inf")
        print("  %-10s %+10.5f %10.5f %10.5f %7.2f %10.5f %10.5f" % (
            mol, dval, iv, ir, ratio, inter, dq))
        if ratio < 1.0:
            short.append((mol, ratio, abs(dval)))
    if short:
        print("  -> the intra-atomic mis-shape is SMALLER than the final class error on %s" %
              ", ".join("%s (%.2f, d(Val) %.5f e)" % kv for kv in short))
        big = [t for t in short if t[2] > 0.05]
        print("     %s" % ("all of those have |d(Val)| under 0.05 e, which is the estimate's own "
                           "resolution - not a refutation on its own" if not big else
                           "and %s has |d(Val)| over 0.05 e, which the estimate cannot explain away"
                           % ", ".join(t[0] for t in big)))
    else:
        print("  -> every molecule's final class error fits inside its intermediate mis-shape")
    #The inter-atomic arm has nowhere to go in the class tables, so it must show up as charge - and
    #it does not, by more than an order of magnitude.  That is the same cancellation that made the
    #NPA charge agreement meaningless, measured one level further up.
    print("  inter-atomic valence weight totals %.5f e, worst atomic charge deviation %.5f e "
          "(%.0fx cancellation)" % (inter_tot, dq_tot, inter_tot / dq_tot if dq_tot else 0.0))


def signed(rows_by_mol, final=None):
    """The same leak with its sign left in - a prediction of d(Val), not a bound on |d(Val)|.

    The unsigned arms answer "is the mis-shape big enough".  They cannot answer "does it push the
    right way", because a mixing weight is positive while the final tables split in two directions:
    native's valence set comes out too SMALL on five molecules and too LARGE on pf5/so2/sf6, and
    that split is the discriminator the whole acceptance gate rests on.  Putting gennbo's own
    occupancies through the mixing gives a signed number:

        pred(a)   = (2l+1) * sum_b w[a][b] * occ_gennbo(b)   what native's shell a would hold if it
                                                             were a weighted average of gennbo's
        d_pred    = sum over native VALENCE shells of pred(a) - gennbo's own valence total
        d_actual  = native's own valence total - gennbo's

    Total electrons are conserved exactly in d_pred (the columns of (2l+1)*w sum to 1 by
    completeness of the native set), so the class sums of d_pred add to zero just as the real ones
    do, and a sign is a statement about direction rather than about normalisation.

    Two gates before the sign is read:

      exact   sum over the class of predfull(a) - gennbo's total must EQUAL d_actual, because
              predfull keeps the off-diagonal coherence and is then an identity, not a model.  If
              this fails the matrices or the density are misread and nothing else here counts.
      resid   d_pred - d_actual is exactly the coherence term dropped by w.  It is reported, not
              hidden: it is the part of the class error that the m-averaged, sign-blind metric
              cannot see, and if it dominates then the unsigned bridge was measuring a proxy.
    """
    print("\nsigned prediction (electrons, native - gennbo; + = native's set is LARGER):")
    print("  %-10s %11s %11s %11s %10s %9s %6s" % (
        "molecule", "d(Val) act", "d(Val) pred", "d(Ryd) pred", "coherence", "|exact|", "sign"))
    agree, total, bad_gate = 0, 0, []
    verdict = {}
    for mol, rows in sorted(rows_by_mol.items()):
        def tot(field, key):
            out = {}
            for r in rows:
                out[r[key]] = out.get(r[key], 0.0) + r[field]
            return out
        g = tot("goccshell", "gcls")
        n = tot("occshell", "cls")
        p = tot("pred", "cls")
        pf = tot("predfull", "cls")
        get = lambda d, k: d.get(k, 0.0)
        act = get(n, "Val") - get(g, "Val")
        pred = get(p, "Val") - get(g, "Val")
        pred_ryd = get(p, "Ryd") - get(g, "Ryd")
        exact = abs((get(pf, "Val") - get(g, "Val")) - act)
        if exact > 1e-6:
            bad_gate.append((mol, exact))
        ok = (act > 0) == (pred > 0)
        agree += int(ok)
        total += 1
        verdict[mol] = (act, pred)
        print("  %-10s %+11.5f %+11.5f %+11.5f %+10.5f %9.1e %6s" % (
            mol, act, pred, pred_ryd, pred - act, exact, "ok" if ok else "WRONG"))
    if bad_gate:
        print("  -> EXACTNESS GATE FAILED on %s: predfull does not reproduce the native class total, "
              "so the signed numbers above are not readable" % ", ".join(
                  "%s (%.1e)" % kv for kv in bad_gate))
        return
    print("  exactness gate passes: predfull reproduces every native valence total to < 1e-6 e")
    print("  -> sign reproduced on %d of %d molecules" % (agree, total))
    #The split is the point, not the count: three molecules go the other way in the final tables and
    #they are the acceptance test.  A metric that gets the majority right by getting the majority
    #sign right has said nothing.
    right = [m for m in ALREADY_RIGHT if m in verdict]
    other = [m for m in verdict if m not in ALREADY_RIGHT]
    if right and other:
        a_r = [verdict[m][0] for m in right]
        p_r = [verdict[m][1] for m in right]
        a_o = [verdict[m][0] for m in other]
        p_o = [verdict[m][1] for m in other]
        split_act = all(x > 0 for x in a_r) and all(x < 0 for x in a_o)
        split_pred = all(x > 0 for x in p_r) and all(x < 0 for x in p_o)
        print("     the split that matters: actual %s on %s / %s on the other %d;  predicted %s / %s"
              % ("+" if all(x > 0 for x in a_r) else "mixed", "/".join(right),
                 "-" if all(x < 0 for x in a_o) else "mixed", len(other),
                 "+" if all(x > 0 for x in p_r) else "mixed",
                 "-" if all(x < 0 for x in p_o) else "mixed"))
        if split_act and split_pred:
            print("     -> the bridge is a PREDICTION: the mixing matrix alone reproduces which "
                  "molecules go which way, so a candidate fix can be screened on it without a run")
        elif split_act:
            print("     -> the signed form does NOT reproduce the split, so the m-averaging is "
                  "throwing away something the final table sees: the sign lives in the coherence "
                  "or in the m resolution, and the unsigned weight stays the only usable screen")
        else:
            print("     -> the actual tables do not split that way in this subset; no split to "
                  "reproduce")


def spectrum(rows_by_mol):
    """Is the same-l channel a wrong SPECTRUM or wrong VECTORS?

    0.42780 e of the intra-atomic valence weight goes into another shell of the SAME (atom, l) -
    the one place where the two sides choose among the same candidate functions, so it is decidable
    without a new gennbo job.  Three things are put side by side:

      |d eigen|    the EIGENVALUES of the density restricted to each side's own block, order-free.
                   This is the spectrum proper.  Agreeing says the two sides hold the same density
                   in that block and only what is done inside it can be wrong - a fix belongs in
                   the diagonalisation.  Disagreeing says the block itself holds different density
                   and a fix belongs UPSTREAM of it.
      self         each side's reported occupancies against its OWN eigenvalues.  This decides at
                   which grouping the question may be asked at all, and it is measured rather than
                   assumed: `block_spectra` computes four candidate levels and only a level where
                   BOTH sides come out self-consistent is quoted.  gennbo is NBO 7 itself, so a
                   level where gennbo fails is a wrong guess about the recipe, not a finding.
      |d sorted| / |d rank|   the two sides' REPORTED shell populations, order-free and rank-paired.
                   `sorted` agreeing while `rank` does not would mean the same populations in a
                   different order - a labelling difference rather than a shape error.
    """
    print("\nspectrum or vectors (per (atom,l) block, electrons):")
    print("  %-10s %7s %11s %11s %11s %9s" % (
        "molecule", "blocks", "|d sorted|", "|d rank|", "same-l off", "vs d(Val)"))
    tag, levels, worst_blocks = {}, {}, []
    for mol, rows in sorted(rows_by_mol.items()):
        blocks = {}
        for r in rows:
            blocks.setdefault((r["atom"], r["l"]), []).append(r)
        dspec = drank = offd = 0.0
        for key, rs in blocks.items():
            nat = sorted((r["occshell"] for r in rs), reverse=True)
            gen = sorted((r["goccshell"] for r in rs), reverse=True)
            s = sum(abs(a - b) for a, b in zip(nat, gen))
            d = sum(abs(r["occshell"] - r["goccshell"]) for r in rs)
            o = sum((2 * r["l"] + 1) * r["occ"] * (r["inblock"] - r["diag"]) for r in rs)
            dspec += s
            drank += d
            offd += o
            for name, vals in rs[0]["spectra"].items():
                cur = levels.setdefault(name, [0.0, 0.0, 0.0, 0])
                for i in range(4):
                    cur[i] += vals[i]
            worst_blocks.append((o, mol, key, s, d, len(rs)))
        ft = final_table(mol)
        dval = abs(ft[0]) if ft else None
        print("  %-10s %7d %11.5f %11.5f %11.5f %9s" % (
            mol, len(blocks), dspec, drank, offd,
            "n/a" if not dval else "%.2f" % (dspec / dval)))
        tag[mol] = (dspec, drank, offd, dval)
    print("  blocks with the most same-l mixing:")
    for o, mol, key, s, d, nsh in sorted(worst_blocks, reverse=True)[:6]:
        print("     %-10s atom %2d l=%d  off-diag %.5f  |d sorted| %.5f  |d rank| %.5f  %d shells"
              % (mol, key[0], key[1], o, s, d, nsh))
    spec = sum(t[0] for t in tag.values())
    rank = sum(t[1] for t in tag.values())
    off = sum(t[2] for t in tag.values())
    print("  reported tables: |d sorted| %.5f   |d rank| %.5f   same-l off-diagonal weight %.5f" % (
        spec, rank, off))
    #Which grouping each side diagonalises, measured.  A level where GENNBO is not self-consistent
    #is a wrong guess about the NAO recipe; a level where only native fails would be a finding in
    #its own right, and it is printed either way rather than being hidden by the choice.
    print("  which level each side reproduces as its own eigenvalues (self ~ 0 = readable):")
    for name, _, _ in LEVELS:
        deig, sn, sg, mism = levels[name]
        print("     %-28s |d eigen| %9.5f   self nat %9.5f   self gen %9.5f%s" % (
            name, deig, sn, sg, "" if not mism else "   (%d class-size mismatches)" % mism))
    readable = [n for n, _, _ in LEVELS if max(levels[n][1], levels[n][2]) < 0.01]
    one_sided = [n for n, _, _ in LEVELS if levels[n][2] < 0.01 <= levels[n][1]]
    if one_sided:
        print("  -> gennbo is self-consistent at %s and native is NOT (self nat %.5f e): native is "
              "not diagonalising what NBO diagonalises, which is a finding about nao.cpp and not "
              "about the metric" % (one_sided[0], levels[one_sided[0]][1]))
    if not readable:
        print("  -> NO level has both sides self-consistent, so no spectrum comparison is readable "
              "here: the grouping each side diagonalises is none of the four, and finding it is the "
              "next step rather than quoting a number")
        return
    name = readable[0]
    eig = levels[name][0]
    print("  -> read at %s (both sides self-consistent to %.1e/%.1e e)" % (
        name, levels[name][1], levels[name][2]))
    #Say the tautology out loud.  At a level where both `self` vanish, each side's eigenvalues ARE
    #its reported diagonal, so |d eigen| is |d sorted| again and must not be quoted as a second,
    #independent measurement.  What IS new is that both sides demonstrably diagonalise this level -
    #and that is what licenses the reading below, because two sides that both diagonalise and hold
    #the same operator cannot report different occupancy multisets.
    if abs(eig - spec) < 0.01 * max(eig, spec, 1e-12):
        print("     (|d eigen| %.5f is |d sorted| %.5f again - self ~ 0 means the eigenvalues ARE "
              "the reported diagonal, so this is one number, not two; what the level table adds is "
              "that both sides DO diagonalise here)" % (eig, spec))
    if off <= 0:
        print("  -> no same-l mixing to explain")
    elif eig > 0.1 * off:
        print("  -> the block EIGENVALUES DIFFER by %.5f e (%.0f%% of the same-l mixing): the two "
              "blocks do not hold the same density, so no choice of vectors inside the block can "
              "repair it and a fix belongs UPSTREAM of the diagonalisation" % (eig, 100 * eig / off))
    #Order before assignment: a swapped pair has agreeing eigenvalues, an agreeing sorted table and
    #an off-diagonal mixing all at once, and reading that as a wrong shape is the mistake this arm
    #exists to avoid.
    elif rank > 2 * spec and rank > 0.01:
        print("  -> the same populations in a DIFFERENT ORDER (rank %.5f vs sorted %.5f): a "
              "labelling and class-assignment difference, not a shape error" % (rank, spec))
    elif spec > 0.1 * off:
        print("  -> the block eigenvalues AGREE (%.5f e) while the reported tables differ by "
              "%.5f e: the block is right and what is done INSIDE it is not, so a fix belongs in "
              "the diagonalisation" % (eig, spec))
    else:
        print("  -> eigenvalues and reported tables both agree while the mixing is not a "
              "permutation (%.5f e): the vectors differ without moving population - a rotation "
              "inside a near-degenerate set, harmless for the class tables" % off)


def demo():
    m = unpack_upper([1.0, 0.5, 2.0], 2)
    assert m[0, 1] == m[1, 0] == 0.5 and m[1, 1] == 2.0
    # The layout test needs a non-trivial metric to be decisive at all: with S = 1 an orthogonal
    # matrix and its transpose are both orthonormal and the test cannot tell them apart, which is
    # why it asserts that the loser fails rather than just taking the winner.
    n = 3
    A = np.arange(1.0, n * n + 1).reshape(n, n)
    S = np.eye(n) + 0.1 * (A + A.T) / A.max()
    w, V = np.linalg.eigh(S)
    Sinv_half = V @ np.diag(w ** -0.5) @ V.T
    Q = np.linalg.qr(A + np.eye(n) * 7)[0]
    C_true = Sinv_half @ Q  # S-orthonormal: C^T S C = 1, and its transpose is not
    vals = C_true.T.reshape(-1)  # column-major flattening
    path = os.path.join(os.environ.get("TEMP", "."), "demo_lfn33.txt")
    with open(path, "w") as f:
        f.write("NAOs in the AO basis:\n" + "\n".join("%.9f" % v for v in vals) + "\n")
    C, layout, err = read_lfn33(path, n, S)
    assert layout == "column-major" and err < 1e-9, (layout, err)
    assert np.abs(C - C_true).max() < 1e-9
    # The same file read as a PNAO matrix must be REFUSED: "NAOs in the AO basis:" is a substring of
    # "PNAOs in the AO basis:", so a plain `in` test would have read one matrix as the other.
    try:
        read_lfn32(path, n, S)
        raise AssertionError("a PNAO read of an NAO file must be refused")
    except AssertionError as e:
        assert "PNAOs in the AO basis:" in str(e), e
    # The PNAO reader's own layout test, on a set that is normalised but NOT orthonormal - which is
    # what a pre-NAO set is, and where read_lfn33's test would fail on both candidates.
    Cp = C_true @ np.array([[1.0, 0.4, 0.0], [0.0, 1.0, 0.3], [0.2, 0.0, 1.0]])
    Cp = Cp / np.sqrt(np.diag(Cp.T @ S @ Cp))          # columns normalised, not orthogonal
    ppath = os.path.join(os.environ.get("TEMP", "."), "demo_lfn32.txt")
    with open(ppath, "w") as f:
        f.write(" PNAOs in the AO basis:\n" + "\n".join("%.9f" % v for v in Cp.T.reshape(-1)))
    Cr, playout, perr = read_lfn32(ppath, n, S)
    assert playout == "column-major" and perr < 1e-8, (playout, perr)  # 9 decimals on disk
    assert np.abs(Cr - Cp).max() < 1e-9
    assert np.abs(Cr.T @ S @ Cr - np.eye(n)).max() > 0.1  # the NAO test would have had nothing here
    # pre_check refuses a set that is not orthonormal inside the block ...
    o, off, _ = pre_check(Cp, [[0], [1], [2]], S, Cp.T @ S @ np.eye(n) @ S @ Cp)
    assert o > 0.1, o
    # ... and the A/B/C verdict says C when it does, whatever the mixing looks like.
    def prow(rank, **kw):
        r = dict(mol="m", atom=0, l=0, rank=rank, nsh=2, cls="Val", gcls="Val", diag=1.0,
                 best=rank, bestval=1.0, inblock=1.0, ortho_n=0.0, off_n=0.0, ortho_g=0.0,
                 off_g=0.0, deig=0.0, occ_n=1.0, occ_g=1.0, degsum=1.0, ndeg=1, bestdeg=False)
        r.update(kw)
        return r
    # A degenerate pair mixed into each other is NOT a disagreement: the eigenspace is reproduced,
    # only its arbitrary basis differs, and the rank "mismatch" inside it is arbitrary on both sides.
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        pre_verdict({"m": [prow(0, diag=0.4, best=1, bestdeg=True, ndeg=2, degsum=1.0),
                           prow(1, diag=0.4, best=0, bestdeg=True, ndeg=2, degsum=1.0)]})
    assert "OUTCOME B" in out.getvalue(), out.getvalue()
    # A failed self-check outranks a perfect mixing: C is checked before A and B, because a metric
    # that fails its own side's gate is not made admissible by coming out flattering.
    for kw, want in ((dict(off_n=0.3), "OUTCOME C"), (dict(inblock=1.5), "OUTCOME C"),
                     (dict(), "OUTCOME B"), (dict(diag=0.4, best=1), "OUTCOME A")):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            pre_verdict({"m": [prow(0, **kw), prow(1, **kw)]})
        assert want in out.getvalue(), (want, out.getvalue())
    # mixing: identical sets give the identity, and a swapped pair shows up off-diagonal
    Cn = np.eye(4)
    Cg = np.eye(4)[:, [1, 0, 2, 3]]
    M, tot, _, _ = mixing(Cn, Cg, np.eye(4), [[0], [1], [2], [3]], [[0], [1], [2], [3]])
    assert abs(M[0, 1] - 1.0) < 1e-12 and abs(M[0, 0]) < 1e-12, M
    assert np.abs(tot - 1.0).max() < 1e-12, tot
    # Two identical p shells stored m-major on the native side and component-major on gennbo's must
    # give M = 1 on the diagonal.  This is the case that read 1/3 while the grouping was wrong, so it
    # is the one the demo has to carry: with the old slicing it is 0.333, not 1.
    I6 = np.eye(6)
    nsh = [[0, 1, 2], [3, 4, 5]]                    # native: shell 1 (x, y, z), then shell 2
    gsh = [[0, 2, 4], [1, 3, 5]]                    # gennbo: x(1, 2), y(1, 2), z(1, 2)
    Cg6 = I6[:, [0, 3, 1, 4, 2, 5]]                 # the same six functions in gennbo's print order
    M6, tot6, _, _ = mixing(I6, Cg6, I6, nsh, gsh)
    assert np.abs(M6 - np.eye(2)).max() < 1e-12, M6
    assert np.abs(tot6 - 1.0).max() < 1e-12, tot6
    os.remove(path)
    # The pairing diagnostic has to call a systematic off-by-one what it is: completeness cannot,
    # because a row sum is the same whichever column carried the weight.
    def row(rank, best, total=1.0, diag=0.9, crossl=0.0, otheratom=0.0, occshell=1.0,
            goccshell=1.0, pred=1.0, predfull=None, cls="Val", l=0, spectra=None,
            degsum=None, ndeg=1):
        return dict(mol="m", atom=1, l=l, rank=rank, cls=cls, gcls=cls, gshell="2s", occ=1.0,
                    gocc=1.0, diag=diag, inblock=1.0, best=best, bestval=0.9, gap=0.5, total=total,
                    degsum=diag if degsum is None else degsum, ndeg=ndeg,
                    crossl=crossl, otheratom=otheratom, occshell=occshell, goccshell=goccshell,
                    pred=pred, predfull=occshell if predfull is None else predfull,
                    spectra=spectra or {n: (0.0, 0.0, 0.0, 0) for n, _, _ in LEVELS})
    # The degeneracy-projected view has to DIFFER from the rank-paired one on a rotation inside a
    # degenerate pair - 1.00000 of rank-paired leak, 0.00000 once projected - and it has to be the
    # SAME number on a shell that is alone, or it is a blanket loosening rather than a gauge fix.
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        by_group([row(0, 1, diag=0.0, degsum=1.0, ndeg=2), row(1, 0, diag=0.0, degsum=1.0, ndeg=2)])
    seg, cur = {}, None
    for line in out.getvalue().splitlines():
        if line.startswith("where the shape error sits ("):
            cur = line.split("(", 1)[1].split(",")[0]
        elif cur and line.startswith("  by class"):
            seg[cur] = line
    assert "Val: 1.00000 (2)" in seg["mean 1-M[a][a]"], seg
    assert "Val: 0.00000 (2)" in seg["mean 1-degsum"], seg
    assert "GAUGE" in out.getvalue() and "2 of 2" in out.getvalue(), out.getvalue()
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        by_group([row(0, 0, diag=0.4), row(1, 1, diag=0.4)])
    seg, cur = {}, None
    for line in out.getvalue().splitlines():
        if line.startswith("where the shape error sits ("):
            cur = line.split("(", 1)[1].split(",")[0]
        elif cur and line.startswith("  by class"):
            seg[cur] = line
    assert seg["mean 1-M[a][a]"] == seg["mean 1-degsum"], seg
    assert "0 of 2" in out.getvalue(), out.getvalue()

    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        zero_check({"shifted": [row(a, a + 1) for a in range(4)]})
    assert "SYSTEMATIC pairing error" in out.getvalue(), out.getvalue()
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        zero_check({"clean": [row(a, a) for a in range(4)]})
    assert "SYSTEMATIC" not in out.getvalue() and "on partner 4/4" in out.getvalue()
    # ... and a row-sum shortfall must be attributed to the read, not to nao.cpp
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        zero_check({"short": [row(a, a, total=0.93) for a in range(4)]})
    assert "MISSING WEIGHT" in out.getvalue(), out.getvalue()
    # The bridge must call a mis-shape too small for the final table a refutation, and it must not
    # count weight on another atom towards an intra-atomic class error.  One shell, occ 1, l = 0,
    # diag 0.9 in-block 1.0 -> intra = 0.1 e; a final d(Val) of 0.5 e cannot come out of that.
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        bridge({"tight": [row(0, 0)]}, final=lambda m: (-0.5, 0.01))
    got = out.getvalue()
    assert "SMALLER than the final class error on tight" in got, got
    assert "cannot explain away" in got, got  # |d(Val)| = 0.5 e, well above the resolution
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        bridge({"loose": [row(0, 0, crossl=0.4)]}, final=lambda m: (-0.2, 0.01))
    assert "SMALLER" not in out.getvalue(), out.getvalue()
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        bridge({"offatom": [row(0, 0, diag=1.0, otheratom=0.4)]}, final=lambda m: (-0.2, 0.01))
    got = out.getvalue()
    #Weight on ANOTHER atom must not be credited to an intra-atomic class error, and it must be
    #visible in the inter-atomic total instead: 0.4 e of it against a 0.01 e charge deviation.
    assert "SMALLER than the final class error on offatom" in got, got
    assert "totals 0.40000 e" in got and "40x cancellation" in got, got
    # The signed arm: predfull is an identity, so a row whose predfull does not reproduce its own
    # occupancy must void the report rather than print a sign.
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        signed({"broken": [row(0, 0, occshell=1.2, predfull=1.0)]})
    assert "EXACTNESS GATE FAILED" in out.getvalue(), out.getvalue()
    # ... and with the gate passing, a native valence total BELOW gennbo's must come out negative,
    # with the coherence residual named.  occshell 0.8 vs goccshell 1.0 -> d(Val) = -0.2; pred 0.85
    # -> -0.15, same sign, residual +0.05.
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        signed({"low": [row(0, 0, occshell=0.8, goccshell=1.0, pred=0.85)]})
    got = out.getvalue()
    assert "-0.20000" in got and "-0.15000" in got and "+0.05000" in got, got
    assert "sign reproduced on 1 of 1" in got, got
    # A wrong sign has to say so: native total ABOVE gennbo's while the prediction is below.
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        signed({"flip": [row(0, 0, occshell=1.2, goccshell=1.0, predfull=1.2, pred=0.9)]})
    assert "WRONG" in out.getvalue() and "sign reproduced on 0 of 1" in out.getvalue()
    # spectrum: two shells of one block holding the same two occupancies in the opposite order is
    # the same spectrum, not a shape error - and that must not be read as vectors-are-wrong.
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        spectrum({"swap": [row(0, 0, occshell=1.9, goccshell=0.1, diag=0.5),
                           row(1, 1, occshell=0.1, goccshell=1.9, diag=0.5)]})
    got = out.getvalue()
    assert "DIFFERENT ORDER" in got, got
    # Eigenvalues AND reported tables agreeing while the mixing is not a permutation is not a
    # finding about the vectors being wrong - it is a rotation that moves no population, and
    # calling it an error is the over-read this branch is full of.
    lvl = lambda **kw: {n: kw.get(n.split(",")[0].replace(" ", "_"), (0.0, 0.0, 0.0, 0))
                        for n, _, _ in LEVELS}
    ok = lvl()
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        spectrum({"vec": [row(0, 0, occshell=1.9, goccshell=1.9, diag=0.5, spectra=ok),
                          row(1, 1, occshell=0.1, goccshell=0.1, diag=0.5, spectra=ok)]})
    assert "without moving population" in out.getvalue(), out.getvalue()
    # The same two occupancies in the opposite order is a labelling difference, not a shape error.
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        spectrum({"swap2": [row(0, 0, occshell=1.9, goccshell=0.1, diag=0.5, spectra=ok),
                            row(1, 1, occshell=0.1, goccshell=1.9, diag=0.5, spectra=ok)]})
    assert "DIFFERENT ORDER" in out.getvalue(), out.getvalue()
    # Eigenvalues that genuinely differ at the readable level point upstream.
    ups = lvl(per_NAO=(0.5, 0.0, 0.0, 0))
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        spectrum({"ups": [row(0, 0, diag=0.9, spectra=ups), row(1, 1, diag=0.9, spectra=ups)]})
    assert "EIGENVALUES DIFFER" in out.getvalue(), out.getvalue()
    # A level where the FIRST candidate is not self-consistent must be skipped, not quoted: here
    # only the m-averaged level is readable, and its 0.0 must drive the verdict, not the 0.9.
    skip = {n: ((0.9, 0.3, 0.3, 0) if not avg else (0.0, 0.0, 0.0, 0)) for n, avg, _ in LEVELS}
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        spectrum({"skip": [row(0, 0, diag=0.5, spectra=skip), row(1, 1, diag=0.5, spectra=skip)]})
    got = out.getvalue()
    assert "read at m-averaged, full block" in got, got
    assert "EIGENVALUES DIFFER" not in got, got
    # And when the readable level's eigenvalues are just the reported table again, say so: quoting
    # it as a second measurement is exactly the double-counting three retired metrics died of.
    tau = {n: ((0.8, 0.0, 0.0, 0) if avg else (0.9, 0.3, 0.3, 0)) for n, avg, _ in LEVELS}
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        spectrum({"tau": [row(0, 0, occshell=1.5, goccshell=1.9, diag=0.9, spectra=tau),
                          row(1, 1, occshell=0.5, goccshell=0.1, diag=0.9, spectra=tau)]})
    assert "is one number, not two" in out.getvalue(), out.getvalue()
    # Native alone failing a level gennbo passes is a statement about nao.cpp, and must be said.
    lop = {n: (0.0, 0.3, 0.0, 0) for n, _, _ in LEVELS}
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        spectrum({"lop": [row(0, 0, spectra=lop), row(1, 1, spectra=lop)]})
    got = out.getvalue()
    assert "native is not diagonalising what NBO diagonalises" in got, got
    assert "NO level has both sides self-consistent" in got, got
    # The cascade arm's four measurements are checked by their own doctests, RUN here - nothing in
    # this repository ran a doctest before, so the examples above this line are prose until someone
    # points a runner at them.  These four are not.
    import doctest
    runner, finder = doctest.DocTestRunner(verbose=False), doctest.DocTestFinder()
    for f in (cascade_transform, sym_shape, sign_frustration, tri_shape,
              cascade_orthogonal, principal_sines, class_mass):
        for t in finder.find(f, f.__name__, globs=globals()):
            runner.run(t)
    assert runner.failures == 0, "%d doctest failure(s) in the cascade arm" % runner.failures
    # And the two shapes are distinguishable, which is the whole premise: a Schmidt-triangular
    # transform and an occupancy-weighted symmetric one must not read the same on either measure.
    c = ["Cor"] * 2 + ["Val"] * 2 + ["Ryd"] * 2
    al = [(0, 0), (1, 0)] * 3
    rng = np.random.default_rng(7)
    M = rng.normal(size=(6, 6)); M = M + M.T
    owso = np.diag([1.0, 2.0, 0.5, 4.0, 0.25, 1.5]) @ M           # W^-1 M
    schmidt = np.triu(rng.normal(size=(6, 6)) + 3.0 * np.eye(6))  # zero below the diagonal
    pairs = [(a, b) for a in range(6) for b in range(a + 1, 6) if c[a] == c[b] and al[a] != al[b]]
    assert sym_shape(owso, pairs)["resid"] < 1e-12, sym_shape(owso, pairs)
    assert sign_frustration(pairs, {(i, j) for i, j in pairs if owso[i, j] * owso[j, i] < 0})[1] == 0
    bo, bs = tri_shape(owso, c, c), tri_shape(schmidt, c, c)
    assert min(bo[("Ryd", "Cor")], bo[("Val", "Cor")]) > 0.05, bo   # OWSO fills both triangles
    assert max(bs[("Ryd", "Cor")], bs[("Val", "Cor")]) == 0.0, bs   # Schmidt empties one
    assert sym_shape(schmidt, pairs)["n"] == 0                      # and its pairs never both clear
    # --space: the instrument is gated on its ability to DISCRIMINATE before any leak number is read.
    # Built on a non-trivial S, because an instrument that only works at S = 1 would not be one.
    rng = np.random.default_rng(11)
    A = rng.normal(size=(8, 8))
    S8 = A @ A.T + 8.0 * np.eye(8)
    X = np.linalg.inv(np.linalg.cholesky(S8)).T
    Q, _ = np.linalg.qr(rng.normal(size=(8, 8)))
    Cn = X @ Q                                                # an S-orthonormal "native" NAO set
    cls = ["Cor"] * 2 + ["Val"] * 3 + ["Ryd"] * 3
    iv, ir = [2, 3, 4], [5, 6, 7]
    # (a) CODE check, not a gate: a space compared with itself gives zero identically, so passing it
    #     says the svd and the clip are right and says nothing whatever about the instrument.  The
    #     tolerance is 1e-07 and not 1e-12 BECAUSE a sine is sqrt(1 - sigma**2): at 1e-15 in the
    #     cosine this comes out near 4e-08, which is where the sqrt floor in the header was found.
    assert principal_sines(Cn[:, iv], Cn[:, iv], S8).max() < 1e-7
    # (b) MUST NOT MOVE under the two gauges that are each side's own free choice - a per-column sign
    #     and the order of columns within a class.  This is what makes the permutation between the
    #     sides irrelevant rather than merely unresolved.
    gauge = Cn.copy() * np.where(rng.random(8) < 0.5, -1.0, 1.0)
    gauge[:, iv] = gauge[:, [4, 2, 3]]
    for c, idx in (("Val", iv), ("Ryd", ir)):
        assert principal_sines(Cn[:, idx], gauge[:, idx], S8).max() < 1e-7, c
    assert abs(class_mass(Cn, gauge, S8, cls, cls)[("Val", "Ryd")]) < 1e-12   # mass is linear, so tight
    # (c) MUST MOVE, by a known amount, under an injected valence -> Rydberg leak: one Givens
    #     rotation of eps between a Val and a Ryd column puts sin(eps) into the largest principal
    #     sine of BOTH classes and exactly sin(eps)**2 dimensions into mass[Val, Ryd].  An exact
    #     analytic expectation, so this calibrates the instrument instead of merely exercising it.
    for eps in (1e-6, 1e-3, 0.1, np.pi / 6):
        leak = Cn.copy()
        leak[:, 4] = np.cos(eps) * Cn[:, 4] + np.sin(eps) * Cn[:, 5]
        leak[:, 5] = -np.sin(eps) * Cn[:, 4] + np.cos(eps) * Cn[:, 5]
        for idx in (iv, ir):
            got = principal_sines(Cn[:, idx], leak[:, idx], S8).max()
            assert abs(got - np.sin(eps)) < 1e-7, (eps, got)   # sqrt floor, see (a)
        m = class_mass(Cn, leak, S8, cls, cls)
        assert abs(m[("Val", "Ryd")] - np.sin(eps) ** 2) < 1e-9, (eps, m[("Val", "Ryd")])
        assert abs(m[("Cor", "Cor")] - 2.0) < 1e-9                # an untouched class stays put
        # (d) and the metric that was proposed as the FIRST gate cannot see any of it: U is
        #     orthogonal, so its singular values are 1 whatever was done to the classes.  At eps =
        #     30 degrees a quarter of a dimension has moved and it still reads 1.000000.
        sv = np.linalg.svd(Cn.T @ S8 @ leak, compute_uv=False)
        assert np.abs(sv - 1.0).max() < 1e-12, sv
    assert abs(np.sin(np.pi / 6) ** 2 - 0.25) < 1e-12             # the "quarter of a dimension"
    print("demo ok")


#The three molecules whose FINAL table is already right (pf5, so2, sf6 - Rydberg ratio 0.93-0.98)
#against the two that are worst (ethane 6.2x, benzene 3.7x).  Not an exact prediction, but the one
#the metric has to satisfy to be measuring the thing it claims to explain.
ALREADY_RIGHT = ("pf5", "so2", "sf6")
WORST = ("ethane", "benzene")


def zero_check(rows_by_mol):
    """The metric's own gates.

    The one free falsification here is EXACT and applies to every molecule: each native shell's
    m-averaged mixing summed over ALL gennbo NAOs is 1, because both sides are S-orthonormal sets
    spanning the same AO space.  A deviation is an error in the matrix read, the pairing or the
    block bookkeeping - in the metric, not in nao.cpp - and nothing else may be quoted until it
    passes.

    It replaces a prediction that does not hold: the out-of-block remainder is NOT expected to be
    zero on LiF or N2.  NAOs are orthonormal over the whole molecule on each side, not per atom, so
    a native NAO on one atom generically has amplitude on another atom's and on other l in gennbo's
    basis - which is the inter-l leak this is here to measure, not an artefact.
    """
    print("\ncompleteness gate (sum_b M[a][b] over ALL gennbo NAOs must be 1, exactly):")
    worst = None
    for mol, rows in sorted(rows_by_mol.items()):
        dev = max(abs(r["total"] - 1.0) for r in rows)
        worst = dev if worst is None else max(worst, dev)
        print("  %-10s max |sum - 1| = %.2e over %d shells" % (mol, dev, len(rows)))
    #Which side a failure is on: completeness holds only if lfn 33 carries the full n x n set, so a
    #SHORTFALL is missing weight - a truncated or pruned print, an error in the read - while an
    #excess would have to be a bookkeeping error in the blocks.  "0.93 instead of 1" reads like a
    #small error and means a missing column.
    if worst < 1e-8:
        print("  -> passes (worst %.2e)" % worst)
    else:
        short = min(r["total"] - 1.0 for rows in rows_by_mol.values() for r in rows)
        over = max(r["total"] - 1.0 for rows in rows_by_mol.values() for r in rows)
        print("  -> METRIC REFUTED: do not quote any leak number (worst %.2e; most negative %.5f, "
              "most positive %.5f)" % (worst, short, over))
        if -short > over:
            print("     a shortfall is MISSING WEIGHT: lfn 33 is not the full n x n set (truncated "
                  "or pruned print), so the fault is on the read side, not in nao.cpp")
        else:
            print("     an excess cannot come from a missing column - suspect the block "
                  "bookkeeping (a column counted in two blocks)")

    #Completeness is blind to pairing, because a row sum does not care which column carried the
    #weight.  This is the pairing-sensitive companion: diagnostic, not pass/fail, since a genuine
    #leak WILL move the argmax off the diagonal - that is the finding.  A systematic off-by-one
    #instead shows up as nearly every shell sitting at the same non-zero offset.
    print("pairing diagnostic (argmax within the (atom, l) block vs the rank partner):")
    for mol, rows in sorted(rows_by_mol.items()):
        offs = {}
        for r in rows:
            offs[r["best"] - r["rank"]] = offs.get(r["best"] - r["rank"], 0) + 1
        on = offs.get(0, 0)
        gaps = [r["gap"] for r in rows if r["best"] != r["rank"]]
        print("  %-10s on partner %d/%d   offsets %s%s" % (
            mol, on, len(rows), " ".join("%+d:%d" % kv for kv in sorted(offs.items())),
            "" if not gaps else "   worst off-partner gap %.5f" % max(gaps)))
    allrows = [r for rows in rows_by_mol.values() for r in rows]
    dominant = max(set(r["best"] - r["rank"] for r in allrows),
                   key=lambda o: sum(1 for r in allrows if r["best"] - r["rank"] == o))
    share = sum(1 for r in allrows if r["best"] - r["rank"] == dominant) / len(allrows)
    if dominant != 0 and share > 0.5:
        print("  -> %.0f %% of all shells sit at offset %+d: that is a SYSTEMATIC pairing error, "
              "not a leak - fix it before reading the numbers above" % (100 * share, dominant))

    def mean_leak(names):
        rows = [r for m in names for r in rows_by_mol.get(m, [])]
        return (sum(1.0 - r["diag"] for r in rows) / len(rows), len(rows)) if rows else (None, 0)

    good, ng = mean_leak(ALREADY_RIGHT)
    bad, nb = mean_leak(WORST)
    if good is None or bad is None:
        print("contrast gate: not both groups present in this set (%d/%d shells)" % (ng, nb))
        return
    print("contrast gate: mean in-block leak %.5f on the already-right %s (%d shells) vs %.5f on "
          "%s (%d shells) -> %s" % (
              good, "/".join(ALREADY_RIGHT), ng, bad, "/".join(WORST), nb,
              "tracks the final-table error" if bad > good else
              "does NOT track it - the metric may be measuring something else"))


def main(argv):
    verbose = "--full" in argv
    root = [a for a in argv[1:] if not a.startswith("--")][0]
    if "--pre" in argv:
        return main_pre(root, verbose)
    if "--cascade" in argv:
        return main_cascade(root)
    if "--space" in argv:
        return main_space(root)
    print("aonao_compare on %s, 1 thread (numpy on matrices of a few hundred), root %s" % (
        socket.gethostname(), root))
    all_rows, rows_by_mol, void, missing = [], {}, [], []
    for mol in sorted(os.listdir(root)):
        d = os.path.join(root, mol)
        if not os.path.isdir(d):
            continue
        if not os.path.isfile(os.path.join(d, mol + ".33")):
            missing.append(mol)
            continue
        try:
            res = compare(mol, d, verbose)
        except AssertionError as e:
            void.append((mol, str(e)))
            print("%-10s VOID: %s" % (mol, e))
            continue
        report(res, verbose)
        rows_by_mol[mol] = res["rows"]
        all_rows += res["rows"]
    #The denominator, printed whether or not it is flattering: a mean over a set that quietly
    #shrank is the failure mode that cost the ELI lane a corpus headline.
    print("\ndenominator: %d molecules compared, %d VOID, %d without an lfn 33" % (
        len(rows_by_mol), len(void), len(missing)))
    for mol, why in void:
        print("  VOID %-10s %s" % (mol, why))
    if missing:
        print("  no lfn 33: %s" % " ".join(missing))
    if all_rows:
        by_group(all_rows)
        zero_check(rows_by_mol)
        bridge(rows_by_mol)
        signed(rows_by_mol)
        spectrum(rows_by_mol)
    return 0


if __name__ == "__main__":
    sys.exit(demo() if "--demo" in sys.argv else main(sys.argv))
