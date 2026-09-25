"""Compare NoSpherA2's AO -> NAO transformation with gennbo 7's own, orbital by orbital.

    py -3.12 aonao_compare.py <dir with one subdirectory per molecule> [--full]
    py -3.12 aonao_compare.py --demo

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
    P = unpack_upper(numbers(section(text, "$DENSITY"))[:packed], n)
    assert len(centers) == n, (len(centers), n)
    return n, S, P


def read_lfn33(path, n, S):
    """gennbo's AO -> NAO matrix, with the layout decided by C^T S C = 1 rather than assumed."""
    body = open(path).read()
    #A VOID, not a traceback: sf6's lfn 33 came back 0 bytes from the staging job.
    assert "NAOs in the AO basis:" in body, "no 'NAOs in the AO basis:' block in %s (%d bytes)" % (
        os.path.basename(path), len(body))
    start = body.index("NAOs in the AO basis:")
    vals = numbers(body[start:])[:n * n]
    assert len(vals) == n * n, "lfn 33 holds %d numbers, need %d" % (len(vals), n * n)
    a = np.array(vals)
    cands = {"column-major": a.reshape(n, n).T, "row-major": a.reshape(n, n)}
    errs = {k: float(np.abs(C.T @ S @ C - np.eye(n)).max()) for k, C in cands.items()}
    best = min(errs, key=errs.get)
    assert errs[best] < 1e-6, "neither layout of lfn 33 is S-orthonormal: %s" % errs
    other = [k for k in errs if k != best][0]
    assert errs[other] > 1e-6, "both layouts pass, the test cannot decide: %s" % errs
    return cands[best], best, errs[best]


def read_naoc(path):
    """The native dump: (labels, C) with labels[i] = (atom, l, m, shell, class, occ)."""
    labels, cols = [], []
    for line in open(path):
        if not line.startswith("NAOC "):
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
        for a in range(M.shape[0]):
            lab = labels[nshells[a][0]]
            g = naos[gshells[a][0]]
            order = np.argsort(M[a])[::-1]
            worst_b = int(order[0])
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
    for what, val in (("mean 1-M[a][a], rank-paired", lambda r: 1.0 - r["diag"]),
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
    #The metric above normalises every shell to 1 whatever its occupancy, so an empty Rydberg shell
    #with a badly wrong shape counts the same as a doubly-occupied valence shell that is nearly right.
    #Populations do not work that way, and the failure being chased is a population failure, so the
    #same three arms are also summed with (2l+1)*occ as the weight.  It is an estimate: M is
    #m-averaged, so this assumes the m components of a shell share its occupancy.
    print("\non the charge scale (sum (2l+1)*occ*leak, electrons):")
    for what, val in (("rank-paired", lambda r: 1.0 - r["diag"]),
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
            goccshell=1.0, pred=1.0, predfull=None, cls="Val", l=0, spectra=None):
        return dict(mol="m", atom=1, l=l, rank=rank, cls=cls, gcls=cls, gshell="2s", occ=1.0,
                    gocc=1.0, diag=diag, inblock=1.0, best=best, bestval=0.9, gap=0.5, total=total,
                    crossl=crossl, otheratom=otheratom, occshell=occshell, goccshell=goccshell,
                    pred=pred, predfull=occshell if predfull is None else predfull,
                    spectra=spectra or {n: (0.0, 0.0, 0.0, 0) for n, _, _ in LEVELS})
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
