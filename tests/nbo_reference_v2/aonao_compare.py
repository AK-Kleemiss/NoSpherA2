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
    """
    flat = [i for sh in rows for i in sh]
    nm = len(rows[0])
    O = Cn[:, flat].T @ S @ Cg  # native block rows against ALL gennbo NAOs
    M = np.zeros((len(rows), len(cols)))
    for a in range(len(rows)):
        ra = O[a * nm:(a + 1) * nm, :]
        for b, cb in enumerate(cols):
            M[a, b] = (ra[:, cb] ** 2).sum() / nm
    total = (O ** 2).sum(axis=1).reshape(len(rows), nm).sum(axis=1) / nm  # 1.0 by completeness
    return M, total


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
    forms = {"C^T S P S C": np.diag(Cg.T @ S @ P @ S @ Cg), "C^T P C": np.diag(Cg.T @ P @ Cg)}
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

    gb = gennbo_blocks(naos)
    rows = []
    for (atom, l), gshells in sorted(gb.items()):
        nm = 2 * l + 1
        nrows = [i for i, lab in enumerate(labels) if lab[1] == atom and lab[2] == l]
        nrows.sort(key=lambda i: (labels[i][4], labels[i][3]))  # shell rank, then m
        assert len(nrows) == nm * len(gshells), "block (%d,%d): %d native vs %d gennbo" % (
            atom, l, len(nrows), nm * len(gshells))
        nshells = [nrows[a * nm:(a + 1) * nm] for a in range(len(gshells))]
        M, total = mixing(Cn, Cg, S, nshells, gshells)
        for a in range(M.shape[0]):
            lab = labels[nshells[a][0]]
            g = naos[gshells[a][0]]
            order = np.argsort(M[a])[::-1]
            worst_b = int(order[0])
            #Gap to the runner-up: a pairing that is merely close is worth knowing about, and a
            #systematic off-by-one has a LARGE gap on the wrong column, not a small one.
            gap = float(M[a, order[0]] - M[a, order[1]]) if M.shape[1] > 1 else float("nan")
            rows.append(dict(mol=mol, atom=atom, l=l, rank=a, cls=CLASS[lab[5]],
                             gcls=g["type"], gshell=g["shell"], occ=lab[6], gocc=g["occupancy"],
                             diag=float(M[a, a]), inblock=float(M[a].sum()),
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
                      ("out-of-block", lambda r: 1.0 - r["inblock"])):
        groups = {}
        for r in all_rows:
            groups.setdefault(r["cls"], []).append((2 * r["l"] + 1) * r["occ"] * val(r))
        print("  %-21s %s   total %.5f" % (
            what, "  ".join("%s: %.5f" % (k, sum(v)) for k, v in sorted(groups.items())),
            sum(sum(v) for v in groups.values())))


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
    M, tot = mixing(Cn, Cg, np.eye(4), [[0], [1], [2], [3]], [[0], [1], [2], [3]])
    assert abs(M[0, 1] - 1.0) < 1e-12 and abs(M[0, 0]) < 1e-12, M
    assert np.abs(tot - 1.0).max() < 1e-12, tot
    # Two identical p shells stored m-major on the native side and component-major on gennbo's must
    # give M = 1 on the diagonal.  This is the case that read 1/3 while the grouping was wrong, so it
    # is the one the demo has to carry: with the old slicing it is 0.333, not 1.
    I6 = np.eye(6)
    nsh = [[0, 1, 2], [3, 4, 5]]                    # native: shell 1 (x, y, z), then shell 2
    gsh = [[0, 2, 4], [1, 3, 5]]                    # gennbo: x(1, 2), y(1, 2), z(1, 2)
    Cg6 = I6[:, [0, 3, 1, 4, 2, 5]]                 # the same six functions in gennbo's print order
    M6, tot6 = mixing(I6, Cg6, I6, nsh, gsh)
    assert np.abs(M6 - np.eye(2)).max() < 1e-12, M6
    assert np.abs(tot6 - 1.0).max() < 1e-12, tot6
    os.remove(path)
    # The pairing diagnostic has to call a systematic off-by-one what it is: completeness cannot,
    # because a row sum is the same whichever column carried the weight.
    def row(rank, best, total=1.0, diag=0.9):
        return dict(mol="m", atom=1, l=0, rank=rank, cls="Val", gcls="Val", gshell="2s", occ=1.0,
                    gocc=1.0, diag=diag, inblock=1.0, best=best, bestval=0.9, gap=0.5, total=total)
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
    return 0


if __name__ == "__main__":
    sys.exit(demo() if "--demo" in sys.argv else main(sys.argv))
