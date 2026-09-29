"""The third code against NoSpherA2's own two cascades, on the same eight molecules.

    py -3.12 janpa_vs_native.py <janpa work dir> [--full]

janpa_third_code.py answers "is gennbo 7 running the published algorithm?" and the answer is no.
That immediately raises the question this script measures, which no reading of a spec can answer:
if the spec's own reference implementation does not reach gennbo either, then how far is NATIVE
from the spec, and is the remaining gap to gennbo the spec's gap or ours?

Three sides, one AO basis (the .47's), one gauge:
  gennbo  unit 33, the arbiter, read with the lane's validated reader.
  JANPA   the published implementation, brought into the .47's order and sign by the constructed,
          asserted gauge of janpa_third_code.py (positive control green at 1e-09 on all eight).
  native  the C++ dumps of BOTH arms, `legacy` (the four-step cascade) and `renat5` (spec steps 5
          and 6 restored, the default since 519cb798), already in the .47's own order and gauge.

THE ARM LABELS ARE VERIFIED BY MEASUREMENT, NOT BY THE DIRECTORY NAME.  A label beside a binary is
not provenance - that trap cost this lane a whole extrapolation column (see the trap note of
25 Sep).  Each arm's Rydberg population error against gennbo's printed table is printed, and the
published numbers for these two arms are 0.3161 e (legacy) and 0.0050 e (renat5); if the two
columns do not reproduce those, the dumps are not the arms the directories claim and every
comparison below is void.

The .47 that gennbo consumed in the matrix run carries extra $NBO keywords, so its md5 differs
from the one in the sines staging; the wavefunction itself is checked element by element instead.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import janpa_third_code as J  # noqa: E402
from aonao_compare import (gennbo_blocks, mixing, principal_sines, read_47,  # noqa: E402
                           read_lfn33, read_naoc)

HERE = os.path.dirname(os.path.abspath(__file__))
SIDES = {"legacy": os.path.join(HERE, "_sinestage", "legacy"),
         "renat5": os.path.join(HERE, "_sinestage", "renat5")}
MOLS = "ammonia benzene ethane lif pf5 sf6 so2 water".split()


def coord_z(path47):
    """Nuclear charges from the .47's own $COORD block (first column of each atom line)."""
    lines = open(path47).read().split("$COORD")[1].split("$END")[0].strip().split("\n")[1:]
    return np.array([int(l.split()[0]) for l in lines])


def npa_from_dump(labels, occ, z):
    """NPA charges the way NPA defines them: Z_A minus the occupancy of the NAOs on A."""
    q = z.astype(float).copy()
    for i, lab in enumerate(labels):
        q[lab[1]] -= occ[i]
    return q


def native_blocks(labels, occ):
    """(atom0, l) -> shells, each the (2l+1) native NAO indices, ranked by m-averaged occupancy.

    The native dump states (atom, l, m, shell, class, occ) per column, so the shell is read from
    its own field - the same rule as JANPA's R number, and for the same reason.
    """
    shells = {}
    for i, lab in enumerate(labels):
        _, atom, l, _, shell, _, _ = lab
        shells.setdefault((atom, l), {}).setdefault(shell, []).append(i)
    out = {}
    for key, by_s in shells.items():
        for s, idx in by_s.items():
            assert len(idx) == 2 * key[1] + 1, "native block %s shell %d has %d" % (
                key, s, len(idx))
        out[key] = [by_s[s] for s in sorted(by_s, key=lambda s: -np.mean(occ[by_s[s]]))]
    return out


def arm_vs(Ca, Cb, S, ba, bb, tg, keys):
    """(worst group metric, worst principal sine, n shells) between two sides."""
    worst, npairs = 0.0, 0
    for key in keys:
        if len(ba[key]) != len(bb[key]):
            continue
        M, _, _, _ = mixing(Ca, Cb, S, ba[key], bb[key])
        npairs += len(ba[key])
        worst = max(worst, J.group_metric(M, tg[key]))
    return worst, npairs


def one(mol, jdir):
    d47 = os.path.join(jdir, mol + ".47")
    n, S, P = read_47(d47)
    SPS = S @ P @ S
    # the same wavefunction on both paths, checked rather than assumed (the md5s differ: the ops
    # .47 carries the W<n> matrix keywords)
    n2, S2, P2 = read_47(os.path.join(SIDES["legacy"], mol, mol + ".47"))
    assert n2 == n
    wf = max(float(np.abs(S - S2).max()), float(np.abs(P - P2).max()))

    _, perm = J.ao_perm(d47)
    Q = np.zeros((n, n))
    Q[perm, np.arange(n)] = 1.0
    Sj = Q @ J.read_janpa(os.path.join(jdir, mol + ".j.S.txt"), n) @ Q.T
    D, resid, _, bad = J.sign_gauge(Sj, S)
    G = Q @ np.diag(D[perm])
    Cj, _, ej, _ = J.read_janpa_c(os.path.join(jdir, mol + ".j.nao2ao.txt"), n, G, S, "nao")
    lab_j = J.janpa_labels(os.path.join(jdir, mol + ".j.nao2ao.txt"))
    occ_j = np.diag(Cj.T @ SPS @ Cj)

    Cg, _, eg = read_lfn33(os.path.join(jdir, mol + ".33"), n, S)
    naos = json.load(open(os.path.join(jdir, mol + ".aopnao.nbo.json")))["nao"]
    gb = gennbo_blocks(naos)
    jb = J.janpa_blocks(lab_j, occ_j)
    tg = {k: J.tie_groups([float(np.mean([naos[i]["occupancy"] for i in sh])) for sh in gb[k]])
          for k in gb}
    keys = sorted(gb)

    z = coord_z(d47)
    nat, nb, npop, nsine, nq = {}, {}, {}, {}, {}
    for arm, root in SIDES.items():
        labs, C = read_naoc(os.path.join(root, mol, mol + ".naoc.txt"))
        assert C.shape == (n, n), (arm, C.shape, n)
        occ = np.diag(C.T @ SPS @ C)
        nat[arm], nb[arm] = C, native_blocks(labs, occ)
        ryd = [i for i, lab in enumerate(labs) if lab[5] == 2]
        npop[arm] = float(occ[ryd].sum())
        # per-class subspace sine against gennbo, the lane's own gauge-free instrument
        s = 0.0
        for c in (0, 1, 2):
            a = [i for i, lab in enumerate(labs) if lab[5] == c]
            b = [e["index"] - 1 for e in naos
                 if {0: "Cor", 1: "Val", 2: "Ryd"}[c] == e["type"]]
            if a and len(a) == len(b):
                s = max(s, float(principal_sines(C[:, a], Cg[:, b], S).max()))
        nsine[arm] = s
        nq[arm] = npa_from_dump(labs, occ, z)
    gryd = sum(e["occupancy"] for e in naos if e["type"] == "Ryd")
    qg = np.array([e["charge"] for e in json.load(
        open(os.path.join(jdir, mol + ".aopnao.nbo.json")))["npa"]])
    qj = np.array([float(x) for x in
                   open(os.path.join(jdir, mol + ".j.npa.txt")).read().split()])
    dq = {"JANPA": float(np.abs(qj - qg).max())}
    for arm in SIDES:
        dq[arm] = float(np.abs(nq[arm] - qg).max())
    dq["legacy-JANPA"] = float(np.abs(nq["legacy"] - qj).max())
    dq["renat5-JANPA"] = float(np.abs(nq["renat5"] - qj).max())

    r = dict(mol=mol, n=n, wf_diff=wf, sign_resid=resid, sign_bad=bad, err_janpa=ej, err_g=eg,
             gennbo_ryd=gryd, pop_err={a: npop[a] - gryd for a in SIDES}, dq=dq,
             sine_vs_gennbo={a: nsine[a] for a in SIDES})
    r["janpa_vs_gennbo"], r["npairs"] = arm_vs(Cj, Cg, S, jb, gb, tg, keys)
    for arm in SIDES:
        r["%s_vs_gennbo" % arm], _ = arm_vs(nat[arm], Cg, S, nb[arm], gb, tg, keys)
        # JANPA as the reference: the tie runs are then JANPA's own, by its own occupancies
        tj = {k: J.tie_groups([float(np.mean(occ_j[sh])) for sh in jb[k]]) for k in jb}
        r["%s_vs_janpa" % arm], _ = arm_vs(nat[arm], Cj, S, nb[arm], jb, tj, keys)
    return r


def main(argv):
    jdir = argv[1]
    rows = [one(m, os.path.join(jdir, m)) for m in MOLS
            if os.path.exists(os.path.join(jdir, m, m + ".j.nao2ao.txt"))]
    assert rows, "no molecule had all three sides - an empty table is not an agreement"
    print("--- the three sides share one wavefunction and one gauge (checked, not assumed) ---")
    print("%-9s %4s %12s %12s %12s %12s" % (
        "mol", "n", "max|dS|,|dP|", "max|DSjD-S|", "janpa CtSC-I", "gennbo CtSC-I"))
    for r in rows:
        print("%-9s %4d %12.3e %12.3e %12.3e %12.3e" % (
            r["mol"], r["n"], r["wf_diff"], r["sign_resid"], r["err_janpa"], r["err_g"]))
    print("\n--- the arm labels, verified by measurement: Rydberg population error vs gennbo's"
          " printed table (published: legacy 0.3161 e, renat5 0.0050 e worst) ---")
    print("%-9s %12s %12s %12s   %12s %12s" % (
        "mol", "gennbo Ryd e", "legacy err", "renat5 err", "legacy sine", "renat5 sine"))
    for r in rows:
        print("%-9s %12.4f %12.4f %12.4f   %12.4f %12.4f" % (
            r["mol"], r["gennbo_ryd"], r["pop_err"]["legacy"], r["pop_err"]["renat5"],
            r["sine_vs_gennbo"]["legacy"], r["sine_vs_gennbo"]["renat5"]))
    print("  worst |pop err|: legacy %.4f e  renat5 %.4f e  (over %d molecules)" % (
        max(abs(r["pop_err"]["legacy"]) for r in rows),
        max(abs(r["pop_err"]["renat5"]) for r in rows), len(rows)))
    print("\n--- worst m-averaged, tie-run block metric, every pair of sides, same metric ---")
    print("%-9s %7s %13s %13s %13s   %13s %13s" % (
        "mol", "shells", "JANPA:gennbo", "legacy:gennbo", "renat5:gennbo",
        "legacy:JANPA", "renat5:JANPA"))
    for r in rows:
        print("%-9s %7d %13.3e %13.3e %13.3e   %13.3e %13.3e" % (
            r["mol"], r["npairs"], r["janpa_vs_gennbo"], r["legacy_vs_gennbo"],
            r["renat5_vs_gennbo"], r["legacy_vs_janpa"], r["renat5_vs_janpa"]))
    print("\n--- NPA charges in ELECTRONS, max |dq| per molecule.  Native's charges are computed"
          " from its own dump (Z_A - occ on A), gennbo's and JANPA's are each code's printed"
          " table; gennbo prints 4 decimals, so its own floor is 5e-05 ---")
    print("%-9s %12s %12s %12s   %12s %12s" % (
        "mol", "JANPA:gennbo", "legacy:gennbo", "renat5:gennbo", "legacy:JANPA", "renat5:JANPA"))
    for r in rows:
        print("%-9s %12.3e %12.3e %12.3e   %12.3e %12.3e" % (
            r["mol"], r["dq"]["JANPA"], r["dq"]["legacy"], r["dq"]["renat5"],
            r["dq"]["legacy-JANPA"], r["dq"]["renat5-JANPA"]))
    print("  worst over %d molecules: JANPA:gennbo %.3e | legacy:gennbo %.3e renat5:gennbo %.3e"
          % (len(rows), max(r["dq"]["JANPA"] for r in rows),
             max(r["dq"]["legacy"] for r in rows), max(r["dq"]["renat5"] for r in rows)))
    med = lambda k: float(np.median([r[k] for r in rows]))  # noqa: E731
    print("  median over %d molecules: JANPA:gennbo %.3e | legacy:gennbo %.3e renat5:gennbo %.3e"
          " | legacy:JANPA %.3e renat5:JANPA %.3e" % (
              len(rows), med("janpa_vs_gennbo"), med("legacy_vs_gennbo"), med("renat5_vs_gennbo"),
              med("legacy_vs_janpa"), med("renat5_vs_janpa")))
    json.dump(rows, open(os.path.join(jdir, "janpa_vs_native.json"), "w"), indent=1, default=str)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
