"""Per-spin NAO comparison for the open-shell held-out molecules.

    py -3.12 spinsplit.py <holdout dir> [mol ...]      # default: all open shells
    py -3.12 spinsplit.py --selfcheck <holdout dir>    # asserts only, no table

Why this exists
---------------
The subspace (principal-sine) instrument refuses open shells, so `cf3` - the
held-out molecule whose d(Ryd) got 2.08x WORSE under renat5 - had no
measurement that could say where the regression sits.  It turns out no new
reference is needed for the POPULATION half of that question: gennbo's own
per-spin NAO tables are already parsed into `nao_alpha` / `nao_beta` in the
committed reference JSONs, and `-nbo_native` emits the same two keys, matched
row for row.  This compares them directly, per spin and per class.

It does NOT close the gate and does not claim agreement with NBO 7: the gate is
an exact basis-function match and native fails it on every molecule of both
sets.  These are signed deltas with their print floors beside them.

Three traps this file is written against
----------------------------------------
1. Read the spin from the KEY (`nao_alpha`/`nao_beta`), never from a column
   position and never inherited from the nearest preceding header.  The same
   defect family has bitten this lane three times (`bf73d301`, `53f11d7f`).
2. `nao[i].energy` on an OPEN shell is gennbo's ALPHA energy, not a total: the
   open-shell composite table prints Occupancy and Spin and has no Energy
   column at all.  Measured: alpha on 124/124 rows of cf3, beta on 0/124.
   So this file never reads `energy` off a composite record.
3. A blank line inside gennbo's NAO table separates ATOMS, it does not end the
   table.  A reader that stops at the first blank returns the first atom only
   (31 of ch3's 49 rows).  This file reads the JSON, not the text, but any
   re-parse of the `.nbo` text must carry the same rule.

Floors: a sum of n printed 5-decimal occupancies has floor n * 0.5e-5.  Each
printed number is quoted with the floor that belongs to it, never a bare 1e-5.
"""
import glob
import json
import os
import sys

DP = 0.5e-5  # half of the last printed decimal place of a 5-dp occupancy


def load(path):
    with open(path) as fh:
        return json.load(fh)


def key(r):
    """Basis-function identity: everything except the class, which is the
    quantity under test.  Never the row's position in the list."""
    return (r["index"], r["element"], r["atom"], r["lang"])


def pair(gen, nat):
    """Match two NAO lists on basis-function identity.  Returns matched pairs
    and the counts that did not match, so the sample size is always stated."""
    gi = {key(r): r for r in gen}
    ni = {key(r): r for r in nat}
    common = [k for k in gi if k in ni]
    return [(gi[k], ni[k]) for k in sorted(common)], len(gi), len(ni)


def ryd_total(rows):
    """Total occupancy in rows THAT LIST CLASS Ryd, and how many there are."""
    sel = [r for r in rows if r["type"] == "Ryd"]
    return sum(r["occupancy"] for r in sel), len(sel)


def one(d, mol, arms):
    gen = load(os.path.join(d, mol, mol + ".gennbo.nbo.json"))
    if not gen.get("open_shell"):
        return None
    out = {"mol": mol, "arms": {}}
    for arm, suffix in arms:
        p = os.path.join(d, mol, mol + suffix)
        if not os.path.exists(p):
            continue
        nat = load(p)
        assert nat.get("open_shell"), f"{mol}/{arm}: native side is not open shell"
        per = {}
        for spin in ("alpha", "beta"):
            k = "nao_" + spin
            assert k in gen and k in nat, f"{mol}: {k} missing"
            pairs, ng, nn = pair(gen[k], nat[k])
            # Rydberg population, each side by its OWN classification
            rg, cg = ryd_total([g for g, _ in pairs])
            rn, cn = ryd_total([n for _, n in pairs])
            # class disagreements, and the occupancy riding on them
            mis = [(g, n) for g, n in pairs if g["type"] != n["type"]]
            # worst single-row occupancy delta, matched row to row
            worst = max((abs(g["occupancy"] - n["occupancy"]), key(g))
                        for g, n in pairs) if pairs else (0.0, None)
            per[spin] = dict(
                n=len(pairs), n_gen=ng, n_nat=nn,
                ryd_gen=rg, ryd_nat=rn, d_ryd=rn - rg,
                cnt_gen=cg, cnt_nat=cn,
                floor=max(cg, cn) * DP,
                mis=len(mis),
                mis_occ=sum(abs(g["occupancy"] - n["occupancy"]) for g, n in mis),
                worst=worst[0], worst_at=worst[1],
                worst_floor=2 * DP,
            )
        # the composite, for contrast with the per-spin split
        pairs, _, _ = pair(gen["nao"], nat["nao"])
        rg, cg = ryd_total([g for g, _ in pairs])
        rn, cn = ryd_total([n for _, n in pairs])
        per["total"] = dict(n=len(pairs), ryd_gen=rg, ryd_nat=rn, d_ryd=rn - rg,
                            cnt_gen=cg, cnt_nat=cn, floor=max(cg, cn) * DP,
                            mis=sum(1 for g, n in pairs if g["type"] != n["type"]))
        out["arms"][arm] = per
    return out


def selfcheck(d):
    """The invariants that would break if the per-spin sets were not what this
    file claims.  Assert-only; prints one line per molecule."""
    mols = sorted(os.path.basename(os.path.dirname(p))
                  for p in glob.glob(os.path.join(d, "*", "*.gennbo.nbo.json")))
    n_open = 0
    for mol in mols:
        gen = load(os.path.join(d, mol, mol + ".gennbo.nbo.json"))
        if not gen.get("open_shell"):
            assert "nao_alpha" not in gen or not gen["nao_alpha"], \
                f"{mol}: closed shell carries a non-empty nao_alpha"
            continue
        n_open += 1
        a, b, t = gen["nao_alpha"], gen["nao_beta"], gen["nao"]
        assert len(a) == len(b) == len(t), f"{mol}: {len(a)}/{len(b)}/{len(t)}"
        # alpha + beta must reproduce the composite OCCUPANCY column
        worst = max(abs(x["occupancy"] + y["occupancy"] - z["occupancy"])
                    for x, y, z in zip(a, b, t))
        assert worst <= 3 * DP, f"{mol}: alpha+beta != composite by {worst}"
        # ... and the composite 'energy' must be the ALPHA energy, not a total,
        # which is trap 2 above stated as an executable claim
        if t[0].get("energy") is not None:
            same_a = sum(abs(x.get("energy", 0) - z.get("energy", 0)) < DP
                         for x, z in zip(a, t))
            assert same_a == len(t), \
                f"{mol}: composite energy is not the alpha energy ({same_a}/{len(t)})"
        print(f"  {mol:10s} open shell, {len(a)} rows/spin, "
              f"alpha+beta==composite within {worst / DP:.1f} floors, "
              f"composite energy == alpha on {len(t)}/{len(t)}")
    assert n_open >= 4, f"expected at least 4 open shells, found {n_open}"
    print(f"selfcheck OK: {n_open} open-shell molecules, {len(mols)} total")


def main(argv):
    if argv and argv[0] == "--selfcheck":
        return selfcheck(argv[1])
    d = argv[0]
    mols = argv[1:] or sorted(os.path.basename(os.path.dirname(p))
                              for p in glob.glob(os.path.join(d, "*", "*.gennbo.nbo.json")))
    arms = [("baseline", ".baseline.native.nbo.json"), ("renat5", ".native.nbo.json")]
    res = [r for r in (one(d, m, arms) for m in mols) if r]
    print("Per-spin Rydberg population, native - gennbo, each side classed by itself.")
    print("NOT an agreement test: the gate is an exact basis-function match and")
    print("native fails it on every molecule of both sets.\n")
    print(f"{'mol':8s} {'arm':9s} {'spin':6s} {'n':>4s} {'nRyd g/n':>9s} "
          f"{'d(Ryd)':>9s} {'floor':>8s} {'xfloor':>7s} {'misclass':>8s} {'worst row':>9s}")
    for r in res:
        for arm in ("baseline", "renat5"):
            if arm not in r["arms"]:
                continue
            for spin in ("alpha", "beta", "total"):
                p = r["arms"][arm][spin]
                xf = abs(p["d_ryd"]) / p["floor"] if p["floor"] else float("nan")
                print(f"{r['mol']:8s} {arm:9s} {spin:6s} {p['n']:4d} "
                      f"{p['cnt_gen']:4d}/{p['cnt_nat']:<4d} {p['d_ryd']:+9.5f} "
                      f"{p['floor']:8.1e} {xf:7.1f} {p['mis']:8d} "
                      f"{p.get('worst', float('nan')):9.5f}")
        print()
    # the question this file was written for
    print("cf3 and the other open shells, arm to arm, per spin:")
    print(f"{'mol':8s} {'spin':6s} {'base d(Ryd)':>12s} {'ren5 d(Ryd)':>12s} "
          f"{'factor':>8s} {'floor':>8s}")
    for r in res:
        if "baseline" not in r["arms"] or "renat5" not in r["arms"]:
            continue
        for spin in ("alpha", "beta", "total"):
            b = r["arms"]["baseline"][spin]
            n = r["arms"]["renat5"][spin]
            f = abs(b["d_ryd"]) / abs(n["d_ryd"]) if n["d_ryd"] else float("inf")
            flag = ""
            if abs(n["d_ryd"]) > abs(b["d_ryd"]):
                flag = "  <-- WORSE under renat5"
            if max(abs(b["d_ryd"]), abs(n["d_ryd"])) < n["floor"]:
                flag += "  (both below floor: no signal)"
            print(f"{r['mol']:8s} {spin:6s} {b['d_ryd']:+12.5f} {n['d_ryd']:+12.5f} "
                  f"{f:8.2f} {n['floor']:8.1e}{flag}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]) or 0)
