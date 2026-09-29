"""The same d(Ryd) summed two legitimate ways, with the floor each one earns.

A sum of k printed 5-decimal values has floor k * 0.5e-5.  The per-ATOM NPA
Rydberg column adds N_atoms values; the per-NAO Ryd rows add n_Ryd values.
Both are the same physical quantity.  The per-NAO sum is the LOOSER instrument
because it adds far more printed numbers - so a regression can be real on one
and invisible on the other, and the convention has to be stated with the number.
"""
import json
import os
import sys

DP = 0.5e-5


def load(d, mol, suffix):
    with open(os.path.join(d, mol, mol + suffix)) as fh:
        return json.load(fh)


def per_atom(j):
    rows = j["npa"]
    return sum(r["rydberg"] for r in rows), len(rows)


def per_nao(j):
    sel = [r for r in j["nao"] if r["type"] == "Ryd"]
    return sum(r["occupancy"] for r in sel), len(sel)


def main(d, mols):
    print(f"{'mol':10s} {'conv':8s} {'k':>4s} {'floor':>9s} "
          f"{'base':>10s} {'ren5':>10s} {'base/fl':>8s} {'ren5/fl':>8s} {'verdict'}")
    flips = []
    for mol in mols:
        g = load(d, mol, ".gennbo.nbo.json")
        b = load(d, mol, ".baseline.native.nbo.json")
        n = load(d, mol, ".native.nbo.json")
        seen, loose = {}, 0.0
        for name, fn in (("per-atom", per_atom), ("per-NAO", per_nao)):
            gv, gk = fn(g)
            bv, bk = fn(b)
            nv, nk = fn(n)
            k = max(gk, bk, nk)
            fl = k * DP
            loose = max(loose, fl)
            db, dn = bv - gv, nv - gv
            if abs(db) < fl and abs(dn) < fl:
                verdict = "both below floor: NO SIGNAL"
            elif abs(dn) > abs(db):
                verdict = "worse under renat5"
            else:
                verdict = "better under renat5"
            seen[name] = (db, dn, verdict)
            print(f"{mol:10s} {name:8s} {k:4d} {fl:9.1e} {db:+10.5f} {dn:+10.5f} "
                  f"{abs(db) / fl:8.1f} {abs(dn) / fl:8.1f} {verdict}")
        # INVARIANT: the two conventions sum the same physical quantity, so the
        # VALUES must agree to within the looser floor.  If they do not, one of
        # them is reading the wrong column and no verdict from either is usable.
        for i in (0, 1):
            spread = abs(seen["per-atom"][i] - seen["per-NAO"][i])
            assert spread <= loose, (
                f"{mol}: the two conventions disagree by {spread:.5f} > "
                f"looser floor {loose:.1e} - one is reading the wrong column")
        if seen["per-atom"][2] != seen["per-NAO"][2]:
            flips.append((mol, seen["per-atom"][2], seen["per-NAO"][2]))
        print()
    if flips:
        print("VERDICT DEPENDS ON THE SUMMATION for:")
        for mol, a, b in flips:
            print(f"  {mol:10s} per-atom says {a!r}, per-NAO says {b!r}")
        print("  Such a molecule is not a robust counter-example: quote the")
        print("  convention with the number, or do not quote the number.")
    else:
        print("every molecule gets the same verdict under both summations")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
