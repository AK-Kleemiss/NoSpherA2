"""Is the open-shell disagreement a factor of two? Report the ratio, per quantity, per pair.

    python ratio_probe.py <mol_dir> <mol>   [--expect 2.0]

"native is exactly twice gennbo" is a claim about every matched pair, and a summary that says
"exactly twice" cannot distinguish 12 pairs at 2.00 from 11 at 2.00 and one at 1.6. So this
prints native/gennbo for every pair it can match and then the spread of those ratios, and it
counts how many are inside 1 % of the expected value and how many are not. A quantity that is
off by something other than two cannot hide inside the summary.

Two ratios have no meaning and are not reported as agreement:

  - a pair where gennbo's value is ~0 (the ratio is unbounded, so the absolute difference is
    printed instead);
  - a quantity only one side printed - that is a `missing`, reported as such, because an empty
    gennbo E2 table against 12 native entries is *consistent* with a factor of two and proves
    nothing. It is why the E2 threshold has to be lowered until gennbo prints its own numbers.

Covalent/ionic shares are computed and printed separately, because a share is a ratio of two
quantities that would both be doubled: it survives the factor and any disagreement in it is a
second, independent finding.
"""
import argparse
import json
import os
import socket
import sys

NEAR_ZERO = 1.0e-6

# gennbo's last printed digit per block, read off its own output format. This is a statement
# about NBO's report, not a tolerance anybody chose: an E2 table printed to two decimals cannot
# resolve a ratio better than a few percent no matter how exact the underlying numbers are.
RESOLUTION = {"E2": 0.01, "NRT valency": 0.001, "NRT bond order": 0.001,
              "NBO occupancy": 0.00001, "NPA": 0.00001}


def ratios(pairs, expect, resolution=0.0, tol_frac=0.01):
    """pairs: list of (label, gennbo, native). Returns rows plus a verdict summary.

    The test is on the residual `native/expect - gennbo`, not on the ratio, because gennbo's
    printed precision sets the floor: its E2 table carries two decimals, so a printed 0.20 is
    anything in [0.195, 0.205) and the ratio against native's 0.39517 is 1.976 whatever the
    true relationship is. Calling that a 1.2 % failure would be measuring NBO's output format.
    `resolution` is gennbo's last printed digit for that block; a pair passes when the residual
    is inside half of it (plus a relative term for the larger values). The ratio is still
    printed, because it is what a reader wants to see.
    """
    rows, good, bad, degenerate = [], 0, 0, 0
    for label, g, n in pairs:
        resid = n / expect - g
        allow = resolution / 2.0 + tol_frac * abs(g)
        inside = abs(resid) <= allow
        if abs(g) < NEAR_ZERO and abs(n) < NEAR_ZERO:
            rows.append((label, g, n, None, resid, "both ~0"))
            degenerate += 1
            continue
        r = n / g if abs(g) >= NEAR_ZERO else None
        good += inside
        bad += not inside
        rows.append((label, g, n, r, resid,
                     "" if inside else "resid %+.4f > %.4f" % (resid, allow)))
    return rows, good, bad, degenerate


def collect(g, n):
    """Every (block, label, gennbo, native) triple whose ratio is worth asking about."""
    out = {}

    def keyed(block, gi, ci, key, fields):
        gd, cd = {}, {}
        for src, dst in ((gi, gd), (ci, cd)):
            for it in src:
                k = key(it)
                dst.setdefault(k, []).append(it)
        pairs, missing = [], []
        for k in sorted(set(gd) | set(cd)):
            gl, cl = gd.get(k, []), cd.get(k, [])
            for i in range(max(len(gl), len(cl))):
                if i >= len(gl) or i >= len(cl):
                    missing.append("%s %s (%s only)" % (k, i, "native" if i < len(cl) else "gennbo"))
                    continue
                for f, short in fields:
                    a, b = gl[i].get(f), cl[i].get(f)
                    if a is None or b is None:
                        continue
                    pairs.append(("%s %s" % (k, short), float(a), float(b)))
        out[block] = {"pairs": pairs, "missing": missing}

    keyed("E2", g.get("e2", []), n.get("e2", []),
          lambda e: "%s -> %s%s" % (e["donor"], e["acceptor"],
                                    (" " + e["spin"]) if e.get("spin") else ""),
          [("energy_kcal", "kcal")])
    gr, nr = g.get("nrt", {}), n.get("nrt", {})
    keyed("NRT valency", gr.get("valencies", []), nr.get("valencies", []),
          lambda v: "atom %d%s" % (v["atom"], (" " + v["spin"]) if v.get("spin") else ""),
          [("valency", "val"), ("covalency", "cov"), ("electrovalency", "ion"),
           ("electron_count", "N")])
    keyed("NRT bond order", gr.get("bond_orders", []), nr.get("bond_orders", []),
          lambda b: "%d-%d%s" % (b["atom1"], b["atom2"], (" " + b["spin"]) if b.get("spin") else ""),
          [("total", "tot"), ("covalent", "cov"), ("ionic", "ion")])
    keyed("NBO occupancy", g.get("nbos", []), n.get("nbos", []),
          lambda o: "%s%s" % (o["description"], (" " + o["spin"]) if o.get("spin") else ""),
          [("occupancy", "occ")])
    keyed("NPA", g.get("npa", []), n.get("npa", []),
          lambda a: "atom %d %s" % (a["atom"], a["element"]),
          [("charge", "q"), ("core", "core"), ("valence", "val")])
    return out


def shares(res):
    """Ionic share of each atom's valency and of each bond order, per spin.

    Both numerator and denominator are per-spin quantities, so a factor of two divides out.
    Anything left here is not the factor-2 story.

    Rows whose spin is "composite" are skipped: NBO prints alpha, beta AND a composite
    alpha+beta table for an open shell, native prints only the two spins, and a composite
    share is not comparable with anything native produces. Until the parser fix that tags
    that third table correctly, it arrives labelled "beta" instead, which is why
    first_wins() below refuses to let a second row for one label replace the first.
    """
    out = []
    nr = res.get("nrt", {})
    for v in nr.get("valencies", []):
        if (v.get("spin") or "") == "composite":
            continue
        tot = v.get("valency") or 0.0
        if abs(tot) > NEAR_ZERO and v.get("electrovalency") is not None:
            out.append(("valency atom %d%s" % (v["atom"], (" " + v["spin"]) if v.get("spin") else ""),
                        100.0 * v["electrovalency"] / tot))
    for b in nr.get("bond_orders", []):
        if (b.get("spin") or "") == "composite":
            continue
        tot = b.get("total") or 0.0
        if abs(tot) > NEAR_ZERO and b.get("ionic") is not None:
            out.append(("bond %d-%d%s" % (b["atom1"], b["atom2"],
                                          (" " + b["spin"]) if b.get("spin") else ""),
                        100.0 * b["ionic"] / tot))
    return out


def first_wins(pairs, side):
    """dict(pairs) that keeps the FIRST value for a label and says so when it drops a later one.

    `dict(shares(g))` kept the last, and that silently replaced ch3's real beta valency share
    (13.61 %) with the composite table's (16.40 %), which put the headline ionic-share gap at
    -11.26 points instead of -8.46. A duplicate label here always means one of the two sides
    printed a table this script does not know how to key, so it is reported, not averaged.
    """
    out = {}
    for k, v in pairs:
        if k in out:
            print("   NOTE %s printed a second row for %-28s (%.2f %%, keeping the first, "
                  "%.2f %%)" % (side, k, v, out[k]))
            continue
        out[k] = v
    return out


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("mol")
    ap.add_argument("--expect", type=float, default=2.0)
    ap.add_argument("--json")
    a = ap.parse_args(argv[1:])

    g = json.load(open(os.path.join(a.dir, a.mol + ".gennbo.nbo.json")))
    n = json.load(open(os.path.join(a.dir, a.mol + ".native.nbo.json")))
    print("ratio_probe on %s, %s, open_shell gennbo=%r native=%r"
          % (socket.gethostname(), a.mol, g.get("open_shell"), n.get("open_shell")))
    print("thresholds gennbo=%r native=%r" % (g.get("thresholds"), n.get("thresholds")))
    print("expecting native/gennbo = %.3f, inside 1 %%\n" % a.expect)

    blocks = collect(g, n)
    summary = {}
    for name in ("E2", "NRT valency", "NRT bond order", "NBO occupancy", "NPA"):
        b = blocks[name]
        rows, good, bad, deg = ratios(b["pairs"], a.expect, RESOLUTION.get(name, 0.0))
        rs = [r[3] for r in rows if r[3] is not None]
        print("== %s: %d matched pairs, %d consistent with %.2f, %d not, %d both-zero, "
              "%d unmatched (gennbo printed to %g)"
              % (name, len(rows), good, a.expect, bad, deg, len(b["missing"]),
                 RESOLUTION.get(name, 0.0)))
        if rs:
            print("   ratio min %.4f  max %.4f  mean %.4f"
                  % (min(rs), max(rs), sum(rs) / len(rs)))
        print("   %-40s %11s %11s %8s %10s %s"
              % ("pair", "gennbo", "native", "ratio", "n/exp-g", "flag"))
        for label, gv, nv, r, resid, flag in rows:
            print("   %-40s %11.5f %11.5f %8s %+10.5f %s"
                  % (label[:40], gv, nv, "%.4f" % r if r is not None else "-", resid, flag))
        for m in b["missing"][:12]:
            print("   UNMATCHED " + m)
        if len(b["missing"]) > 12:
            print("   ... and %d more unmatched" % (len(b["missing"]) - 12))
        resids = [abs(r[4]) for r in rows]
        summary[name] = {"pairs": len(rows), "at_expect": good, "not_at_expect": bad,
                         "both_zero": deg, "unmatched": len(b["missing"]),
                         "ratio_min": min(rs) if rs else None,
                         "ratio_max": max(rs) if rs else None,
                         "worst_residual": max(resids) if resids else None,
                         "gennbo_resolution": RESOLUTION.get(name, 0.0)}
        print()

    print("== ionic share, which a factor of two cannot explain (both parts are per-spin)")
    gs, ns = first_wins(shares(g), "gennbo"), first_wins(shares(n), "native")
    print("   %-36s %10s %10s %10s" % ("quantity", "gennbo %", "native %", "diff"))
    worst, worst_label = 0.0, None
    for k in sorted(set(gs) | set(ns)):
        if k not in gs or k not in ns:
            print("   %-36s %10s %10s   one side only"
                  % (k[:36], "%.2f" % gs[k] if k in gs else "-",
                     "%.2f" % ns[k] if k in ns else "-"))
            continue
        d = ns[k] - gs[k]
        if abs(d) > abs(worst):
            worst, worst_label = d, k
        print("   %-36s %10.2f %10.2f %+10.2f" % (k[:36], gs[k], ns[k], d))
    print("   worst ionic-share difference %+.2f points at %s" % (worst, worst_label))
    summary["ionic_share"] = {"worst_points": worst, "worst_at": worst_label}

    if a.json:
        with open(a.json, "w") as f:
            json.dump({"host": socket.gethostname(), "molecule": a.mol,
                       "expect": a.expect, "summary": summary}, f, indent=1)
        print("wrote " + a.json)
    return 0


def demo():
    """One runnable check: the ratio logic must catch an outlier hiding among good pairs."""
    rows, good, bad, deg = ratios([("a", 1.0, 2.0), ("b", 0.5, 1.0), ("c", 1.0, 1.6),
                                   ("d", 0.0, 0.0)], 2.0)
    assert (good, bad, deg) == (2, 1, 1), (good, bad, deg)
    assert rows[2][5].startswith("resid") and rows[3][3] is None
    # gennbo's printed 0.20 against native's 0.39517 is 2.00 at gennbo's own resolution, and
    # calling it a 1.2 % failure would be measuring the report format, not the method.
    rows, good, bad, _ = ratios([("bd", 0.20, 0.39517)], 2.0, resolution=0.01)
    assert (good, bad) == (1, 0) and abs(rows[0][3] - 1.976) < 1e-3
    # but a genuinely different number must still fail at that same resolution
    assert ratios([("bd", 0.20, 0.30)], 2.0, resolution=0.01)[2] == 1
    # and a share must be insensitive to the doubling
    per_spin = {"nrt": {"valencies": [{"atom": 1, "spin": "alpha", "valency": 1.5,
                                       "electrovalency": 0.3}]}}
    doubled = {"nrt": {"valencies": [{"atom": 1, "spin": "alpha", "valency": 3.0,
                                      "electrovalency": 0.6}]}}
    assert abs(shares(per_spin)[0][1] - shares(doubled)[0][1]) < 1e-9
    print("demo ok on %s" % socket.gethostname())


if __name__ == "__main__":
    sys.exit(demo() if "--demo" in sys.argv else main(sys.argv))
