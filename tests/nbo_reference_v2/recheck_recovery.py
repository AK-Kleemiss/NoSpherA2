"""Re-run control 1 (recovery of a planted w) with the SAME multistart budget as the treatment.

In fitw.py's main() the recovery control is fitted from ONE start while the real fit gets STARTS,
so a recovery that misses is not yet evidence that the solver is incompetent on that arm - the two
were not given the same budget.  This re-runs the identical planted target (same SEED, same draw
order, so the same wp) at starts=STARTS and prints both numbers side by side.

    py -3.12 recheck_recovery.py <ops dir> <mol>,<arm> [<mol>,<arm> ...]

An arm can only be CERTIFIED by this, never condemned: if the recovery still misses at the same
budget the treatment had, the arm's MISS stays unclaimable.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fitw import SEED, STARTS, compose, fit_one, norm_w, prepare, verdict  # noqa: E402

d = sys.argv[1]
rows = {}
for p in sorted(__import__("glob").glob(os.path.join(d, "fitw_result*.json"))):
    for r in json.load(open(p)):
        rows[r["mol"]] = r

outp = os.path.join(d, "fitw_recheck.json")
out = json.load(open(outp)) if os.path.exists(outp) else {}

for spec in sys.argv[2:]:
    mol, arm = spec.split(",")
    renat = arm == "renat"
    g = prepare(mol, os.path.join(d, "ops_" + mol))
    # replay main()'s rng exactly: one draw per arm, in BACKBONES order, after the two WRONG verdicts
    rng = np.random.default_rng(SEED)
    for name in ("base", "renat"):
        rng.uniform(0.5, 2.0, g["n"])                       # the `moves` draw
        wp = norm_w(10.0 ** rng.uniform(-4, 0, g["n"]), g["cls"])
        if name == arm:
            break
    Ct = compose(wp, g, renat)
    rc = fit_one(g, renat, Ct, np.log(np.full(g["n"], 0.1)), starts=STARTS)
    worst = verdict(compose(rc["w"], g, renat), Ct, g)["worst"]
    v = rows[mol]["arms"][arm]
    old = v["recover"]["worst"]
    fit = v["vfit"]["worst"]
    print("%-8s %-6s floor %.2e | recovery worst 1 start %.4e -> %d starts %.4e | fit %.4e | "
          "margin fit/max(floor,rec) %.1e | %s" % (
              mol, arm, v["floor"], old, STARTS, worst, fit,
              fit / max(v["floor"], worst),
              "CERTIFIED" if fit / max(v["floor"], worst) >= 10.0 else "STILL UNCLAIMABLE"))
    sys.stdout.flush()
    out.setdefault(mol, {})[arm] = dict(worst=float(worst), res=rc["res"], starts=rc["starts"],
                                        worst_1start=old, secs=rc["secs"])
    json.dump(out, open(outp, "w"), indent=1)
print("wrote", outp)
