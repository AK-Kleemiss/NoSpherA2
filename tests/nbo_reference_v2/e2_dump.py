"""Both full E2 tables side by side, grouped by spin, for one molecule directory."""
import json
import os
import socket
import sys

d, mol = sys.argv[1], sys.argv[2]
g = json.load(open(os.path.join(d, mol + ".gennbo.nbo.json")))
n = json.load(open(os.path.join(d, mol + ".native.nbo.json")))
print("e2_dump on %s, %s" % (socket.gethostname(), mol))
for side, r in (("gennbo", g), ("native", n)):
    e2 = r.get("e2", [])
    print("\n== %s: %d entries, thresholds %r" % (side, len(e2), r.get("thresholds")))
    byspin = {}
    for e in e2:
        byspin.setdefault(e.get("spin") or "(closed)", []).append(e)
    for spin in sorted(byspin):
        rows = sorted(byspin[spin], key=lambda e: -e["energy_kcal"])
        print("   %s: %d entries, sum %.4f kcal" % (spin, len(rows),
                                                   sum(e["energy_kcal"] for e in rows)))
        for e in rows:
            print("      %-26s -> %-26s %8.4f  dE=%s F=%s"
                  % (e["donor"], e["acceptor"], e["energy_kcal"],
                     e.get("energy_diff"), e.get("fock")))
