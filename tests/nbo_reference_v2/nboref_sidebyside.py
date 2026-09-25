"""Side by side, gennbo 7 against -nbo_native, for one molecule.

    python nboref_sidebyside.py <dir> <mol>

Printed so the pattern is visible rather than summarised: where the two codes agree and
where they do not, per orbital type, because a per-quantity max deviation cannot tell a
whole-quantity disagreement from a Rydberg tail.
"""
import json
import os
import sys

d, mol = sys.argv[1], sys.argv[2]
g = json.load(open(os.path.join(d, mol + ".gennbo.nbo.json")))
n = json.load(open(os.path.join(d, mol + ".native.nbo.json")))

print("== top-level keys")
print("gennbo:", sorted(g))
print("native:", sorted(n))

print("\n== NPA")
print("%-10s %10s %10s   %10s %10s   %10s %10s" % ("atom", "q gennbo", "q native",
                                                   "core g", "core n", "val g", "val n"))
for a, b in zip(g["npa"], n["npa"]):
    print("%-10s %10.5f %10.5f   %10.5f %10.5f   %10.5f %10.5f"
          % ("%d %s" % (a["atom"], a["element"]), a["charge"], b["charge"],
             a.get("core", 0), b.get("core", 0), a.get("valence", 0), b.get("valence", 0)))

print("\n== NAO (first 20)")
print("%4s %-12s %10s %10s %12s %12s" % ("i", "label", "occ g", "occ n", "E g", "E n"))
for a, b in list(zip(g["nao"], n["nao"]))[:20]:
    print("%4d %-12s %10.5f %10.5f %12.5f %12.5f"
          % (a["index"], "%s %s" % (a.get("element", ""), a.get("type", "")),
             a["occupancy"], b["occupancy"], a.get("energy", 0), b.get("energy", 0)))

print("\n== NBO, grouped by type")
by = {}
for a, b in zip(g["nbos"], n["nbos"]):
    t = a["description"].split()[0]
    by.setdefault(t, []).append((a, b))
for t in sorted(by):
    dev_o = max(abs(a["occupancy"] - b["occupancy"]) for a, b in by[t])
    dev_e = max(abs(a.get("energy", 0) - b.get("energy", 0)) for a, b in by[t])
    print("%-6s %3d orbitals   max |d occ| %9.6f   max |d E| %12.6f" % (t, len(by[t]), dev_o, dev_e))
print("\n   first 12, in file order")
print("%-26s %9s %9s %10s %10s" % ("description", "occ g", "occ n", "E g", "E n"))
for a, b in list(zip(g["nbos"], n["nbos"]))[:12]:
    print("%-26s %9.5f %9.5f %10.5f %10.5f"
          % (a["description"] + (" " + a["spin"] if a.get("spin") else ""),
             a["occupancy"], b["occupancy"], a.get("energy", 0), b.get("energy", 0)))

print("\n== E2")
for side, r in (("gennbo", g), ("native", n)):
    print(side, len(r.get("e2", [])), "entries")
    for e in r.get("e2", [])[:8]:
        print("   %-24s -> %-24s %8.2f" % (e["donor"], e["acceptor"], e["energy_kcal"]))

print("\n== NRT")
for side, r in (("gennbo", g), ("native", n)):
    nr = r.get("nrt", {})
    print(side, "weights", [(w.get("weight_percent"), w.get("weight_fraction"),
                             w.get("changes", "")[:30], w.get("spin", "")) for w in nr.get("weights", [])][:6])
    print(side, "valencies", [(v["atom"], v.get("spin", ""), v.get("valency"), v.get("covalency"),
                               v.get("electrovalency"), v.get("electron_count"))
                              for v in nr.get("valencies", [])][:8])
    print(side, "bond orders", [(b["atom1"], b["atom2"], b.get("spin", ""), b.get("total"),
                                 b.get("covalent"), b.get("ionic"))
                                for b in nr.get("bond_orders", [])][:8])
