"""Print the NPA Rydberg column of gennbo and of two native JSONs side by side.

    python3 whichryd.py <mol> <dir_with_gennbo_and_old_native> <dir_with_new_native>

The point is to find out WHERE a 0.31 e change in native's Rydberg population came from when
the only thing that changed was which NoSpherA2 binary ran: the cascade, or what the run
reports.  Every field is read by its own key, never by position.
"""
import json
import sys


def npa(path):
    d = json.load(open(path))
    rows = d.get("npa_charges") or d.get("npa") or []
    out = {}
    for r in rows:
        if not isinstance(r, dict):
            continue
        key = (r.get("atom"), r.get("element") or r.get("symbol"))
        out[key] = {k: r.get(k) for k in ("charge", "core", "valence", "rydberg", "total")}
    return d, out


mol, dold, dnew = sys.argv[1], sys.argv[2], sys.argv[3]
dg, g = npa("%s/%s.gennbo.nbo.json" % (dold, mol))
do, o = npa("%s/%s.native.nbo.json" % (dold, mol))
dn, n = npa("%s/%s.native.nbo.json" % (dnew, mol))
print("parser_version  gennbo=%s old_native=%s new_native=%s" % (
    dg.get("parser_version"), do.get("parser_version"), dn.get("parser_version")))
for k in ("nao_count", "n_electrons", "electrons", "spin_resolved", "open_shell"):
    print("  %-14s gennbo=%s old=%s new=%s" % (k, dg.get(k), do.get(k), dn.get(k)))
print("%-10s %10s %10s %10s   %10s %10s %10s" % ("atom", "Ryd_gen", "Ryd_old", "Ryd_new",
                                                 "Val_gen", "Val_old", "Val_new"))
keys = sorted(set(g) | set(o) | set(n), key=lambda t: (t[0] is None, t[0]))
sg = so = sn = 0.0
for k in keys:
    rg, ro, rn = g.get(k, {}), o.get(k, {}), n.get(k, {})
    for src, tot in ((rg, "g"), (ro, "o"), (rn, "n")):
        pass
    sg += rg.get("rydberg") or 0.0
    so += ro.get("rydberg") or 0.0
    sn += rn.get("rydberg") or 0.0
    print("%-10s %10s %10s %10s   %10s %10s %10s" % (
        k, rg.get("rydberg"), ro.get("rydberg"), rn.get("rydberg"),
        rg.get("valence"), ro.get("valence"), rn.get("valence")))
print("%-10s %10.5f %10.5f %10.5f" % ("SUM Ryd", sg, so, sn))
print("d(Ryd) old-gennbo %+.5f   new-gennbo %+.5f   new-old %+.5f" % (so - sg, sn - sg, sn - so))
