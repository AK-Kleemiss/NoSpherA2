"""Discrimination gate: exit 0 only if two native JSONs differ in a NUMERIC field.

    python3 armdisc.py <a.native.nbo.json> <b.native.nbo.json>

Timings are excluded on purpose - two runs of the SAME binary differ in `timings` and in
`nrt.seconds` and in nothing else, so a differ-anywhere test would pass on a knob that does
nothing.  What is compared is the NPA table (charge/core/valence/rydberg/total, per atom, read
by its own key) and the NAO occupancies if present.  Prints what it compared and how many
values it had, because an instrument that does not state its own sample size can pass on zero
rows.
"""
import json
import sys

FIELDS = ("charge", "core", "valence", "rydberg", "total")


def numbers(path):
    d = json.load(open(path))
    out = {}
    for r in d.get("npa") or []:
        if not isinstance(r, dict):
            continue
        key = (r.get("atom"), r.get("element"))
        for f in FIELDS:
            v = r.get(f)
            if isinstance(v, (int, float)):
                out[("npa", key, f)] = float(v)
    for i, r in enumerate(d.get("nao") or []):
        if not isinstance(r, dict):
            continue
        for f in ("occupancy", "energy"):
            v = r.get(f)
            if isinstance(v, (int, float)):
                out[("nao", r.get("nao", i), r.get("atom"), f)] = float(v)
    return out


a, b = numbers(sys.argv[1]), numbers(sys.argv[2])
common = sorted(set(a) & set(b), key=repr)
if not common:
    print("armdisc: NO comparable numeric fields (a=%d b=%d) - cannot discriminate" % (len(a), len(b)))
    sys.exit(4)
worst, where = 0.0, None
for k in common:
    d = abs(a[k] - b[k])
    if d > worst:
        worst, where = d, k
print("armdisc: compared %d numeric fields (a=%d b=%d), worst |delta| %.3e at %s"
      % (len(common), len(a), len(b), worst, where))
# 1e-5 is gennbo's print floor; the arms must differ by more than the floor to be distinguishable
sys.exit(0 if worst > 1e-5 else 1)
