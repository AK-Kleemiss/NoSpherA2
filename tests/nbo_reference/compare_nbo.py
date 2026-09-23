"""Per-quantity agreement between a stored NBO reference and a candidate result.

    python compare_nbo.py water.nbo.json candidate.json
    python compare_nbo.py --all <candidate_dir>     # every molecule of the dataset

Both files are the JSON that `NoSpherA2 -nbo` writes (see nbo_run.h for the structures);
an in-house implementation only has to emit the same keys. The reference drives: every
reference entry must have a counterpart in the candidate, extra candidate entries are
ignored. Exit code 0 means every quantity agreed within tolerance.

The C++ side of the same check is compare_nbo_results() in Src/core/nbo_run.h, used by
tests/src/Nbo47Tests.cpp; this script is for comparing against the stored dataset, which
lives next to it as <molecule>.nbo.json with index.json describing the set.
"""
import json
import os
import sys

TOL = {  # mirrors NboTolerances in Src/core/nbo_run.h
    "charge": 2.0e-3,
    "occupancy": 2.0e-3,
    "energy": 2.0e-3,        # a.u.
    "hybrid_percent": 0.5,   # percentage points
    "e2_kcal": 0.1,
    "bond_order": 5.0e-3,
    "weight_percent": 0.5,
}


class Quantity:
    def __init__(self, name, tol):
        self.name, self.tol = name, tol
        self.compared = self.missing = self.failed = 0
        self.max_dev, self.worst = 0.0, ""

    def check(self, label, ref, cand):
        self.compared += 1
        dev = abs(float(ref) - float(cand))
        if dev > self.max_dev:
            self.max_dev, self.worst = dev, label
        if dev > self.tol:
            self.failed += 1

    def ok(self):
        return self.failed == 0 and self.missing == 0

    def line(self):
        return "  %-18s %4d compared %3d missing %3d failed  max dev %10.5f  %s" % (
            self.name, self.compared, self.missing, self.failed, self.max_dev, self.worst)


def compare_keyed(q, ref_items, cand_items, key, fields):
    """fields: list of (json key, label suffix). Entries are matched by key()."""
    index = {}
    for it in cand_items:
        index.setdefault(key(it), it)
    for it in ref_items:
        k = key(it)
        other = index.get(k)
        if other is None:
            q.missing += 1
            continue
        for f, suffix in fields:
            if f in it:
                q.check("%s %s" % (k, suffix), it[f], other.get(f, float("nan")))


def compare(reference, candidate):
    qs = []

    def add(name, tol):
        qs.append(Quantity(name, TOL[tol]))
        return qs[-1]

    compare_keyed(add("npa charge", "charge"), reference["npa"], candidate["npa"],
                  lambda a: "atom %d %s" % (a["atom"], a["element"]),
                  [("charge", "charge"), ("core", "core"), ("valence", "val"),
                   ("rydberg", "ryd"), ("spin_density", "spin")])
    compare_keyed(add("nao occupancy", "occupancy"), reference["nao"], candidate["nao"],
                  lambda a: "nao %d" % a["index"], [("occupancy", "occ")])
    compare_keyed(add("nao energy", "energy"), reference["nao"], candidate["nao"],
                  lambda a: "nao %d" % a["index"], [("energy", "E")])

    key_nbo = lambda o: "%s%s" % (o["description"], (" " + o["spin"]) if o["spin"] else "")
    compare_keyed(add("nbo occupancy", "occupancy"), reference["nbos"], candidate["nbos"],
                  key_nbo, [("occupancy", "occ")])
    compare_keyed(add("nbo energy", "energy"), reference["nbos"], candidate["nbos"],
                  key_nbo, [("energy", "E")])

    # hybridisations: flattened to one entry per (nbo, hybrid position)
    def hybrids(res):
        out = []
        for o in res["nbos"]:
            for n, h in enumerate(o.get("hybrids", [])):
                e = dict(h)
                e["_key"] = "%s h%d" % (key_nbo(o), n)
                out.append(e)
        return out

    compare_keyed(add("hybrid %s/%p", "hybrid_percent"), hybrids(reference), hybrids(candidate),
                  lambda h: h["_key"],
                  [("s", "%s"), ("p", "%p"), ("d", "%d"), ("f", "%f"),
                   ("weight_percent", "pol")])

    compare_keyed(add("E2", "e2_kcal"), reference["e2"], candidate["e2"],
                  lambda e: "%s -> %s%s" % (e["donor"], e["acceptor"],
                                            (" " + e["spin"]) if e["spin"] else ""),
                  [("energy_kcal", "E(2)")])

    rn, cn = reference.get("nrt", {}), candidate.get("nrt", {})
    # Structures are matched by what changed relative to the reference structure, not by their
    # number: starting the search from a different Lewis structure renumbers the list and swaps
    # the descriptions of equivalent structures without moving a single weight.
    def wkey(w):
        return "%s%s" % (" ".join(w.get("changes", "").split()) or "(leading)",
                         (" " + w["spin"]) if w["spin"] else "")

    compare_keyed(add("NRT weight", "weight_percent"), rn.get("weights", []), cn.get("weights", []),
                  wkey, [("weight_percent", "%")])

    # And by rank, which survives that relabelling. This is the one that has to hold; the
    # by-description comparison above tells you whether the labels also match.
    def ranked(n):
        out = []
        for spin in sorted({w["spin"] for w in n.get("weights", [])}):
            ws = sorted((w for w in n["weights"] if w["spin"] == spin),
                        key=lambda w: -w["weight_percent"])
            for i, w in enumerate(ws):
                out.append({"_key": "rank %d%s" % (i + 1, (" " + spin) if spin else ""),
                            "weight_percent": w["weight_percent"]})
        return out

    compare_keyed(add("NRT weight by rank", "weight_percent"), ranked(rn), ranked(cn),
                  lambda w: w["_key"], [("weight_percent", "%")])
    compare_keyed(add("NRT valency", "bond_order"), rn.get("valencies", []), cn.get("valencies", []),
                  lambda v: "atom %d%s" % (v["atom"], (" " + v["spin"]) if v["spin"] else ""),
                  [("valency", "val"), ("covalency", "cov"), ("electrovalency", "ion"),
                   ("electron_count", "N")])
    compare_keyed(add("NRT bond order", "bond_order"), rn.get("bond_orders", []),
                  cn.get("bond_orders", []),
                  lambda b: "%d-%d%s" % (b["atom1"], b["atom2"],
                                         (" " + b["spin"]) if b["spin"] else ""),
                  [("total", "tot"), ("covalent", "cov"), ("ionic", "ion")])
    return qs


def report(name, qs, notes):
    ok = all(q.ok() for q in qs)
    print("%-14s %s" % (name, "PASS" if ok else "FAIL"))
    for n in notes:
        print("  note: " + n)
    for q in qs:
        if q.compared or q.missing:
            print(q.line())
    return ok


def run(ref_path, cand_path):
    with open(ref_path) as f:
        reference = json.load(f)
    with open(cand_path) as f:
        candidate = json.load(f)
    notes = []
    if reference.get("open_shell") != candidate.get("open_shell"):
        notes.append("open_shell differs: reference %s, candidate %s"
                     % (reference.get("open_shell"), candidate.get("open_shell")))
    rt, ct = reference.get("thresholds", {}), candidate.get("thresholds", {})
    for k in sorted(rt):
        if k in ct and rt[k] != ct[k]:
            notes.append("threshold %s: reference %s, candidate %s" % (k, rt[k], ct[k]))
    return report(reference.get("name", ref_path), compare(reference, candidate), notes)


def main(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    if len(argv) == 3 and argv[1] == "--all":
        names = [f[:-len(".nbo.json")] for f in sorted(os.listdir(here))
                 if f.endswith(".nbo.json")]
        ok, skipped = True, []
        for n in names:
            cand = os.path.join(argv[2], n + ".nbo.json")
            if not os.path.exists(cand):
                skipped.append(n)
                continue
            ok &= run(os.path.join(here, n + ".nbo.json"), cand)
        if skipped:
            print("no candidate for: " + " ".join(skipped))
        return 0 if ok else 1
    if len(argv) != 3:
        print(__doc__)
        return 2
    return 0 if run(argv[1], argv[2]) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
