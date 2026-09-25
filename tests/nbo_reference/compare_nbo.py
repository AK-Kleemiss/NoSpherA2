"""Per-quantity agreement between a stored NBO reference and a candidate result.

    python compare_nbo.py water.nbo.json candidate.json
    python compare_nbo.py --all <candidate_dir>     # every molecule of the dataset

<candidate_dir> is the directory holding the candidate <molecule>.nbo.json files; the
references are always read from the directory this script lives in. In --all mode a reference
molecule with no candidate file is a failure, so pointing the gate at the wrong directory
cannot pass by comparing nothing - and <candidate_dir> may not BE the reference directory,
because `--all .` run from tests/nbo_reference compared the dataset with itself and reported
22 of 22 at max dev 0.00000. A check that passes when handed its own answers cannot see the
failure mode it exists for, so it is refused rather than reported as a pass.

Both files are the JSON that `NoSpherA2 -nbo` writes (see nbo_run.h for the structures);
an in-house implementation only has to emit the same keys. The reference drives, but the
comparison is symmetric on membership: a reference entry with no counterpart is `missing`, a
candidate entry with no counterpart is `surplus`, and either fails. Exit code 0 means every
quantity agreed within tolerance and neither side printed anything the other did not.

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
    "weight_fraction": 5.0e-3,   # the same 0.5 percentage points, on the 5-decimal value
}


class Quantity:
    def __init__(self, name, tol):
        self.name, self.tol = name, tol
        self.compared = self.missing = self.surplus = self.failed = 0
        self.max_dev, self.worst = 0.0, ""

    def check(self, label, ref, cand):
        self.compared += 1
        dev = abs(float(ref) - float(cand))
        if dev > self.max_dev:
            self.max_dev, self.worst = dev, label
        if dev > self.tol:
            self.failed += 1

    def ok(self):
        return self.failed == 0 and self.missing == 0 and self.surplus == 0

    def line(self):
        return ("  %-18s %4d compared %3d missing %3d surplus %3d failed  max dev %10.5f  %s"
                % (self.name, self.compared, self.missing, self.surplus, self.failed,
                   self.max_dev, self.worst))


def compare_keyed(q, ref_items, cand_items, key, fields):
    """fields: list of (json key, label suffix). Entries are matched by key(). A key is not
    unique - two resonance structures can carry the same Added(Removed) description, and
    acetylene has four such pairs - so entries sharing a key pair up in order of appearance
    instead of all comparing against the first one. Without that a file fails against itself.

    Unmatched entries are counted on BOTH sides. The loop walks the reference, so for a long
    time an entry the candidate printed and the reference lacks was neither `compared` nor
    `missing` and left no trace at all - twelve extra native E2 rows on ch3 went unreported by
    the very check that exists to find them. A candidate that invents entries disagrees with
    the reference exactly as much as one that omits them, so the leftovers are `surplus` and
    `ok()` fails on them."""
    index = {}
    for it in cand_items:
        index.setdefault(key(it), []).append(it)
    taken = {}
    for it in ref_items:
        k = key(it)
        queue, n = index.get(k, []), taken.get(k, 0)
        if n >= len(queue):
            q.missing += 1
            continue
        taken[k], other = n + 1, queue[n]
        for f, suffix in fields:
            if f in it:
                q.check("%s %s" % (k, suffix), it[f], other.get(f, float("nan")))
    q.surplus += sum(len(queue) - taken.get(k, 0) for k, queue in index.items())


def compare(reference, candidate):
    qs = []

    def add(name, tol):
        qs.append(Quantity(name, TOL[tol]))
        return qs[-1]

    compare_keyed(add("npa charge", "charge"), reference["npa"], candidate["npa"],
                  lambda a: "atom %d %s" % (a["atom"], a["element"]),
                  [("charge", "charge"), ("core", "core"), ("valence", "val"),
                   ("rydberg", "ryd"), ("spin_density", "spin")])
    # NBO prints every Rydberg NAO the basis can form, and the empty ones - 0.00000 occupancy to
    # all five printed decimals - carry no information: nothing downstream of the NAO table reads
    # them (collect_reference.py says so at its NAO section), and their energies are basis tails
    # that two codes have no reason to place identically. They were 89 % of this gate's failures.
    def occupied(res):
        return [n for n in res["nao"]
                if not (n["type"] == "Ryd" and abs(n["occupancy"]) < 1.0e-5)]

    compare_keyed(add("nao occupancy", "occupancy"), occupied(reference), occupied(candidate),
                  lambda a: "nao %d" % a["index"], [("occupancy", "occ")])
    compare_keyed(add("nao energy", "energy"), occupied(reference), occupied(candidate),
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

    # Compare the 5-decimal weight vector, not the 2-decimal table column: it is the same number
    # with three more digits, and it is the one NRT actually minimised.
    def weights(n):
        out = []
        for w in n.get("weights", []):
            e = dict(w)
            e["weight"] = w["weight_fraction"] or w["weight_percent"] / 100.0
            out.append(e)
        return out

    compare_keyed(add("NRT weight", "weight_fraction"), weights(rn), weights(cn),
                  wkey, [("weight", "w")])

    # And by rank, which survives that relabelling. This is the one that has to hold; the
    # by-description comparison above tells you whether the labels also match.
    def ranked(n):
        out = []
        all_w = weights(n)
        for spin in sorted({w["spin"] for w in all_w}):
            ws = sorted((w for w in all_w if w["spin"] == spin), key=lambda w: -w["weight"])
            for i, w in enumerate(ws):
                out.append({"_key": "rank %d%s" % (i + 1, (" " + spin) if spin else ""),
                            "weight": w["weight"]})
        return out

    compare_keyed(add("NRT weight by rank", "weight_fraction"), ranked(rn), ranked(cn),
                  lambda w: w["_key"], [("weight", "w")])
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
        if q.compared or q.missing or q.surplus:
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


def selftest():
    """`python compare_nbo.py --selftest` - the membership check, both directions.

    The interesting case is the one that used to pass: an identical file plus one extra
    candidate entry. Nothing is out of tolerance, nothing is missing, and before Quantity
    gained `surplus` the whole comparison reported PASS."""
    base = {"name": "t", "npa": [{"atom": 1, "element": "H", "charge": 0.1}],
            "nao": [], "nbos": [], "nrt": {},
            "e2": [{"donor": "BD 1", "acceptor": "BD* 2", "spin": "alpha",
                    "energy_kcal": 1.0}]}
    same = json.loads(json.dumps(base))
    qs = {q.name: q for q in compare(base, same)}
    assert all(q.ok() for q in qs.values()), "a file must agree with itself"

    extra = json.loads(json.dumps(base))
    extra["e2"].append({"donor": "BD 1", "acceptor": "BD* 3", "spin": "alpha",
                        "energy_kcal": 0.4})
    qs = {q.name: q for q in compare(base, extra)}
    assert qs["E2"].surplus == 1, qs["E2"].line()
    assert qs["E2"].missing == 0 and qs["E2"].failed == 0, qs["E2"].line()
    assert not qs["E2"].ok(), "a surplus entry must fail, not be a footnote"

    short = json.loads(json.dumps(base))
    short["e2"] = []
    qs = {q.name: q for q in compare(base, short)}
    assert (qs["E2"].missing, qs["E2"].surplus) == (1, 0), qs["E2"].line()

    # Two entries sharing a key pair up in order of appearance, so duplicates are not a surplus.
    dup = {"name": "t", "npa": [], "nao": [], "nbos": [], "nrt": {},
           "e2": [{"donor": "BD 1", "acceptor": "BD* 2", "spin": "", "energy_kcal": 1.0},
                  {"donor": "BD 1", "acceptor": "BD* 2", "spin": "", "energy_kcal": 2.0}]}
    qs = {q.name: q for q in compare(dup, json.loads(json.dumps(dup)))}
    assert qs["E2"].ok() and qs["E2"].surplus == 0, qs["E2"].line()
    import socket
    print("compare_nbo selftest OK on %s" % socket.gethostname())
    return 0


def main(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    if len(argv) == 2 and argv[1] == "--selftest":
        return selftest()
    if len(argv) == 3 and argv[1] == "--all":
        if os.path.realpath(argv[2]) == os.path.realpath(here):
            print("refusing --all %s: that is the reference directory itself, so every molecule\n"
                  "would be compared with its own file and pass at max dev 0.00000. Point it at\n"
                  "the directory holding the candidate results." % argv[2])
            return 2
        names = [f[:-len(".nbo.json")] for f in sorted(os.listdir(here))
                 if f.endswith(".nbo.json")]
        passed, failed, missing = 0, 0, []
        for n in names:
            cand = os.path.join(argv[2], n + ".nbo.json")
            if not os.path.exists(cand):
                missing.append(n)
                continue
            if run(os.path.join(here, n + ".nbo.json"), cand):
                passed += 1
            else:
                failed += 1
        if missing:
            print("no candidate in %s for: %s" % (argv[2], " ".join(missing)))
        print("%d of %d reference molecules PASS, %d FAIL, %d without a candidate file"
              % (passed, len(names), failed, len(missing)))
        return 0 if passed and not failed and not missing else 1
    if len(argv) != 3:
        print(__doc__)
        return 2
    return 0 if run(argv[1], argv[2]) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
