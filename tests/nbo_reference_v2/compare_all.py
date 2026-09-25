"""Native NoSpherA2 against gennbo 7, on the same wavefunction and the same .47 archive.

    python compare_all.py <root>  [--json summary.json] [--only mol,mol]

<root> holds one directory per molecule, each with the two JSONs stage 2 wrote:

    <mol>.gennbo.nbo.json   NoSpherA2 -nbo_parse of gennbo 7's own output
    <mol>.native.nbo.json   NoSpherA2 -nbo_native on the same .47

The two sides have independent provenance - one is the external program's output, the other
is our own analysis - which is the point. `compare_nbo.py --all .` run inside
tests/nbo_reference compares that dataset with itself and reports 22/22 at 0.00000; this
does not.

The comparison engine is tests/nbo_reference/compare_nbo.py, imported from there rather than
copied here. gennbo is the reference side, so a quantity gennbo printed and we did not is a
`missing` and one we printed and gennbo did not is a `surplus` - neither is a silent pass.
NRT weights are read off the by-rank comparison, never by label: starting the search from a
different Lewis structure renumbers the list without moving a weight.

Every row carries the hostname and thread count the native side ran with, read back from
provenance_nbo.json, because a number without its thread count and node is not a
measurement anyone can reuse.

Two views, always both printed
------------------------------
The engine already drops empty Rydberg *NAOs* and says why: nothing downstream reads them, and
their energies are basis tails two codes have no reason to place identically. The same argument
applies verbatim to empty Rydberg *NBOs* and to E2 rows accepting into one, so this adds the
matching gate - as a second column, never as a replacement. A gate that is only ever shown
after it has removed most of the failures cannot be told apart from a gate that was tuned until
the test passed, so `ungated` and `gated` are printed side by side with the number of points
each dropped, and both verdicts are reported.
"""
import argparse
import json
import os
import socket
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# The engine lives in tests/nbo_reference and is imported from there. This directory used to
# hold a byte-identical duplicate of it, and "make the engine symmetric" is not a change you can
# make in one of two copies and trust.
sys.path.insert(0, os.path.join(HERE, os.pardir, "nbo_reference"))
import compare_nbo  # noqa: E402  the shared engine
from config import load_nbo  # refuses a reference an older parser wrote

EMPTY = 1.0e-5  # the engine's own threshold for "this orbital holds nothing"


def gate_empty_rydberg(res):
    """A shallow copy of one side without its empty Rydberg NBOs, and without the E2 rows that
    accept into a Rydberg orbital.

    Comparing an empty orbital by label compares two different orbitals. Any unitary rotation
    inside the empty space is an equally valid NBO set, so its labels, its ordering and its
    one-electron energies are that code's arbitrary choice - ch3 showed gennbo's first five RY
    energies (3.50, 26.26, 1.35, 12.33, 1.99) against native's (2.95, 1.56, 5.68, 1.40, 2.21),
    a permutation-like set rather than a disagreement. The occupancies, which are physical,
    agreed to 9.8e-4 across all 36.
    """
    out = dict(res)
    out["nbos"] = [o for o in res.get("nbos", [])
                   if not (o["description"].split()[0].startswith("RY")
                           and abs(o["occupancy"]) < EMPTY)]
    # E2 labels spell the same orbital differently from the NBO list ("RY ( 1) H 2" against
    # "RY (3) H 4 s"), so the gate keys on the acceptor's type token rather than trying to
    # match a label: an acceptor whose occupancy is what makes the interaction meaningful.
    out["e2"] = [e for e in res.get("e2", [])
                 if not e["acceptor"].split()[0].startswith("RY")]
    return out


def view(reference, candidate):
    qs = compare_nbo.compare(reference, candidate)
    return {"status": "PASS" if all(q.ok() for q in qs) else "FAIL",
            "quantities": [{"name": q.name, "tol": q.tol, "compared": q.compared,
                            "missing": q.missing, "surplus": q.surplus, "failed": q.failed,
                            "max_dev": q.max_dev, "worst": q.worst}
                           for q in qs if q.compared or q.missing or q.surplus],
            "n_points": sum(q.compared for q in qs),
            "n_failed": sum(q.failed for q in qs),
            "n_missing": sum(q.missing for q in qs),
            "n_surplus": sum(q.surplus for q in qs)}


def one(root, mol):
    d = os.path.join(root, mol)
    ref = os.path.join(d, mol + ".gennbo.nbo.json")
    cand = os.path.join(d, mol + ".native.nbo.json")
    out = {"molecule": mol, "gennbo_json": ref, "native_json": cand}
    for p in (ref, cand):
        if not os.path.exists(p):
            out["status"] = "missing file: " + os.path.basename(p)
            return out
    reference = load_nbo(ref)
    candidate = load_nbo(cand)

    prov = os.path.join(d, "provenance_nbo.json")
    if os.path.exists(prov):
        with open(prov) as f:
            p = json.load(f)
        out["provenance"] = {k: p.get(k) for k in
                             ("threads", "node_cores_total", "slurm_job_id",
                              "nospher_a2_commit", "nbo_keywords", "native_flags",
                              "seconds", "exit_codes")}
        # env.sh writes "hostname", the probe writes "host"; a row with no host is a row nobody
        # can reuse, so take whichever is there rather than printing None.
        out["provenance"]["host"] = p.get("hostname") or p.get("host")

    out["ungated"] = view(reference, candidate)
    out["gated"] = view(gate_empty_rydberg(reference), gate_empty_rydberg(candidate))
    # Surplus entries - printed by native, absent from gennbo - now come from the engine itself
    # (Quantity.surplus). They used to be invisible: the engine walked the reference side only,
    # so native's 12 ch3 E2 rows against gennbo's empty table left no trace at all.
    out["surplus"] = {"total": out["ungated"]["n_surplus"],
                      "by_quantity": {q["name"]: q["surplus"]
                                      for q in out["ungated"]["quantities"] if q["surplus"]}}
    out["dropped"] = {
        "nbos_gennbo": len(reference.get("nbos", [])) - len(gate_empty_rydberg(reference)["nbos"]),
        "nbos_native": len(candidate.get("nbos", [])) - len(gate_empty_rydberg(candidate)["nbos"]),
        "e2_gennbo": len(reference.get("e2", [])) - len(gate_empty_rydberg(reference)["e2"]),
        "e2_native": len(candidate.get("e2", [])) - len(gate_empty_rydberg(candidate)["e2"]),
        "points": out["ungated"]["n_points"] - out["gated"]["n_points"]}
    # The ungated verdict is the molecule's verdict. The gate is an explanation of *which*
    # comparisons failed, not permission to stop counting them.
    out["status"] = out["ungated"]["status"]
    for k in ("quantities", "n_points", "n_failed", "n_missing", "n_surplus"):
        out[k] = out["ungated"][k]
    # open_shell and the printing thresholds are part of "the same question was asked":
    # a native run at a different E2 cut-off compares a shorter table, not a worse one.
    for k in ("open_shell", "thresholds"):
        if reference.get(k) != candidate.get(k):
            out.setdefault("notes", []).append(
                "%s differs: gennbo %r, native %r" % (k, reference.get(k), candidate.get(k)))
    return out


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--json")
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv[1:])

    mols = ([m for m in a.only.split(",") if m] or
            sorted(d for d in os.listdir(a.root)
                   if os.path.isdir(os.path.join(a.root, d))))
    rows = [one(a.root, m) for m in mols]

    print("compare_all on %s, python %s" % (socket.gethostname(), sys.version.split()[0]))
    print("gennbo 7 (reference side) vs -nbo_native (candidate side), same .47\n")
    print("two views: `ungated` is every printed orbital, `gated` drops empty Rydberg NBOs and")
    print("E2 rows accepting into one. Both are shown; the ungated verdict is the verdict.\n")
    hdr = ("%-14s | %-6s %6s %6s %6s %6s | %-6s %6s %6s %6s %6s | %6s  %s"
           % ("molecule", "ungat", "points", "fail", "miss", "surp",
              "gated", "points", "fail", "miss", "surp", "-pts", "native host/threads"))
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        pv = r.get("provenance") or {}
        where = "%s/%s" % (pv.get("host", "?"), pv.get("threads", "?"))
        u, g = r.get("ungated"), r.get("gated")
        if not u:
            print("%-14s | %s" % (r["molecule"], r["status"]))
            continue
        print("%-14s | %-6s %6d %6d %6d %6d | %-6s %6d %6d %6d %6d | %6d  %s"
              % (r["molecule"], u["status"], u["n_points"], u["n_failed"], u["n_missing"],
                 u["n_surplus"],
                 g["status"], g["n_points"], g["n_failed"], g["n_missing"], g["n_surplus"],
                 r["dropped"]["points"], where))
        dd = r["dropped"]
        print("    gate dropped %d/%d empty Rydberg NBOs (gennbo/native) and %d/%d Rydberg-acceptor"
              " E2 rows, %d comparisons in total"
              % (dd["nbos_gennbo"], dd["nbos_native"], dd["e2_gennbo"], dd["e2_native"],
                 dd["points"]))
        sp = r.get("surplus") or {}
        if sp.get("total"):
            print("    native printed %d entries gennbo did not: %s"
                  % (sp["total"], ", ".join("%s %d" % kv
                                            for kv in sorted(sp["by_quantity"].items()))))
        for n in r.get("notes", []):
            print("    note: " + n)
        gq = {q["name"]: q for q in g["quantities"]}
        for q in u["quantities"]:
            flag = "FAIL" if (q["failed"] or q["missing"] or q["surplus"]) else "ok  "
            h = gq.get(q["name"])
            gtxt = ("gated %4d pts %3d fail %3d surp  max dev %12.6f" %
                    (h["compared"], h["failed"], h["surplus"], h["max_dev"])) if h else "gated (gone)"
            print("    %-20s %s %5d pts %3d fail %3d miss %3d surp  max dev %12.6f | %s"
                  % (q["name"], flag, q["compared"], q["failed"], q["missing"], q["surplus"],
                     q["max_dev"], gtxt))
            print("        worst ungated: %s" % q["worst"])

    ok = [r for r in rows if r["status"] == "PASS"]
    bad = [r for r in rows if r["status"] == "FAIL"]
    absent = [r for r in rows if r["status"] not in ("PASS", "FAIL")]
    pts = sum(r.get("n_points", 0) for r in rows)
    gpts = sum((r.get("gated") or {}).get("n_points", 0) for r in rows)
    gok = [r for r in rows if (r.get("gated") or {}).get("status") == "PASS"]
    print("\n%d PASS, %d FAIL, %d without both JSONs; %d data points compared"
          % (len(ok), len(bad), len(absent), pts))
    print("gated view for the same molecules: %d PASS, %d FAIL, %d data points (%d dropped)"
          % (len(gok), len(rows) - len(gok) - len(absent), gpts, pts - gpts))
    if absent:
        print("no pair for: " + " ".join("%s (%s)" % (r["molecule"], r["status"])
                                         for r in absent))

    if a.json:
        with open(a.json, "w") as f:
            json.dump({"host": socket.gethostname(), "rows": rows,
                       "n_pass": len(ok), "n_fail": len(bad), "n_points": pts,
                       "gated": {"n_pass": len(gok), "n_points": gpts,
                                 "n_dropped": pts - gpts}}, f, indent=1)
        print("wrote " + a.json)
    return 0 if ok and not bad and not absent else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
