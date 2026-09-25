"""Native NoSpherA2 against gennbo 7, on the same wavefunction and the same .47 archive.

    python compare_all.py <root>  [--json summary.json] [--only mol,mol]

<root> holds one directory per molecule, each with the two JSONs stage 2 wrote:

    <mol>.gennbo.nbo.json   NoSpherA2 -nbo_parse of gennbo 7's own output
    <mol>.native.nbo.json   NoSpherA2 -nbo_native on the same .47

The two sides have independent provenance - one is the external program's output, the other
is our own analysis - which is the point. `compare_nbo.py --all .` run inside
tests/nbo_reference compares that dataset with itself and reports 22/22 at 0.00000; this
does not.

The comparison engine is tests/nbo_reference/compare_nbo.py, imported unchanged (a copy
sits next to this file on the cluster). gennbo is the reference side, so a quantity gennbo
printed and we did not is a `missing`, not a silent pass. NRT weights are read off the
by-rank comparison, never by label: starting the search from a different Lewis structure
renumbers the list without moving a weight.

Every row carries the hostname and thread count the native side ran with, read back from
provenance_nbo.json, because a number without its thread count and node is not a
measurement anyone can reuse.
"""
import argparse
import json
import os
import socket
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_nbo  # noqa: E402  the engine, unchanged


def one(root, mol):
    d = os.path.join(root, mol)
    ref = os.path.join(d, mol + ".gennbo.nbo.json")
    cand = os.path.join(d, mol + ".native.nbo.json")
    out = {"molecule": mol, "gennbo_json": ref, "native_json": cand}
    for p in (ref, cand):
        if not os.path.exists(p):
            out["status"] = "missing file: " + os.path.basename(p)
            return out
    with open(ref) as f:
        reference = json.load(f)
    with open(cand) as f:
        candidate = json.load(f)

    prov = os.path.join(d, "provenance_nbo.json")
    if os.path.exists(prov):
        with open(prov) as f:
            p = json.load(f)
        out["provenance"] = {k: p.get(k) for k in
                             ("host", "threads", "node_cores_total", "slurm_job_id",
                              "nospher_a2_commit", "nbo_keywords", "native_flags",
                              "seconds", "exit_codes")}

    qs = compare_nbo.compare(reference, candidate)
    out["status"] = "PASS" if all(q.ok() for q in qs) else "FAIL"
    out["quantities"] = [{"name": q.name, "tol": q.tol, "compared": q.compared,
                          "missing": q.missing, "failed": q.failed,
                          "max_dev": q.max_dev, "worst": q.worst}
                         for q in qs if q.compared or q.missing]
    out["n_points"] = sum(q.compared for q in qs)
    out["n_failed"] = sum(q.failed for q in qs)
    out["n_missing"] = sum(q.missing for q in qs)
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
    hdr = "%-14s %-6s %7s %7s %7s  %s" % ("molecule", "state", "points", "failed",
                                          "missing", "native host/threads")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        pv = r.get("provenance") or {}
        where = "%s/%s" % (pv.get("host", "?"), pv.get("threads", "?"))
        print("%-14s %-6s %7s %7s %7s  %s"
              % (r["molecule"], r["status"][:6], r.get("n_points", "-"),
                 r.get("n_failed", "-"), r.get("n_missing", "-"), where))
        for n in r.get("notes", []):
            print("    note: " + n)
        for q in r.get("quantities", []):
            flag = "FAIL" if (q["failed"] or q["missing"]) else "ok  "
            print("    %-20s %s %5d pts %3d fail %3d miss  max dev %12.6f  %s"
                  % (q["name"], flag, q["compared"], q["failed"], q["missing"],
                     q["max_dev"], q["worst"]))

    ok = [r for r in rows if r["status"] == "PASS"]
    bad = [r for r in rows if r["status"] == "FAIL"]
    absent = [r for r in rows if r["status"] not in ("PASS", "FAIL")]
    pts = sum(r.get("n_points", 0) for r in rows)
    print("\n%d PASS, %d FAIL, %d without both JSONs; %d data points compared"
          % (len(ok), len(bad), len(absent), pts))
    if absent:
        print("no pair for: " + " ".join("%s (%s)" % (r["molecule"], r["status"])
                                         for r in absent))

    if a.json:
        with open(a.json, "w") as f:
            json.dump({"host": socket.gethostname(), "rows": rows,
                       "n_pass": len(ok), "n_fail": len(bad), "n_points": pts}, f, indent=1)
        print("wrote " + a.json)
    return 0 if ok and not bad and not absent else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
