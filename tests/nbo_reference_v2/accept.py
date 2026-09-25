"""Accept or reject each regenerated wavefunction as a stand-in for the lost original.

    python accept.py <run_dir> [--json accepted.json]

<run_dir> holds one directory per molecule with the ORCA output in it.  For each molecule
this reads back from that output exactly what collect_reference.py read - basis dimension,
charge, multiplicity, HF type, the final single point energy, whether the optimisation
converged - and compares it against what tests/nbo_reference/index.json recorded for the
run whose wavefunction no longer exists.

A molecule is ACCEPTED only if all of:

  * the optimisation converged,
  * charge, multiplicity and HF type match,
  * `basis_functions` matches EXACTLY.  This is the identity test.  def2-TZVP on a different
    molecule, or on the same atoms in a different composition, gives a different count; the
    `water` reference of the old set (43 basis functions, OHH) against the `water` test
    wavefunction that still exists in the tree (73, OHHHe) would have been caught by it at
    once.  One function of difference is a different molecule, not a tolerance.
  * |final_energy - recorded| <= the per-molecule tolerance below.

TOLERANCE.  Both runs are ORCA 6.1.1, B3LYP/def2-TZVP, TightSCF, default Opt convergence.
What can differ is only where inside the same convergence basin the optimiser stopped.
ORCA's default geometry criteria are 5e-6 Eh on the energy change and 3e-4 Eh/a0 on the
maximum gradient; at a stationary point the energy is quadratic in the displacement, so a
residual gradient g with curvature k leaves E above the minimum by about g^2/2k.  For the
stiff bonds in this set (k of order 0.5 Eh/a0^2) that is under 1e-7 Eh, and the observed
spread is dominated instead by the DFT integration grid moving with the nuclei.  The default
here is therefore 2.0e-5 Eh, which is four times ORCA's own per-step energy criterion.

Two groups get a looser number, declared here rather than widened after seeing the result:

  * ni_co_4 and ticl4 (1.0e-4 Eh): a transition metal with soft M-L bending modes; the
    gradient tolerance buys much less energy resolution on a flat surface.
  * nitromethane (1.0e-4 Eh): an almost free methyl rotor, so the optimiser can stop at a
    different rotamer of an essentially flat torsion.

The deviation actually observed is printed for every molecule whether it passes or not, so
a tolerance that turns out to be too generous is visible instead of hidden.
"""
import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REF_DIR = os.path.join(HERE, os.pardir, "nbo_reference")

DEFAULT_TOL_HA = 2.0e-5
TOL_HA = {
    "ni_co_4": 1.0e-4,
    "ticl4": 1.0e-4,
    "nitromethane": 1.0e-4,
}

# the same regexes collect_reference.py uses, so both sides read the same fields
FIELDS = {
    "version":         (r"Program Version (\S+)", str),
    "keyword_line":    (r"^\|\s*1>\s*(!.*?)\s*$", str),
    "basis":           (r"Your calculation utilizes the basis:\s*(\S+)", str),
    "hf_type":         (r"Hartree-Fock type\s+HFTyp\s+\.+\s*(\S+)", str),
    "charge":          (r"Total Charge\s+Charge\s+\.+\s*(-?\d+)", int),
    "multiplicity":    (r"Multiplicity\s+Mult\s+\.+\s*(\d+)", int),
    "basis_functions": (r"Basis Dimension\s+Dim\s+\.+\s*(\d+)", int),
}


def read_orca(path):
    out = {}
    pats = {k: re.compile(p, re.M) for k, (p, _) in FIELDS.items()}
    energy, converged = None, False
    with open(path, errors="replace") as f:
        for line in f:
            for k, pat in pats.items():
                if k not in out:
                    m = pat.search(line)
                    if m:
                        out[k] = FIELDS[k][1](m.group(1))
            if line.startswith("FINAL SINGLE POINT ENERGY"):
                energy = float(line.split()[-1])
            elif "THE OPTIMIZATION HAS CONVERGED" in line:
                converged = True
    out["final_energy_hartree"] = energy
    out["optimisation_converged"] = converged
    return out


def judge(name, ref, got):
    """(accepted, tol, dev, list of reasons it failed)"""
    tol = TOL_HA.get(name, DEFAULT_TOL_HA)
    bad = []
    if not got.get("optimisation_converged"):
        bad.append("optimisation did not converge")
    for k in ("charge", "multiplicity", "hf_type", "basis"):
        if got.get(k) != ref.get(k):
            bad.append("%s %r != recorded %r" % (k, got.get(k), ref.get(k)))
    if got.get("basis_functions") != ref.get("basis_functions"):
        bad.append("basis_functions %r != recorded %r - DIFFERENT MOLECULE"
                   % (got.get("basis_functions"), ref.get("basis_functions")))
    if got.get("version") != ref.get("version"):
        bad.append("ORCA version %r != recorded %r" % (got.get("version"), ref.get("version")))
    dev = None
    if got.get("final_energy_hartree") is None:
        bad.append("no FINAL SINGLE POINT ENERGY in the output")
    else:
        dev = abs(got["final_energy_hartree"] - ref["final_energy_hartree"])
        if dev > tol:
            bad.append("energy deviation %.3e Ha > tolerance %.1e Ha" % (dev, tol))
    return (not bad), tol, dev, bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()

    idx = json.load(open(os.path.join(REF_DIR, "index.json")))
    rec = {m["name"]: m for m in idx["molecules"]}
    result, n_ok = {}, 0
    print("%-13s %5s %5s  %-18s %-18s %10s %9s  %s"
          % ("molecule", "nbf", "ref", "E regenerated", "E recorded", "dev/Ha", "tol/Ha", "verdict"))
    for name in sorted(rec):
        out = os.path.join(a.run_dir, name, name + ".out")
        if not os.path.exists(out):
            print("%-13s %s" % (name, "MISSING " + out))
            result[name] = {"accepted": False, "reasons": ["output missing"]}
            continue
        got = read_orca(out)
        ok, tol, dev, bad = judge(name, rec[name]["orca"], got)
        prov = os.path.join(a.run_dir, name, "provenance_orca.json")
        p = json.load(open(prov)) if os.path.exists(prov) else {}
        print("%-13s %5s %5d  %18.9f %18.9f %10s %9.1e  %s"
              % (name, got.get("basis_functions"), rec[name]["orca"]["basis_functions"],
                 got.get("final_energy_hartree") or float("nan"),
                 rec[name]["orca"]["final_energy_hartree"],
                 ("%.3e" % dev) if dev is not None else "-", tol,
                 "ACCEPT" if ok else "REJECT: " + "; ".join(bad)))
        n_ok += ok
        result[name] = {
            "accepted": ok, "tolerance_hartree": tol, "energy_deviation_hartree": dev,
            "reasons": bad, "orca": got,
            "recorded": rec[name]["orca"],
            "host": p.get("hostname"), "threads": p.get("threads"),
            "node_cores_total": p.get("node_cores_total"),
            "wall_seconds": p.get("wall_seconds"),
        }
    print("\n%d of %d accepted" % (n_ok, len(rec)))
    if a.json:
        json.dump(result, open(a.json, "w"), indent=1, sort_keys=True)
        print("wrote " + a.json)
    return 0 if n_ok == len(rec) else 1


if __name__ == "__main__":
    sys.exit(main())
