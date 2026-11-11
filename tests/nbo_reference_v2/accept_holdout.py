"""Accept or reject each HELD-OUT wavefunction, with the weaker criterion that is available.

    python accept_holdout.py <run_dir> [--json accepted_holdout.json] [--demo]

WHY NOT accept.py.  accept.py's gate is an identity test: basis_functions equal to
index.json's recorded value EXACTLY, energy within 2.0e-5 Eh of the recorded value.  That is
the right gate for the 22, because for them a recorded value exists.  None of the 13 molecules
in make_holdout.py was in the lost original run, so there is nothing to compare against and
accept.py - which iterates index.json - cannot judge them at all.  This is the same situation
allyl, hco and no2 are in, and data/accepted_radicals.json states it the same way.

WHAT IS CHECKED INSTEAD, all of it necessary and none of it an identity test:

  1. THE OPTIMIZATION HAS CONVERGED appears in the output.
  2. The last SCF of the run says SCF CONVERGED AFTER n CYCLES.
  3. Charge and multiplicity are the ones make_holdout.py asked for - read back out of ORCA's
     own echo, not out of the input file, so a mis-written input is visible.
  4. Basis and keyword line are the ones every other molecule in this tree used.
  5. The last <S**2> is within tolerance of S(S+1) for the multiplicity asked for: 0.03 for a
     doublet, which is the tolerance data/accepted_radicals.json already uses (allyl came in
     0.028 away), and the same 4 % of the exact value for a triplet, i.e. 0.08 of 2.00.  An
     open shell that converged to a spin-contaminated state is a different wavefunction from
     the one asked for, and the purpose of this set is to stress the SPIN path among others.
  6. basis_functions and the final energy are RECORDED, so the next regeneration of this set
     does have the exact identity test this run could not have.

WHAT THIS CANNOT DO, stated because the difference matters: it cannot prove the wavefunction
reproduces somebody else's earlier wavefunction.  It can only prove the wavefunction is the
one this tree asked for.  Every comparison the set then feeds is still two-sided - ORCA plus
gennbo 7 on one side, -nbo_native on the other, both reading the same .47 - which is the
property arbitration needs.

The electron count and any ECP line are recorded for every molecule.  For hi that is the
measurement rather than a detail: def2-TZVP on iodine is an ECP basis, 28 core electrons are
not in the wavefunction at all, and whether the rest of the pipeline survives that is the
question hi was added to answer.
"""
import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import accept                       # read_orca and its FIELDS regexes, unchanged
import make_holdout as mh
import make_inputs as mi

BASIS = "def2-TZVP"
# S(S+1) for the multiplicity, and the tolerance: 4 % of the exact value, which IS the 0.03 of
# 0.75 that data/accepted_radicals.json already accepted allyl at.
S2_TOL_FRACTION = 0.04

S2_LAST = re.compile(r"Expectation value of <S\*\*2>\s*:\s*(-?\d+\.\d+)")
NEL_LAST = re.compile(r"Number of Electrons\s+NEL\s+\.+\s*(\d+)")
SCF_OK = re.compile(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES")
ECP_LINE = re.compile(r"^\s*(.*ECP.*)$")


def s2_exact(mult):
    s = (mult - 1) / 2.0
    return s * (s + 1.0)


def read_extra(path):
    """The fields accept.read_orca does not read, each taken from the LAST line that carries
    its own label - never from a position or from the nearest preceding header."""
    s2, nel, scf, ecp = None, None, None, []
    with open(path, errors="replace") as f:
        for line in f:
            m = S2_LAST.search(line)
            if m:
                s2 = float(m.group(1))
            m = NEL_LAST.search(line)
            if m:
                nel = int(m.group(1))
            m = SCF_OK.search(line)
            if m:
                scf = int(m.group(1))
            if "ECP" in line and len(ecp) < 6:
                t = line.strip()
                if t and t not in ecp:
                    ecp.append(t)
    return {"s_squared": s2, "electrons": nel, "scf_cycles_last": scf, "ecp_lines": ecp}


def judge(name, got, extra):
    q, mult = mh.CHARGE_MULT[name]
    bad = []
    if not got.get("optimisation_converged"):
        bad.append("optimisation did not converge")
    if extra["scf_cycles_last"] is None:
        bad.append("no SCF CONVERGED AFTER line - the last SCF did not converge")
    if got.get("charge") != q:
        bad.append("charge %r != asked %r" % (got.get("charge"), q))
    if got.get("multiplicity") != mult:
        bad.append("multiplicity %r != asked %r" % (got.get("multiplicity"), mult))
    if got.get("basis") != BASIS:
        bad.append("basis %r != %s" % (got.get("basis"), BASIS))
    if got.get("keyword_line") != mi.KEYWORD_LINE:
        bad.append("keyword line %r != %r" % (got.get("keyword_line"), mi.KEYWORD_LINE))
    if got.get("final_energy_hartree") is None:
        bad.append("no FINAL SINGLE POINT ENERGY in the output")
    if got.get("basis_functions") is None:
        bad.append("no Basis Dimension in the output")
    dev2 = None
    if mult > 1:
        want = s2_exact(mult)
        tol = S2_TOL_FRACTION * want
        if extra["s_squared"] is None:
            bad.append("multiplicity %d but no <S**2> printed" % mult)
        else:
            dev2 = abs(extra["s_squared"] - want)
            if dev2 > tol:
                bad.append("<S**2> %.6f is %.4f from %.2f, tolerance %.4f"
                           % (extra["s_squared"], dev2, want, tol))
    return (not bad), dev2, bad


def demo():
    """The gate must reject the failures it exists to catch, on synthetic inputs only."""
    assert abs(s2_exact(2) - 0.75) < 1e-12
    assert abs(s2_exact(3) - 2.00) < 1e-12
    assert abs(s2_exact(1)) < 1e-12
    good = {"optimisation_converged": True, "charge": 0, "multiplicity": 2,
            "basis": BASIS, "keyword_line": mi.KEYWORD_LINE,
            "final_energy_hartree": -1.0, "basis_functions": 42}
    ex = {"s_squared": 0.7534, "electrons": 9, "scf_cycles_last": 6, "ecp_lines": []}
    ok, dev, bad = judge("hs", good, ex)
    assert ok, bad
    assert abs(dev - 0.0034) < 1e-9, dev
    ok, _, bad = judge("hs", dict(good, optimisation_converged=False), ex)
    assert not ok and "did not converge" in bad[0], bad
    ok, _, bad = judge("hs", good, dict(ex, scf_cycles_last=None))
    assert not ok and "SCF" in bad[0], bad
    ok, _, bad = judge("hs", good, dict(ex, s_squared=1.10))   # spin contamination
    assert not ok and "<S**2>" in bad[0], bad
    ok, _, bad = judge("hs", dict(good, multiplicity=4), ex)   # wrong state
    assert not ok, bad
    ok, _, bad = judge("hs", dict(good, basis="def2-SVP"), ex)
    assert not ok, bad
    # a closed shell must not be failed for having no <S**2>, and a triplet must be
    ok, _, bad = judge("cl2", dict(good, multiplicity=1), dict(ex, s_squared=None))
    assert ok, bad
    ok, _, bad = judge("ch2", dict(good, multiplicity=3), dict(ex, s_squared=None))
    assert not ok, bad
    ok, _, bad = judge("ch2", dict(good, multiplicity=3), dict(ex, s_squared=2.05))
    assert ok, bad                                             # 0.05 < 0.08
    print("demo OK, %d molecules in the held-out set" % len(mh.CHARGE_MULT))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", nargs="?")
    ap.add_argument("--json", default=None)
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    if a.demo:
        demo()
        return
    if not a.run_dir:
        ap.error("run_dir is required unless --demo")
    demo()
    result, n_ok, n_seen = {}, 0, 0
    print("%-10s %4s %4s %5s %-20s %9s %9s  %s"
          % ("molecule", "nbf", "nel", "mult", "E / Eh", "<S**2>", "dev", "verdict"))
    for name in sorted(mh.CHARGE_MULT):
        out = os.path.join(a.run_dir, name, name + ".out")
        if not os.path.exists(out):
            print("%-10s MISSING %s" % (name, out))
            result[name] = {"accepted": False, "reasons": ["output missing: " + out]}
            continue
        n_seen += 1
        got = accept.read_orca(out)
        extra = read_extra(out)
        ok, dev2, bad = judge(name, got, extra)
        n_ok += ok
        result[name] = {
            "accepted": bool(ok),
            "basis_functions": got.get("basis_functions"),
            "final_single_point_energy": got.get("final_energy_hartree"),
            "charge": got.get("charge"),
            "multiplicity": got.get("multiplicity"),
            "electrons": extra["electrons"],
            "s_squared": extra["s_squared"],
            "s_squared_exact": s2_exact(got.get("multiplicity") or 1),
            "optimisation_converged": bool(got.get("optimisation_converged")),
            "scf_cycles_last": extra["scf_cycles_last"],
            "orca_version": got.get("version"),
            "ecp_lines": extra["ecp_lines"],
            "reasons": bad or ["converged, charge/multiplicity/basis as asked"],
        }
        print("%-10s %4s %4s %5s %-20s %9s %9s  %s"
              % (name, got.get("basis_functions"), extra["electrons"],
                 got.get("multiplicity"),
                 "%.8f" % got["final_energy_hartree"] if got.get("final_energy_hartree") else "-",
                 "%.6f" % extra["s_squared"] if extra["s_squared"] is not None else "-",
                 "%.4f" % dev2 if dev2 is not None else "-",
                 "ACCEPT" if ok else "REJECT: " + "; ".join(bad)))
    print("\n%d of %d present accepted (%d of %d in the set)"
          % (n_ok, n_seen, n_ok, len(mh.CHARGE_MULT)))
    if a.json:
        doc = {
            "_note": [
                "The held-out set of make_holdout.py. None of these molecules is in",
                "tests/nbo_reference/index.json, so accept.py's exact basis-function identity",
                "test is not available and this weaker criterion is what was checked instead.",
                "basis_functions and the energy are RECORDED here so the next regeneration has",
                "the exact test.",
            ],
            "_criterion": ("converged optimisation + converged last SCF + charge/multiplicity/"
                           "basis/keyword line as asked + <S**2> within 4 %% of S(S+1); "
                           "basis functions recorded for next time"),
            "_keyword_line": mi.KEYWORD_LINE,
        }
        doc.update(result)
        with open(a.json, "w") as f:
            json.dump(doc, f, indent=1, sort_keys=False)
        print("wrote " + a.json)
    sys.exit(0 if n_ok == len(mh.CHARGE_MULT) else 1)


if __name__ == "__main__":
    main()
