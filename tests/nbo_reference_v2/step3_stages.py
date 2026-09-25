"""Where the Rydberg population comes from: pre-NAO, after step 3, final, gennbo 7.

Step 4 was exonerated by measurement (NAO_CLASS_SPLIT's arm is 5.6x worse and the remaining leak
sits in l blocks holding no valence shell, which a unitary intra-block rotation cannot feed).  That
leaves step 3 - the Schmidt projection of each class out of the ones above it, then the OWSO within
the class - as the only stage that can move population between l values.  This script measures it
instead of arguing about it, over four stages of the same run:

  pre_occ      the pre-NAO occupancies of step 1: one generalised eigenproblem per (atom, l),
               nothing inter-atomic has happened yet - and, as the run below showed, a gross
               population that does not sum to N, so it is kept only as the retired metric it is
  step3        the same shells after the projection + OWSO, m-averaged, as step 4 inherits them
  final        the NAOs NoSpherA2 reports
  gennbo       NBO 7's own NAO table for the same wavefunction

PRE-REGISTERED, written before the run (and deliberately not as an exhaustive list of causes - the
last two-way prediction in this lane was answered by a third branch):

  The metric is the molecule's total Rydberg-class population at each stage, against gennbo's
  final Rydberg total, plus the inter-l split the previous arm built: population in (atom, l)
  blocks that contain no valence shell.  The number to move is +1.15335 e.

  Expectation: pre_occ's Rydberg total is CLOSE to gennbo's - the pre-NAO stage is shared between
  the two codes - and step3's is much larger.  That would put the leak in the projection and the
  renormalisation that follows it, and would also explain why step 4 recovers most but not all of
  it: it can only recover what sits in a block that also holds a valence shell.

  What would refute it: a pre_occ Rydberg total already several times gennbo's.  Step 1 is
  per-(atom, l) and cannot move population between l values either, so that outcome would mean the
  excess is not created by any inter-l mixing at all and the whole localisation argument is wrong.

OUTCOME, job 586429, all 25 molecules: the pre-NAO Rydberg total is 235.54142 e over 22 molecules
against gennbo's 2.74986 - 86x, not "several times".  By the letter of the pre-registration that
refutes the hypothesis.  It does not, and the reason is the THIRD retired metric in this lane:

  pre_occ is the eigenvalue w of (S^A)^-1/2 (SPS)^A (S^A)^-1/2.  That is a gross population of
  functions that are not orthogonal between atoms, and the `pre_all` column below shows its total
  over all classes is far above the molecule's electron count.  A number whose total is not N
  cannot be compared with a post-orthogonalisation population that sums to N: the excess is not
  Rydberg character, it is the double counting that orthogonalisation removes.

  So the pre-NAO comparison is retired, next to the NPA charge metric, and for the same kind of
  reason: it answered, and the answer meant nothing.  Stage 3's own column is unaffected - after
  step 3 the set IS S-orthonormal, and the class totals do sum to N (checked below), which is why
  the step-3 and final columns can be compared with gennbo at all.

  What it costs: there is now no arbitrated intermediate.  NBO 7 prints no pre-NAO occupancy
  table, so no stage before its final NAOs can be compared with ours at all, and a change to step 3
  can only be judged by the gate below.

Acceptance gate for ANY later change to step 3, fixed here so it cannot be relaxed afterwards:
pf5, so2 and sf6 must keep a final Rydberg total within 0.90-1.05 of gennbo's (they are at
0.93, 0.94 and 0.98 now and they are the only three whose leak has the opposite sign); no molecule
may get worse by more than 0.01 e; and the no-valence-block excess must fall.  A change that
improves the mean while flattening those three has found a fudge factor.

    python3 step3_stages.py <root> [mol ...]

<root> holds one directory per molecule with NoSpherA2.log from a NAO_DUMP_STEP3=1 run and both
JSONs.  Prints one row per molecule and the totals.
"""
import json
import os
import socket
import sys
from collections import Counter

from config import load_nbo  # refuses a reference an older parser wrote
from ryd_excess_by_l import by_block  # the (atom, l) keying the previous arm used

CLASSES = {0: "Cor", 1: "Val", 2: "Ryd"}


def step3_rows(log_path):
    """(atom, l, shell, class, step3_occ, pre_occ) per shell, occupancies as whole-shell sums."""
    out = []
    for line in open(log_path, errors="replace"):
        if not line.startswith("STEP3 ") or line.split()[1] == "atom":
            continue
        f = line.split()
        if len(f) < 6:
            continue
        atom, l, sh, cls = int(f[1]), int(f[2]), int(f[3]), int(f[4])
        nm = 2 * l + 1
        occ = float(f[5]) * nm            #the dump divides by nm; put the components back
        pre = float(f[6]) * nm if len(f) > 6 else float("nan")
        out.append((atom, l, sh, CLASSES[cls], occ, pre))
    return out


def totals(rows):
    s3, pre = Counter(), Counter()
    for _, _, _, cls, occ, p in rows:
        s3[cls] += occ
        pre[cls] += p
    return s3, pre


def json_totals(path):
    t = Counter()
    for r in load_nbo(path)["nao"]:
        t[r["type"]] += r["occupancy"]
    return t


def no_valence_excess(rows, gennbo_nao):
    """Rydberg population in (atom, l) blocks that hold no valence shell, ours minus gennbo's.

    Ours is keyed by (atom, l) directly; gennbo's rows carry the printed lang letter instead, so
    they go through ryd_excess_by_l.by_block - the same mapping the previous arm measured with.
    """
    mine = {}
    for atom, l, _sh, cls, occ, _pre in rows:
        mine.setdefault((atom, l), []).append((cls, occ))
    theirs = {k: [(r["type"], r["occupancy"]) for r in v]
              for k, v in by_block(gennbo_nao).items()}
    tot = 0.0
    for key, block in mine.items():
        if any(c == "Val" for c, _ in block):
            continue
        tot += sum(o for c, o in block if c == "Ryd")
        tot -= sum(o for c, o in theirs.get(key, []) if c == "Ryd")
    return tot


def demo():
    """A p shell's three components must be summed back, and the pre column must be optional."""
    import tempfile
    p = os.path.join(tempfile.mkdtemp(), "NoSpherA2.log")
    open(p, "w").write("STEP3 atom l shell class occ pre\n"
                       "STEP3 1 1 0 2 0.10000000 0.01000000\n"
                       "STEP3 1 0 0 1 0.50000000 0.60000000\n")
    rows = step3_rows(p)
    s3, pre = totals(rows)
    assert abs(s3["Ryd"] - 0.3) < 1e-12 and abs(pre["Ryd"] - 0.03) < 1e-12, (s3, pre)
    assert abs(s3["Val"] - 0.5) < 1e-12 and abs(pre["Val"] - 0.6) < 1e-12, (s3, pre)
    print("demo OK on %s" % socket.gethostname())


def main(argv):
    if argv[:1] == ["--demo"]:
        demo()
        return 0
    root = argv[0]
    mols = argv[1:] or sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d)))
    print("step3_stages on %s, 1 thread, root %s" % (socket.gethostname(), root))
    print("%-13s %9s %9s %9s %9s %8s %9s %9s %9s" %
          ("molecule", "pre_Ryd", "s3_Ryd", "fin_Ryd", "gen_Ryd", "fin/gen", "noVal_d",
           "pre_all/N", "s3_all-N"))
    sums = Counter()
    n = 0
    for mol in mols:
        d = os.path.join(root, mol)
        log = os.path.join(d, "NoSpherA2.log")
        gp = os.path.join(d, "%s.gennbo.nbo.json" % mol)
        np_ = os.path.join(d, "%s.native.nbo.json" % mol)
        if not (os.path.exists(log) and os.path.exists(gp) and os.path.exists(np_)):
            continue
        rows = step3_rows(log)
        if not rows:
            print("%-13s no STEP3 lines - was NAO_DUMP_STEP3=1 set?" % mol)
            continue
        s3, pre = totals(rows)
        fin = json_totals(np_)
        gen = json_totals(gp)
        nv = no_valence_excess(rows, load_nbo(gp)["nao"])
        ratio = fin["Ryd"] / gen["Ryd"] if gen["Ryd"] else float("nan")
        n_el = sum(gen.values())
        print("%-13s %9.5f %9.5f %9.5f %9.5f %8.2f %+9.5f %9.2f %+9.2e" %
              (mol, pre["Ryd"], s3["Ryd"], fin["Ryd"], gen["Ryd"], ratio, nv,
               sum(pre.values()) / n_el, sum(s3.values()) - n_el))
        for k, v in (("pre", pre["Ryd"]), ("s3", s3["Ryd"]), ("fin", fin["Ryd"]),
                     ("gen", gen["Ryd"]), ("nv", nv),
                     ("pre_all", sum(pre.values())), ("s3_all", sum(s3.values())),
                     ("n", sum(gen.values()))):
            sums[k] += v
        n += 1
    if not n:
        print("no molecule had all three files")
        return 1
    print("\n%d molecules.  Rydberg totals summed: pre-NAO %.5f, after step 3 %.5f, final %.5f, "
          "gennbo %.5f" % (n, sums["pre"], sums["s3"], sums["fin"], sums["gen"]))
    print("Excess over gennbo: pre-NAO %+.5f, after step 3 %+.5f, final %+.5f" %
          (sums["pre"] - sums["gen"], sums["s3"] - sums["gen"], sums["fin"] - sums["gen"]))
    print("Rydberg population in l blocks with no valence shell, ours minus gennbo's, after "
          "step 3: %+.5f e" % sums["nv"])
    print("Trace check: summed over all classes, pre-NAO %.5f, after step 3 %.5f, N = %.5f" %
          (sums["pre_all"], sums["s3_all"], sums["n"]))
    print("pre-NAO's total is %.1fx N, so it is a gross non-orthogonal population and NOT a "
          "yardstick - see the OUTCOME paragraph above.  Step 3 sums to N and is." %
          (sums["pre_all"] / sums["n"]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
