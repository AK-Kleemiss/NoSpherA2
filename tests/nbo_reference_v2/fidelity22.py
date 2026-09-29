"""Port fidelity: the C++ NAO matrix against spec_steps.py's own cascade, arm by arm.

    py -3.12 fidelity22.py <sides dir> <fid dir>

<sides dir>/<mol>/ holds the .47, .32, .33 and .aonao.nbo.json the numpy side reads.
<fid dir>/{renat5,legacy}/<mol>/ holds the .naoc.txt/.naocpre.txt the NEW binary wrote, default
and with NAO_LEGACY_CASCADE=1.

Four comparisons per molecule, all with the SAME instrument (step4_core_block.fidelity: worst
column angle, and worst degeneracy-projected subspace angle, both in radians):

    legacy C++  vs numpy base      the port must not have moved the legacy path
    renat5 C++  vs numpy renat5    the claim - the port implements the arm that was measured
    renat5 C++  vs numpy base      DISCRIMINATION, must be large or the instrument is blind
    legacy C++  vs renat5 C++      the two dumps differ at all

Each row carries the molecule's own floor (sqrt(2 max(|C'SC-1|, e33)) from load_sides): a
deviation below the floor is the floor, not an agreement.  The subspace cell prints the number of
degenerate columns it was taken over, because fidelity() returns 0.0 for an EMPTY maximum and a
0.00e+00 with 0 degenerate columns is a sample size of zero, not a perfect agreement.
"""
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from step3_span import load_sides                     # noqa: E402
from step4_core_block import fidelity, step4          # noqa: E402
from spec_steps import cascade, ARMS                  # noqa: E402

SIDES, FID = sys.argv[1], sys.argv[2]
MOLS = "ammonia benzene ethane lif pf5 sf6 so2 water".split()
KW = dict(ARMS)


def staged(mol, arm, tmp):
    """A sides directory for `mol` whose two dumps come from `arm` and nothing else changes."""
    d = os.path.join(tmp, arm, mol)
    os.makedirs(d, exist_ok=True)
    for ext in (".47", ".32", ".33", ".aonao.nbo.json"):
        shutil.copyfile(os.path.join(SIDES, mol, mol + ext), os.path.join(d, mol + ext))
    for ext in (".naoc.txt", ".naocpre.txt"):
        shutil.copyfile(os.path.join(FID, arm, mol, mol + ext), os.path.join(d, mol + ext))
    return d


tmp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_fidstage")
COLS = ("legacy C++ vs numpy base", "renat5 C++ vs numpy renat5",
        "renat5 C++ vs numpy base", "legacy C++ vs renat5 C++")
hdr = "%-9s %-9s | %-31s | %-31s | %-31s | %s" % (("molecule", "floor") + COLS)
print(hdr)
print("-" * len(hdr))
worst = {k: 0.0 for k in ("lb", "rr", "rb", "lr")}
wcol = {k: 0.0 for k in ("lb", "rr", "rb", "lr")}
for mol in MOLS:
    r = {}
    for arm in ("legacy", "renat5"):
        r[arm] = load_sides(mol, staged(mol, arm, tmp), core_own_block=True)
    x = r["legacy"]
    rep = {}
    for name in ("base", "renat5"):
        C = cascade(x["Cnpre"], x["pre_occ"], x["ncls"], x["S"], x["SPS"], x["lbl"], **KW[name])
        rep[name] = step4(C, x["SPS"], x["lbl"], True)
    cells = []
    for tag, cnat, crep in (("lb", r["legacy"]["Cnat"], rep["base"]),
                            ("rr", r["renat5"]["Cnat"], rep["renat5"]),
                            ("rb", r["renat5"]["Cnat"], rep["base"]),
                            ("lr", r["legacy"]["Cnat"], r["renat5"]["Cnat"])):
        wc, ws, nd, ns = fidelity(cnat, crep, x["S"], x["SPS"], x["lbl"])
        worst[tag] = max(worst[tag], ws)
        wcol[tag] = max(wcol[tag], wc)
        cells.append("col %.2e sub %.2e %2d/%-3d" % (wc, ws, nd, ns))
    print("%-9s %.2e | %-31s | %-31s | %-31s | %s" % ((mol, x["floor"]) + tuple(cells)))
print("-" * len(hdr))
print("sub cell = worst subspace angle, then degenerate columns / shells in the block set.")
print("0.00e+00 over 0 degenerate columns is an EMPTY maximum, not an agreement.")
for tag, name in zip(("lb", "rr", "rb", "lr"), COLS):
    print("%-27s worst over the 8: column %.3e   subspace %.3e" % (name, wcol[tag], worst[tag]))
