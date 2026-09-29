"""Per-class subspace sine against gennbo, and the Rydberg population, from the C++ dumps.

    py -3.12 sines22.py <sides dir> <fid dir>

The gauge-free measurement the lane uses, taken on the binary's own final NAO matrix rather than
on a numpy replica: for each class, sin of the largest principal angle between the C++ class
subspace and gennbo's class subspace of the same size, read out of the arbiter's own .33.  The
Rydberg population column is sum_j c_j' S P S c_j over the Rydberg columns, against gennbo's
printed Rydberg occupancy sum, so the "population error" is a number in electrons.

Nothing here is a pass/fail: the acceptance gate is compare_all.py's, and it fails on all 22
before and after this change.  Each row carries the molecule's own instrument floor.
"""
import json
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aonao_compare import principal_sines               # noqa: E402
from step3_span import class_cols, load_sides, orthonormal  # noqa: E402

SIDES, FID = sys.argv[1], sys.argv[2]
MOLS = "ammonia benzene ethane lif pf5 sf6 so2 water".split()
CLS = {0: "Cor", 1: "Val", 2: "Ryd"}


def staged(mol, arm, tmp):
    d = os.path.join(tmp, arm, mol)
    os.makedirs(d, exist_ok=True)
    for ext in (".47", ".32", ".33", ".aonao.nbo.json"):
        shutil.copyfile(os.path.join(SIDES, mol, mol + ext), os.path.join(d, mol + ext))
    for ext in (".naoc.txt", ".naocpre.txt"):
        shutil.copyfile(os.path.join(FID, arm, mol, mol + ext), os.path.join(d, mol + ext))
    return d


def measure(x):
    C, S, Cg = x["Cnat"], x["S"], x["Cg"]
    out = {}
    for c in (0, 1, 2):
        a, b = class_cols(x["ncls"], c), class_cols(x["gcls"], c)
        if a and len(a) == len(b):
            out[CLS[c]] = float(principal_sines(orthonormal(C[:, a], S)[0],
                                               orthonormal(Cg[:, b], S)[0], S).max())
    ryd = class_cols(x["ncls"], 2)
    out["pop"] = float(np.einsum("ij,ij->j", C[:, ryd], x["SPS"] @ C[:, ryd]).sum())
    return out


tmp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_sinestage")
hdr = ("%-9s %-9s | %-23s | %-23s | %-23s | %s" %
       ("molecule", "floor", "Cor sine legacy->renat5", "Val sine legacy->renat5",
        "Ryd sine legacy->renat5", "Ryd pop e: gennbo / legacy / renat5"))
print(hdr)
print("-" * len(hdr))
n_better = {k: 0 for k in ("Cor", "Val", "Ryd", "pop")}
n_worse = dict(n_better)
n_flat = dict(n_better)
at_floor = 0
worst = {"legacy": 0.0, "renat5": 0.0}
wpop = {"legacy": 0.0, "renat5": 0.0}
for mol in MOLS:
    r, m = {}, {}
    for arm in ("legacy", "renat5"):
        r[arm] = load_sides(mol, staged(mol, arm, tmp), core_own_block=True)
        m[arm] = measure(r[arm])
    floor = r["legacy"]["floor"]
    # gennbo's own Rydberg population, read by FIELD NAME from its printed NAO table (5 decimals,
    # so this column's own print floor is 1e-5 per entry, not zero).
    gnao = json.load(open(os.path.join(SIDES, mol, mol + ".aonao.nbo.json")))["nao"]
    gpop = sum(e["occupancy"] for e in gnao if e["type"] == "Ryd")
    cells = []
    for k in ("Cor", "Val", "Ryd"):
        a, b = m["legacy"].get(k), m["renat5"].get(k)
        if a is None or b is None:
            cells.append("%-23s" % "n/a (class sizes differ)")
            continue
        d = b - a
        n_better[k] += d < -1e-12
        n_worse[k] += d > 1e-12
        n_flat[k] += abs(d) <= 1e-12
        if k == "Ryd":
            worst["legacy"] = max(worst["legacy"], a)
            worst["renat5"] = max(worst["renat5"], b)
            at_floor += b <= 10.0 * floor
        cells.append("%.4f -> %.4f" % (a, b))
    da = abs(m["legacy"]["pop"] - gpop)
    db = abs(m["renat5"]["pop"] - gpop)
    wpop["legacy"] = max(wpop["legacy"], da)
    wpop["renat5"] = max(wpop["renat5"], db)
    n_better["pop"] += db < da - 1e-12
    n_worse["pop"] += db > da + 1e-12
    n_flat["pop"] += abs(db - da) <= 1e-12
    print("%-9s %.2e | %-23s | %-23s | %-23s | %8.4f / %8.4f / %8.4f  |dpop| %.4f -> %.4f" % (
        mol, floor, cells[0], cells[1], cells[2], gpop, m["legacy"]["pop"], m["renat5"]["pop"],
        da, db))
print("-" * len(hdr))
for k in ("Cor", "Val", "Ryd", "pop"):
    print("%-4s closer %d / further %d / flat %d of %d" % (
        k, n_better[k], n_worse[k], n_flat[k], len(MOLS)))
print("worst Ryd sine over the 8: %.4f -> %.4f" % (worst["legacy"], worst["renat5"]))
print("worst |Ryd pop - gennbo| over the 8: %.4f e -> %.4f e" % (wpop["legacy"], wpop["renat5"]))
print("renat5 Ryd sines within 10x their own floor: %d of %d" % (at_floor, len(MOLS)))
print("The Val and Ryd sine columns come out EQUAL on all 8 and that is an identity, not a\ncoincidence: with the core exact, Val and Ryd are complementary subspaces of one fixed span,\nand complementary subspaces of the same space share their principal angles.  Read the two\ncolumns as ONE measurement, reported twice.")
