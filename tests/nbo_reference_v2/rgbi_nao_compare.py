"""Compare RGBI's molecular-local fallback operator with NBO 7 pre-NAOs."""

import os
import sys

import numpy as np

from dmpnao_step4 import read_47
from fitw import MOLS, prepare


def compare(mol, directory):
    g = prepare(mol, directory)
    n, S, C, lbl = g["n"], g["S"], g["C32"], g["lbl"]
    _, _, P = read_47(os.path.join(directory, mol + ".47"))
    SPS = g["SPS"]
    local = np.zeros((n, n))
    for atom in sorted({t[1] for t in lbl}):
        cols = [i for i, t in enumerate(lbl) if t[1] == atom]
        rows = np.flatnonzero(np.max(np.abs(C[:, cols]), axis=1) > 1e-8)
        assert len(rows) == len(cols), (mol, atom, len(rows), len(cols))
        assert np.max(np.abs(C[np.ix_(np.setdiff1d(np.arange(n), rows), cols)])) < 1e-8
        ca = C[np.ix_(rows, cols)]
        sa = S[np.ix_(rows, rows)]
        pa = P[np.ix_(rows, rows)]
        local[np.ix_(cols, cols)] = ca.T @ sa @ pa @ sa @ ca
    global_density = C.T @ SPS @ C
    delta = np.diag(local - global_density)
    off = []
    for atom, l in sorted({t[1:3] for t in lbl}):
        shells = sorted({t[4] for t in lbl if t[1:3] == (atom, l)})
        for s1 in shells:
            for s2 in shells:
                if s1 >= s2:
                    continue
                ij = [(i, j) for i, a in enumerate(lbl) for j, b in enumerate(lbl)
                      if a[1:3] == b[1:3] == (atom, l) and a[4] == s1 and b[4] == s2
                      and a[3] == b[3]]
                off.append(abs(np.mean([local[i, j] for i, j in ij])))
    return np.max(np.abs(delta)), np.linalg.norm(delta), max(off, default=0.0)


if __name__ == "__main__":
    root = sys.argv[1]
    for name in MOLS:
        result = compare(name, os.path.join(root, "ops_" + name))
        print("%-9s local/global pre-occ max %.6f norm %.6f; local radial offdiag %.6f" %
              (name, *result))
