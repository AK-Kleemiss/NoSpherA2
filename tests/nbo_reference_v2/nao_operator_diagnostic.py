"""Check NBO 7 final NAO density blocks against the native final operators.

    py -3.12 nao_operator_diagnostic.py <directory containing ops_<molecule>>
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from fit144 import rediag
from fitw import MOLS, prepare, verdict
from spec_steps import al_groups, cascade, m_average


def offdiagonal(g, key):
    C = g[key]
    D = C.T @ g["SPS"] @ C
    same = cross = 0.0
    for (_, l), cols in al_groups(g["lbl"], range(g["n"])).items():
        nm = 2 * l + 1
        B = m_average(D, cols, nm)
        for s in range(len(B)):
            for t in range(s):
                if g["cls"][cols[s * nm]] == g["cls"][cols[t * nm]]:
                    same = max(same, abs(B[s, t]))
                else:
                    cross = max(cross, abs(B[s, t]))
    return same, cross


def main(root):
    for mol in MOLS:
        g = prepare(mol, os.path.join(root, "ops_" + mol))
        C = cascade(g["C32"], g["pre_occ"], g["cls"], g["S"], g["SPS"], g["lbl"],
                    renat=True, ryd_w="step5")
        results = {}
        for key in ("core", "class"):
            A = rediag(C, g["SPS"], g["lbl"], key)
            pop = np.diag(A.T @ g["SPS"] @ A)
            results[key] = (verdict(A, g["C33"], g)["worst"],
                            sum(pop[i] for i, c in enumerate(g["cls"]) if c == 2))
        p32 = offdiagonal(g, "C32")
        p33 = offdiagonal(g, "C33")
        print("%-8s PNAO same %.3e | NAO same %.3e cross %.3e | "
              "native core sine %.4f Ryd %.5f | class sine %.4f Ryd %.5f" %
              (mol, p32[0], p33[0], p33[1], *results["core"], *results["class"]))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
