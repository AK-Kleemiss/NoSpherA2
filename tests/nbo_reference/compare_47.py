"""Block-by-block comparison of two NBO FILE.47 archives.

    python compare_47.py reference.47 candidate.47

Compares $OVERLAP, $DENSITY, $FOCK and $LCAOMO element by element and reports the
largest deviation per block. Both files must describe the same system in the same
basis and the same atom order; only the numbers are compared.

This is the check that caught the open-shell writer: an archive that is missing its
second spin block, or that writes the same block twice, has the right element count
nowhere but $OVERLAP, and an archive with a subtly wrong $DENSITY produces plausible
NBO output that is wrong in every downstream number. The reference for that check was
the FILE.47 ORCA writes itself for the same wavefunction (`! ... NBO`), which is
independent of anything in this repository.

Exit code 0 means every block matched within tolerance.
"""
import sys

# $OVERLAP and $DENSITY come from the same integrals and the same density matrix on both
# sides, so they have to agree to print precision. $FOCK does not: a writer that rebuilds
# it from orbital energies and coefficients carries the SCF's residual, which is the size
# of the convergence threshold and only moves the energies of empty orbitals.
TOL = {"OVERLAP": 1e-10, "DENSITY": 1e-8, "FOCK": 1e-3, "LCAOMO": 1e-7}


def block(path, name):
    vals, inside = [], False
    with open(path, errors="replace") as f:
        for line in f:
            s = line.strip()
            if s.startswith("$" + name):
                inside = True
                continue
            if inside:
                if s.startswith("$END"):
                    break
                vals += [float(x) for x in s.split()]
    return vals


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    ok = True
    for name, tol in TOL.items():
        a, b = block(argv[1], name), block(argv[2], name)
        if not a and not b:
            print("%-8s absent from both" % name)
            continue
        if len(a) != len(b):
            # The usual cause on an open-shell system: one side wrote one spin block.
            print("%-8s LENGTH MISMATCH: reference %d, candidate %d" % (name, len(a), len(b)))
            ok = False
            continue
        dev = max((abs(x - y) for x, y in zip(a, b)), default=0.0)
        print("%-8s %6d values  max |diff| %.3e  tol %.0e  %s"
              % (name, len(a), dev, tol, "ok" if dev <= tol else "FAIL"))
        ok &= dev <= tol
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
