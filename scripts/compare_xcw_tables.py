#!/usr/bin/env python3
"""Compare the XCW lambda-scan tables of two NoSpherA2 logs.

The table is what the fit actually produces, so it is the right thing to diff when
checking whether a faster I tensor moved the physics. Reports the worst absolute and
relative difference per column, and flags a row-count mismatch, which means the scan
itself diverged rather than merely shifting in the last digits.
"""
import io
import re
import sys

COLS = ["lambda", "criterion", "GooF", "total energy", "perturbation", "target", "A^2 (halt)"]


def table(path):
    """Scan rows are 6 columns, or 7 with -xcw_gaussian_halt (the A^2 halting statistic)."""
    rows = []
    for line in io.open(path, encoding="utf-8", errors="replace"):
        parts = line.split()
        if len(parts) not in (6, 7):
            continue
        try:
            vals = [float(p) for p in parts]
        except ValueError:
            continue
        # lambda in [0,1] and an energy of order -1e3 is the signature of the scan rows
        if 0.0 <= vals[0] <= 1.0 and abs(vals[3]) > 1.0:
            rows.append(vals)
    return rows


def halting_lambdas(path):
    """The halting recommendation is a decision, not a number - worth diffing separately."""
    out = []
    for line in io.open(path, encoding="utf-8", errors="replace"):
        m = re.search(r"halting lambda\* = ([0-9.]+) \(A\^2=([0-9.]+)\)", line)
        if m:
            out.append((float(m.group(1)), float(m.group(2))))
    return out


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("usage: compare_xcw_tables.py REFERENCE.log OTHER.log\n")
        return 2
    a, b = table(sys.argv[1]), table(sys.argv[2])
    print("rows: %d reference, %d other" % (len(a), len(b)))
    if len(a) != len(b):
        print("ROW COUNT DIFFERS - the scan diverged, not just the digits")
    n = min(len(a), len(b))
    if n == 0:
        print("no table rows parsed")
        return 1
    worst = 0.0
    ncol = min(len(a[0]), len(b[0]))
    for c in range(ncol):
        da = max(abs(a[i][c] - b[i][c]) for i in range(n))
        base = max(abs(a[i][c]) for i in range(n)) or 1.0
        rel = da / base
        worst = max(worst, rel)
        print("  %-14s max|abs| %.3e   max rel %.3e" % (COLS[c], da, rel))
    print("worst relative difference across the table: %.3e" % worst)

    ha, hb = halting_lambdas(sys.argv[1]), halting_lambdas(sys.argv[2])
    if ha or hb:
        print("halting recommendations: %d reference, %d other" % (len(ha), len(hb)))
        if ha == hb:
            print("  IDENTICAL - the halting decision is unchanged: %s" % (ha,))
        else:
            print("  DIFFER - fp32 changed the halting decision")
            print("    reference: %s" % (ha,))
            print("    other    : %s" % (hb,))
    return 0


if __name__ == "__main__":
    sys.exit(main())
