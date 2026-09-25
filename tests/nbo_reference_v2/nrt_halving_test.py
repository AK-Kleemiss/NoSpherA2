"""Would halving the open-shell NRT numbers make them right?

nrt_ratio_stats.py showed a wide spread of native/gennbo ratios, but a raw ratio on a small number
is mostly a statement about gennbo's printed precision: its NRT tables carry three decimals, so a
printed 0.051 is anything in [0.0505, 0.0515) and the ratio against native's 0.00001 says little.
ratio_probe.py already settled how to ask this honestly - test the RESIDUAL native/2 - gennbo
against half of gennbo's last printed digit plus a relative term - and this applies that same
criterion to every spin-resolved NRT number, not to the subset one probe matched.

    python nrt_halving_test.py <root> [<root> ...]

Reported per molecule: how many numbers halving repairs, how many it does not, and the worst
offenders with their residuals.  "16 of 16 ch3 pairs at exactly 2.0000" is only the whole story
if the fails column is zero everywhere.
"""
import os
import socket
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nrt_ratio_stats import pairs   # the pairing, unchanged, so both tests see the same numbers

RESOLUTION = 0.001   # gennbo's last printed digit in the NRT tables
TOL_FRAC = 0.01      # ratio_probe.py's relative term


def verdict(g, n, expect=2.0):
    """True when dividing native by `expect` reproduces gennbo inside its printed precision."""
    return abs(n / expect - g) <= RESOLUTION / 2.0 + TOL_FRAC * abs(g)


def demo():
    """Halving must repair an exactly-doubled number and must not repair a tripled one."""
    assert verdict(0.404, 0.808)
    assert not verdict(0.404, 1.212)
    assert verdict(0.500, 1.0004)        # inside the printed-precision band
    print("demo OK on %s" % socket.gethostname())


if __name__ == "__main__":
    if sys.argv[1:2] == ["--demo"]:
        demo()
        raise SystemExit(0)
    print("nrt_halving_test on %s, 1 thread, residual |native/2 - gennbo| <= %.4f + %.2f|gennbo|"
          % (socket.gethostname(), RESOLUTION / 2.0, TOL_FRAC))
    print("%-12s %5s %6s %6s   %s" % ("molecule", "n", "repair", "fail", "worst three residuals"))
    for root in sys.argv[1:]:
        for mol in sorted(os.listdir(root)):
            d = os.path.join(root, mol)
            if not (os.path.isfile(os.path.join(d, "%s.native.nbo.json" % mol))
                    and os.path.isfile(os.path.join(d, "%s.gennbo.nbo.json" % mol))):
                continue
            rows = list(pairs(d, mol))
            if not rows:
                continue
            bad = [(abs(n / 2.0 - g), lbl, g, n) for lbl, g, n in rows if not verdict(g, n)]
            bad.sort(reverse=True)
            worst = "; ".join("%s g=%.4f n=%.4f r=%+.4f" % (l, g, n, n / 2.0 - g)
                              for _, l, g, n in bad[:3])
            print("%-12s %5d %6d %6d   %s" % (mol, len(rows), len(rows) - len(bad), len(bad), worst))
    print("\nEvery failing row is a number that halving leaves wrong: a second error under the two.")
