"""The check on read_47's open-shell guard: accept one packed triangle, refuse two.

    py -3.12 check_read47.py <a CLOSED-shell .47>

The .47 archives are not in git (derived, one command each - see data/README.md), so this takes
one as an argument. Any closed-shell member does; water is the cheapest.

The refusal case is built here by duplicating that archive's own $DENSITY payload, which is the
layout measured on real OPEN archives: ch3 1225 -> 2450 and o2 1953 -> 3906 values against
packed = n(n+1)/2, ratio exactly 2.0000, while water and benzene sit at 1.0000. read_47 used to
slice [:packed], so on those two it silently returned the ALPHA density wearing a total-density
label - and no layer downstream could catch it, because an open-shell .33's alpha block passes
read_lfn33's own C^T S C = 1 test on its own.
"""
import os
import re
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aonao_compare as ac  # noqa: E402

if len(sys.argv) != 2:
    raise SystemExit(__doc__)
src = sys.argv[1]

n, S, P = ac.read_47(src)
assert P.shape == (n, n), (n, P.shape)
ne = float((S @ P).trace())
assert ne > 0.5, ne  # a density this reader mangled would not trace to an electron count
print("closed shell accepted: n=%d, trace(S P) = %.5f electrons" % (n, ne))

t = open(src).read()
m = re.search(r"(\$DENSITY)(.*?)(\$END)", t, re.S)
assert m, "no $DENSITY section in " + src
doubled = t[:m.start()] + m.group(1) + m.group(2) + m.group(2) + m.group(3) + t[m.end():]
p2 = os.path.join(tempfile.mkdtemp(), "doubled_density.47")
open(p2, "w", newline="\n").write(doubled)

try:
    ac.read_47(p2)
except AssertionError as e:
    msg = str(e)
    assert "ratio 2.0000" in msg, msg          # it must report what it actually found
    assert "spin" in msg, msg                  # and say whose job the choice is
    print("doubled $DENSITY refused, with the reason in the message:")
    print("   ", msg.split(": ", 1)[1][:150])
else:
    raise SystemExit("FAIL: read_47 accepted two packed triangles and kept one of them silently")
print("OK")
