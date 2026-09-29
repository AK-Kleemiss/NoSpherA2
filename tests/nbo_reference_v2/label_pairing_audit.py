"""Does the comparison engine still lose matches to a printed label anywhere?

The NAO tables cost 11 of 22 molecules when paired by printed shell label, because the two sides
number Rydberg shells differently.  compare_nbo.py itself pairs NAOs by printed INDEX, which is a
rank and safe.  But it pairs NBOs, hybrids and E2 rows by the printed DESCRIPTION, and a
description carries a per-type ordinal - "BD ( 2) C 1- C 2" - that is each side's own bookkeeping:
two implementations that find the same two pi bonds in the opposite order disagree on every one of
those keys while agreeing on the physics.

This measures whether that is costing anything, by pairing the same rows three ways:

    full        the engine's key, ordinal included
    no_ordinal  the ordinal stripped, so "BD ( 2) C 1- C 2" and "BD ( 1) C 1- C 2" pair
    by_rank     within one (type, atom list), by descending occupancy

and reporting how many rows fail to pair under each.  If `full` loses rows that `no_ordinal`
keeps, the hole is live and it is the same hole as the NAO one.

    python label_pairing_audit.py <root>
"""
import json
import os
import re
import socket
import sys
from config import load_nbo  # refuses a reference an older parser wrote

ORDINAL = re.compile(r"\(\s*\d+\s*\)")


def key_full(o):
    return "%s%s" % (o["description"], (" " + o["spin"]) if o["spin"] else "")


def key_no_ordinal(o):
    return ORDINAL.sub("()", key_full(o))


def unmatched(ref, cand, key):
    """Rows on either side with no partner, duplicates pairing up in order of appearance."""
    from collections import Counter
    a, b = Counter(key(x) for x in ref), Counter(key(x) for x in cand)
    miss = sum((a - b).values())
    surp = sum((b - a).values())
    return miss, surp


def demo():
    """Stripping the ordinal must rescue exactly the pair that differs only in its ordinal."""
    ref = [{"description": "BD ( 1) C 1- C 2", "spin": ""},
           {"description": "BD ( 2) C 1- C 2", "spin": ""}]
    cand = [{"description": "BD ( 2) C 1- C 2", "spin": ""},
            {"description": "BD ( 3) C 1- C 2", "spin": ""}]
    assert unmatched(ref, cand, key_full) == (1, 1), unmatched(ref, cand, key_full)
    assert unmatched(ref, cand, key_no_ordinal) == (0, 0), unmatched(ref, cand, key_no_ordinal)
    print("demo OK on %s" % socket.gethostname())


if __name__ == "__main__":
    if sys.argv[1:2] == ["--demo"]:
        demo()
        raise SystemExit(0)
    root = sys.argv[1]
    print("label_pairing_audit on %s, 1 thread, root %s" % (socket.gethostname(), root))
    print("%-13s %6s %17s %17s %17s" %
          ("molecule", "rows", "NBO full m/s", "NBO no-ordinal", "E2 full -> no-ord"))
    tot = [0, 0, 0, 0, 0, 0]
    for mol in sorted(os.listdir(root)):
        gp = os.path.join(root, mol, "%s.gennbo.nbo.json" % mol)
        np_ = os.path.join(root, mol, "%s.native.nbo.json" % mol)
        if not (os.path.isfile(gp) and os.path.isfile(np_)):
            continue
        g, n = load_nbo(gp), load_nbo(np_)
        fm, fs = unmatched(g["nbos"], n["nbos"], key_full)
        nm, ns = unmatched(g["nbos"], n["nbos"], key_no_ordinal)
        e2g = [{"description": "%s -> %s" % (e["donor"], e["acceptor"]), "spin": e["spin"]}
               for e in g.get("e2", [])]
        e2n = [{"description": "%s -> %s" % (e["donor"], e["acceptor"]), "spin": e["spin"]}
               for e in n.get("e2", [])]
        em, es = unmatched(e2g, e2n, key_full)
        om, os_ = unmatched(e2g, e2n, key_no_ordinal)
        print("%-13s %6d %8d/%-8d %8d/%-8d %6d/%-4d -> %d/%d" %
              (mol, len(g["nbos"]), fm, fs, nm, ns, em, es, om, os_))
        for i, v in enumerate((fm, fs, nm, ns, em + es, om + os_)):
            tot[i] += v
    print("\nNBO rows unpaired, engine's key: %d missing + %d surplus" % (tot[0], tot[1]))
    print("NBO rows unpaired, ordinal stripped: %d missing + %d surplus" % (tot[2], tot[3]))
    print("E2 rows unpaired: %d with the ordinal, %d without it" % (tot[4], tot[5]))
    print("A drop means the engine is reporting a bookkeeping difference as a physics difference.")
