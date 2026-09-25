"""Prove the in-place re-parse reproduces the fixed parse, field for field.

nbo_ref_v2_fixed/ was the correct copy; nbo_ref_v2/ now holds an in-place re-parse of the same
kept .nbo text.  If the two agree everywhere except the new parser_version stamp (and the timings,
which are wall clock), there is one truth and the split copy can be retired rather than trusted.
A difference anywhere else is a finding, not a formatting detail, so print it.

    python3 verify_one_truth.py
"""
import json
import os
import socket

NEW = "/work/akkleemiss/florian/nbo_ref_v2"
FIXED = "/work/akkleemiss/florian/nbo_ref_v2_fixed"
IGNORE = {"parser_version", "timings", "source"}


def diffs(a, b, path=""):
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            if k in IGNORE:
                continue
            if k not in a or k not in b:
                out.append("%s/%s only in %s" % (path, k, "new" if k in a else "fixed"))
            else:
                out += diffs(a[k], b[k], path + "/" + k)
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return ["%s length %d vs %d" % (path, len(a), len(b))]
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += diffs(x, y, "%s[%d]" % (path, i))
        return out
    if isinstance(a, float) or isinstance(b, float):
        try:
            if abs(float(a) - float(b)) <= 1e-12:
                return []
        except (TypeError, ValueError):
            pass
    return [] if a == b else ["%s %r vs %r" % (path, a, b)]


def demo():
    """An ignored key must not count, a changed number must."""
    assert diffs({"parser_version": 2, "x": 1.0}, {"x": 1.0}) == []
    assert diffs({"x": 1.0}, {"x": 1.5})
    print("demo OK on %s" % socket.gethostname())


if __name__ == "__main__":
    demo()
    print("verify_one_truth on %s, 1 thread" % socket.gethostname())
    total = clean = 0
    for mol in sorted(os.listdir(FIXED)):
        f = os.path.join(FIXED, mol, "%s.gennbo.nbo.json" % mol)
        n = os.path.join(NEW, mol, "%s.gennbo.nbo.json" % mol)
        if not (os.path.exists(f) and os.path.exists(n)):
            continue
        total += 1
        d = diffs(json.load(open(n)), json.load(open(f)))
        if d:
            print("%-14s %d differences, first 3: %s" % (mol, len(d), "; ".join(d[:3])))
        else:
            clean += 1
    print("=== %d of %d molecules identical to the fixed copy outside %s" % (clean, total, sorted(IGNORE)))
