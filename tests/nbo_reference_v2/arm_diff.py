"""Compare two arms' native JSONs NUMERICALLY, so a wall clock cannot be mistaken for a result.

    python3 arm_diff.py <arm A dir> <arm B dir> [--tol 0]

`cmp` called all 22 different because every run writes its own timings; that is not a finding.
This walks both trees, skips paths under a timing dict, and reports the largest absolute
difference per molecule together with where it sits.  A path present on one side only, or a list
of different length, is reported as STRUCTURAL and is never folded into a max.
"""
import json
import os
import sys

MOLS = ("acetylene ammonia benzene ch3 ethane ethene formaldehyde formate hcn lif n2 ni_co_4 "
        "nitromethane no o2 ozone pf5 pyridine sf6 so2 ticl4 water").split()


def timing_path(path):
    #Deliberately explicit: a path is a clock if some key on it is a timing container or ends in
    #seconds.  Nothing else is excluded, so a field that happens to move is a result, not noise.
    return any(k in ("timings", "seconds", "wall_seconds", "qp_seconds", "profile") or
               (isinstance(k, str) and k.endswith("seconds")) for k in path)


def walk(a, b, path, out, struct):
    if timing_path(path):
        return
    if type(a) is not type(b) and not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
        struct.append(("type", path, type(a).__name__, type(b).__name__))
        return
    if isinstance(a, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                if not timing_path(path + (k,)):
                    struct.append(("key", path + (k,), k in a, k in b))
                continue
            walk(a[k], b[k], path + (k,), out, struct)
    elif isinstance(a, list):
        if len(a) != len(b):
            struct.append(("len", path, len(a), len(b)))
            return
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, path + (i,), out, struct)
    elif isinstance(a, bool) or isinstance(a, str) or a is None:
        if a != b:
            struct.append(("value", path, a, b))
    elif isinstance(a, (int, float)):
        out.append((abs(float(a) - float(b)), path, float(a), float(b)))


def main(argv):
    A, B = argv[1], argv[2]
    tol = float(argv[argv.index("--tol") + 1]) if "--tol" in argv else 0.0
    print("A = %s\nB = %s\ntolerance %g, timing paths excluded" % (A, B, tol))
    nmoved = nstruct = 0
    worst_all = (0.0, None, None)
    for mol in MOLS:
        pa = os.path.join(A, mol, mol + ".native.nbo.json")
        pb = os.path.join(B, mol, mol + ".native.nbo.json")
        if not (os.path.exists(pa) and os.path.exists(pb)):
            print("%-13s MISSING" % mol)
            nstruct += 1
            continue
        out, struct = [], []
        walk(json.load(open(pa)), json.load(open(pb)), (), out, struct)
        over = [o for o in out if o[0] > tol]
        w = max(out, default=(0.0, (), 0.0, 0.0))
        if w[0] > worst_all[0]:
            worst_all = (w[0], mol, w[1])
        if struct:
            nstruct += 1
            print("%-13s STRUCTURAL %d: %s" % (mol, len(struct), struct[:3]))
        if over:
            nmoved += 1
            top = sorted(over, reverse=True)[:3]
            print("%-13s %5d fields > tol, worst %.3e at %s (%.8g vs %.8g)" % (
                mol, len(over), top[0][0], ".".join(map(str, top[0][1])), top[0][2], top[0][3]))
            for d, path, x, y in top[1:]:
                print("%-13s        %.3e at %s (%.8g vs %.8g)" % ("", d, ".".join(map(str, path)), x, y))
        elif not struct:
            print("%-13s identical to %g (worst numeric difference %.3e)" % (mol, tol, w[0]))
    print("\n%d of %d molecules moved beyond %g, %d structural, worst overall %.3e on %s at %s" % (
        nmoved, len(MOLS), tol, nstruct, worst_all[0], worst_all[1],
        ".".join(map(str, worst_all[2])) if worst_all[2] else "-"))
    return 0


def demo():
    out, struct = [], []
    walk({"timings": {"x": 1.0}, "nao": [{"occupancy": 1.0}], "n": 3},
         {"timings": {"x": 9.0}, "nao": [{"occupancy": 1.5}], "n": 3}, (), out, struct)
    #the clock must be invisible and the occupancy must not be
    assert not struct, struct
    assert len(out) == 2, out                                   # occupancy and n
    assert max(out)[0] == 0.5 and max(out)[1] == ("nao", 0, "occupancy"), out
    out, struct = [], []
    walk({"nao": [1.0, 2.0]}, {"nao": [1.0]}, (), out, struct)
    assert struct and struct[0][0] == "len", struct             # a length change is never a max
    print("demo ok")
    return 0


if __name__ == "__main__":
    sys.exit(demo() if "--demo" in sys.argv else main(sys.argv))
