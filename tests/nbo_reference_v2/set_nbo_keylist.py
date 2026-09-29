"""Add or replace keywords on a .47 archive's $NBO line, leaving every other byte alone.

    python set_nbo_keylist.py <file.47> "E2PERT=0.02" [more keywords ...]

A threshold probe must change the threshold and nothing else: rewriting the whole keylist would
silently drop NRTSYM=off or an NRT request that stage 2 put there, and the comparison would then
be measuring two different analyses. So a keyword given here replaces the existing one with the
same name (matched on the part before '=') and is otherwise appended.
"""
import sys


def new_keylist(old, additions):
    """The $NBO line's keyword list with `additions` applied.

    >>> new_keylist("$NBO NRT NRTSYM=off E2PERT=0.5 $END", ["E2PERT=0.02"])
    '$NBO NRT NRTSYM=off E2PERT=0.02 $END'
    """
    toks = old.split()
    assert toks[0].upper() == "$NBO", toks[:1]
    end = toks[-1].upper() == "$END"
    body = toks[1:-1] if end else toks[1:]
    for add in additions:
        name = add.split("=")[0].upper()
        for i, t in enumerate(body):
            if t.split("=")[0].upper() == name:
                body[i] = add
                break
        else:
            body.append(add)
    return " ".join(["$NBO"] + body + (["$END"] if end else []))


def main(argv):
    path, additions = argv[1], [w for a in argv[2:] for w in a.split()]
    lines = open(path).read().splitlines(True)
    for i, line in enumerate(lines):
        if line.strip().upper().startswith("$NBO"):
            lines[i] = " " + new_keylist(line.strip(), additions) + "\n"
            break
    else:
        raise SystemExit("no $NBO line in " + path)
    open(path, "w").writelines(lines)
    print("keylist now: " + lines[i].strip())
    return 0


def demo():
    assert new_keylist("$NBO NRT E2PERT=0.5 $END", ["E2PERT=0.02"]) == \
        "$NBO NRT E2PERT=0.02 $END"
    assert new_keylist("$NBO NRTSYM=off $END", ["E2PERT=0.02"]) == \
        "$NBO NRTSYM=off E2PERT=0.02 $END"
    # the existing keylist must survive: this is the whole point of not rewriting the line
    assert "NRTSYM=off" in new_keylist("$NBO NRT NRTSYM=off NRTLST=0.1 $END", ["E2PERT=0.02"])
    assert new_keylist("$NBO", ["NRT"]) == "$NBO NRT"
    print("demo ok")


if __name__ == "__main__":
    sys.exit(demo() if "--demo" in sys.argv else main(sys.argv))
