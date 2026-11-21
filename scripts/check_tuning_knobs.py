#!/usr/bin/env python3
"""Check that docs/wiki/Tuning-Knobs.md lists every -tune knob the sources read, and nothing else.

A knob is a literal name passed to tuning() or to one of its wrappers (tuning_double in basin_integration.cpp,
read_threshold in rgbi.cpp, nao_knob in nao.cpp). A knob read through a variable or a
new wrapper is invisible here: add the wrapper to CALL. It also checks that the table -tune list prints
(documented_knobs in Src/core/base/tuning.cpp) holds the page's first two columns, without the backticks.
Exit status 1 when the sets differ.
"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAGE = os.path.join(ROOT, "docs", "wiki", "Tuning-Knobs.md")
CALL = re.compile(r'\b(?:tuning|tuning_double|read_threshold|nao_knob)\(\s*"([A-Z][A-Z0-9_]+)"')
#A knob name in the first column of a table row, with or without =VALUE
ROW = re.compile(r"^\| `([A-Z][A-Z0-9_]+)[`=]", re.M)
REGISTRY = os.path.join(ROOT, "Src", "core", "base", "tuning.cpp")
ENTRY = re.compile(r'^\t\{"([^"]+)", "([^"]*)"\},\r?$', re.M)


def source_knobs():
    found = {}
    for top in ("Src", "app"):
        for d, _, files in os.walk(os.path.join(ROOT, top)):
            for n in files:
                if not n.endswith((".cpp", ".h", ".cu", ".cuh", ".hip")):
                    continue
                p = os.path.join(d, n)
                for k in CALL.findall(io.open(p, encoding="utf-8", errors="replace").read()):
                    found.setdefault(k, os.path.relpath(p, ROOT).replace(os.sep, "/"))
    return found


def main():
    src = source_knobs()
    #the environment variables at the end are not knobs
    text = io.open(PAGE, encoding="utf-8").read().split("\n## Not knobs")[0]
    page = set(ROW.findall(text))
    missing = sorted(set(src) - page)
    stale = sorted(page - set(src))
    for k in missing:
        print("not on the page: %s (%s)" % (k, src[k]))
    for k in stale:
        print("on the page, read by no source: %s" % k)
    rows = set()
    for line in text.split("\n"):
        if line.startswith("| `"):
            c = [x.strip().replace("`", "").replace("\\|", "|") for x in line.strip().strip("|").split(" | ")]
            rows.add((c[0], c[1]))
    table = set(ENTRY.findall(io.open(REGISTRY, encoding="utf-8").read()))
    for r in sorted(rows - table):
        print("on the page, not in documented_knobs: %s = %s" % r)
    for r in sorted(table - rows):
        print("in documented_knobs, not on the page: %s = %s" % r)
    differ = missing or stale or rows != table
    if not differ:
        print("%d knobs, page, sources and -tune list agree" % len(src))
    return 1 if differ else 0


if __name__ == "__main__":
    sys.exit(main())
