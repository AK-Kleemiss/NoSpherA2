#!/usr/bin/env python3
"""Check that docs/wiki/Tuning-Knobs.md lists every -tune knob the sources read, and nothing else.

A knob is a literal name passed to tuning() or to one of its wrappers (tuning_double in b2c.cpp,
read_threshold in bondwise_analysis.cpp, nao_knob in nao.cpp). A knob read through a variable or a
new wrapper is invisible here: add the wrapper to CALL. Exit status 1 when the two sets differ.
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
    if not missing and not stale:
        print("%d knobs, page and sources agree" % len(src))
    return 1 if missing or stale else 0


if __name__ == "__main__":
    sys.exit(main())
