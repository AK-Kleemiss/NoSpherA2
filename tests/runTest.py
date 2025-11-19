#!/usr/bin/python
import os.path
import difflib
from subprocess import run
import sys
import re

if __name__ == "__main__":
    print("argv:", sys.argv)
    run(sys.argv[3], shell=True)
    name = sys.argv[1]
    good = sys.argv[2]
    with open(good) as f:
        good = f.read().strip()
        good = re.sub(r'[ \t]+', ' ', good)
    with open(name) as f:
        log = f.read().strip()
        log = re.sub(r'[ \t]+', ' ', log)

    if good != log:
        raise Exception(f"Test failed:\nDiff: {''.join(difflib.unified_diff(good.splitlines(keepends=True), log.splitlines(keepends=True)))}")
