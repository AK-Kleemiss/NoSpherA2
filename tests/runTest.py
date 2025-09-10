#!/usr/bin/python
import os.path
from subprocess import run
import sys

if __name__ == "__main__":
    run(sys.argv[2], shell=True)
    name = sys.argv[1]
    with open(f"{name}") as f:
        good = f.read().strip()

    with open("NoSpherA2.log") as f:
        log = f.read().strip()

    if good != log:
        raise Exception("Test failed")
