#!/usr/bin/python
import os.path
from subprocess import run
import sys

if __name__ == "__main__":
    run(sys.argv[1:])
    name = os.path.abspath(os.path.basename("."))
    with open(f"{name}.good") as f:
        good = f.read()

    with open(f"{name}.log") as f:
        log = f.read()

    if good != log:
        raise Exception("Test failed")
