#!/usr/bin/python
import os.path
from subprocess import run
import sys

if __name__ == "__main__":
    print("argv:", sys.argv)
    run(sys.argv[3], shell=True)
    name = sys.argv[1]
    good = sys.argv[2]
    with open(good) as f:
        good = f.read().strip()

    with open(name) as f:
        log = f.read().strip()

    if good != log:
        raise Exception(f"Test failed:\nLog:\n{log}")
