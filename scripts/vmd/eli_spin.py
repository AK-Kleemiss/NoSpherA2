#!/usr/bin/env python3
"""One command from a wavefunction to a VMD scene of its spin-resolved ELI-D.

    python scripts/vmd/eli_spin.py mol.fchk [--resolution 0.1] [--radius 3.5] [--rho] [--vmd]
                                            [make_eli_spin_vmd.py options ...]

Runs ``NoSpherA2 -wfn mol.fchk -eli_spin -resolution R -radius A`` (with ``-rho`` when --rho or
--rho-clip is given), writes ``mol_eli_spin.vmd`` with make_eli_spin_vmd.py, and with --vmd starts
``vmd -e mol_eli_spin.vmd`` when vmd is on PATH. NoSpherA2 is taken from --nospherA2, $NOS_EXE, PATH,
or the build trees of this checkout, in that order.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_eli_spin_vmd  # noqa: E402


def find_nospherA2(given: str | None) -> str:
    for c in (given, os.environ.get("NOS_EXE"), shutil.which("NoSpherA2")):
        if c:
            return c
    for p in sorted((HERE.parent.parent / "build").glob("*/bin/NoSpherA2*")):
        if p.is_file() and os.access(p, os.X_OK):
            return str(p)
    raise SystemExit("NoSpherA2 not found: pass --nospherA2 or set NOS_EXE.")


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("wfn", type=Path)
    ap.add_argument("--resolution", default="0.1")
    ap.add_argument("--radius", default="3.5")
    ap.add_argument("--rho", action="store_true", help="also write <stem>_rho.cube")
    ap.add_argument("--vmd", action="store_true", help="open the scene in VMD")
    ap.add_argument("--nospherA2")
    a, rest = ap.parse_known_args(argv)
    wfn = a.wfn.resolve()
    cmd = [find_nospherA2(a.nospherA2), "-wfn", wfn.name, "-eli_spin", "-resolution", a.resolution, "-radius", a.radius]
    if a.rho or "--rho-clip" in rest:
        cmd.append("-rho")
    print(" ".join(cmd))
    subprocess.run(cmd, cwd=wfn.parent, check=True)
    make_eli_spin_vmd.main([str(wfn.parent / wfn.stem)] + rest)
    script = wfn.parent / f"{wfn.stem}_eli_spin.vmd"
    if a.vmd:
        vmd = shutil.which("vmd")
        if vmd is None:
            print("vmd is not on PATH; open it with:  vmd -e " + str(script))
        else:
            subprocess.Popen([vmd, "-e", str(script)], cwd=wfn.parent)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
