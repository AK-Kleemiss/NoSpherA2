"""Re-run gennbo on an existing .47 with a different $NBO keylist / extra keylists.

Used for (a) the uniform final-keyword pass over the reference set and (b) the
start-dependence spread study, where the only thing that changes between runs is
the keylist section of the archive.

  py -3.12 nrt_runs.py --src bench/benzene/benzene.47 --outdir spread/benzene \
      --tag nbi_off --nbo "NRT E2PERT NRTLST=0.1 NRTNBI=off"

--extra takes a file whose content is inserted verbatim after the $NBO line
(a $CHOOSE or $NRTSTR keylist).
"""
import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


def wsl_path(p: Path) -> str:
    s = str(p.resolve()).replace("\\", "/")
    if len(s) > 1 and s[1] == ":":
        s = "/mnt/" + s[0].lower() + s[2:]
    return s


def run(src: Path, outdir: Path, tag: str, nbo: str, extra: str = "") -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    dst = outdir / (tag + ".47")
    text = src.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines(True)
    for i, line in enumerate(lines):
        if line.lstrip().startswith("$NBO"):
            lines[i] = " $NBO " + nbo + " $END\n"
            break
    else:
        raise SystemExit("no $NBO line in " + str(src))
    #An extra keylist ($CHOOSE, $NRTSTR) goes after the data blocks. Put right behind $NBO it
    #fails with "$CHOOSE error: 'END' is not an acceptable orbital type"; at the end it is read.
    if extra:
        lines.append(extra.rstrip("\n") + "\n")
    dst.write_text("".join(lines), encoding="utf-8", newline="")

    #gennbo's script deletes the old output with an interactive rm; left in place it asks and the
    #run silently keeps the previous .nbo.
    old = outdir / (tag + ".nbo")
    if old.exists():
        old.unlink()

    t0 = time.time()
    cmd = ["wsl", "bash", "-lc", "cd '%s' && ~/nbo7/gennbo %s" % (wsl_path(outdir), tag)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.time() - t0
    log = outdir / (tag + ".gennbo.log")
    log.write_text((proc.stdout or "") + (proc.stderr or ""), encoding="utf-8", errors="replace")

    out = outdir / (tag + ".nbo")
    info = {"tag": tag, "keywords": nbo, "extra": extra, "wrapper_wall_seconds": round(wall, 2),
            "exit": proc.returncode, "output": str(out), "bytes": out.stat().st_size if out.exists() else 0}
    if out.exists():
        txt = out.read_text(encoding="utf-8", errors="replace")
        m = re.search(r"QPNRT\((\d+)/(\d+)\): D\(0\)=([\d.]+); D\(w\)=([\d.]+)", txt)
        if m:
            info.update(retained=int(m.group(1)), candidates=int(m.group(2)),
                        d_0=float(m.group(3)), d_w=float(m.group(4)))
        m = re.search(r"completed in ([\d.]+) CPU seconds \((\d+) wall seconds\)", txt)
        if m:
            info.update(nbo_cpu_seconds=float(m.group(1)), nbo_wall_seconds=int(m.group(2)))
    return info


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--nbo", required=True)
    ap.add_argument("--extra", default="")
    a = ap.parse_args()
    extra = Path(a.extra).read_text(encoding="utf-8") if a.extra else ""
    info = run(Path(a.src), Path(a.outdir), a.tag, a.nbo, extra)
    print(json.dumps(info, indent=1))
    return 0 if info["exit"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
