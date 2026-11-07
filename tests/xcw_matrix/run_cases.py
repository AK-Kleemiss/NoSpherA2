"""Run XCW cases on the NH3 fixture, time them, and print R1 / wR2 / GooF / scale / extinction.

    python run_cases.py --exe D:/git/NoSpherA2/spaces/3f316b01-r1sym/NoSpherA2.exe [--cases timing|timing_H2|eri|reference|matrix|matrix_sigma] [--cpus 1]

Each case runs in its own folder under --work. wR2 is not printed by NoSpherA2, so it is computed here from
nh3.hkl and the last lambda's NA2_*_Fcalc.txt (-print_F): wR2 = sqrt(sum w (Fo^2 - Fc^2)^2 / sum w Fo^4) over all
reflections, w = 1/sigma(Fo^2)^2. Fc is the file's scale*sqrt(y)*|F_calc| column. A case's 'weight' entry is passed
through verbatim into the settings, for an exe that knows the keyword (Briefs 1-2).
"""
import argparse
import itertools
import math
import pathlib
import re
import shutil
import subprocess
import time

import psutil

HERE = pathlib.Path(__file__).resolve().parent
CIF = {"none": "nh3.cif", "C": "nh3_C.cif", "CD": "nh3_CD.cif"}
CONV = "normal normal_conv"  # --conv overrides, e.g. 'sloppy normal_conv'
EXT = {"off": "", "fixed": "extinction shelx fixed 0.2", "refined": "extinction shelx"}


def settings(c):
    s = (f"{c.get('conv', CONV)} params 22 basis_set {c.get('basis', 'sto-3g')} {c['FF']} charge 0 mult 1 rhf "
         f"start 0 step_size 0.01 end {0.01 * (c.get('nlam', 1) - 1):.2f}")
    # extinction last only for older exes, whose bare `extinction shelx` read the next keyword as its value
    for extra in ("grown" if c["grown"] else "", "weighted" if c.get("target") == "H2" else "", c.get("weight", ""), c.get("extra", ""), EXT[c["ext"]]):
        if extra:
            s += " " + extra
    return s + "\n"


def run(exe, c, work, cpus):
    d = work / c["name"]
    if d.exists() and not c.get("keep"):
        shutil.rmtree(d)
    d.mkdir(parents=True, exist_ok=True)
    for f in ("nh3.hkl", CIF[c["anh"]], "nh3_grown.xyz"):
        shutil.copy(HERE / f, d / f)
    (d / "s.txt").write_text(settings(c))
    cmd = [exe, "-do_XCW", "-no_gpu", "-no_date", "-print_F", "-cpus", str(cpus), "-cif", CIF[c["anh"]],
           "-hkl", "nh3.hkl", "-XCW_settings", "s.txt"] + (["-xyz", "nh3_grown.xyz"] if c["grown"] else [])
    # cpu = the child's user+sys seconds, last sample before exit (10 ms polling); wall is load-sensitive, cpu much less
    t0, cpu = time.perf_counter(), None
    with open(d / "stdout.txt", "w") as out:
        p = psutil.Popen(cmd, cwd=d, stdout=out, stderr=subprocess.STDOUT)
        while p.poll() is None:
            try:
                ct = p.cpu_times()
                cpu = ct.user + ct.system
            except psutil.Error:
                pass
            time.sleep(0.01)
    wall = time.perf_counter() - t0
    res = {"name": c["name"], "wall": wall, "cpu": cpu, "rc": p.returncode}
    log = (d / "NoSpherA2.log").read_text(errors="replace") if (d / "NoSpherA2.log").exists() else ""
    # columns by the tab-separated header ("Total Energy" is one column of two words), values by name from the
    # "XCW final model:" line
    head = next((l for l in log.splitlines() if l.lstrip().startswith("Lambda\t")), None)
    rows = [l.split() for l in log.splitlines() if re.match(r"^\s*\d+\.\d{5}\s", l)]
    if head and rows:
        row = dict(zip(re.split(r"\t+", head.strip()), rows[-1]))
        res.update(lam=float(row["Lambda"]), crit=float(row["Criterion"]))
    m = re.search(r"^XCW final model: (.*)$", log, re.M)
    if m:
        # pairs only up to the scale: "extinction(SHELX (Zachariasen)) p" and the weight report are not key/value
        kv = dict(zip(*[iter(m.group(1).split(" scale ")[0].split())] * 2))
        res.update(lam=float(kv["lambda"]), r1=float(kv["R1"]), wr2=float(kv["wR2"]), r1all=float(kv["R1(all)"]),
                   wr2all=float(kv["wR2(all)"]), goof=float(kv["GooF"]), wscheme=m.group(1).split(" weight ")[1].split()[0])
    m = re.search(r"Generated total of (\d+) grid", log)
    res["grid"] = int(m.group(1)) if m else None
    scf = (d / "SCF.log").read_text(errors="replace") if (d / "SCF.log").exists() else ""
    m = re.search(r"(\d+) grown,", scf)
    res["grown_atoms"] = int(m.group(1)) if m else None
    ext = re.findall(r"lambda [\d.]+: extinction.*\s([-\d.eE+]+)\s*$", scf, re.M)
    res["ext"] = float(ext[-1]) if ext else None
    fc = sorted(d.glob("NA2_*_Fcalc.txt"))
    if fc:
        lines = fc[-1].read_text().splitlines()
        res["scale"] = float(lines[0].split("scale =")[1])
        fcalc = {tuple(map(int, l.split()[:3])): float(l.split()[5]) for l in lines[1:] if l.strip()}
        num = den = 0.0
        for l in (HERE / "nh3.hkl").read_text().splitlines():
            h, k, ll, i_o, s = int(l[0:4]), int(l[4:8]), int(l[8:12]), float(l[12:20]), float(l[20:28])
            fc2 = fcalc[(h, k, ll)] ** 2
            num += (i_o - fc2) ** 2 / s ** 2
            den += i_o ** 2 / s ** 2
        # independent cross-check of the printed wR2(all): sigma weights, so it only has to agree on sigma rows
        res["wr2s"] = math.sqrt(num / den)
        if res.get("wscheme") == "sigma" and res.get("wr2all") is not None:
            assert abs(res["wr2s"] - res["wr2all"]) < 2e-5, (c["name"], res["wr2s"], res["wr2all"])
    eri = list(d.glob("*.eri*")) + list(d.glob("eri*"))
    res["eri_bytes"] = sum(f.stat().st_size for f in eri) if eri else None
    if p.returncode:
        res["err"] = (d / "stdout.txt").read_text(errors="replace").strip().splitlines()[-3:]
    return res


def timing_cases(targets=("normal", "H2")):
    out = []
    for tgt, grown, FF, anh, ext, nlam in itertools.product(targets, (False, True), ("F", "F2"), ("none", "C", "CD"),
                                                             ("off", "fixed", "refined"), (1, 3)):
        out.append(dict(name=f"{tgt}_{'g' if grown else 'u'}_{FF}_{anh}_{ext}_{nlam}", target=tgt, grown=grown, FF=FF,
                        anh=anh, ext=ext, nlam=nlam))
    return out


def eri_cases():
    base = dict(grown=False, FF="F", anh="CD", ext="refined", nlam=3)
    gbase = dict(base, grown=True)
    return [dict(base, name="eri_save", extra="save_eri nh3.eri"),
            dict(base, name="eri_save", extra="read_eri nh3.eri", keep=True),
            dict(gbase, name="eri_gsave", extra="save_eri nh3.eri"),
            dict(gbase, name="eri_gsave", extra="read_eri nh3.eri", keep=True)]


def reference_cases():
    out = []
    for tgt, grown, FF, anh, ext in itertools.product(("normal", "H2"), (False, True), ("F", "F2"), ("none", "CD"),
                                                      ("off", "refined")):
        out.append(dict(name=f"ref_{tgt}_{'g' if grown else 'u'}_{FF}_{anh}_{ext}", target=tgt, grown=grown, FF=FF,
                        anh=anh, ext=ext, nlam=3))
    return out


# Pairwise covering design over weight(8) x anh(3) x ext(3) x grown(2) x F/F2(2) x target(2): every pair of levels of
# any two factors appears in at least one of rows 1-24 (24 = 8 x 3, the minimum); 'D' is the all-default row and
# 'ON' the all-on row. Every row runs 3 lambda steps so the *_opt and robust schemes get an update between steps.
MATRIX = [  # id, weight, anh, ext, grown, F/F2, target
    ("D", "default", "none", "off", "u", "F", "normal"),
    ("1", "sigma", "none", "refined", "u", "F", "H2"),
    ("2", "sigma", "C", "fixed", "g", "F2", "normal"),
    ("3", "sigma", "CD", "off", "g", "F2", "H2"),
    ("4", "unit", "none", "off", "g", "F", "normal"),
    ("5", "unit", "C", "refined", "u", "F2", "normal"),
    ("6", "unit", "CD", "fixed", "u", "F", "H2"),
    ("7", "stl", "none", "fixed", "u", "F2", "normal"),
    ("8", "stl", "C", "off", "g", "F", "H2"),
    ("9", "stl", "CD", "refined", "g", "F", "normal"),
    ("10", "shelx", "none", "off", "u", "F2", "H2"),
    ("11", "shelx", "C", "fixed", "g", "F", "normal"),
    ("12", "shelx", "CD", "refined", "g", "F", "normal"),
    ("13", "shelx_opt", "none", "off", "g", "F", "H2"),
    ("14", "shelx_opt", "C", "fixed", "u", "F2", "normal"),
    ("15", "shelx_opt", "CD", "refined", "g", "F", "H2"),
    ("16", "chebyshev", "none", "refined", "g", "F", "normal"),
    ("17", "chebyshev", "C", "off", "u", "F2", "H2"),
    ("18", "chebyshev", "CD", "fixed", "g", "F", "H2"),
    ("19", "chebyshev_opt", "none", "refined", "g", "F", "H2"),
    ("20", "chebyshev_opt", "C", "off", "u", "F2", "normal"),
    ("21", "chebyshev_opt", "CD", "fixed", "g", "F2", "H2"),
    ("22", "robust", "none", "refined", "u", "F2", "normal"),
    ("23", "robust", "C", "fixed", "g", "F", "H2"),
    ("24", "robust", "CD", "off", "g", "F", "normal"),
    ("ON", "chebyshev_opt", "CD", "refined", "g", "F2", "H2"),
]
# Placeholder parameters until Briefs 1-2 fix the syntax. sigma is the default, so its rows pass no keyword and run
# on today's exe; 'default' never passes one.
WEIGHT_KW = {"default": "", "sigma": "", "unit": "weight unit", "stl": "weight stl 1",
             "shelx": "weight shelx 0.1 0", "shelx_opt": "weight shelx_opt 0.1 0",
             "chebyshev": "weight chebyshev 1 0.5 0.2", "chebyshev_opt": "weight chebyshev_opt 3",
             "robust": "weight robust shelx 0.1 0 6"}


# normal: with the TRAH-DIIS stall fix every row is under 3 s; sloppy moved the chebyshev_opt rows (m19 wR2 0.143 vs 0.172)
MATRIX_CONV = "normal normal_conv"


def check_pairwise():
    """Assert rows 1-24 of MATRIX cover every pair of levels of every two factors."""
    rows = [r[1:] for r in MATRIX if r[0].isdigit()]
    levels = [sorted({r[i] for r in rows}) for i in range(6)]
    assert [len(l) for l in levels] == [8, 3, 3, 2, 2, 2], levels
    for i, j in itertools.combinations(range(6), 2):
        missing = set(itertools.product(levels[i], levels[j])) - {(r[i], r[j]) for r in rows}
        assert not missing, (i, j, missing)


def matrix_cases(only=None):
    check_pairwise()
    return [dict(name=f"m{rid}_{w}_{anh}_{ext}_{g}_{FF}_{tgt}", weight=WEIGHT_KW[w], anh=anh, ext=ext, grown=g == "g",
                 FF=FF, target=tgt, nlam=3, conv=MATRIX_CONV) for rid, w, anh, ext, g, FF, tgt in MATRIX if only is None or w in only]


def fmt(r):
    keys = ("wall", "cpu", "rc", "lam", "r1", "r1all", "wr2", "goof", "crit", "scale", "ext", "grid", "grown_atoms",
            "eri_bytes")
    return r["name"] + "  " + "  ".join(f"{k}={r[k]:.4g}" if isinstance(r.get(k), float) else f"{k}={r.get(k)}"
                                        for k in keys) + ("  ERR " + " | ".join(r["err"]) if "err" in r else "")


def main():
    global CONV
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--work", default="E:/Temp/xcw_matrix_scratch/cases")
    ap.add_argument("--cases", default="timing", choices=("timing", "timing_H2", "eri", "reference", "matrix", "matrix_sigma"))
    ap.add_argument("--cpus", type=int, default=1)
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--conv", default=CONV)
    ap.add_argument("--golden", help="with --cases matrix: write the ctest reference (XcwMatrixTests) to this file")
    a = ap.parse_args()
    CONV = a.conv
    cases = {"timing": timing_cases, "timing_H2": lambda: timing_cases(("H2",)), "eri": eri_cases, "reference": reference_cases,
             "matrix": matrix_cases, "matrix_sigma": lambda: matrix_cases(("default", "sigma"))}[a.cases]()
    total = 0.0
    golden = ["# XcwMatrixTests reference, written by run_cases.py --cases matrix --golden (no timings, by design)",
              "# id\tcif\tgrown\tlambda\tR1\twR2\tGooF\tscale\textinction\tsettings"]
    for c in cases:
        best = None
        for _ in range(a.repeat):
            r = run(a.exe, c, pathlib.Path(a.work), a.cpus)
            best = r if best is None or r["wall"] < best["wall"] else best
        total += best["wall"]
        print(fmt(best), flush=True)
        assert not a.golden or best["rc"] == 0, best
        golden.append("\t".join([c["name"].split("_")[0][1:], CIF[c["anh"]], str(int(c["grown"]))] +
                                [repr(best[k]) for k in ("lam", "r1", "wr2", "goof", "scale")] +
                                [repr(best["ext"]) if best["ext"] is not None else "-", settings(c).strip()]))
    print(f"total {total:.1f} s over {len(cases)} cases")
    if a.golden:
        pathlib.Path(a.golden).write_text("\n".join(golden) + "\n")


if __name__ == "__main__":
    main()
