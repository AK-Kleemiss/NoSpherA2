import os
import shutil
import platform
import subprocess
import difflib
from pathlib import Path

import pytest

class NosTest:
    def __init__(self, name, _dir=None, _good=None, _actual="NoSpherA2.log", _full=False, _skip_mac=False, **kwargs):
        self.name = name
        # Mirrors CMake logic: if no directory is specified, assume it's the test name
        self.dir = _dir or name
        self.good = _good or f"{name}.good"
        self.actual = _actual
        self.full = _full
        self.skip_mac = _skip_mac
        self.cli_args = kwargs

    def build_cmd(self, exe_path):
        cmd = [exe_path]
        for k, v in self.cli_args.items():
            flag = f"-{k}"
            if isinstance(v, bool):
                if v: cmd.append(flag)
            elif isinstance(v, (list, tuple)):
                cmd.append(flag)
                cmd.extend(str(x) for x in v)
            else:
                cmd.extend([flag, str(v)])
        cmd.append("-all_charges")
        return cmd

TESTS = [
    NosTest("alanine_occ", cif="alanine.cif", dmin=0.5, wfn="alanine.owf.fchk", acc=1, no_date=True),
    NosTest("alanine_integrated_occ", occ="alanine.toml", cif="alanine.cif", dmin=0.5, acc=1, no_date=True),
    NosTest("disorder_THPP", _dir="disorder", _good="disorder_THPP.good",
            cif="thpp.cif", hkl="thpp.hkl", acc=1, no_date=True,
            mtc=["olex2/Wfn_job/Part_1/thpp.wfx", 0.1, "olex2/Wfn_job/Part_2/thpp.wfx", 0.2],
            mtc_mult=[1, 1], mtc_ECP=[0, 0], mtc_charge=[0, 0]),
    NosTest("fourier_transform", _dir="SALTED", test_analytical=True, no_date=True),
    NosTest("fractal", _dir="sucrose_fchk_SF", _good="fractal.good", _actual="sucrose_diff.cube_fractal_plot",
            fractal="sucrose_diff.cube", no_date=True),
    NosTest("grown_water", _dir="grown", cif="water.cif", hkl="water.hkl", wfn="water.wfx", acc=1, no_date=True),
    NosTest("Hybrid_mode", _dir="Hybrid", cif="ZP2.cif", dmin=0.9, acc=1, no_date=True,
            mtc=["ZP2_part1.gbw", 0.1, "ZP2_part2.gbw", 0.2], mtc_mult=[1, 1], mtc_charge=[0, 0], mtc_ECP=[0, 0]),
    NosTest("malbac_SF_ECP", _dir="ECP_SF", cif="malbac.cif", hkl="malbac.hkl", wfn="malbac.gbw", acc=1, ECP=1, no_date=True),
    NosTest("openBLAS", _dir="OpenBLAS", blastest=True, no_date=True),
    NosTest("properties", _dir="sucrose_fchk_SF", wfn="olex2/Wfn_job/sucrose.wfx", cif="sucrose.cif", lap=True, eli=True, rdg=True, DEF=True, resolution=1.5, no_date=True),
    NosTest("reading_SALTED", _dir="SALTED", test_reading_SALTED_binary=True),
    NosTest("ri_fit", _dir="epoxide_gbw", _good="ri_fit.good", _skip_mac=True,
            wfn="epoxide.gbw", cif="epoxide.cif", dmin=0.4, ri_fit="combo_basis_fit", no_date=True),
    NosTest("rubredoxin_cmtc", _dir="rubredoxin_cmtc", hkl="1yk4_h.hkl", acc=1, no_date=True,
            cmtc=["residues/1.gbw", "residues/1.cif", 0, "residues/2.gbw", "residues/2.cif", 0.1, "residues/3.gbw", "residues/3.cif", 0.2, "residues/4.gbw", "residues/4.cif", 0],
            mtc_mult=[1, 1, 1, 1], mtc_charge=[0, 0, 0, 0], mtc_ECP=[0, 0, 0, 0]),
    NosTest("SALTED", _dir="SALTED", SALTED="Model", cif="test_cysteine.cif", wfn="test_cysteine.xyz", dmin=0.73, no_date=True),
    NosTest("sucrose_IAM", _dir="sucrose_IAM_SF", cif="sucrose.cif", hkl="sucrose.hkl", xyz="sucrose.xyz", IAM=True, no_date=True),
    NosTest("sucrose_ptb", _dir="sucrose_IAM_SF", cif="sucrose.cif", dmin=0.8, wfn="wfn.xtb", mult=0, charge=0, acc=1, ECP=3, no_date=True),
    NosTest("sucrose_SF", _dir="sucrose_fchk_SF", cif="sucrose.cif", hkl="olex2/Wfn_job/sucrose.hkl", wfn="olex2/Wfn_job/sucrose.wfx", acc=1, no_date=True),
    NosTest("sucrose_twin", _dir="sucrose_fchk_SF", cif="sucrose.cif", hkl="olex2/Wfn_job/sucrose.hkl", wfn="olex2/Wfn_job/sucrose.wfx", acc=1, twin=[1, 0, 0, 0, -1, 0, 0, 1, -2], no_date=True),
    NosTest("wfn_reading", _dir="wfn_reading", hkl="test.hkl", acc=1, cif="test.cif", wfn="test.wfn", no_date=True),
    NosTest("fchk_conversion", _dir="NiP3_fchk", _good="good.fchk", _full=True, b="dev2-TZVP", d="./", wfn="in.ffn", no_date=True),
    NosTest("fourier_transform_full", _dir="SALTED", _full=True, test_analytical="full", no_date=True)
]

@pytest.fixture(scope="session")
def exe_path():
    exe = os.environ.get("NOS_EXE", "NoSpherA2")
    if not shutil.which(exe) and not os.path.exists(exe):
        pytest.fail("Executable not found. Pass via NOS_EXE environment variable.", pytrace=False)
    return os.path.abspath(exe)

@pytest.fixture(autouse=True)
def setup_env():
    if "OCC_DATA_PATH" not in os.environ:
        os.environ["OCC_DATA_PATH"] = os.path.abspath("occ/share")


def check_differences(good_path, actual_path):
    with open(good_path, 'r') as fg, open(actual_path, 'r') as fa:
        expected = [line.strip() for line in fg if line.strip()]
        actual = [line.strip() for line in fa if line.strip()]

    if expected != actual:
        diff = list(difflib.unified_diff(
            expected, actual,
            fromfile=f'Expected ({os.path.basename(good_path)})',
            tofile=f'Actual ({os.path.basename(actual_path)})',
            lineterm=''
        ))

        colored_diff = []
        for line in diff:
            if line.startswith('+') and not line.startswith('+++'):
                colored_diff.append(f"\033[92m{line}\033[0m") # Green
            elif line.startswith('-') and not line.startswith('---'):
                colored_diff.append(f"\033[91m{line}\033[0m") # Red
            else:
                colored_diff.append(line)

        pytest.fail("Output Mismatch!\n" + "\n".join(colored_diff), pytrace=False)

@pytest.mark.parametrize("test", TESTS, ids=lambda t: t.name)
def test_nos(test, exe_path, tmp_path):
    if test.full and not os.environ.get("RUN_FULL_TEST"):
        pytest.skip("RUN_FULL_TEST not set")
    if test.skip_mac and platform.system() == "Darwin":
        pytest.skip("Skipped on macOS")

    # Isolate test execution: Copy source files to pytest's temporary directory
    src_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), test.dir)
    work_dir = tmp_path / test.dir

    if os.path.exists(src_dir):
        shutil.copytree(src_dir, work_dir)
    else:
        pytest.fail(f"Source directory missing: {src_dir}", pytrace=False)

    cmd_args = test.build_cmd(exe_path)
    print(" ".join(cmd_args))
    try:
        subprocess.run(cmd_args, cwd=work_dir, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        pytest.fail(f"Execution failed.\nStderr: {e.stderr}", pytrace=False)

    good_path = os.path.join(work_dir, test.good)
    actual_gen = os.path.join(work_dir, test.actual)
    extra_logs = work_dir.glob("*.log")
    failed_dir = Path(os.path.dirname(os.path.abspath(__file__))) / "failed_logs" / test.name
    failed_dir.mkdir(parents=True, exist_ok=True)
    for log in extra_logs:
        log.copy(failed_dir / log.name)

    if test.actual == "NoSpherA2.log":
        actual_path = os.path.join(work_dir, test.good.replace("good", "log"))
        if os.path.exists(actual_gen):
            shutil.move(actual_gen, actual_path)
    else:
        actual_path = actual_gen

    if not os.path.exists(actual_path):
        pytest.fail(f"Log {actual_path} missing.", pytrace=False)

    if not os.path.exists(good_path):
        pytest.fail(f"Expected log {good_path} missing.", pytrace=False)


    check_differences(good_path, actual_path)
    
    if test.name != "properties":
        return
    
    #Check the extra cubefiles generated from properties calculation
    cube_dir = os.path.join(work_dir, "olex2", "Wfn_job")
    cubes = ["sucrose_def.cube", "sucrose_eli.cube", "sucrose_lap.cube", "sucrose_rdg.cube", "sucrose_rho.cube", "sucrose_signed_rho.cube"]
    for cube in cubes:
        good_cube = os.path.join(cube_dir, f"{cube}.good")
        actual_cube = os.path.join(cube_dir, f"{cube}")
        
        if not os.path.exists(actual_cube):
            pytest.fail(f"Log {actual_path} missing.", pytrace=False)

        if not os.path.exists(good_cube):
            pytest.fail(f"Expected log {good_path} missing.", pytrace=False)
        
        check_differences(good_cube, actual_cube)
    
