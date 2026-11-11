"""The shared library the way Olex2 loads it (ctypes, nosphera2_run): an error comes back as a return code
instead of ending the host, and a repeated call in the same process gives the same result (no static that
only the first run clears). Run with pytest after a build; NOS_DLL names the library, else build/*/bin is
searched, and the test skips without one."""
import ctypes
import glob
import os
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
NAME = {"win32": "NoSpherA2.dll", "darwin": "libNoSpherA2.dylib"}.get(sys.platform, "libNoSpherA2.so")


def find_dll():
    if os.environ.get("NOS_DLL"):
        return os.environ["NOS_DLL"]
    hits = sorted(glob.glob(str(ROOT / "build" / "*" / "bin" / NAME)), key=os.path.getmtime)
    return hits[-1] if hits else None


@pytest.fixture(scope="module")
def dll():
    path = find_dll()
    if not path:
        pytest.skip(f"no {NAME} built")
    lib = ctypes.CDLL(os.path.abspath(path))
    lib.nosphera2_run.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)]
    lib.nosphera2_run.restype = ctypes.c_int
    return lib


def run(lib, *args):
    argv = [a.encode("mbcs" if os.name == "nt" else "utf-8") for a in ("NoSpherA2",) + args]
    return lib.nosphera2_run(len(argv), (ctypes.c_char_p * len(argv))(*argv))


def write_cube(path: Path, n: int = 12) -> None:
    lines = ["fractal test", "", "    1    0.000000    0.000000    0.000000"]
    lines += [f"{n:5d}" + "".join(f"{0.3 if i == j else 0.0:12.6f}" for j in range(3)) for i in range(3)]
    lines.append("    6    6.000000    0.000000    0.000000    0.000000")
    vals = [((x - n / 2) ** 2 + (y - n / 2) ** 2 + (z - n / 2) ** 2) for x in range(n) for y in range(n) for z in range(n)]
    lines += [" ".join(f"{2.0 ** (-v / 8):.5E}" for v in vals[k:k + 6]) for k in range(0, len(vals), 6)]
    path.write_text("\n".join(lines) + "\n")


def test_error_returns_code_and_host_survives(dll):
    with tempfile.TemporaryDirectory() as t:
        missing = str(Path(t) / "missing.wfn")
        first = run(dll, "-wfn", missing, "-cif", str(Path(t) / "missing.cif"), "-out", str(Path(t) / "a.log"))
        second = run(dll, "-wfn", missing, "-cif", str(Path(t) / "missing.cif"), "-out", str(Path(t) / "b.log"))
        assert first != 0 and first == second


def test_fractal_twice_in_one_process(dll):
    with tempfile.TemporaryDirectory() as t:
        cube = Path(t) / "map.cube"
        write_cube(cube)
        plots = []
        for k in range(2):
            assert run(dll, "-fractal", str(cube), "-out", str(Path(t) / f"f{k}.log")) == 0
            plot = Path(str(cube) + "_fractal_plot")
            assert plot.exists()
            plots.append(plot.read_text())
            plot.unlink()
        assert plots[0] == plots[1] and plots[0].strip()
