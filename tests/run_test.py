import glob
import os
import shutil
import platform
import subprocess
import difflib
from pathlib import Path

try:
    import tomllib
except ImportError:
    import tomli as tomllib

import pytest


class NosTest:
    def __init__(
        self,
        name,
        _dir=None,
        _good=None,
        _actual="NoSpherA2.log",
        _full=False,
        _skip_mac=False,
        **kwargs,
    ):
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
                if v:
                    cmd.append(flag)
            elif v == "":
                cmd.append(flag)
            elif isinstance(v, (list, tuple)):
                cmd.append(flag)
                cmd.extend(str(x) for x in v)
            else:
                cmd.extend([flag, str(v)])
        return cmd


_TOML_PATH = Path(__file__).parent / "tests.toml"


def _load_tests():
    with open(_TOML_PATH, "rb") as f:
        data = tomllib.load(f)
    defaults = data.pop("defaults", {})
    tests = []
    for name, entry in data.items():
        args = {**defaults, **entry.get("args", {})}
        tests.append(
            NosTest(
                name,
                _dir=entry.get("directory"),
                _good=entry.get("good"),
                _actual=entry.get("actual", "NoSpherA2.log"),
                _full=entry.get("full", False),
                _skip_mac=entry.get("skip_mac", False),
                **args,
            )
        )
    return tests


TESTS = _load_tests()


@pytest.fixture(scope="session")
def exe_path():
    exe = os.environ.get("NOS_EXE", "NoSpherA2")
    if not shutil.which(exe) and not os.path.exists(exe):
        pytest.fail(
            "Executable not found. Pass via NOS_EXE environment variable.",
            pytrace=False,
        )
    return os.path.abspath(exe)

@pytest.mark.parametrize("nos_test", TESTS, ids=lambda t: t.name)
def test(nos_test, exe_path, tmp_path, request):
    if nos_test.full and not os.environ.get("RUN_FULL_TEST"):
        pytest.skip("RUN_FULL_TEST not set")
    if nos_test.skip_mac and platform.system() == "Darwin":
        pytest.skip("Skipped on macOS")

    # Isolate test execution: Copy source files to pytest's temporary directory
    src_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), nos_test.dir)
    work_dir = tmp_path / nos_test.dir

    if os.path.exists(src_dir):
        shutil.copytree(src_dir, work_dir)
    else:
        pytest.fail(f"Source directory missing: {src_dir}", pytrace=False)

    cmd_args = nos_test.build_cmd(exe_path)
    print(" ".join(cmd_args))
    OCC_DATA_PATH = str(request.config.rootpath / "occ" / "share")
    try:
        subprocess.run(
            cmd_args, cwd=work_dir, check=True, capture_output=True, text=True,env={"OCC_DATA_PATH":OCC_DATA_PATH}
        )
    except subprocess.CalledProcessError as e:
        pytest.fail(f"Execution failed.\nStderr: {e.stderr}", pytrace=False)

    good_path = os.path.join(work_dir, nos_test.good)
    actual_gen = os.path.join(work_dir, nos_test.actual)
    extra_logs = work_dir.glob("*.log")
    failed_dir = (
        Path(os.path.dirname(os.path.abspath(__file__))) / "failed_logs" / nos_test.name
    )
    failed_dir.mkdir(parents=True, exist_ok=True)
    for log in extra_logs:
        log.copy(failed_dir / log.name)

    if nos_test.actual == "NoSpherA2.log":
        actual_path = os.path.join(work_dir, nos_test.good.replace("good", "log"))
        if os.path.exists(actual_gen):
            shutil.move(actual_gen, actual_path)
    else:
        actual_path = actual_gen

    if not os.path.exists(actual_path):
        pytest.fail(f"Log {actual_path} missing.", pytrace=False)

    if not os.path.exists(good_path):
        pytest.fail(f"Expected log {good_path} missing.", pytrace=False)

    with open(good_path, "r") as fg, open(actual_path, "r") as fa:
        expected = [line.strip() for line in fg if line.strip()]
        actual = [line.strip() for line in fa if line.strip()]

    if expected != actual:
        diff = list(
            difflib.unified_diff(
                expected,
                actual,
                fromfile=f"Expected ({nos_test.good})",
                tofile=f"Actual ({os.path.basename(actual_path)})",
                lineterm="",
            )
        )

        colored_diff = []
        for line in diff:
            if line.startswith("+") and not line.startswith("+++"):
                colored_diff.append(f"\033[92m{line}\033[0m")  # Green
            elif line.startswith("-") and not line.startswith("---"):
                colored_diff.append(f"\033[91m{line}\033[0m")  # Red
            else:
                colored_diff.append(line)

        pytest.fail("Output Mismatch!\n" + "\n".join(colored_diff), pytrace=False)
