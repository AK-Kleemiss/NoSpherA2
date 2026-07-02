#!/usr/bin/env python3
"""Update test .good baselines from actual test output.

Two sources of "actual" output are recognized when --run is not used, and
whichever exists is used (pytest checked first, then ctest):

- pytest: tests/failed_logs/<test_name>/<file>, written by the Python golden-file
  harness when a comparison fails.
- ctest: tests/<directory>/<actual>, written directly by the in-process C++
  integration test runner (tests/src/IntegrationTests.cpp), which cd's into the
  test directory and leaves its output right next to the .good file.

Pass --run to have this script execute each test itself (mirroring the argv
construction in tests/src/IntegrationTests.cpp: -all_charges -no_date plus
[test.args]) right before diffing, instead of trusting whatever actual output
happens to be lying around. This avoids "oopsies" from stale or cross-test
leftovers: with --run, a test's .good is only ever updated from output that
was just produced by that exact test, under the exe you're currently pointing
at, with no dependency on artifacts from some earlier, possibly different, run.

By default this script shows a diff for each test whose actual output differs
from the current .good file and asks for confirmation before overwriting.
Pass --yes to accept every change without prompting.
"""

from __future__ import annotations

import argparse
import difflib
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

try:
    import tomllib  # Python 3.11+
except ImportError:  # pragma: no cover
    import tomli as tomllib  # type: ignore


@dataclass
class TestEntry:
    name: str
    directory: str
    good: str
    actual: str
    args: dict = field(default_factory=dict)


@dataclass
class PendingUpdate:
    test: TestEntry
    target: Path
    source: Path


def load_tests(toml_path: Path) -> list[TestEntry]:
    with toml_path.open("rb") as fd:
        data = tomllib.load(fd)

    data.pop("defaults", None)
    entries: list[TestEntry] = []
    for name, entry in data.items():
        directory = entry.get("directory", name)
        good = entry.get("good", f"{name}.good")
        actual = entry.get("actual", "NoSpherA2.log")
        test_args = entry.get("args", {})
        entries.append(TestEntry(name=name, directory=directory, good=good, actual=actual, args=test_args))
    return entries


def pytest_source_candidates(test: TestEntry) -> list[str]:
    if test.actual == "NoSpherA2.log":
        # File is stored under the good-stem name: <good-stem>.log
        good_name = Path(test.good).name
        return [f"{Path(good_name).stem}.log"]
    else:
        # File is stored as-is (actual filename, may not be a .log)
        return [Path(test.actual).name]


def resolve_source(root: Path, failed_logs_dir: Path, test: TestEntry) -> Path | None:
    """Prefer a pytest failed_logs artifact if present, else fall back to the
    ctest in-process runner's actual-output file (written directly next to
    the .good file in tests/<directory>/)."""
    test_failed_dir = failed_logs_dir / test.name
    for candidate in pytest_source_candidates(test):
        candidate_path = test_failed_dir / candidate
        if candidate_path.exists():
            return candidate_path

    ctest_actual = root / "tests" / test.directory / test.actual
    if ctest_actual.exists():
        return ctest_actual

    return None


def ctest_actual_path(root: Path, test: TestEntry) -> Path:
    return root / "tests" / test.directory / test.actual


def build_argv(exe: str, args_table: dict) -> list[str]:
    """Mirrors the argv construction in tests/src/IntegrationTests.cpp
    (run_inprocess_test): -all_charges -no_date defaults, then "-key" per
    [test.args] entry followed by its value token(s). A bare bool (true or
    false) contributes no value token -- only the "-key" flag itself.

    exe may be a multi-token command (e.g. "py -3 fake_exe.py") and is split
    with shlex; use --exe/NOS_EXE with a single path for the normal case."""
    argv = shlex.split(exe, posix=(os.name != "nt")) + ["-all_charges", "-no_date"]
    for key, value in args_table.items():
        argv.append(f"-{key}")
        values = value if isinstance(value, list) else [value]
        for v in values:
            if isinstance(v, bool):
                continue
            argv.append(str(v))
    return argv


def default_executable_candidates(root: Path) -> list[Path]:
    """Known build output locations, most-preferred first: the legacy Visual
    Studio solution build (Windows only) ahead of the CMake preset build
    trees (see CMakePresets.json), release presets ahead of debug ones."""
    exe_name = "NoSpherA2.exe" if os.name == "nt" else "NoSpherA2"
    candidates: list[Path] = []
    if os.name == "nt":
        candidates.append(root / "Windows" / "x64" / "Release" / exe_name)
    for preset in (
        "release-windows", "release-linux", "release-macos-arm64", "release-macos-x86_64",
        "debug-windows", "debug-linux", "debug-macos-arm64", "debug-macos-x86_64",
    ):
        candidates.append(root / "build" / preset / "bin" / exe_name)
    return candidates


def resolve_executable(root: Path, explicit: str | None) -> str:
    """Priority: --exe > $NOS_EXE > first existing known build output
    (VS build preferred over CMake, see default_executable_candidates) >
    bare 'NoSpherA2' resolved via PATH as a last resort."""
    if explicit:
        return explicit
    if os.environ.get("NOS_EXE"):
        return os.environ["NOS_EXE"]
    for candidate in default_executable_candidates(root):
        if candidate.is_file():
            return str(candidate)
    return "NoSpherA2"


def run_test(
    root: Path,
    test: TestEntry,
    exe: str,
    occ_data_path: str | None,
    timeout: float,
) -> tuple[bool, str]:
    test_dir = root / "tests" / test.directory
    if not test_dir.exists():
        return False, f"SKIP {test.name}: test directory missing ({test_dir})"

    argv = build_argv(exe, test.args)

    env = os.environ.copy()
    env.setdefault("OCC_DATA_PATH", occ_data_path or str(root / "Lib" / "occ" / "share" / "occ"))

    try:
        proc = subprocess.run(
            argv,
            cwd=test_dir,
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except FileNotFoundError:
        return False, f"SKIP {test.name}: executable not found ({exe}); pass --exe or set NOS_EXE"
    except subprocess.TimeoutExpired:
        return False, f"SKIP {test.name}: timed out after {timeout}s ({' '.join(argv)})"

    if proc.returncode != 0:
        tail = "\n".join((proc.stdout + proc.stderr).splitlines()[-20:])
        return False, f"SKIP {test.name}: run failed with exit code {proc.returncode} ({' '.join(argv)})\n{tail}"

    return True, f"RAN  {test.name}: {' '.join(argv)}"


def read_lines(path: Path) -> list[str]:
    return path.read_text(errors="replace").splitlines()


def files_equal(a: Path, b: Path) -> bool:
    return read_lines(a) == read_lines(b)


def print_diff(name: str, good: Path, src: Path, max_lines: int = 200) -> None:
    diff_lines = list(
        difflib.unified_diff(
            read_lines(good),
            read_lines(src),
            fromfile=str(good),
            tofile=str(src),
            lineterm="",
        )
    )
    print(f"\n--- Diff for {name} ---")
    if not diff_lines:
        print("(no textual differences)")
    else:
        for line in diff_lines[:max_lines]:
            print(line)
        if len(diff_lines) > max_lines:
            print(f"... ({len(diff_lines) - max_lines} more diff lines omitted)")
    print("--- end diff ---")


def plan_one(
    root: Path,
    failed_logs_dir: Path,
    test: TestEntry,
    force_ctest: bool = False,
) -> tuple[PendingUpdate | None, str]:
    target = root / "tests" / test.directory / test.good

    tests_root = (root / "tests").resolve()
    try:
        target.resolve().relative_to(tests_root)
    except ValueError:
        return None, f"SKIP {test.name}: target escapes tests/ directory ({target})"

    if not target.exists():
        return None, f"SKIP {test.name}: target missing ({target})"

    if force_ctest:
        src = ctest_actual_path(root, test)
        if not src.exists():
            return None, f"SKIP {test.name}: no actual output found after running ({src})"
    else:
        src = resolve_source(root, failed_logs_dir, test)
        if src is None:
            ctest_actual = ctest_actual_path(root, test)
            tried = ", ".join(str(failed_logs_dir / test.name / c) for c in pytest_source_candidates(test))
            return None, f"SKIP {test.name}: no actual output found (tried {tried}; and {ctest_actual})"

    if files_equal(target, src):
        return None, f"OK   {test.name}: already matches ({src})"

    return PendingUpdate(test, target, src), f"DIFF {test.name}: {src} -> {target}"


def apply_update(update: PendingUpdate, backup: bool, dry_run: bool) -> str:
    if dry_run:
        return f"DRYRUN {update.test.name}: {update.source} -> {update.target}"

    if backup:
        backup_path = update.target.with_suffix(update.target.suffix + ".bak")
        shutil.copy2(update.target, backup_path)

    shutil.copy2(update.source, update.target)
    return f"UPDATED {update.test.name}: {update.source} -> {update.target}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Update .good files from pytest failed_logs, ctest actual outputs, or a fresh run."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Repository root (default: script directory)",
    )
    parser.add_argument(
        "--tests-file",
        type=Path,
        default=Path("tests/tests.toml"),
        help="Path to tests TOML relative to --root or absolute",
    )
    parser.add_argument(
        "--failed-logs",
        type=Path,
        default=Path("tests/failed_logs"),
        help="Path to pytest failed logs directory relative to --root or absolute (ignored with --run)",
    )
    parser.add_argument(
        "--test",
        action="append",
        default=[],
        help="Only update the given test name (can be passed multiple times)",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Execute each test itself before diffing, instead of trusting leftover output "
        "from a previous/unrelated run. Only that test's own fresh output is ever used "
        "to update its .good file.",
    )
    parser.add_argument(
        "--exe",
        default=None,
        help="Path to the NoSpherA2 executable for --run (default: $NOS_EXE, else the first existing "
        "known build output -- VS Release build preferred over CMake preset builds -- else "
        "'NoSpherA2' resolved via PATH; see default_executable_candidates())",
    )
    parser.add_argument(
        "--occ-data-path",
        default=None,
        help="OCC_DATA_PATH for --run (default: <root>/Lib/occ/share/occ if not already set)",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=600.0,
        help="Per-test timeout in seconds for --run (default: 600)",
    )
    parser.add_argument(
        "-y", "--yes",
        action="store_true",
        help="Accept all differing tests without prompting for confirmation",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be updated without writing files (implies --yes)",
    )
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Do not create .good.bak backups before overwrite",
    )
    return parser.parse_args()


def prompt_for_confirmation(name: str) -> str:
    """Returns one of 'yes', 'no', 'all', 'quit'."""
    while True:
        answer = input(f"Update {name} good file from the diff above? [y]es/[n]o/[a]ll/[q]uit: ").strip().lower()
        if answer in ("y", "yes"):
            return "yes"
        if answer in ("n", "no", ""):
            return "no"
        if answer in ("a", "all"):
            return "all"
        if answer in ("q", "quit"):
            return "quit"
        print("Please answer y, n, a, or q.")


def main() -> int:
    args = parse_args()

    root = args.root.resolve()
    tests_file = args.tests_file if args.tests_file.is_absolute() else (root / args.tests_file)
    failed_logs = args.failed_logs if args.failed_logs.is_absolute() else (root / args.failed_logs)
    exe = resolve_executable(root, args.exe)
    if args.run:
        print(f"Using executable: {exe}")

    if not tests_file.exists():
        print(f"ERROR: tests file not found: {tests_file}", file=sys.stderr)
        return 2

    entries = load_tests(tests_file)

    requested = set(args.test)
    if requested:
        entries = [entry for entry in entries if entry.name in requested]
        missing = sorted(requested - {entry.name for entry in entries})
        for name in missing:
            print(f"WARN: unknown test name: {name}")

    updated = 0
    skipped = 0
    accept_all = args.yes or args.dry_run
    quit_requested = False

    for entry in entries:
        if quit_requested:
            print(f"SKIP {entry.name}: aborted by user")
            skipped += 1
            continue

        if args.run:
            ran_ok, run_message = run_test(root, entry, exe, args.occ_data_path, args.timeout)
            print(run_message)
            if not ran_ok:
                skipped += 1
                continue

        pending, message = plan_one(root, failed_logs, entry, force_ctest=args.run)
        if pending is None:
            print(message)
            skipped += 1
            continue

        proceed = accept_all
        if not proceed:
            print_diff(entry.name, pending.target, pending.source)
            answer = prompt_for_confirmation(entry.name)
            if answer == "yes":
                proceed = True
            elif answer == "all":
                proceed = True
                accept_all = True
            elif answer == "quit":
                quit_requested = True

        if not proceed:
            print(f"SKIP {entry.name}: rejected by user")
            skipped += 1
            continue

        print(apply_update(pending, backup=not args.no_backup, dry_run=args.dry_run))
        updated += 1

    print(f"SUMMARY: updated={updated}, skipped={skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
