#!/usr/bin/env python3
"""Update test .good baselines from pytest failed_logs outputs.

This script reads tests/tests.toml to map test names to their expected baseline files,
then copies matching *.log files from tests/failed_logs/<test_name>/ into the
corresponding tests/<directory>/*.good files.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from dataclasses import dataclass
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


def load_tests(toml_path: Path) -> list[TestEntry]:
    with toml_path.open("rb") as fd:
        data = tomllib.load(fd)

    data.pop("defaults", None)
    entries: list[TestEntry] = []
    for name, entry in data.items():
        directory = entry.get("directory", name)
        good = entry.get("good", f"{name}.good")
        actual = entry.get("actual", "NoSpherA2.log")
        entries.append(TestEntry(name=name, directory=directory, good=good, actual=actual))
    return entries


def source_candidates(test: TestEntry) -> list[str]:
    if test.actual == "NoSpherA2.log":
        # File is stored under the good-stem name: <good-stem>.log
        good_name = Path(test.good).name
        return [f"{Path(good_name).stem}.log"]
    else:
        # File is stored as-is (actual filename, may not be a .log)
        return [Path(test.actual).name]


def update_one(
    root: Path,
    failed_logs_dir: Path,
    test: TestEntry,
    dry_run: bool,
    backup: bool,
) -> tuple[bool, str]:
    target = root / "tests" / test.directory / test.good

    if target.suffix != ".good":
        return False, f"SKIP {test.name}: target is not a .good file ({target})"

    if not target.exists():
        return False, f"SKIP {test.name}: target missing ({target})"

    test_failed_dir = failed_logs_dir / test.name
    if not test_failed_dir.exists():
        return False, f"SKIP {test.name}: no failed logs directory ({test_failed_dir})"

    src = None
    for candidate in source_candidates(test):
        candidate_path = test_failed_dir / candidate
        if candidate_path.exists():
            src = candidate_path
            break

    if src is None:
        tried = ", ".join(source_candidates(test))
        return False, f"SKIP {test.name}: no matching source log in {test_failed_dir} (tried: {tried})"

    if dry_run:
        return True, f"DRYRUN {test.name}: {src} -> {target}"

    if backup:
        backup_path = target.with_suffix(target.suffix + ".bak")
        shutil.copy2(target, backup_path)

    shutil.copy2(src, target)
    return True, f"UPDATED {test.name}: {src} -> {target}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Update .good files from pytest failed_logs outputs."
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
        help="Path to failed logs directory relative to --root or absolute",
    )
    parser.add_argument(
        "--test",
        action="append",
        default=[],
        help="Only update the given test name (can be passed multiple times)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be updated without writing files",
    )
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Do not create .good.bak backups before overwrite",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    root = args.root.resolve()
    tests_file = args.tests_file if args.tests_file.is_absolute() else (root / args.tests_file)
    failed_logs = args.failed_logs if args.failed_logs.is_absolute() else (root / args.failed_logs)

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
    for entry in entries:
        ok, message = update_one(
            root=root,
            failed_logs_dir=failed_logs,
            test=entry,
            dry_run=args.dry_run,
            backup=not args.no_backup,
        )
        print(message)
        if ok:
            updated += 1
        else:
            skipped += 1

    print(f"SUMMARY: updated={updated}, skipped={skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
