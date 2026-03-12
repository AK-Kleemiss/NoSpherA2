# NoSpherA2 Copilot Working Notes

- Repo: mixed C++ core + Python pytest harness.
- Main test entry in CI: `uv run pytest` from repo root.
- Python test config: `pyproject.toml`, tests under `tests/`.
- Test harness expects executable via `NOS_EXE`; default fallback is `NoSpherA2`.
- Set `OCC_DATA_PATH` to `occ/share` when running tests locally if missing.
- Full/slow tests gated by env var `RUN_FULL_TEST`.

## Build/Test Ground Truth

- Linux/macOS quick build: `make` (produces `NoSpherA2` at repo root).
- Windows quick build: `make.exe` in VS developer shell (produces `NoSpherA2.exe`).
- CMake presets exist for OCC-focused workflows, e.g. `linux-occ-gcc`, `macos-release-*`, `windows-clang-cl`.

## Prompting Rules For This Repo

- Ask Copilot to preserve existing numerics/output formatting in golden-file tests.
- Ask Copilot to prefer minimal diffs and avoid broad refactors in `tests/run_test.py`.
- Ask Copilot to validate env vars (`NOS_EXE`, `OCC_DATA_PATH`, `RUN_FULL_TEST`) before debugging failures.
- Ask Copilot to report exact file/line touched and why behavior is safe.
- For CI-only failures, ask Copilot to compare workflow steps in `.github/workflows/c-cpp_all.yml` with local commands.

## High-Value Prompt Templates

- "Debug this failing pytest by tracing from failing test case to parser/helper and patch minimally; keep existing thresholds and test intent."
- "Before editing, list assumptions about env vars, executable path, and platform skips; then implement and run targeted validation."
- "When fixing golden-output comparison, support scientific notation and multi-value lines; avoid changing unrelated test logic."
- "For GitHub Actions failures, align local reproduction with workflow install/run commands and artifact paths first."

## Known Pitfalls

- Diff parsers that call `float(line)` can fail when one line has multiple values.
- Boolean-like workflow/action inputs are stringly typed; avoid truthiness mistakes in conditions.
- On Windows, architecture/toolchain context matters (x64 dev shell for clang-cl preset flows).
