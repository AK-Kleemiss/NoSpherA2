# NoSpherA2 CI Triage Prompts

- CI baseline for Python tests: install with `uv sync --locked --all-extras --dev`, then run `uv run pytest`.
- CI test harness relies on `NOS_EXE` and `OCC_DATA_PATH`.
- Optional heavy tests are gated by `RUN_FULL_TEST`.

## Prompt Templates

- "Reproduce this CI failure by matching workflow steps exactly, including env vars and executable path assumptions, before proposing code changes."
- "Diff local vs CI execution for install, binary location, env vars, and platform skips. Then propose smallest reliable fix."
- "Focus only on this failing matrix job; identify whether root cause is workflow config, environment, artifact path, or test logic."
- "When CI-only failure appears, inspect `.github/workflows/c-cpp_all.yml` and the corresponding test helper path first."

## CI Debug Checklist

- Is `NOS_EXE` set to the built/downloaded binary path?
- Is `OCC_DATA_PATH` valid for the job workspace?
- Is executable permission needed on Linux/macOS?
- Are platform skip conditions expected (e.g., macOS skip flags)?
- Are diff/parsing helpers robust to current output format?

## Good Follow-up Prompt

- "Propose one immediate hotfix and one hardening improvement (tests/workflow guard) to prevent recurrence."
