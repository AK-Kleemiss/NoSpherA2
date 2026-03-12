# NoSpherA2 Review Prompts (Bug-Risk First)

- Preferred review output order: findings first by severity, then assumptions/questions, then short summary.
- Require concrete file/line evidence for each finding.
- Emphasize behavior regressions over style.

## Prompt Templates

- "Review these changes with a bug-risk focus. List findings by severity with exact file/line references, then open questions, then brief summary."
- "Check for regressions in numerical stability, parser robustness, test determinism, and cross-platform behavior."
- "Identify missing tests for each behavior change and propose the smallest additional test coverage."
- "Ignore pure formatting unless it affects readability or hides logic changes."

## Review Checklist

- Parsing assumptions: single-value vs multi-value lines.
- Floating-point comparison strategy and threshold appropriateness.
- Golden-file comparison behavior and false-positive risk.
- Environment assumptions (`NOS_EXE`, `OCC_DATA_PATH`, optional full-test gates).
- Workflow condition correctness for string boolean inputs.

## Good Follow-up Prompt

- "From your findings, generate a minimal patch plan with ordered steps and validation commands."
