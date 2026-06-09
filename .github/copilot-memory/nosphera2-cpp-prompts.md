# NoSpherA2 C++ Change Prompts

- Default request style: minimal diff, no broad refactor, preserve CLI/output semantics.
- Prefer touching one subsystem at a time (`Src/`, `tests/`, then build glue).
- Ask for explicit risk check on numerics, parsing, and platform-specific behavior.

## Prompt Templates

- "Implement a minimal C++ fix in the failing code path only. Do not reformat unrelated files. Keep public behavior and CLI flags unchanged."
- "Before editing, list assumptions about input format, numeric precision, and expected output files; then patch only the parser/logic that fails."
- "For this C++ change, show impacted files/functions, then compile or run targeted tests only for the affected area."
- "If changing output formatting or parsing, preserve scientific notation compatibility and multi-value line handling."

## Safety Checklist For Copilot

- Confirm where behavior is defined (source + tests).
- Keep thresholds/tolerances unless test intent requires a change.
- Avoid altering existing `.good` files unless explicitly requested.
- Call out cross-platform implications (Windows/Linux/macOS).

## Good Follow-up Prompt

- "Now review your own patch for regressions in file parsing, floating-point comparisons, and backward compatibility with existing test data."
