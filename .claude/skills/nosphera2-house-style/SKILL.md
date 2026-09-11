---
name: nosphera2-house-style
description: Use before committing any C++ change to NoSpherA2, and when writing new code in Src/core or app/. Carries the house layout conventions (terse, comment-sparse, index-for, OpenMP at column 0, project typedefs) and the pre-commit self-review checklist that catches code which reads foreign to this codebase.
---

# NoSpherA2 house code style

This codebase is **terse, dense and comment-sparse**. Code that is correct but
laid out differently still gets rejected by human reviewers, and has been: the
August 2026 edits were called "nonsense" by readers who had to scroll them.

There is no formatter. No `.clang-format`, no `.editorconfig` outside the
vendored `occ/`. The `cpplint` entry in `.pre-commit-config.yaml` has no config
and applies Google defaults that contradict nearly every convention here —
it is vestigial, not authoritative.

**Style is enforced by imitation. Match the file you are editing, not a global
rule.** Where this skill and the surrounding 50 lines disagree, the file wins.

## Run this before every commit that touches C++

Re-read `git diff --cached` for **layout alone**, ignoring whether the code is
correct. Nine questions:

1. **Blank lines.** Did you insert any inside a function body to "group"
   statements? Core files run 0–9 % blank lines; `fchk.cpp` is 0 %. Remove them.
2. **Comment ratio.** Count comment lines against added lines in the diff.
   Baseline commits here run **0–13 %**. Above ~15 %, justify each block or cut
   it. Rationale belongs in the commit message.
3. **Markdown in comments.** `**bold**`, bullet lists, headings, ASCII tables
   inside `//` or `/* */`. These appear nowhere else in `Src/`. Delete them.
4. **Benchmark numbers in source.** Measurements go in the commit message and
   the vault, never in a comment above a default value or a function.
5. **Defensive prologues.** Any validation added that has not fired on a real
   structure? Delete it. The house guard style is one or two compact statements
   at the top of the function, using `err_checkf`.
6. **Loops.** Index-`for` with `int` and a short name (`i`, `j`, `a`, `s`,
   `x`/`y`/`z`)? Not range-`for` in an index-`for` file, not `std::accumulate`
   over a hot numeric loop.
7. **Braces.** Did you brace a single-statement body in a file that leaves them
   off? 40 % of `if` and 24 % of `for` bodies here are braceless.
8. **Scaffolding.** New RAII guard, `member_` trailing underscores, deleted copy
   constructors, condition variables — in a file that has none? Justify or drop.
9. **Whitespace.** Tabs in a tab file, spaces in a space file. Did the diff
   convert any? `scattering_factors.cpp`, `XCW.cpp`, `basis_set.cpp`,
   `convenience.h`, `XCW.h`, `cell.cpp` are tab-indented.

If a hunk fails any of these, fix it before committing. Do not reformat
untouched code to compensate.

## Conventions

### Density
- Blank lines 0–9 % of a file. No blank-line grouping inside functions.
- Function length: median 24 lines, p90 191. Long functions are accepted.
- Multiple declarations per line are idiomatic:
  `double temp_occ = -1.0, temp_ener = 0.0, last_ener = -DBL_MAX;`
- Lines are short: median 32–46 chars, 1–3 % over 120.

### Comments
- Bare `//` above the line, often with no space after the slashes.
- Doxygen only at API boundaries (`scattering_factors.h`, `cell.h`).
- Explain **mathematics**, not decisions. `nos_math.cpp` is 39 % comments
  because the maths needs it; `properties.cpp` is 1 %.
- **Leave commented-out code in place.** It is not debt here.

### Loops
- `for (int i = 0; i < n; i++)` — 2155 of 2620 loops. `int` even against
  `.size()`. Postfix `i++` dominates 7:1. Zero iterator loops in the codebase.
- STL algorithms are for one-liners only; `std::find_if` appears once in total.

### Types
- **Use the project typedefs** from `convenience.h`: `vec`, `vec2`, `vec3`,
  `ivec*`, `cvec*`, `svec`, `hkl_list`, `dMatrix1..4`. Never spell out
  `std::vector<double>` (392 vs 35 in favour of the alias).
- `auto` sparingly, usually `const auto`. Write the type for ordinary locals.
- `resize` + index-assign over `reserve` + `push_back`.
- Header one-liner accessors carry a redundant trailing `;` after the closing
  brace. That is house style.
- `using namespace std;` at **function** scope, not file scope.

### Parallelism
- **OpenMP only.** No TBB in `Src/core`.
- **`#pragma omp` goes at column 0**, always, regardless of the loop's nesting
  depth. All 201 directives do this.
- Pragma sits immediately above the `for`, no blank line between.

### Errors and logging
- `err_checkf(cond, "msg", file)` — 513 uses, the universal check. Use this
  spelling, not the legacy `err_chkf` / `err_chekf` aliases.
- `err_not_impl_f("msg", file)` for unsupported branches.
- Do not convert `err_checkf` sites to exceptions.
- Log to the injected `std::ostream& file` parameter, not `std::cout`, whenever
  the function has one. Gate verbose output behind `if (debug)`.

### Files
- `#include "pch.h"` is line 1 of every `.cpp` in `Src/core`.
- `#pragma once` in headers; zero include guards.
- Includes: project-local quoted first, rough dependency order, no grouping and
  no alphabetization.
- `app/` is deliberately trivial — `main.cpp` is 9 lines. Logic goes in
  `Src/core`.

## Do not invent uniformity that isn't there

These are genuinely split. Match the file; never "fix" them repo-wide:

- **Braces**: 47/53 K&R vs Allman overall. Decided per file (`XCW.cpp` 99 % K&R,
  `SALTED_equicomb.cpp` 100 % Allman), but a coin flip inside `fchk.cpp`,
  `cube.cpp`, `basis_set.cpp`. Function-definition braces sit on their own line
  even in K&R files.
- **Tabs vs spaces**: six large files are tab-indented.
- **Function naming**: `snake_case` (dominant), `camelCase` (`computeELF`),
  `Calc_Pascal_Snake` (`Calc_Rho`). New free functions → `snake_case`; new
  members of an existing family → follow the family.
- **Members**: no trailing underscore in the old core (`WFN`, `atom`, `cube`,
  `cell`). Confined to seven newer headers. Do not spread it.
- `&` placement, `NULL` vs `nullptr`, `i++` vs `++i` — split by file age, no
  migration in progress.

## Commit messages

They drifted too — ~30 lines each across the August 2026 AI-co-authored commits,
against one-liners before. Keep the subject to one line that names the change.
Put rationale in the body only when a reviewer needs it, and keep it prose, not
a report.

## Background

Full evidence, the measured percentages and the four divergence tells:
`Software-Notes\NoSpherA2-Codebase\NoSpherA2-House-Code-Style-27-Aug-V1.0.md`
in the Obsidian vault.
