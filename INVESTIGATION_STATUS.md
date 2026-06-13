# Investigation Status: Resolved

Last updated: 2026-06-13

The previous `alanine_integrated_occ`, Visual Studio Debug, MSVC dependency-build,
and fchk conversion investigations are resolved on the current branch.

## Current Verified Results

- Python pytest Release full: `21 passed` in 81.88s, with one tolerated
  `fchk_conversion` numeric warning at line `261148`
- Python pytest Debug full: `21 passed` in 1103.46s
- Visual Studio native Debug full: `21 passed` in 24.4986 minutes
- Visual Studio native Release full: `21 passed` in 1.3985 minutes

`Release+Copy` was intentionally not rerun in the latest validation pass.

## Resolved Items

- `alanine_integrated_occ` heap corruption is fixed in both Release and Debug.
- Visual Studio tests run in-process through `NoSpherA2_DLL.dll`; the subprocess path
  has been removed and must not be reintroduced.
- oneTBB is built and deployed as the shared core runtime. `tbbmalloc_proxy` and
  `tbbmalloc` are not linked.
- MSVC Debug dependency builds have a matching `windows-msvc-debug` workflow preset.
- `fchk_conversion` now generates and compares `log.fchk`, and the fchk writer handles
  restricted/beta output, coefficient bounds, virtual orbitals, and G shells correctly.

Use `HANDOFF.md` for the detailed current handoff, changed files, rules, and
validation commands.
