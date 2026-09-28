# renat5 plus open-shell NRT unit correction, 2026-09-28

Source: `nao_renat5/tree` on AKL007, the stored renat5 source with `Src/core/nrt.cpp` and `tests/src/NrtTests.cpp` from local `nbo_external_reference` commit `418a3eb3`. Original remote copies have `.pre_unit_20260928` suffix. This is a focused build, not the exact full local branch.

Binary: `/work/akkleemiss/florian/nao_renat5/tree/build/release-linux/bin/NoSpherA2`, MD5 `eae35c17d924584319e95b2b88228ef9`.

Input: the 22 stored `<mol>_native.47` and `<mol>.gbw` pairs under `/work/akkleemiss/florian/nbo_ref_v2/<mol>/`, with flags from that root's `bin/config.py --native-flags <mol>`. Eight threads, `OMP_NUM_THREADS=8`, `nice -n 15`. All 22 analyses exited 0. Native JSON and logs: `/work/akkleemiss/florian/nrt_renat5_20260928/`.

Reference: gennbo 7 JSON in `tests/nbo_reference_v2/data/<mol>/`. Comparison: `python compare_all.py <paired-root> --only <22 molecule names> --json summary.json`. `comparison.txt` and `summary.json` are the output. The gated view excludes empty Rydberg NBOs and selected E2 rows; the ungated verdict is the verdict.

`NrtTests.*` on the same build: 7/7 passed. The stored renat5 baseline has NPA 63/452, NRT valency 320/476 and NRT bond order 326/702 failed entries; this combined run has 63/452, 284/476 and 305/702, respectively. Both are 0/22 overall. The three additional open-shell references in the 25-molecule local set were not present in this cluster root and were not included.