"""scripts/vmd/make_eli_spin_vmd.py on a synthetic -eli_spin cube set: the .vmd loads the molecule once with three volumes
(one mol new, two mol addfile) and draws three isosurface representations; a closed-shell set (aa only) draws one.
No VMD needed. Run with pytest, or as `python tests/test_make_eli_spin_vmd.py`."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts" / "vmd"))
import make_eli_spin_vmd  # noqa: E402


def write_cube(path: Path, peak: float) -> None:
    lines = ["test cube", "", "    2    0.000000    0.000000    0.000000",
             "    2    1.000000    0.000000    0.000000", "    2    0.000000    1.000000    0.000000",
             "    2    0.000000    0.000000    1.000000",
             "    6    6.000000    0.000000    0.000000    0.000000", "    1    1.000000    2.000000    0.000000    0.000000"]
    lines += [f" {0.1:.5E} {peak:.5E} {0.2:.5E} {0.3:.5E} {0.4:.5E} {0.5:.5E}", f" {0.6:.5E} {0.7:.5E}"]
    path.write_text("\n".join(lines) + "\n")


def make_set(d: Path, members) -> None:
    for suffix in members:
        write_cube(d / f"mol{suffix}", 2.0)


def test_three_members_three_volumes_three_surfaces():
    with tempfile.TemporaryDirectory() as t:
        d = Path(t)
        make_set(d, ["_eli_aa.cube", "_eli_bb.cube", "_eli_triplet.cube"])
        assert make_eli_spin_vmd.main([str(d)]) == 0
        text = (d / "mol_eli_spin.vmd").read_text()
        assert text.count("mol new ") == 1
        assert text.count("mol addfile ") == 2
        assert text.count("mol representation Isosurface 1.600000") == 3
        for vol in range(3):
            assert f"Isosurface 1.600000 {vol} 0 0 1 1" in text
        assert make_eli_spin_vmd.main([str(d), "--iso-bb", "1.2", "--side-by-side", "-o", str(d / "s.vmd")]) == 0
        side = (d / "s.vmd").read_text()
        assert side.count("mol new ") == 3 and side.count("Isosurface 1.200000 0") == 1
        assert side.count("transoffset") == 3


def test_closed_shell_set_draws_aa_only():
    with tempfile.TemporaryDirectory() as t:
        d = Path(t)
        make_set(d, ["_eli_aa.cube"])
        assert make_eli_spin_vmd.main([str(d / "mol_eli_aa.cube")]) == 0
        text = (d / "mol_eli_spin.vmd").read_text()
        assert text.count("mol addfile ") == 0 and text.count("Isosurface") == 1


if __name__ == "__main__":
    test_three_members_three_volumes_three_surfaces()
    test_closed_shell_set_draws_aa_only()
    print("make_eli_spin_vmd checks passed")
