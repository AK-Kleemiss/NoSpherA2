#!/usr/bin/env python3
"""Generate a VMD script that shows NoSpherA2's spin-resolved ELI-D cubes.

NoSpherA2 run with ``-eli_spin`` writes, next to the wavefunction::

    <stem>_eli_aa.cube        ELI-D(aa), alpha-alpha (same-spin) pairs
    <stem>_eli_bb.cube        ELI-D(bb), beta-beta pairs          (open shell only)
    <stem>_eli_triplet.cube   ELI-D(t),  triplet-coupled pairs    (open shell only)

and, with ``-eli`` / ``-s_rho`` / ``-rho`` on the same command, the spin-summed
``<stem>_eli.cube``, ``<stem>_s_rho.cube`` and ``<stem>_rho.cube`` on the same grid.

This script writes ``<stem>_eli_spin.vmd``: the molecule is loaded once, every
member is a volume on it, and each member gets its own isosurface
representation (aa blue, bb red, triplet green) so it can be toggled in
Graphics > Representations. ``--side-by-side`` instead loads one copy of the
molecule per member, shifted along x, one member each.

Usage
-----
    python make_eli_spin_vmd.py <directory | stem | one of the cubes> [options]
    vmd -e <stem>_eli_spin.vmd

Options
-------
    -o, --output PATH     output .vmd (default <stem>_eli_spin.vmd)
    --iso V               ELI-D level for every member (default 1.6)
    --iso-aa/--iso-bb/--iso-triplet V   per-member levels
    --rho-clip [V]        colour the surfaces by the density cube (needs <stem>_rho.cube, i.e. -rho);
                          with V the colour scale starts at rho = V, so lower-density parts fade out
    --side-by-side        three molecule copies, translated, one member each
    --shift A             translation between the copies (default: 1.5 x the molecule's width)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# (suffix, key, VMD ColorID, description); the order is the volume order on the molecule
MEMBERS = [
    ("_eli_aa.cube", "aa", 0, "ELI-D(aa)  alpha-alpha pairs"),
    ("_eli_bb.cube", "bb", 1, "ELI-D(bb)  beta-beta pairs"),
    ("_eli_triplet.cube", "triplet", 7, "ELI-D(t)   triplet-coupled pairs"),
]
RHO_SUFFIX = "_rho.cube"
DEFAULT_ISO = 1.6


def resolve_stem(target: Path) -> tuple[Path, str]:
    """(directory, stem) from a directory, a stem or one of the cubes."""
    target = target.expanduser().resolve()
    if target.is_dir():
        matches = sorted(target.glob("*_eli_aa.cube"))
        if not matches:
            raise SystemExit(
                f"No spin-resolved ELI-D cubes in {target}.\n"
                "Run NoSpherA2 with -eli_spin first, e.g.\n"
                "  NoSpherA2 -wfn mol.fchk -eli_spin -resolution 0.1 -radius 3.5"
            )
        if len(matches) > 1:
            raise SystemExit(
                f"Several ELI-D cube sets in {target}: " + ", ".join(m.name for m in matches)
                + "\nPass the stem instead of the directory."
            )
        return target, matches[0].name[: -len("_eli_aa.cube")]
    for suffix, *_ in MEMBERS:
        if target.name.endswith(suffix):
            return target.parent, target.name[: -len(suffix)]
    return target.parent, target.name


def vmd_path(path: Path) -> str:
    return str(path).replace("\\", "/")


def cube_extent(path: Path) -> tuple[float, float, float]:
    """(min, max, x width in Angstrom of the atoms) of a cube; the width sets the side-by-side shift."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    nat_field = int(lines[2].split()[0])
    nat = abs(nat_field)
    header = 6 + nat + (1 if nat_field < 0 else 0)
    xs = [float(lines[6 + a].split()[2]) * 0.529177 for a in range(nat)]
    lo, hi = float("inf"), float("-inf")
    for line in lines[header:]:
        for token in line.split():
            try:
                v = float(token)
            except ValueError:
                continue
            lo, hi = min(lo, v), max(hi, v)
    width = (max(xs) - min(xs)) if xs else 0.0
    return (lo if lo != float("inf") else 0.0), (hi if hi != float("-inf") else 0.0), width


def build_script(directory: Path, stem: str, iso: dict[str, float], rho_clip: float | None,
                 use_rho: bool, side_by_side: bool, shift: float | None) -> str:
    present = [(directory / f"{stem}{s}", s, k, c, d) for s, k, c, d in MEMBERS if (directory / f"{stem}{s}").exists()]
    if not present:
        raise SystemExit(f"No {stem}_eli_*.cube in {directory}.")
    rho = directory / f"{stem}{RHO_SUFFIX}"
    use_rho = use_rho and rho.exists()
    if (rho_clip is not None or use_rho) and not rho.exists():
        print(f"No {rho.name}: --rho-clip ignored (add -rho to the NoSpherA2 run).")
    width = 0.0
    print("Cube value ranges:")
    for path, s, k, *_ in present:
        lo, hi, w = cube_extent(path)
        width = max(width, w)
        print(f"  {path.name:<34} min {lo:.4f}  max {hi:.4f}")
        if iso[k] > hi:
            print(f"  WARNING: level {iso[k]} lies above the maximum of {path.name}; that surface will be empty.")
    if len(present) == 1:
        print("Only ELI-D(aa): a closed-shell wavefunction, where bb equals aa and the triplet member is a multiple of it.")
    step = shift if shift is not None else max(1.5 * width, 4.0)

    out: list[str] = []
    add = out.append
    add("# VMD scene of NoSpherA2 spin-resolved ELI-D, generated by make_eli_spin_vmd.py")
    add(f"# Load with:  vmd -e {stem}_eli_spin.vmd")
    for i, (_, _, k, _, d) in enumerate(present):
        add(f"#   {d}: level {iso[k]}")
    add("display projection Orthographic")
    add("display depthcue off")
    add("axes location Off")
    add("color Display Background white")
    add("")

    def load(tag: str, members: list) -> None:
        first = members[0][0]
        add(f'mol new "{vmd_path(first)}" type cube first 0 last -1 step 1 waitfor all')
        add(f"set {tag} [molinfo top]")
        for path, *_ in members[1:]:
            add(f'mol addfile "{vmd_path(path)}" type cube first 0 last -1 step 1 waitfor 1 volsets {{0 }} ${tag}')
        if use_rho:
            add(f'mol addfile "{vmd_path(rho)}" type cube first 0 last -1 step 1 waitfor 1 volsets {{0 }} ${tag}')
        add(f"while {{[molinfo ${tag} get numreps] > 0}} {{ mol delrep 0 ${tag} }}")
        add("mol representation Licorice 0.100000 12.000000 12.000000")
        add("mol color Name")
        add("mol selection all")
        add("mol material Opaque")
        add(f"mol addrep ${tag}")

    def surface(tag: str, vol: int, rho_vol: int, member) -> None:
        _, _, k, colour, d = member
        add(f"# {d}")
        add(f"mol representation Isosurface {iso[k]:.6f} {vol} 0 0 1 1")
        add(f"mol color Volume {rho_vol}" if use_rho else f"mol color ColorID {colour}")
        add("mol selection all")
        add("mol material Transparent")
        add(f"mol addrep ${tag}")
        if use_rho and rho_clip is not None:
            add(f"mol scaleminmax ${tag} [expr {{[molinfo ${tag} get numreps] - 1}}] {rho_clip:.6f} [lindex [voltool info minmax -mol ${tag} -vol {rho_vol}] 1]")

    if not side_by_side:
        load("eli", present)
        for vol, m in enumerate(present):
            surface("eli", vol, len(present), m)
        add(f'mol rename $eli "{stem} ELI-D spin"')
    else:
        for i, m in enumerate(present):
            tag = f"eli_{m[2]}"
            load(tag, [m])
            surface(tag, 0, 1, m)
            add(f'mol rename ${tag} "{stem} {m[4].split()[0]}"')
    add("display resetview")
    if side_by_side:
        #after resetview, which puts every molecule on one centre: the copies are shifted through their centre matrix,
        #which moves atoms and volumes together (moveby would leave the isosurface behind)
        for i, m in enumerate(present):
            tag = f"eli_{m[2]}"
            add(f"molinfo ${tag} set center_matrix [list [transmult [lindex [molinfo ${tag} get center_matrix] 0] [transoffset {{{i * step:.3f} 0 0}}]]]")
    add("")
    return "\n".join(out)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="VMD script for NoSpherA2 -eli_spin cubes.")
    ap.add_argument("target", type=Path, help="directory, shared stem, or one of the cubes")
    ap.add_argument("-o", "--output", type=Path)
    ap.add_argument("--iso", type=float, default=DEFAULT_ISO, help="ELI-D level for every member (default 1.6)")
    ap.add_argument("--iso-aa", type=float)
    ap.add_argument("--iso-bb", type=float)
    ap.add_argument("--iso-triplet", type=float)
    ap.add_argument("--rho-clip", type=float, nargs="?", const=-1.0, default=None,
                    help="colour by the density cube; with a value, start the colour scale there")
    ap.add_argument("--side-by-side", action="store_true")
    ap.add_argument("--shift", type=float)
    a = ap.parse_args(argv)
    iso = {"aa": a.iso_aa if a.iso_aa is not None else a.iso,
           "bb": a.iso_bb if a.iso_bb is not None else a.iso,
           "triplet": a.iso_triplet if a.iso_triplet is not None else a.iso}
    directory, stem = resolve_stem(a.target)
    clip = None if a.rho_clip is None or a.rho_clip < 0 else a.rho_clip
    script = build_script(directory, stem, iso, clip, a.rho_clip is not None, a.side_by_side, a.shift)
    output = a.output if a.output is not None else directory / f"{stem}_eli_spin.vmd"
    output.write_text(script, encoding="utf-8")
    print(f"Wrote {output}")
    print(f"Load it with:  vmd -e {vmd_path(output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
