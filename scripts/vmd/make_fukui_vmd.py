#!/usr/bin/env python3
"""Generate a VMD script that visualises NoSpherA2 Fukui-function cubes.

NoSpherA2 run with ``-fukui`` writes four cube files next to the wavefunction::

    <stem>_fukui_plus.cube        f+(r) = |psi_LUMO|^2
    <stem>_fukui_minus.cube       f-(r) = |psi_HOMO|^2
    <stem>_fukui_zero.cube        f0(r) = (f+ + f-) / 2
    <stem>_dual_descriptor.cube   df(r) = f+ - f-

This script writes a ``.vmd`` file that loads them and sets up sensible
representations. VMD needs absolute paths, which is the main reason to generate
the script rather than ship a fixed one.

Reading the pictures
--------------------
f+ is large where the molecule ACCEPTS electron density, i.e. where it is
attacked by a nucleophile - the electrophilic sites. f- is large where it gives
density up, i.e. where it is attacked by an electrophile - the nucleophilic
sites. Note the wording is the opposite way round from what the names suggest,
and it is the single easiest thing to get backwards.

The dual descriptor folds both into one signed field and is usually the more
readable picture: positive (blue by default) marks electrophilic regions,
negative (red) nucleophilic ones.

Usage
-----
    python make_fukui_vmd.py <stem-or-cube-or-directory> [options]

    # after a NoSpherA2 -fukui run in the current directory
    python make_fukui_vmd.py .
    vmd -e fukui.vmd

Options
-------
    -o, --output PATH     Output .vmd path (default: <stem>_fukui.vmd)
    --iso VALUE           Isovalue for f+/f-/f0 surfaces (default: 0.01)
    --dual-iso VALUE      Isovalue magnitude for the dual descriptor (default: 0.01)
    --style STYLE         Atom representation: Licorice, CPK, Lines (default: Licorice)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Suffixes NoSpherA2 writes for -fukui, in the order they are loaded as VMD
# volume IDs 0, 1, 2, 3.
FUKUI_SUFFIXES = [
    ("_fukui_plus.cube", "f+  (LUMO density; nucleophilic attack -> electrophilic site)"),
    ("_fukui_minus.cube", "f-  (HOMO density; electrophilic attack -> nucleophilic site)"),
    ("_fukui_zero.cube", "f0  (radical attack)"),
    ("_dual_descriptor.cube", "df  (dual descriptor; + electrophilic / - nucleophilic)"),
]

# VMD ColorID values.
COLOR_BLUE = 0
COLOR_RED = 1
COLOR_ORANGE = 3
COLOR_GREEN = 7


def resolve_stem(target: Path) -> tuple[Path, str]:
    """Return (directory, stem) for a directory, a stem, or one of the cubes."""
    target = target.expanduser().resolve()

    if target.is_dir():
        matches = sorted(target.glob("*_dual_descriptor.cube"))
        if not matches:
            matches = sorted(target.glob("*_fukui_plus.cube"))
        if not matches:
            raise SystemExit(
                f"No Fukui cube files found in {target}.\n"
                "Run NoSpherA2 with -fukui first, e.g.\n"
                "  NoSpherA2 -wfn mol.wfx -fukui -resolution 0.1 -radius 3.0"
            )
        if len(matches) > 1:
            names = ", ".join(m.name for m in matches)
            raise SystemExit(
                f"Several Fukui cube sets found in {target}: {names}\n"
                "Pass the stem explicitly instead of the directory."
            )
        name = matches[0].name
        for suffix, _ in FUKUI_SUFFIXES:
            if name.endswith(suffix):
                return target, name[: -len(suffix)]

    name = target.name
    for suffix, _ in FUKUI_SUFFIXES:
        if name.endswith(suffix):
            return target.parent, name[: -len(suffix)]

    # Treat it as a bare stem.
    return target.parent, name


def vmd_path(path: Path) -> str:
    """VMD/Tcl wants forward slashes even on Windows."""
    return str(path).replace("\\", "/")


def cube_extremes(path: Path) -> tuple[float, float]:
    """Return (min, max) of the volumetric data in a Gaussian cube file.

    Needed because an isovalue above the data range renders an empty screen,
    which is by far the most common way one of these plots 'does not work'.
    The frontier-orbital densities here peak around 0.01-0.1, and a value
    carried over from a rho or ELI plot will silently show nothing at all.
    """
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    # Line 3 holds the atom count (negative when an MO-index block follows).
    natoms_field = int(lines[2].split()[0])
    natoms = abs(natoms_field)
    header = 6 + natoms
    if natoms_field < 0:
        header += 1  # the extra orbital-index line

    lo = float("inf")
    hi = float("-inf")
    for line in lines[header:]:
        for token in line.split():
            try:
                value = float(token)
            except ValueError:
                continue
            lo = min(lo, value)
            hi = max(hi, value)

    if lo == float("inf"):
        return 0.0, 0.0
    return lo, hi


def build_script(
    directory: Path,
    stem: str,
    iso: float | None,
    dual_iso: float | None,
    style: str,
) -> str:
    present: list[tuple[Path, str, str]] = []
    for suffix, description in FUKUI_SUFFIXES:
        candidate = directory / f"{stem}{suffix}"
        if candidate.exists():
            present.append((candidate, suffix, description))

    if not present:
        raise SystemExit(
            f"None of the expected Fukui cubes exist for stem '{stem}' in {directory}.\n"
            "Expected at least one of: " + ", ".join(s for s, _ in FUKUI_SUFFIXES)
        )

    # Scan the data ranges so the default isovalue actually lands inside them.
    extremes = {suffix: cube_extremes(path) for path, suffix, _ in present}
    print("Cube value ranges:")
    for path, suffix, _ in present:
        lo, hi = extremes[suffix]
        print(f"  {path.name:<34} min {lo:+.6f}   max {hi:+.6f}")

    # A surface at ~30% of the peak is a readable default for these densities.
    #
    # Deliberately the SMALLEST peak across f+/f-/f0, not the largest. f+ and f-
    # routinely differ several-fold in height (here the HOMO peaks at 0.078 and
    # the LUMO at only 0.014, the LUMO being the more diffuse orbital), so 30% of
    # the largest would sit above the f+ maximum entirely and render that surface
    # empty. Keeping ONE shared isovalue also matters scientifically: f+ and f-
    # are only visually comparable when drawn at the same level, so per-cube
    # autoscaling would make the more diffuse orbital look artificially as
    # concentrated as the tighter one.
    if iso is None:
        peaks = [
            extremes[s][1]
            for s in ("_fukui_plus.cube", "_fukui_minus.cube", "_fukui_zero.cube")
            if s in extremes and extremes[s][1] > 0
        ]
        iso = 0.3 * min(peaks) if peaks else 0.01
        print(
            f"Chose isovalue {iso:.6f} for f+/f-/f0 "
            "(30% of the smallest peak, so every surface is visible and all are "
            "drawn at the same level)."
        )
    if dual_iso is None:
        if "_dual_descriptor.cube" in extremes:
            lo, hi = extremes["_dual_descriptor.cube"]
            span = max(abs(lo), abs(hi))
            dual_iso = 0.3 * span if span > 0 else 0.01
        else:
            dual_iso = iso
        print(f"Chose isovalue {dual_iso:.6f} for the dual descriptor.")

    for path, suffix, _ in present:
        lo, hi = extremes[suffix]
        target = dual_iso if suffix == "_dual_descriptor.cube" else iso
        if abs(target) > max(abs(lo), abs(hi)):
            print(
                f"  WARNING: isovalue {target:.6f} lies outside the range of "
                f"{path.name} - that surface will be empty."
            )

    lines: list[str] = []
    add = lines.append

    add("# VMD visualisation of NoSpherA2 Fukui functions and the dual descriptor.")
    add(f"# Generated by make_fukui_vmd.py for stem '{stem}'.")
    add("#")
    add("# Load with:  vmd -e " + f"{stem}_fukui.vmd")
    add("#")
    add("# Volumes loaded into the single molecule, in this order:")
    for index, (_, _, description) in enumerate(present):
        add(f"#   vol {index}: {description}")
    add("")
    add("display projection Orthographic")
    add("display depthcue off")
    add("axes location Off")
    add("color Display Background white")
    # Deliberately no "display rendermode GLSL": it is purely cosmetic and it is
    # an error in text mode (vmd -dispdev text), which is exactly how you would
    # run this script to check it non-interactively.
    add("")

    first_path, _, _ = present[0]
    add("# The cube file carries the atom positions as well as the volumetric data,")
    add("# so no separate coordinate file is needed.")
    add(f'mol new "{vmd_path(first_path)}" type cube first 0 last -1 step 1 waitfor all')
    add("set fukui [molinfo top]")
    for path, _, _ in present[1:]:
        add(
            f'mol addfile "{vmd_path(path)}" type cube first 0 last -1 step 1 '
            "waitfor 1 volsets {0 } $fukui"
        )
    add("")
    add("# Drop the default representation.")
    add("while {[molinfo $fukui get numreps] > 0} {")
    add("    mol delrep 0 $fukui")
    add("}")
    add("")

    add("# --- the molecule itself -------------------------------------------------")
    if style.lower() == "licorice":
        add("mol representation Licorice 0.100000 12.000000 12.000000")
    elif style.lower() == "cpk":
        add("mol representation CPK 1.000000 0.300000 12.000000 12.000000")
    else:
        add("mol representation Lines 1.000000")
    add("mol color Name")
    add("mol selection all")
    add("mol material Opaque")
    add("mol addrep $fukui")
    add("")

    suffix_to_index = {suffix: index for index, (_, suffix, _) in enumerate(present)}

    def isosurface_rep(vol_id: int, value: float, color_id: int, comment: str) -> None:
        add(f"# {comment}")
        # Isosurface <isovalue> <volID> <show> <draw> <step> <size>
        #   show 0 = isosurface only, draw 0 = solid surface
        add(f"mol representation Isosurface {value:.6f} {vol_id} 0 0 1 1")
        add(f"mol color ColorID {color_id}")
        add("mol selection all")
        add("mol material Transparent")
        add("mol addrep $fukui")
        add("")

    if "_fukui_plus.cube" in suffix_to_index:
        isosurface_rep(
            suffix_to_index["_fukui_plus.cube"],
            iso,
            COLOR_BLUE,
            "f+ : where the molecule accepts density (attacked by a nucleophile).",
        )
    if "_fukui_minus.cube" in suffix_to_index:
        isosurface_rep(
            suffix_to_index["_fukui_minus.cube"],
            iso,
            COLOR_RED,
            "f- : where the molecule gives density up (attacked by an electrophile).",
        )
    if "_fukui_zero.cube" in suffix_to_index:
        isosurface_rep(
            suffix_to_index["_fukui_zero.cube"],
            iso,
            COLOR_GREEN,
            "f0 : radical attack. Off by default - enable in the GUI if wanted.",
        )
    if "_dual_descriptor.cube" in suffix_to_index:
        dual_id = suffix_to_index["_dual_descriptor.cube"]
        add("# --- dual descriptor -----------------------------------------------------")
        add("# Signed field, so it needs two surfaces: one at +iso and one at -iso.")
        isosurface_rep(
            dual_id,
            abs(dual_iso),
            COLOR_BLUE,
            "df > 0 : ELECTROPHILIC region (susceptible to nucleophilic attack).",
        )
        isosurface_rep(
            dual_id,
            -abs(dual_iso),
            COLOR_RED,
            "df < 0 : NUCLEOPHILIC region (susceptible to electrophilic attack).",
        )

    # Turn f0 off by default - it is rarely the interesting picture and having
    # every surface visible at once is unreadable.
    if "_fukui_zero.cube" in suffix_to_index:
        f0_rep = 1 + list(suffix_to_index).index("_fukui_zero.cube")
        add(f"# f0 is hidden by default; 'mol showrep $fukui {f0_rep} 1' re-enables it.")
        add(f"mol showrep $fukui {f0_rep} 0")
        add("")

    add("mol rename $fukui \"" + stem + " Fukui\"")
    add("display resetview")
    add("")
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Generate a VMD script for NoSpherA2 Fukui cube files.",
    )
    parser.add_argument(
        "target",
        type=Path,
        help="Directory containing the cubes, the shared stem, or one of the cube files.",
    )
    parser.add_argument("-o", "--output", type=Path, default=None, help="Output .vmd path.")
    parser.add_argument(
        "--iso",
        type=float,
        default=None,
        help="Isovalue for f+/f-/f0. Default: 30%% of the largest peak in the data.",
    )
    parser.add_argument(
        "--dual-iso",
        type=float,
        default=None,
        help="Isovalue magnitude for the dual descriptor. Default: 30%% of its peak.",
    )
    parser.add_argument(
        "--style",
        default="Licorice",
        choices=["Licorice", "CPK", "Lines", "licorice", "cpk", "lines"],
        help="Atom representation style (default: Licorice).",
    )
    args = parser.parse_args(argv)

    directory, stem = resolve_stem(args.target)
    script = build_script(directory, stem, args.iso, args.dual_iso, args.style)

    output = args.output if args.output is not None else directory / f"{stem}_fukui.vmd"
    output.write_text(script, encoding="utf-8")

    print(f"Wrote {output}")
    print(f"Load it with:  vmd -e {vmd_path(output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
