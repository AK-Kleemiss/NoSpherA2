#!/usr/bin/env python3
"""
NoSpherA2-distro.py

Repackages downloaded/extracted GitHub Actions artifacts into the
hart-XX.zip distribution zips Olex2 fetches on update, merging in the
pTB executable, the occ/share data folders and extra external directories.

Expected layout (defaults, all overridable via CLI flags):

    artifacts/
        NoSpherA2-windows-x64/NoSpherA2-windows-x64_zip/   <- platform build
        NoSpherA2-macos-universal/NoSpherA2-macos-universal_tar/
        NoSpherA2-linux-x86_64/NoSpherA2-linux-x86_64_tar/
    hart-win64/         <- sibling of artifacts/, extra files for windows
    hart-mac64/         <- sibling of artifacts/, extra files for macos
    hart-lin64/         <- sibling of artifacts/, extra files for linux
    basis_sets/         <- sibling of artifacts/
    etc/                <- default: the etc/ next to this script (in the
                           NoSpherA2 repo, so the model is updated with it);
                           models Olex2 reads from its own etc/
    occ/share/          <- sibling of artifacts/ (fallback for --occ-share-dir)
    ptb-cache/<tag>/    <- created here by the pTB release download

The zip root is the Olex2 base directory (index.ind extracts hart-XX.zip
into it), and Olex2 looks executables up in the base directory first, so
NoSpherA2 and ptb both sit at the zip root.

Each output zip is built by flattening these into the zip root:
    - the platform build directory (NoSpherA2[.exe] + runtime libs)
    - the pTB executable of the platform from a GitHub release of
      AK-Kleemiss/ptb (ptb.exe / ptb_linux_static / ptb_macos -> ptb.exe / ptb)
    - the matching extra hart-XX directory
and adding
    - the basis_sets directory as a "basis_sets/" folder
    - the etc directory as "etc/" (merges into the Olex2 etc/ on extraction)
    - the occ/share directory as "occ/share/" (basis/, methods/, solvent/,
      sgdata.json; dftd3/, dftd4/ and xtb/ are left out unless
      --occ-share-full is given - the D4 data is compiled into NoSpherA2).

Any file named NoSpherA2_Tests (with or without an extension, e.g.
NoSpherA2_Tests.exe) found in the platform build or extras directory is
skipped by default - pass --include-tests to keep them.

Resulting zip structure:
    basis_sets/...
    etc/geometry_aid_model.npz
    occ/share/basis/...  occ/share/methods/...  occ/share/solvent/...  occ/share/sgdata.json
    <platform build files>   (NoSpherA2_Tests[.exe] excluded by default)
    ptb[.exe]                (mode 0755 in the unix zips)
    <hart-XX extra files>

pTB runtime libraries: the Windows ptb.exe needs only libiomp5md.dll, which
the NoSpherA2 build already ships; the Linux asset is static. From
v3.10-NoSpherA2 on the macOS asset links only Accelerate and libSystem; the
v3.9 and older ptb_macos needed Homebrew's libgfortran/libgomp/libquadmath -
check an older asset with `otool -L ptb_macos` and put any non-system
libraries into --ptb-libs-dir-mac (they land at the zip root).

Unix modes: the artifact tarballs lose their modes when extracted on Windows,
so files without an extension in the platform build (NoSpherA2) and the pTB
executable are stored as 0755. Olex2 itself chmods NoSpherA2 and ptb at
startup on Linux/macOS (initpy_funcs.py), so this only matters for zips that
are unpacked by hand.

Usage:
    python NoSpherA2-distro.py                       # ptb from the latest release
    python NoSpherA2-distro.py --ptb-tag v3.9-NoSpherA2_with_Lanthanides
    python NoSpherA2-distro.py --ptb-dir ./my-ptb    # local ptb.exe / ptb_linux_static / ptb_macos
    python NoSpherA2-distro.py --no-ptb --occ-share-dir D:/git/NoSpherA2/occ/share
    python NoSpherA2-distro.py --include-tests
    python NoSpherA2-distro.py --artifacts-dir ./artifacts --extras-dir . --out-dir ./dist

GITHUB_TOKEN is optional (the ptb releases are public) but raises the API
rate limit when set.
"""

import argparse
import json
import os
import stat
import sys
import urllib.request
import zipfile

# name -> (platform build dir relative to artifacts-dir, extra dir name relative to extras-dir)
PLATFORMS = {
    "hart-win64.zip": (
        os.path.join("NoSpherA2-windows-x64", "NoSpherA2-windows-x64_zip"),
        "hart-win64",
    ),
    "hart-mac64.zip": (
        os.path.join("NoSpherA2-macos-universal", "NoSpherA2-macos-universal_tar"),
        "hart-mac64",
    ),
    "hart-lin64.zip": (
        os.path.join("NoSpherA2-linux-x86_64", "NoSpherA2-linux-x86_64_tar"),
        "hart-lin64",
    ),
}

PTB_REPO = "AK-Kleemiss/ptb"
# zip name -> (release asset name, name inside the zip)
PTB_ASSETS = {
    "hart-win64.zip": ("ptb.exe", "ptb.exe"),
    "hart-lin64.zip": ("ptb_linux_static", "ptb"),
    "hart-mac64.zip": ("ptb_macos", "ptb"),
}
# zip name -> --ptb-libs-dir-XX option that supplies extra runtime libraries
PTB_LIBS_OPTION = {
    "hart-win64.zip": "ptb_libs_dir_win",
    "hart-lin64.zip": "ptb_libs_dir_lin",
    "hart-mac64.zip": "ptb_libs_dir_mac",
}

# occ/share folders left out unless --occ-share-full: D4 data is compiled in,
# D3 and xtb are not used by NoSpherA2
OCC_SHARE_SKIP = ("dftd3", "dftd4", "xtb")
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# this script lives in <repo>/scripts, the occ submodule and the shipped models beside it
OCC_SHARE_DEFAULT_SOURCE = os.path.join(SCRIPT_DIR, "..", "occ", "share")
ETC_DEFAULT_SOURCE = os.path.join(SCRIPT_DIR, "etc")

TEST_FILE_STEM = "nospherA2_tests".lower()


def _is_test_file(filename: str) -> bool:
    stem = os.path.splitext(filename)[0]
    return stem.lower() == TEST_FILE_STEM


def add_file(zf: zipfile.ZipFile, path: str, arcname: str, executable: bool = False) -> None:
    """Add one file; executable=True stores unix mode 0755 so the file runs after
    extraction on Linux/macOS even when the source (a download on Windows) has no mode."""
    info = zipfile.ZipInfo.from_file(path, arcname)
    info.compress_type = zipfile.ZIP_DEFLATED
    if executable:
        info.external_attr = (stat.S_IFREG | 0o755) << 16
    with open(path, "rb") as fh:
        zf.writestr(info, fh.read())


def add_dir_to_zip(zf: zipfile.ZipFile, src_dir: str, arc_prefix: str = "", skip_tests: bool = False,
                   skip_top_dirs=(), exec_no_ext: bool = False) -> int:
    """Add every file under src_dir into zf. Returns number of files added.
    skip_top_dirs: names of direct sub-directories of src_dir to leave out.
    exec_no_ext: store files without an extension (unix executables) as 0755."""
    count = 0
    for root, dirs, files in os.walk(src_dir):
        if root == src_dir:
            dirs[:] = [d for d in dirs if d not in skip_top_dirs]
        for fname in files:
            if skip_tests and _is_test_file(fname):
                continue
            full_path = os.path.join(root, fname)
            rel_path = os.path.relpath(full_path, src_dir)
            arcname = os.path.join(arc_prefix, rel_path) if arc_prefix else rel_path
            if exec_no_ext and not os.path.splitext(fname)[1]:
                add_file(zf, full_path, arcname, executable=True)
            else:
                zf.write(full_path, arcname)
            count += 1
    return count


def _github_get(url: str, accept: str = "application/vnd.github+json") -> bytes:
    headers = {"Accept": accept, "User-Agent": "NoSpherA2-distro"}
    token = os.environ.get("GITHUB_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req) as resp:
        return resp.read()


def fetch_ptb(tag: str, cache_dir: str) -> str:
    """Download the pTB release assets of PTB_ASSETS into cache_dir/<tag>/ and
    return that directory. Files already present are kept."""
    if tag == "latest":
        url = f"https://api.github.com/repos/{PTB_REPO}/releases/latest"
    else:
        url = f"https://api.github.com/repos/{PTB_REPO}/releases/tags/{tag}"
    print(f"pTB release: {url}")
    release = json.loads(_github_get(url).decode("utf-8"))
    tag = release["tag_name"]
    dest = os.path.join(cache_dir, tag)
    os.makedirs(dest, exist_ok=True)
    print(f"  tag {tag} ({release.get('published_at', '?')}), cache {dest}")
    assets = {a["name"]: a for a in release.get("assets", [])}
    wanted = {asset for asset, _ in PTB_ASSETS.values()}
    for name in sorted(wanted):
        target = os.path.join(dest, name)
        if os.path.isfile(target) and os.path.getsize(target) > 0:
            print(f"  = {name} (cached, {os.path.getsize(target)} bytes)")
            continue
        if name not in assets:
            print(f"  Warning: release {tag} has no asset {name}", file=sys.stderr)
            continue
        print(f"  v {name} ({assets[name]['size']} bytes) ...")
        data = _github_get(assets[name]["url"], accept="application/octet-stream")
        with open(target, "wb") as fh:
            fh.write(data)
    return dest


def build_zip(
    out_path: str,
    basis_sets_dir: str,
    etc_dir: str,
    occ_share_dir: str,
    occ_share_full: bool,
    platform_dir: str,
    extra_dir: str,
    ptb_dir: str,
    ptb_libs_dir: str,
    include_tests: bool,
):
    zip_name = os.path.basename(out_path)
    print(f"Building {out_path} ...")
    manifest = []
    with zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED) as zf:
        if os.path.isdir(basis_sets_dir):
            n = add_dir_to_zip(zf, basis_sets_dir, arc_prefix="basis_sets")
            print(f"  + basis_sets/  ({n} files from {basis_sets_dir})")
            manifest.append(f"basis_sets: {n} files")
        else:
            print(f"  Warning: basis_sets dir not found: {basis_sets_dir}", file=sys.stderr)

        if os.path.isdir(etc_dir):
            n = add_dir_to_zip(zf, etc_dir, arc_prefix="etc")
            print(f"  + etc/  ({n} files from {etc_dir})")
            manifest.append(f"etc: {n} files")
        else:
            print(f"  Warning: etc dir not found (no geometry_aid_model.npz): {etc_dir}", file=sys.stderr)

        if os.path.isdir(occ_share_dir):
            skip = () if occ_share_full else OCC_SHARE_SKIP
            n = add_dir_to_zip(zf, occ_share_dir, arc_prefix=os.path.join("occ", "share"), skip_top_dirs=skip)
            print(f"  + occ/share/  ({n} files from {occ_share_dir}"
                  f"{'' if occ_share_full else ', without ' + '/'.join(OCC_SHARE_SKIP)})")
            manifest.append(f"occ/share: {n} files")
        else:
            print(f"  Warning: occ/share dir not found: {occ_share_dir}", file=sys.stderr)

        if os.path.isdir(platform_dir):
            n = add_dir_to_zip(zf, platform_dir, skip_tests=not include_tests, exec_no_ext=True)
            print(f"  + platform build  ({n} files from {platform_dir})")
            manifest.append(f"build: {n} files")
        else:
            print(f"  Warning: platform build dir not found: {platform_dir}", file=sys.stderr)

        if ptb_dir and zip_name in PTB_ASSETS:
            asset, arcname = PTB_ASSETS[zip_name]
            src = os.path.join(ptb_dir, asset)
            if os.path.isfile(src):
                add_file(zf, src, arcname, executable=not arcname.endswith(".exe"))
                print(f"  + {arcname}  ({os.path.getsize(src)} bytes from {src})")
                manifest.append(f"ptb: {arcname}")
            else:
                print(f"  Warning: pTB asset not found: {src}", file=sys.stderr)
            if ptb_libs_dir:
                if os.path.isdir(ptb_libs_dir):
                    n = add_dir_to_zip(zf, ptb_libs_dir, exec_no_ext=True)
                    print(f"  + pTB libraries  ({n} files from {ptb_libs_dir})")
                    manifest.append(f"ptb libs: {n} files")
                else:
                    print(f"  Warning: pTB libs dir not found: {ptb_libs_dir}", file=sys.stderr)
            elif zip_name == "hart-mac64.zip":
                print("  Note: ptb_macos before v3.10-NoSpherA2 needs Homebrew's libgfortran/libgomp/"
                      "libquadmath - check `otool -L` and pass --ptb-libs-dir-mac for such an asset")

        if os.path.isdir(extra_dir):
            n = add_dir_to_zip(zf, extra_dir, skip_tests=not include_tests)
            print(f"  + extras  ({n} files from {extra_dir})")
            manifest.append(f"extras: {n} files")
        else:
            print(f"  Warning: extras dir not found: {extra_dir}", file=sys.stderr)

        root_files = sorted(i.filename for i in zf.infolist() if "/" not in i.filename)
    print(f"  root: {', '.join(root_files)}")
    print(f"Done: {out_path}  [{'; '.join(manifest)}]")


def main():
    parser = argparse.ArgumentParser(description="Repackage artifacts into distribution zips")
    parser.add_argument("--artifacts-dir", default="./artifacts", help="downloaded/extracted artifacts directory")
    parser.add_argument("--extras-dir", default=None,
                        help="directory containing hart-winXX/hart-macXX/hart-linXX and basis_sets "
                             "(default: parent directory of --artifacts-dir)")
    parser.add_argument("--include-tests", action="store_true",
                        help="keep NoSpherA2_Tests / NoSpherA2_Tests.exe files (skipped by default)")
    parser.add_argument("--out-dir", default="./dist", help="output directory for the distribution zips")
    parser.add_argument("--ptb-tag", default="latest",
                        help=f"release tag of {PTB_REPO} to take the pTB executables from (default: latest)")
    parser.add_argument("--ptb-dir", default=None,
                        help="directory holding ptb.exe / ptb_linux_static / ptb_macos instead of downloading")
    parser.add_argument("--no-ptb", action="store_true", help="do not add pTB")
    parser.add_argument("--ptb-libs-dir-win", default=None, help="extra runtime libraries for ptb.exe (zip root)")
    parser.add_argument("--ptb-libs-dir-lin", default=None, help="extra runtime libraries for the Linux ptb")
    parser.add_argument("--ptb-libs-dir-mac", default=None, help="extra runtime libraries for the macOS ptb")
    parser.add_argument("--etc-dir", default=None,
                        help="files to ship as etc/ (default: the etc/ next to this script)")
    parser.add_argument("--occ-share-dir", default=None,
                        help=f"occ/share source (default: {OCC_SHARE_DEFAULT_SOURCE} if present, "
                             "else <extras-dir>/occ/share)")
    parser.add_argument("--occ-share-full", action="store_true",
                        help=f"also ship {', '.join(OCC_SHARE_SKIP)} from occ/share")
    args = parser.parse_args()

    artifacts_dir = os.path.abspath(args.artifacts_dir)
    extras_dir = os.path.abspath(args.extras_dir) if args.extras_dir else os.path.dirname(artifacts_dir)
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    basis_sets_dir = os.path.join(extras_dir, "basis_sets")
    etc_dir = os.path.abspath(args.etc_dir) if args.etc_dir else ETC_DEFAULT_SOURCE

    if args.occ_share_dir:
        occ_share_dir = os.path.abspath(args.occ_share_dir)
    elif os.path.isdir(OCC_SHARE_DEFAULT_SOURCE):
        occ_share_dir = OCC_SHARE_DEFAULT_SOURCE
    else:
        occ_share_dir = os.path.join(extras_dir, "occ", "share")

    ptb_dir = None
    if not args.no_ptb:
        if args.ptb_dir:
            ptb_dir = os.path.abspath(args.ptb_dir)
        else:
            try:
                ptb_dir = fetch_ptb(args.ptb_tag, os.path.join(extras_dir, "ptb-cache"))
            except Exception as e:  # network or API failure: build without pTB rather than abort
                print(f"Warning: pTB download failed ({e}); building without pTB", file=sys.stderr)

    for zip_name, (platform_rel, extra_name) in PLATFORMS.items():
        platform_dir = os.path.join(artifacts_dir, platform_rel)
        extra_dir = os.path.join(extras_dir, extra_name)
        out_path = os.path.join(out_dir, zip_name)
        ptb_libs_dir = getattr(args, PTB_LIBS_OPTION[zip_name])
        build_zip(out_path, basis_sets_dir, etc_dir, occ_share_dir, args.occ_share_full, platform_dir, extra_dir,
                  ptb_dir, os.path.abspath(ptb_libs_dir) if ptb_libs_dir else None, args.include_tests)

    print(f"\nAll done. Distribution zips are in: {out_dir}")


if __name__ == "__main__":
    main()
