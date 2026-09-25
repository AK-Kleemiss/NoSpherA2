"""Per-molecule NBO settings, and why each one differs from the uniform case.

    python config.py --nbo-keywords <mol>     # the $NBO keylist for gennbo
    python config.py --native-flags <mol>     # the matching -nbo_native flags
    python config.py --nrtstr <mol>           # the $NRTSTR keylist, empty unless so2
    python config.py --list

The uniform keylist is the one the surviving ch3 archive carries:

    $NBO NRT E2PERT NRTLST=0.1 NRTDTL $END

Four molecules deviate, all of them for reasons collect_reference.py's NOTES dict recorded
against the original run, not for reasons discovered here:

  sf6      NRTSYM=off.  With NRT's symmetry detection on, NBO aborts the whole analysis
           after 0.08 s ("generated 48 symmetry operator(s) for Th but expected 24") and
           prints no NAO, NBO, E2 or NRT at all.
  ni_co_4  NRTE2=5.  At the default delocalisation threshold the search does not finish.
  ticl4    NRTE2=20.  Same.
  so2      a hand-written $NRTSTR.  NBO localises SO2 as S-O single plus S=O triple and NRT
           then rejects both symmetry-related parents ("Unable to find reasonable resonance
           structure").  The four structures supplied below are the ones the original archive
           carried: O=S=O, the two O(-)-S(+)=O forms and the doubly ionic O(-)-S(2+)-O(-).
           Each is 9 valence pairs, which is what NBO checks ("Structure n of the $NRTSTR
           keylist has 1 too many electron pair(s)"), so a wrong one fails loudly.
           Atom numbering is the one make_inputs.py writes: 1 = S, 2 and 3 = O.

The native flags are chosen to put the in-house run on the SAME thresholds NBO reported for
the original run, because below a printing threshold a missing entry and a small one look
alike and the comparison would then be measuring the threshold instead of the method:

  -nbo_e2min   NBO's E2 printing threshold, 0.5 kcal closed shell and 0.25 open shell
               (index.json thresholds_reported_by_nbo.e2_kcal).  NoSpherA2's default is 0.5
               for both, so ch3, no and o2 need it passed.
  -nrt_e2      NBO's NRT delocalisation threshold, 1 kcal for the main-group cases and the
               5 / 20 above for the two transition metals.  NoSpherA2's default is 2.0, so
               every molecule needs it passed.
  -nrt_no_symmetry  for sf6 only, mirroring NRTSYM=off.
"""
import argparse
import json
import os
import sys

BASE_KEYWORDS = "NRT E2PERT NRTLST=0.1 NRTDTL"

EXTRA_KEYWORDS = {
    "sf6": "NRTSYM=off",
    "ni_co_4": "NRTE2=5",
    "ticl4": "NRTE2=20",
}

# NBO's own reported NRT delocalisation threshold per molecule, read off index.json at
# import time rather than duplicated here.
HERE = os.path.dirname(os.path.abspath(__file__))
INDEX = os.path.join(HERE, os.pardir, "nbo_reference", "index.json")

NRTSTR_SO2 = """ $NRTSTR
  STR              ! O=S=O
   LONE 1 1 2 2 3 2 END
   BOND D 1 2 D 1 3 END
  END
  STR              ! O(-)-S(+)=O
   LONE 1 1 2 3 3 2 END
   BOND S 1 2 D 1 3 END
  END
  STR              ! O=S(+)-O(-)
   LONE 1 1 2 2 3 3 END
   BOND D 1 2 S 1 3 END
  END
  STR              ! O(-)-S(2+)-O(-)
   LONE 1 1 2 3 3 3 END
   BOND S 1 2 S 1 3 END
  END
 $END
"""


def thresholds():
    idx = json.load(open(INDEX))
    return {m["name"]: (m["thresholds_reported_by_nbo"], m["orca"]["multiplicity"])
            for m in idx["molecules"]}


def nbo_keywords(mol):
    kw = BASE_KEYWORDS
    if mol in EXTRA_KEYWORDS:
        kw += " " + EXTRA_KEYWORDS[mol]
    return kw


def native_flags(mol):
    th, mult = thresholds()[mol]
    flags = ["-nrt",
             "-nbo_e2min", "%g" % th["e2_kcal"],
             "-nrt_e2", "%g" % th["nrt_deloc_kcal"]]
    if mol == "sf6":
        flags.append("-nrt_no_symmetry")
    return " ".join(flags)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nbo-keywords")
    ap.add_argument("--native-flags")
    ap.add_argument("--nrtstr")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.nbo_keywords:
        print(nbo_keywords(a.nbo_keywords))
    elif a.native_flags:
        print(native_flags(a.native_flags))
    elif a.nrtstr:
        sys.stdout.write(NRTSTR_SO2 if a.nrtstr == "so2" else "")
    elif a.list:
        for m in sorted(thresholds()):
            print("%-13s | %-40s | %s" % (m, nbo_keywords(m), native_flags(m)))
    else:
        ap.error("pick one")


if __name__ == "__main__":
    main()
