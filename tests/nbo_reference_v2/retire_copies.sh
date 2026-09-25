#!/bin/bash
# Leave exactly one readable copy of every reference, and make the others say why they are not it.
set -u
V2=/work/akkleemiss/florian/nbo_ref_v2
FIXED=/work/akkleemiss/florian/nbo_ref_v2_fixed
echo "host $(hostname)"

# 1. The archive directory's "native" JSON was computed from ch3.gbw, our own molecule, while
#    ch3_original.nbo is a DIFFERENT molecule with the same name - the pairing that made the old
#    dataset unusable, sitting under a filename that reads as a legitimate pair.  The archive's own
#    wavefunction is gone, so it cannot be regenerated; take the name away from it instead.
cd $V2/ch3_original || exit 1
if [ -f ch3_original.native.nbo.json ]; then
  mv -n ch3_original.native.nbo.json ch3_original.native-from-ch3-gbw.NOT-A-PAIR.json
  echo "renamed the archive's native JSON: $(ls *NOT-A-PAIR* 2>/dev/null)"
fi
cat > README-ch3_original.txt <<'TXT'
ch3_original.nbo is the ONE surviving archive of the lost original run.  Its wavefunction is gone,
which is why this lane exists, and the molecule it describes is not the ch3 regenerated here.

There is therefore no native side for it and there cannot be one.  The file that used to be called
ch3_original.native.nbo.json was produced from ch3.gbw - a different molecule - and is kept only as
evidence, under a name that cannot be mistaken for a pair.  Do not compare it with anything.

ch3_original.gennbo.nbo.json is a current re-parse of the kept .nbo text and carries
parser_version 2.
TXT
echo "wrote $PWD/README-ch3_original.txt"

# 2. nbo_ref_v2_fixed was the correct copy while nbo_ref_v2 held the old parse.  nbo_ref_v2 is now
#    re-parsed in place and verified identical to it field for field, so the mirror has no job left.
#    Its files are deliberately left WITHOUT the parser_version stamp: config.load_nbo refuses them,
#    which is the behaviour wanted from a second copy of a reference.
cat > $FIXED/README-SUPERSEDED.txt <<'TXT'
SUPERSEDED on 25 Sep 2026.  These JSONs were the correct re-parse at a time when
/work/akkleemiss/florian/nbo_ref_v2 still held the old one.  nbo_ref_v2 has since been re-parsed in
place from the same kept .nbo text and verify_one_truth.py confirmed all 23 molecules identical
outside parser_version, source and timings.

Read nbo_ref_v2.  The files here carry no parser_version, so config.load_nbo refuses them by
design: a second copy of a reference is how the wrong one gets measured, and that already happened
once in this lane.
TXT
echo "wrote $FIXED/README-SUPERSEDED.txt"
echo "unstamped files left as fossils: $(grep -L parser_version $FIXED/*/*.json 2>/dev/null | wc -l)"
