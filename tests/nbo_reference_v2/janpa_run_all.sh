#!/bin/bash
# JANPA 2.02 (the Java original, a THIRD CODE) on the eight molecules that carry gennbo unit
# 32/33 dumps, with a matrix exported after EVERY stage of the published cascade.
#
# The read: orca_2mkl's molden has ORCA's contraction coefficients, which are not MOLDEN-
# normalised - raw, JANPA reports a normalisation problem, non-orthogonal MOs and "input data
# seems to be improper", and returns Li +15.5 e on a 12-electron molecule.  `-fromorca3bf` is
# molden2molden's documented inverse of exactly that convention and takes the warning count to 0.
# That is the arm used here; the read is then ASSERTED downstream against the .47's own $OVERLAP
# and against gennbo's printed NPA charges, not against the absence of warnings.
# janpa_java/ (janpa.jar md5 fc7b3805c866117ea2ee51e5f7d00fd0, molden2molden.jar md5
# 622ffbaf67aa6f2732503dfa52fc851b) and janpa_work/<mol>/ are expected BESIDE this script; the
# jars are not redistributed here.  janpa_work/<mol>/ needs <mol>.molden.input and <mol>.47.
set -u
JAVA="/c/Program Files/Angry IP Scanner/jre/bin/java.exe"
HERE=$(cd "$(dirname "$0")" && pwd)
JAR=$HERE/janpa_java
W=$HERE/janpa_work
MOLS=${MOLS:-"lif water ammonia ethane benzene pf5 so2 sf6"}

echo "host=$(hostname) date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
"$JAVA" -version 2>&1 | head -1
for f in "$JAR/janpa.jar" "$JAR/molden2molden.jar" "/c/Program Files/Angry IP Scanner/jre/bin/java.exe"; do
    printf 'BIN %-20s %9s  %s\n' "$(basename "$f")" "$(stat -c %s "$f")" "$(md5sum "$f" | cut -c1-32)"
done

for MOL in $MOLS; do
    d=$W/$MOL
    [ -d "$d" ] || { echo "$MOL: no work dir"; continue; }
    cd "$d"
    printf 'IN  %-10s molden %9s %s   .47 %9s %s\n' "$MOL" \
        "$(stat -c %s $MOL.molden.input)" "$(md5sum $MOL.molden.input | cut -c1-32)" \
        "$(stat -c %s $MOL.47)" "$(md5sum $MOL.47 | cut -c1-32)"
    "$JAVA" -jar $JAR/molden2molden.jar -i $MOL.molden.input -o $MOL.fob.molden \
        -fromorca3bf -orca3signs > m2m.log 2>&1
    T0=$(date +%s)
    "$JAVA" -Xmx4g -jar $JAR/janpa.jar -i $MOL.fob.molden \
        -MatrixFloatNumberFormat "%.16e" \
        -npacharges     $MOL.j.npa.txt \
        -S_Matrix_File  $MOL.j.S.txt \
        -D_Matrix_File  $MOL.j.D.txt \
        -PNAO2AO_File   $MOL.j.pnao2ao.txt \
        -NAO2AO_File    $MOL.j.nao2ao.txt \
        -SDS_NAO_File   $MOL.j.sds_nao.txt \
        -PNAO_OverlapMatrix_File      $MOL.j.s1_pnao_ovl.txt \
        -PNAO_SDS_Matrix_File         $MOL.j.s1_pnao_sds.txt \
        -NMB_old_Overlap_Matrix_File  $MOL.j.s3_nmb_ovl.txt \
        -NMB_old_SDS_Matrix_File      $MOL.j.s3_nmb_sds.txt \
        -NRB_old_Overlap_Matrix_File  $MOL.j.s4_nrb_old_ovl.txt \
        -NRB_new_Overlap_Matrix_File  $MOL.j.s4_nrb_new_ovl.txt \
        -S_Matrix_after_ON2_File      $MOL.j.s5_ovl.txt \
        -SDS_Matrix_after_ON2_File    $MOL.j.s5_sds.txt \
        -NRB_Overlap_after_Schmidt2_File   $MOL.j.s4b_nrb_schmidt.txt \
        -NRB_Overlap_after_OW_heavy_File   $MOL.j.s6_nrb_owheavy.txt \
        -NRB_Overlap_after_OW2_final_File  $MOL.j.s6_nrb_ow2.txt \
        -OW2_File                     $MOL.j.s6_ow2.txt \
        -S_Matrix_after_OW2_File      $MOL.j.s6_ovl.txt \
        -SDS_Matrix_after_OW2_File    $MOL.j.s6_sds.txt \
        -verboseprint > janpa.log 2>&1
    RC=$?; T1=$(date +%s)
    printf 'RUN %-10s rc=%d %4ds warnings=%s  npa=%s\n' "$MOL" "$RC" "$((T1-T0))" \
        "$(grep -c WARNING janpa.log)" "$(tr '\n' ' ' < $MOL.j.npa.txt | cut -c1-60)"
    grep -i WARNING janpa.log | sed 's/^/      /' | head -4
done
echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
