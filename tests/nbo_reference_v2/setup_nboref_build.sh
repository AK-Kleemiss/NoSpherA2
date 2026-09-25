#!/bin/bash
# Set up MY OWN NoSpherA2 checkout on the cluster at the branch that carries the parser fix,
# without touching any sibling session's tree. Source is copied from nos_base (read-only use,
# for the five submodules and to skip a fresh clone); build/ is deliberately NOT copied, because
# a CMake cache bakes absolute paths to the tree it was configured in.
set -eu
SRC=/work/akkleemiss/florian/nos_base
DST=/work/akkleemiss/florian/nos_nboref
BRANCH=nbo_external_reference

if [ -d "$DST/.git" ]; then
    echo "already present, updating"
else
    rsync -a --exclude build/ --exclude .git/modules/*/objects/pack/tmp* "$SRC/" "$DST/"
fi

cd "$DST"
git fetch origin "$BRANCH"
git checkout -B "$BRANCH" FETCH_HEAD
git submodule update --init --recursive
echo "=== $DST now at:"
git log --oneline -3
echo "=== the fix is present if the next line prints valency_spin:"
grep -n "valency_spin" Src/core/nbo_run.cpp | head -4
echo "=== host $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores"
