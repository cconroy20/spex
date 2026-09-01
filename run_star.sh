#!/bin/bash
# Every SYNTHE run one star needs:
#   eos/      a 0.2 nm run with more_output=yes, purely for the .mol number
#             densities the line-index build reads
#   chunks/   the full 355-1000 nm spectrum, in 25 nm pieces (memory)
#   species/  one run per species, whole band, each with a private lines.list
#
#   ./run_star.sh <tag> <model.atm>
#
# Re-running is cheap: anything with a non-empty .spec already on disk is
# skipped, so adding species to run/species.list only costs the new ones.
set -e
A=/Users/cconroy/kurucz/atlas12
BASE=$(cd "$(dirname "$0")" && pwd)
TAG=$1; MODEL=$2
[ -f "$MODEL" ] || { echo "no such model: $MODEL"; exit 1; }
# every run below cds into its own directory first, so the model path must be
# absolute or SYNTHE will not find it
MODEL=$(cd "$(dirname "$MODEL")" && pwd)/$(basename "$MODEL")
STEM=$(basename "$MODEL"); STEM=${STEM%.*}
OUT=$BASE/run/$TAG
mkdir -p $OUT/chunks $OUT/species
echo "=== $TAG  <- $MODEL  (spec stem $STEM) ==="

if [ ! -s "$OUT/eos/$STEM.mol" ]; then
  d=$OUT/eos; rm -rf $d; mkdir -p $d
  ( cd $d && ATLAS12=$A $A/bin/synthe.exe "$MODEL" wlbeg=500 wlend=500.2 \
      resolu=300000 more_output=yes > log 2>&1 )
  echo "eos: $(ls $d | tr '\n' ' ')"
fi

for (( w=355; w<1000; w+=25 )); do
  d=$OUT/chunks/w${w}
  [ -s "$d/$STEM.spec" ] && continue
  rm -rf $d; mkdir -p $d
  ( cd $d && ATLAS12=$A $A/bin/synthe.exe "$MODEL" wlbeg=$w wlend=$((w+25)) \
      resolu=300000 > log 2>&1 )
  echo "  chunk ${w}: $(wc -l < $d/$STEM.spec) pts $(grep -m1 'Lines:' $d/log | tr -s ' ')"
done

# The environment directory is per star, not per species: several stars run at
# once, and a shared one gets deleted out from under whichever is reading it.
while IFS=$'\t' read -r sp kind rest; do
  [ -z "$sp" ] && continue
  d=$OUT/species/$sp
  [ -s "$d/$STEM.spec" ] && { echo "  $sp: have it"; continue; }
  D=$OUT/env; rm -rf $D; mkdir -p $D/data $D/bin
  ln -sf $A/bin/synthe.exe $D/bin/synthe.exe
  for f in $A/data/*; do ln -sf "$f" $D/data/ 2>/dev/null || true; done
  rm -f $D/data/lines.list; : > $D/data/lines.list
  for f in $rest; do
    if [ "$kind" = gfall ]; then
      ln -sf $BASE/cache/species/$f $D/data/$f
      printf 'gfall\t%s\n' "$f" >> $D/data/lines.list
    else
      printf '%s\tmol/%s\n' "$kind" "$f" >> $D/data/lines.list
    fi
  done
  rm -rf $d; mkdir -p $d
  ( cd $d && ATLAS12=$D $D/bin/synthe.exe "$MODEL" wlbeg=355 wlend=1000 \
      resolu=300000 > log 2>&1 )
  echo "  $sp: $(grep -m1 'Lines:' $d/log | tr -s ' ')"
done < $BASE/run/species.list
rm -rf $OUT/env
echo "DONE $TAG"
