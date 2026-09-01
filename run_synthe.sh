#!/bin/bash
# Full-optical solar synthesis with SYNTHE, in memory-safe wavelength chunks.
export ATLAS12=/Users/cconroy/kurucz/atlas12
MODEL=/Users/cconroy/kurucz/grids/THESUN/atm/ap00t5777g4.44at12.dat
OUT=$(cd "$(dirname "$0")" && pwd)/run/chunks
mkdir -p $OUT
STEP=25
for (( w=380; w<1000; w+=STEP )); do
  w2=$(( w + STEP ))
  d=$OUT/w${w}
  rm -rf $d; mkdir -p $d
  ( cd $d && $ATLAS12/bin/synthe.exe $MODEL wlbeg=$w wlend=$w2 resolu=300000 \
      > synthe.log 2>&1 )
  n=$(wc -l < $d/ap00t5777g4.44at12.spec 2>/dev/null || echo 0)
  echo "chunk ${w}-${w2} nm: $n points  $(grep -m1 'Lines:' $d/synthe.log)"
done
echo "ALL CHUNKS DONE"
