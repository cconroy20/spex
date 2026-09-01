#!/bin/bash
export ATLAS12=/Users/cconroy/kurucz/atlas12
MODEL=/Users/cconroy/kurucz/grids/THESUN/atm/ap00t5777g4.44at12.dat
OUT=$(cd "$(dirname "$0")" && pwd)/run/chunks
for w in 355 380; do :; done
d=$OUT/w355; rm -rf $d; mkdir -p $d
( cd $d && $ATLAS12/bin/synthe.exe $MODEL wlbeg=355 wlend=380 resolu=300000 > synthe.log 2>&1 )
echo "chunk 355-380 nm: $(wc -l < $d/ap00t5777g4.44at12.spec) points  $(grep -m1 'Lines:' $d/synthe.log)"
echo BLUE_DONE
