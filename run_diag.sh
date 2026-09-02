#!/bin/bash
A=/Users/cconroy/kurucz/atlas12
MODEL=/Users/cconroy/kurucz/grids/THESUN/atm/ap00t5777g4.44at12.dat
RUN=$(cd "$(dirname "$0")" && pwd)/run
BASE=$RUN/diag
mkdir -p $BASE
for win in "430 435" "500 505" "660 665"; do
  set -- $win; w1=$1; w2=$2
  for cfg in full gfall gfpred; do
    if [ "$cfg" = "full" ]; then export ATLAS12=$A
    else export ATLAS12=$RUN/alt_$cfg; fi
    d=$BASE/${w1}_${cfg}; rm -rf $d; mkdir -p $d
    ( cd $d && $ATLAS12/bin/synthe.exe $MODEL wlbeg=$w1 wlend=$w2 resolu=300000 > log 2>&1 )
    echo "$w1-$w2 $cfg: $(grep -m1 'Lines:' $d/log)"
  done
done
echo DIAG_DONE
