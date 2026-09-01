#!/bin/bash
# One SYNTHE run per species, over the whole optical band.  Each species gets a
# lines.list of its own, so the profile physics is SYNTHE's throughout --
# hydrogen included, which needs the Holtsmark quasi-static Stark treatment a
# Voigt cannot reproduce.
set -e
A=/Users/cconroy/kurucz/atlas12
BASE=$PWD
MODEL=/Users/cconroy/kurucz/grids/THESUN/atm/ap00t5777g4.44at12.dat
OUT=$BASE/run/species; mkdir -p $OUT

mkdir_env () {   # $1 = tag, rest = lines.list body on stdin
  D=$BASE/run/env_$1; rm -rf $D; mkdir -p $D/data $D/bin
  ln -sf $A/bin/synthe.exe $D/bin/synthe.exe
  for f in $A/data/*; do ln -sf "$f" $D/data/ 2>/dev/null || true; done
  rm -f $D/data/lines.list
  cat > $D/data/lines.list
  echo $D
}

run_one () {   # $1 = tag, $2 = env dir
  d=$OUT/$1; rm -rf $d; mkdir -p $d
  ( cd $d && ATLAS12=$2 $2/bin/synthe.exe $MODEL wlbeg=355 wlend=1000 resolu=300000 > log 2>&1 )
  echo "  $1: $(grep -m1 'Lines:' $d/log | tr -s ' ') $(grep -m1 Runtime $d/log | tr -s ' ')"
}

for code in 1.00 6.00 11.00 12.00 13.00 14.00 19.00 20.00 20.01 21.01 \
            22.00 22.01 23.00 24.00 24.01 25.00 26.00 26.01 27.00 28.00; do
  f=$BASE/cache/species/gf_${code}.dat
  D=$(printf 'gfall\tgf_%s.dat\n' $code | mkdir_env a$code)
  ln -sf $f $D/data/gf_${code}.dat
  run_one a$code $D
done

# molecules: already one file per species in Kurucz's data directory
mol_group () {  # $1 = name, rest = files
  tag=$1; shift
  D=$( for f in "$@"; do printf 'mol\tmol/%s\n' $f; done | mkdir_env m$tag )
  run_one m$tag $D
}
mol_group CH  chjorg.dat
mol_group CN  cnaxbrooke.dat cnbxbrooke.dat cnxx12brooke.dat
mol_group MgH mgh24_owens22.dat mgh25_owens22.dat mgh26_owens22.dat mghax.dat mghbx.dat
mol_group C2  c2ax.dat c2ba.dat c2dabrooke.dat c2ea.dat
mol_group SiH sihaxsightly.asc sihxxsightly.asc
echo SPECIES_DONE
