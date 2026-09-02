"""Split gfallvac into one file per species code, for the per-species runs.

    python3 split_gfall.py 3.00 7.00 ...        (default: everything we ship)

SYNTHE has no species filter, but it reads whatever file `lines.list` names,
so giving each species its own copy of gfall is the whole trick.
"""
import sys
from pathlib import Path

GFALL = Path.home() / 'kurucz' / 'atlas12' / 'data' / 'gfallvac08oct17.dat'
OUT = Path(__file__).resolve().parent / 'cache' / 'species'
OUT.mkdir(parents=True, exist_ok=True)

LAM_LO, LAM_HI = 3540.0, 25010.0        # the band build_lines.py indexes

codes = sys.argv[1:]
want = {f'{float(c):6.2f}': c for c in codes}
fh = {c: open(OUT / f'gf_{c}.dat', 'w') for c in codes}
n = {c: 0 for c in codes}
inband = {c: 0 for c in codes}
for line in open(GFALL, errors='replace'):
    k = line[18:24]
    c = want.get(k)
    if c is None:
        continue
    fh[c].write(line)
    n[c] += 1
    try:
        w = float(line[0:11]) * 10.0
    except ValueError:
        continue
    if LAM_LO <= w <= LAM_HI:
        inband[c] += 1
for f in fh.values():
    f.close()
for c in codes:
    print(f'  {c:>6s}  {n[c]:9,d} records  {inband[c]:8,d} in '
          f'{LAM_LO:.0f}-{LAM_HI:.0f} A')
