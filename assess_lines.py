"""How much the line index grows once it has to serve several stars.

    python3 assess_lines.py tag [tag ...]

The catalogue half of the index (wavelength, species, log gf, chi, level
labels) is a property of the atom and is identical for every star; only the
two depth columns differ.  So the question is how much bigger the shared
catalogue gets when it has to hold every line that is visible in ANY star,
against shipping a separate index per star.
"""
import sys
import zlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CUT = 0.002

tags = sys.argv[1:]
keys, depth = None, {}
for t in tags:
    d = np.load(HERE / 'cache' / f'{t}_lineindex.npz', allow_pickle=True)
    names = [str(x) for x in d['names']]
    lam, dep, sp = [], [], []
    for i, nm in enumerate(names):
        lam.append(d[f'{i}_lam']); dep.append(d[f'{i}_depth'])
        sp.append(np.full(len(d[f'{i}_lam']), i))
    lam = np.concatenate(lam); dep = np.concatenate(dep); sp = np.concatenate(sp)
    k = np.rint(lam * 1e4).astype(np.int64) * 100 + sp      # line identity
    if keys is None:
        keys = k
    elif len(k) != len(keys) or not np.array_equal(k, keys):
        raise SystemExit(f'{t}: line list differs from the first star')
    depth[t] = dep
    print(f'  {t:10s} {int((dep > CUT).sum()):8,d} lines above {CUT:.1%}')

alive = np.zeros(len(keys), bool)
for t in tags:
    alive |= depth[t] > CUT
print(f'\n  union     {int(alive.sum()):8,d} lines'
      f'   ({alive.sum()/max((depth[t] > CUT).sum() for t in tags):.2f}x the largest star)')

# real columns, deflated, split into the half that is shared and the half
# that is not
d0 = np.load(HERE / 'cache' / f'{tags[0]}_lineindex.npz', allow_pickle=True)
names = [str(x) for x in d0['names']]
cat = {k: np.concatenate([d0[f'{i}_{k}'] for i in range(len(names))])[alive]
       for k in ('lam', 'loggf', 'elo', 'lo', 'up', 'jlo', 'jup')}
sp = np.concatenate([np.full(len(d0[f'{i}_lam']), i)
                     for i in range(len(names))])[alive]
o = np.argsort(cat['lam'])
labs = sorted(set(cat['lo']) | set(cat['up']))
li = {v: k for k, v in enumerate(labs)}
cols = [np.diff(np.rint(cat['lam'][o] * 1e4).astype(np.int64), prepend=0).astype('<u4'),
        sp[o].astype('u1'),
        np.rint(cat['loggf'][o] * 1000).clip(-32768, 32767).astype('<i2'),
        np.rint(cat['elo'][o] / 8065.54429 * 1000).clip(0, 65535).astype('<u2'),
        np.array([li[v] for v in cat['lo'][o]], '<u2'),
        np.array([li[v] for v in cat['up'][o]], '<u2'),
        np.rint(cat['jlo'][o] * 2).clip(0, 255).astype('u1'),
        np.rint(cat['jup'][o] * 2).clip(0, 255).astype('u1')]
shared = sum(len(zlib.compress(c.tobytes(), 9)) for c in cols) \
    + len(zlib.compress('\n'.join(labs).encode(), 9))
print(f'\n  shared catalogue {shared/1e6:.2f} MB deflated, shipped once')
for t in tags:
    dep = depth[t][alive][o]
    per = len(zlib.compress(np.rint(dep * 255).clip(0, 255).astype('u1').tobytes(), 9))
    print(f'    + {t:10s} depths {per/1e6:.2f} MB')
