"""Collect one star's per-species SYNTHE runs onto a single grid.

    python3 assemble_species.py [tag]        (default: sun)
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from starcfg import star, species                                             # noqa: E402

HERE = Path(__file__).resolve().parent
NAMES = {t: n for t, n, *_ in species()}

S = star(sys.argv[1] if len(sys.argv) > 1 else 'sun')
out, lam = {}, None
for tag, name in NAMES.items():
    f = S['run'] / 'species' / tag / f'{S["stem"]}.spec'
    if not f.exists() or not f.stat().st_size:
        print(f'  MISSING {tag} ({name})')
        continue
    d = np.loadtxt(f)
    if lam is None:
        lam = d[:, 0]
    elif len(d) != len(lam) or abs(d[0, 0] - lam[0]) > 1e-6:
        print(f'  regridding {name}')
        out[name] = np.interp(lam, d[:, 0], d[:, 1] / d[:, 2])
        continue
    out[name] = d[:, 1] / d[:, 2]

s = np.load(HERE / 'cache' / f'{S["tag"]}_synthe.npz')
out['all'] = np.interp(lam, s['lam_vac'], s['norm'])

np.savez_compressed(HERE / 'cache' / f'{S["tag"]}_species.npz', lam_vac=lam,
                    **{k.replace(' ', '_'): v for k, v in out.items()})
print(f'\n{S["name"]}: {len(out)} species on {len(lam)} points '
      f'-> cache/{S["tag"]}_species.npz')
print(f'\nmean line absorption over {lam[0]:.0f}-{lam[-1]:.0f} A:')
for v, k in sorted(((1 - out[k].mean(), k) for k in out), reverse=True):
    print(f'   {k:9s} {100*v:6.3f}%')
