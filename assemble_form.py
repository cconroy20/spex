"""Stitch the per-chunk formation depths into one array per star.

    python3 assemble_form.py [tag]

Two quantities at the layer where tau_lambda = 1: log tau_5000, which is the
standard reference scale and compares across stars, and the temperature there,
which is what the line actually responds to.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from starcfg import star                                             # noqa: E402

S = star(sys.argv[1] if len(sys.argv) > 1 else 'sun')
lam, lt, tf = [], [], []
missing = 0
for d in sorted((S['run'] / 'chunks').glob('w*'), key=lambda p: int(p.name[1:])):
    f = d / 'form.npz'
    if not f.exists():
        missing += 1
        continue
    z = np.load(f)
    lam.append(z['lam_vac']); lt.append(z['logtau5000']); tf.append(z['tform'])
if missing:
    print(f'  WARNING: {missing} chunks have no form.npz')
lam = np.concatenate(lam); lt = np.concatenate(lt); tf = np.concatenate(tf)
o = np.argsort(lam)
lam, lt, tf = lam[o], lt[o], tf[o]
keep = np.r_[True, np.diff(lam) > 1e-6]
lam, lt, tf = lam[keep], lt[keep], tf[keep]

np.savez_compressed(HERE / 'cache' / f'{S["tag"]}_form.npz',
                    lam_vac=lam, logtau5000=lt, tform=tf)
print(f'{S["name"]}: {len(lam):,} wavelengths -> cache/{S["tag"]}_form.npz')
print(f'  log tau5000  median {np.median(lt):+.2f}   5th {np.percentile(lt,5):+.2f}'
      f'   min {lt.min():+.2f}')
print(f'  T(tau_lam=1) median {np.median(tf):.0f} K  5th {np.percentile(tf,5):.0f} K'
      f'  min {tf.min():.0f} K')
