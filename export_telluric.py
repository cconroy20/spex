"""Telluric transmission for the explorer.

Source: ESO SkyCalc (Noll et al. 2012, A&A 543, A92; Jones et al. 2013, A&A
560, A91), the line-by-line molecular transmission table shipped with alf.

Two things make this the right table rather than SkyCalc's `skytable` product:
it is molecular absorption ONLY -- Rayleigh optical depth is 0.0016 at 0.34 um
against 0.51 for the full-transmission product -- which is what "telluric
lines" means and what can be multiplied into a stellar spectrum without
imposing a spurious blue slope; and it is sampled at R = 60,000 rather than
20,000, so the lines are resolved rather than smoothed.

The file is at altitude 30 degrees, i.e. airmass 2 (its `A30` label), and is
rescaled to the zenith by T -> T^(1/2), which is exact for pure absorption.
That airmass reading was checked against SkyCalc's own documented airmass-1
model: the O2 A-band equivalent width then agrees to 1.1%.  Weak bands differ
by up to 39% between the two SkyCalc products, so the absolute depths carry
that much uncertainty; line positions are unaffected.

The grid is identical for every star, so this ships once at data/telluric.bin.
"""
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path.home() / 'memos' / 'resline'))
sys.path.insert(0, str(HERE))
import resline_lib as R                                              # noqa: E402
import binfmt                                                        # noqa: E402
from starcfg import star                                             # noqa: E402

SRC = Path.home() / 'sps' / 'alf' / 'sky' / 'SkyCalc' / 'LBL_A30_s0_w025_R0060000_T.fits'
AIRMASS_FILE = 2.0
AIRMASS_OUT = 1.0
PWV_MM = 2.5

d = fits.open(SRC)[1].data
lam_vac = d['lam'] * 1e4                      # micron -> Angstrom, vacuum
trans = np.clip(d['trans'], 0.0, 1.0) ** (AIRMASS_OUT / AIRMASS_FILE)

S = star('sun')
spec = np.load(HERE / 'cache' / f'{S["tag"]}_species.npz')
grid = spec['lam_vac']
t = np.interp(grid, lam_vac, trans, left=1.0, right=1.0)

OUT = HERE / 'web' / 'data'
OUT.mkdir(parents=True, exist_ok=True)
q = np.clip(t * 65535.0, 0, 65535).astype('<u2')
n_ov = len(grid) // 52
binfmt.write_series(OUT / 'telluric.bin', q)
b = np.concatenate([x.reshape(-1, 52).min(axis=1) for x in (q[:n_ov * 52],)]
                   + [q[:n_ov * 52].reshape(-1, 52).max(axis=1),
                      q[:n_ov * 52].reshape(-1, 52).mean(axis=1)]).astype('<u2')
binfmt.write_series(OUT / 'telluric_ov.bin', b, chunk=n_ov)

air = R.vac_to_air(grid)
print(f'telluric: {len(grid):,} points on the shared grid, '
      f'airmass {AIRMASS_OUT:.1f}, PWV {PWV_MM} mm')
print(f'  mean transmission {t.mean():.4f}, deepest {t.min():.4f}')
for lo, hi, lab in [(6860, 6960, 'O2 B'), (7580, 7720, 'O2 A'),
                    (9300, 9600, 'H2O 0.94'), (13000, 15000, 'H2O 1.4'),
                    (18000, 20000, 'H2O 1.9')]:
    m = (air > lo) & (air < hi)
    print(f'  {lab:9s} {lo/1e4:.2f}-{hi/1e4:.2f} um: mean {t[m].mean():.3f}  min {t[m].min():.3f}')
print(f'  files: {(OUT/"telluric.bin").stat().st_size/1e3:.0f} kB full, '
      f'{(OUT/"telluric_ov.bin").stat().st_size/1e3:.0f} kB overview')
