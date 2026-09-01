"""Per-species solar spectra over 3800-10000 A."""
import sys, time
from pathlib import Path
import numpy as np
sys.path.append(str(Path.home() / 'memos' / 'resline'))
import sunlib, specflux as SF

GROUPS = ['Fe I', 'Fe II', 'Ca I', 'Ca II', 'Mg I', 'Mg II', 'Na I', 'Si I',
          'Al I', 'Ti I', 'Ti II', 'Cr I', 'Cr II', 'Mn I', 'Ni I', 'Co I',
          'V I', 'Sc II', 'K I', 'C I', 'Ba II', 'Sr II', 'Y II', 'Zr II',
          'H I', 'CH', 'CN', 'MgH', 'C2', 'SiH', 'NH', 'OH']

phot = sunlib.Photosphere()
LO = sunlib.LineOpacity(phot)
L, nU, labels = SF.build_lines(phot, LO.potion)
eng = SF.Engine(phot, L, nU, LO.lam)
d = np.load('cache/cont_cal.npz'); kc_lg, kc_cal = d['lg'], d['kc']
w = np.load('cache/lines_hw.npz'); ix, hw, keep = w['ix'], w['hw'], w['keep']

names = [g for g in GROUPS if g in labels] + ['other']
gof = {g: i for i, g in enumerate(names)}
sid2g = np.array([gof.get(lab, gof['other']) for lab in labels])
gid = sid2g[L['sidx']]
print(f'{len(names)} groups, {keep.sum()} lines above the floor')

lo, hi = (float(sys.argv[1]), float(sys.argv[2])) if len(sys.argv) > 2 else (0., 1e9)
m = (LO.lam >= lo) & (LO.lam <= hi)
i0, i1 = int(np.argmax(m)), int(np.argmax(m) + m.sum())
sub = LO.lam[i0:i1]
eng.lam = sub
t0 = time.time()
res = SF.run(phot, eng, kc_lg, kc_cal, ix - i0, hw, keep, gid, len(names), names,
             chunk=int(sys.argv[3]) if len(sys.argv) > 3 else 6000)
print(f'done in {time.time()-t0:.0f}s')
out = sys.argv[4] if len(sys.argv) > 4 else 'cache/species_flux.npz'
np.savez_compressed(out, lam_vac=sub, **{k.replace(' ', '_'): v for k, v in res.items()})
print('wrote', out)
