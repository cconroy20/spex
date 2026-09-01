"""Stitch one star's SYNTHE chunks into a single spectrum.

    python3 assemble.py [tag]        (default: sun)

The Sun additionally gets checked against the IAG flux atlas; for the others
there is nothing to check against here, so only the internal diagnostics run.
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path.home() / 'memos' / 'resline'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import resline_lib as R                                              # noqa: E402
from starcfg import star                                             # noqa: E402

HERE = Path(__file__).resolve().parent
C_A = 2.99792458e18          # speed of light [A/s]
SIGMA_SB = 5.670374419e-5

S = star(sys.argv[1] if len(sys.argv) > 1 else 'sun')
rows = []
for d in sorted((S['run'] / 'chunks').glob('w*'), key=lambda p: int(p.name[1:])):
    f = d / f'{S["stem"]}.spec'
    if f.exists() and f.stat().st_size:
        rows.append(np.loadtxt(f))
    else:
        print(f'  MISSING {d.name}')
d = np.vstack(rows)
lam_vac, hnu, hnu_c = d[:, 0], d[:, 1], d[:, 2]

o = np.argsort(lam_vac)
lam_vac, hnu, hnu_c = lam_vac[o], hnu[o], hnu_c[o]
keep = np.r_[True, np.diff(lam_vac) > 1e-6]
print(f'{S["name"]}: {len(lam_vac)} points, {np.sum(~keep)} duplicates at chunk joins')
lam_vac, hnu, hnu_c = lam_vac[keep], hnu[keep], hnu_c[keep]

step = np.diff(np.log(lam_vac))
print(f'log-lambda step: median {np.median(step):.3e}  max {step.max():.3e} '
      f'(R = {1/np.median(step):.0f}); largest gap {step.max()/np.median(step):.2f} x median')

lam_air = R.vac_to_air(lam_vac)
flux = 4.0 * np.pi * hnu * C_A / lam_vac**2          # erg/s/cm^2/A at the surface
cont = 4.0 * np.pi * hnu_c * C_A / lam_vac**2
norm = flux / cont

np.savez_compressed(HERE / 'cache' / f'{S["tag"]}_synthe.npz', lam_air=lam_air,
                    lam_vac=lam_vac, flux=flux, cont=cont, norm=norm)

print('\n' + '=' * 70 + '\nVALIDATION\n' + '=' * 70)
band = np.trapz(flux, lam_air)
tot = SIGMA_SB * float(S['teff'])**4
print(f'  integral of F_lam over {lam_air[0]:.0f}-{lam_air[-1]:.0f} A = {band:.4e} erg/s/cm^2')
print(f'  sigma Teff^4                        = {tot:.4e}   -> band holds {100*band/tot:.1f}%')
print(f'  continuum-only integral             = {np.trapz(cont, lam_air):.4e}'
      f'   -> line blanketing removes {100*(1-band/np.trapz(cont,lam_air)):.1f}%')
i5 = np.argmin(abs(lam_air - 5000.0))
print(f'  F_lam(5000 A) continuum = {cont[i5]:.3e} erg/s/cm^2/A')
print(f'  mean normalized flux    = {norm.mean():.4f}   (line absorption {100*(1-norm.mean()):.1f}%)')

if S['tag'] == 'sun':
    print()
    for lo, hi in [(4300, 4320), (5160, 5190), (5885, 5900), (6555, 6570), (8490, 8510)]:
        wo, fo = R.read_iag(lo - 3, hi + 3)
        if len(wo) < 10:
            continue
        m = (lam_air > lo) & (lam_air < hi)
        fs = np.interp(wo, lam_air[m], norm[m])
        print(f'  {lo}-{hi} A  observed mean {fo.mean():.4f}   synthetic mean {fs.mean():.4f}'
              f'   ratio {fs.mean()/fo.mean():.3f}')
