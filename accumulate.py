import time
import numpy as np
import sunlib

t0 = time.time()
g = sunlib.parse_gfall_full()
p = sunlib.Photosphere()
LO = sunlib.LineOpacity(p)
L = LO.prepare(g)

TOP = ['Ca II', 'Fe I', 'Mg I', 'H I', 'Na I', 'Al I', 'Ca I', 'Si I',
       'Fe II', 'Cr I', 'Ti II', 'K I', 'Ni I', 'Mn I', 'Ti I', 'C I',
       'Sc II', 'Co I', 'Cr II', 'Mg II', 'Ba II', 'Sr II', 'V I', 'Y II',
       'Zr II', 'S I', 'O I', 'Zn I', 'Cu I', 'Sr I']
groups = {}
for c in np.unique(L['code']):
    nm = sunlib.species_name(c)
    if nm in TOP:
        groups[c] = nm

out, nused, nalive, ntot = LO.accumulate(L, groups)
print(f'accumulated {nused} lines ({nalive} above the floor of {ntot} usable), '
      f'{time.time()-t0:.0f}s')

np.savez_compressed('cache/opacity_species.npz', lam_vac=LO.lam,
                    kcont=LO.kcont, T=p.T, ne=p.ne, nH=p.nH, jdepth=p.j,
                    **{k.replace(' ', '_'): v for k, v in out.items()})

tot = sum(out.values())
print()
print(f'reference depth: tau_5000 = 1, layer {p.j}, T = {p.T:.0f} K')
print(f'{"species":9s} {"peak k_line/k_cont":>19s} {"frac of grid with k_l>k_c":>27s}')
order = sorted(out, key=lambda k: -np.trapz(np.minimum(out[k] / LO.kcont, 1e3), LO.lam))
for k in order[:16]:
    r = out[k] / LO.kcont
    print(f'{k:9s} {r.max():19.3e} {100*np.mean(r > 1):26.2f}%')
print(f'{"TOTAL":9s} {(tot/LO.kcont).max():19.3e} {100*np.mean(tot/LO.kcont > 1):26.2f}%')
