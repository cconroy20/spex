"""Line identification index for the web app.

Every line the app can show comes from the same primary data the synthesis
used -- Kurucz's gfallvac for atoms, his molecular line files for molecules.
The only thing added here is a predicted central depth, computed at one
reference depth with the machinery in sunlib, which is what decides whether
a line is worth carrying to the browser at all.

The depth actually reported to the user is measured off the synthesized
per-species spectrum; the predicted depth is used to rank blends and to cull.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sunlib as S                                                   # noqa: E402
import resline_lib as R                                              # noqa: E402
from starcfg import star, species as _species                                             # noqa: E402

ST = star(sys.argv[1] if len(sys.argv) > 1 else 'sun')
_SP = _species()

ATOMS = [t[1:] for t, *_ in _SP if t.startswith('a')]
# TiO ships as a 4.2 GB binary with no text form, so it is synthesized but not
# line-indexed: naming one rotational line out of 97 million is not the point
MOLS = {n: fs for t, n, g, k, fs in _SP
        if k == 'mol' and not any(f.endswith('.bin') for f in fs)}
MOLDIR = Path.home() / 'kurucz' / 'atlas12' / 'data' / 'mol'

LAM_LO, LAM_HI = 3540.0, 10010.0          # vacuum A, a little outside the grid


def read_atom(path):
    """gfall records with the level labels kept."""
    def f(s, d=0.0):
        try:
            return float(s)
        except ValueError:
            return d

    cols = {k: [] for k in ('wl_nm', 'loggf', 'code', 'elo', 'eup',
                            'loggr', 'loggs', 'loggw', 'x1', 'x2', 'aut')}
    lo_lab, up_lab, lo_j, up_j = [], [], [], []
    for line in open(path, errors='replace'):
        c = f(line[18:24], -1.0)
        if c <= 0:
            continue
        w = f(line[0:11]) * 10.0
        if not (LAM_LO <= w <= LAM_HI):
            continue
        flag = line[141:144] if len(line) > 144 else '   '
        if flag == 'COR':
            continue
        e1, e2 = f(line[24:36]), f(line[52:64])
        lower_first = abs(e1) <= abs(e2)
        cols['wl_nm'].append(f(line[0:11]))
        cols['loggf'].append(f(line[11:18], -99.0))
        cols['code'].append(c)
        cols['elo'].append(min(abs(e1), abs(e2)))
        cols['eup'].append(max(abs(e1), abs(e2)))
        cols['loggr'].append(f(line[80:86]))
        cols['loggs'].append(f(line[86:92]))
        cols['loggw'].append(f(line[92:98]))
        cols['x1'].append(f(line[109:115]))
        cols['x2'].append(f(line[118:124]))
        cols['aut'].append(1.0 if flag == 'AUT' else 0.0)
        a, b = line[41:52].strip(), line[69:80].strip()
        ja, jb = f(line[36:41]), f(line[64:69])
        lo_lab.append(a if lower_first else b)
        up_lab.append(b if lower_first else a)
        lo_j.append(ja if lower_first else jb)
        up_j.append(jb if lower_first else ja)
    g = {k: np.array(v) for k, v in cols.items()}
    return g, np.array(lo_lab), np.array(up_lab), np.array(lo_j), np.array(up_j)


def read_mol(paths):
    def f(s, d=0.0):
        try:
            return float(s)
        except ValueError:
            return d

    cols = {k: [] for k in ('wl_nm', 'loggf', 'elo', 'eup', 'jlo', 'icode',
                            'iso', 'loggr', 'upperX')}
    lo_lab, up_lab = [], []
    for p in paths:
        for line in open(p, errors='replace'):
            e, ep = f(line[22:32]), f(line[37:48])
            lo_e, up_e = min(abs(e), abs(ep)), max(abs(e), abs(ep))
            dE = up_e - lo_e
            if dE <= 0:
                continue
            w = 1.0e8 / dE                                   # vacuum A
            if not (LAM_LO <= w <= LAM_HI):
                continue
            xj, xjp = f(line[17:22]), f(line[32:37])
            cols['wl_nm'].append(w / 10.0)
            cols['loggf'].append(f(line[10:17], -99.0))
            cols['elo'].append(lo_e)
            cols['eup'].append(up_e)
            cols['jlo'].append(xj if abs(e) <= abs(ep) else xjp)
            # fehfx.dat writes 156 where the species code belongs; SYNTHE's
            # own reader remaps it (read_molec_ascii), so do the same here
            ic = f(line[48:52])
            cols['icode'].append(126.0 if ic == 156.0 else ic)
            cols['iso'].append(f(line[68:70]))
            cols['loggr'].append(f(line[70:74]))
            cols['upperX'].append(1.0 if line[60:61] == 'X' else 0.0)
            a, b = line[52:60].strip(), line[60:68].strip()
            lower_first = abs(e) <= abs(ep)
            lo_lab.append(a if lower_first else b)
            up_lab.append(b if lower_first else a)
    return ({k: np.array(v) for k, v in cols.items()},
            np.array(lo_lab), np.array(up_lab))


# Two reference depths.  tau = 1 is where most atomic lines are set, but the
# molecules live in the cool upper photosphere and are all but dissociated
# there, so ranking on tau = 1 alone throws away CN and C2 lines that SYNTHE
# puts at 10% depth.  A line is kept if it is strong at either depth.
TAUS = (1.0, 0.1)


def main():
    phots = [S.Photosphere(tau_ref=t, atm_path=ST['model']) for t in TAUS]
    Ls = [S.LineOpacity(p, lam_lo_A=LAM_LO, lam_hi_A=LAM_HI, mol_path=ST['mol'])
          for p in phots]
    print(f'{ST["name"]}: {ST["model"].name}')
    for t, p in zip(TAUS, phots):
        print(f'reference depth tau={t}: j={p.j}  T={p.T:.0f} K  '
              f'ne={p.ne:.3e}  vturb={p.vturb/1e5:.2f} km/s')
    L = Ls[0]

    rows = []
    for code in ATOMS:
        g, lo_lab, up_lab, lo_j, up_j = read_atom(HERE / 'cache' / 'species' /
                                                  f'gf_{code}.dat')
        if not len(g['wl_nm']):
            continue
        keep = g['aut'] == 0
        g = {k: v[keep] for k, v in g.items()}
        lo_lab, up_lab = lo_lab[keep], up_lab[keep]
        lo_j, up_j = lo_j[keep], up_j[keep]
        Ps = [x.prepare(g) for x in Ls]
        assert len(Ps[0]['lam']) == len(lo_lab), (code, len(Ps[0]['lam']), len(lo_lab))
        rows.append(dict(name=S.species_name(float(code)), lam=Ps[0]['lam'],
                         peaks=[P['peak'] for P in Ps],
                         loggf=g['loggf'] + g['x1'] + g['x2'],
                         elo=g['elo'], eup=g['eup'], lo=lo_lab, up=up_lab,
                         jlo=lo_j, jup=up_j))
        print(f'  {S.species_name(float(code)):6s} {len(Ps[0]["lam"]):8d} lines')

    for name, files in MOLS.items():
        m, lo_lab, up_lab = read_mol([MOLDIR / f for f in files])
        Ps = [x.prepare_molecules(m) for x in Ls]
        P = Ps[0]
        # prepare_molecules drops records its dispatch table cannot place;
        # rebuild the same mask so the labels stay aligned
        import moldisp
        disp = moldisp.parse_dispatch()
        Tq, Q = S.read_mol_partition()
        ic, iso = m['icode'].astype(int), m['iso'].astype(int)
        ok = np.array([c in S.MOL_NAME and S.MOL_NAME[c] in Q for c in ic])
        ok &= np.array([(i, c) in disp or (i, None) in disp
                        for i, c in zip(iso, ic)])
        assert len(P['lam']) == int(ok.sum())
        rows.append(dict(name=name, lam=P['lam'],
                         peaks=[q['peak'] for q in Ps],
                         loggf=m['loggf'][ok], elo=m['elo'][ok],
                         eup=m['eup'][ok], lo=lo_lab[ok], up=up_lab[ok],
                         jlo=m['jlo'][ok], jup=np.zeros(int(ok.sum()))))
        print(f'  {name:6s} {len(P["lam"]):8d} lines')

    # predicted central depth, weak-line against the continuum at this depth
    for r in rows:
        ix = np.clip(np.rint(np.log(r['lam'] / L.lam[0]) / L.dln).astype(int),
                     0, L.n - 1)
        r['depth'] = np.max([pk / x.kcont[ix] / (1.0 + pk / x.kcont[ix])
                             for pk, x in zip(r['peaks'], Ls)], axis=0)

    alld = np.concatenate([r['depth'] for r in rows])
    print(f'\ntotal {len(alld):,} lines in band')
    for cut in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05):
        print(f'  depth > {cut*100:5.1f}% : {int((alld > cut).sum()):8,d}')
    np.savez(HERE / 'cache' / f'{ST["tag"]}_lineindex.npz',
             **{f'{i}_{k}': r[k] for i, r in enumerate(rows)
                for k in ('lam', 'depth', 'loggf', 'elo', 'eup', 'lo', 'up',
                          'jlo', 'jup')},
             names=np.array([r['name'] for r in rows]))
    print(f'wrote cache/{ST["tag"]}_lineindex.npz')


main()
