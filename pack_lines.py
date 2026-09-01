"""Pack the line index: one shared catalogue, plus two depth bytes per star.

    python3 pack_lines.py tag [tag ...]

Wavelength, species, log gf, excitation potential and the level labels are
properties of the atom, identical for every star; only the depths differ.  So
the catalogue is built once over the union of every star's visible lines and
shipped at web/data/lines/, and each star adds only
web/data/<tag>/lines/<block>.bin -- a predicted and a measured depth per line.

The union turns out to be barely larger than the richest single star: the
strong lines of an F dwarf and a metal-poor giant are close to a subset of a
cool giant's, so sharing costs almost nothing.
"""
import json
import re
import sys
import zlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path.home() / 'memos' / 'resline'))
sys.path.insert(0, str(HERE))
import resline_lib as R                                              # noqa: E402
from starcfg import star, species                                    # noqa: E402

CUT = 0.002            # predicted central depth; below this a line is invisible
BLOCK = 500.0          # A of air per block
MERGE_A = 0.01         # hyperfine / isotopic components closer than this merge
CM_EV = 1.0 / 8065.54429
LOGCUT = np.log10(CUT)
MOLECULES = {n for t, n, *_ in species() if t.startswith('m')}

# state letter with the vibrational level after it.  Most files pad it to two
# digits (A00e1); CO right-justifies it instead, so v=2 is written "X 2" and
# v=23 is "X23".  Allow either, or CO falls back to an unparsed label.
LEVEL = re.compile(r'^(\d*)([A-Za-z])\s*(\d+)')
RYD_H = 109677.58
SERIES = {1: 'Lyman', 2: 'Balmer', 3: 'Paschen', 4: 'Brackett', 5: 'Pfund'}
GREEK = {1: 'α', 2: 'β', 3: 'γ', 4: 'δ'}


def band(lo, up):
    """Molecular level codes -> band designation, A00e1/X01f1 -> A-X (0,1)."""
    a, b = LEVEL.match(up), LEVEL.match(lo)
    if not a or not b:
        return f'{up} – {lo}'.strip()
    return f'{a.group(2)}–{b.group(2)} ({int(a.group(3))},{int(b.group(3))}) band'


def hydrogen(elo, eup):
    """Kurucz's hydrogen records carry the text AVERAGE ENERGIES where the
    level labels go, so name the transition from the energies instead."""
    def n(e):
        return int(round(np.sqrt(RYD_H / max(RYD_H - e, 1e-6))))
    a, b = n(elo), n(eup)
    if not (1 <= a < b <= 40):
        return ''
    g = GREEK.get(b - a, '')
    return f'{SERIES.get(a, f"n={a}")}{" " + g if g else ""} (n = {a} → {b})'


def qdepth(d):
    """Predicted depth on a log scale: it is only used to rank, and ranking
    cares about ratios, so 1/255 of a decade beats 1/255 of the range."""
    q = np.rint(255.0 * (np.log10(np.maximum(d, 1e-9)) - LOGCUT) / (-LOGCUT))
    return np.where(d > CUT, np.clip(q, 1, 255), 0).astype('u1')


tags = sys.argv[1:] or ['sun']
raw = {t: np.load(HERE / 'cache' / f'{t}_lineindex.npz', allow_pickle=True)
       for t in tags}
spec = {t: np.load(HERE / 'cache' / f'{t}_species.npz') for t in tags}
names = [str(x) for x in raw[tags[0]]['names']]
lam_vac = spec[tags[0]]['lam_vac']
lam0 = float(lam_vac[0])
dln = float(np.median(np.diff(np.log(lam_vac))))
nsp = len(lam_vac)

recs = []
for i, name in enumerate(names):
    d0 = raw[tags[0]]
    lam = d0[f'{i}_lam']
    dep = {t: raw[t][f'{i}_depth'] for t in tags}
    keep = np.zeros(len(lam), bool)
    for t in tags:
        keep |= dep[t] > CUT
    if not keep.any():
        continue
    lam = lam[keep]
    gf, elo, eup = d0[f'{i}_loggf'][keep], d0[f'{i}_elo'][keep], d0[f'{i}_eup'][keep]
    lo, up = d0[f'{i}_lo'][keep], d0[f'{i}_up'][keep]
    jl, ju = d0[f'{i}_jlo'][keep], d0[f'{i}_jup'][keep]
    dep = {t: v[keep] for t, v in dep.items()}

    o = np.argsort(lam)
    lam, gf, elo, eup = lam[o], gf[o], elo[o], eup[o]
    lo, up, jl, ju = lo[o], up[o], jl[o], ju[o]
    dep = {t: v[o] for t, v in dep.items()}

    newgrp = np.ones(len(lam), bool)
    newgrp[1:] = ~((np.diff(lam) < MERGE_A) & (np.abs(np.diff(elo)) < 0.05)
                   & (lo[1:] == lo[:-1]) & (up[1:] == up[:-1]))
    g = np.cumsum(newgrp) - 1
    w = 10.0**gf
    tw = np.bincount(g, weights=w)
    mlam = np.bincount(g, weights=w * lam) / tw
    first = np.flatnonzero(newgrp)
    mdep = {t: np.maximum.reduceat(v, first) for t, v in dep.items()}

    L, U, J1, J2 = lo[first], up[first], jl[first], ju[first]
    if name in MOLECULES:
        U = np.array([band(a, b) for a, b in zip(L, U)])
        L = np.full(len(U), '')
        J1 = J2 = np.zeros(len(U))
    elif name == 'H I':
        U = np.array([hydrogen(a, b) for a, b in zip(elo[first], eup[first])])
        L = np.full(len(U), '')
        J1 = J2 = np.zeros(len(U))

    air = R.vac_to_air(mlam)
    ix = np.rint(np.log(mlam / lam0) / dln).astype(int)
    ok = (ix >= 2) & (ix < nsp - 2)
    key = name.replace(' ', '_')
    dsyn = {}
    for t in tags:
        v = np.zeros(len(mlam))
        if key in spec[t].files:
            f = spec[t][key]
            w5 = np.stack([f[np.clip(ix + k, 0, nsp - 1)] for k in range(-2, 3)])
            v = np.where(ok, 1.0 - w5.min(axis=0), 0.0)
        dsyn[t] = v
    recs.append(dict(name=name, lam=air, gf=np.log10(tw), elo=elo[first],
                     lo=L, up=U, jl=J1, ju=J2, dep=mdep, dsyn=dsyn))
    print(f'  {name:6s} {len(mlam):7,d} lines in the union')

allsp = [r['name'] for r in recs]
lo_e = min(r['lam'].min() for r in recs)
hi_e = max(r['lam'].max() for r in recs)
edges = np.arange(np.floor(lo_e / BLOCK) * BLOCK,
                  np.ceil(hi_e / BLOCK) * BLOCK + 1, BLOCK)

CAT = HERE / 'web' / 'data' / 'lines'
CAT.mkdir(parents=True, exist_ok=True)
DEP = {t: HERE / 'web' / 'data' / t / 'lines' for t in tags}
for p in DEP.values():
    p.mkdir(parents=True, exist_ok=True)

index, tot, cat_b, dep_b = [], 0, 0, {t: 0 for t in tags}
for b0, b1 in zip(edges[:-1], edges[1:]):
    cols = {k: [] for k in ('lam', 'sp', 'gf', 'chi', 'll', 'ul', 'jl', 'ju')}
    dp = {t: [] for t in tags}
    ds = {t: [] for t in tags}
    labs, li = [], {}
    for si, r in enumerate(recs):
        m = (r['lam'] >= b0) & (r['lam'] < b1)
        if not m.any():
            continue
        cols['lam'].append(r['lam'][m]); cols['sp'].append(np.full(int(m.sum()), si))
        cols['gf'].append(r['gf'][m]); cols['chi'].append(r['elo'][m] * CM_EV)
        cols['jl'].append(r['jl'][m]); cols['ju'].append(r['ju'][m])
        cols['ll'].append(r['lo'][m]); cols['ul'].append(r['up'][m])
        for t in tags:
            dp[t].append(r['dep'][t][m]); ds[t].append(r['dsyn'][t][m])
    if not cols['lam']:
        continue
    lam = np.concatenate(cols['lam']); o = np.argsort(lam)
    lam = lam[o]
    sp = np.concatenate(cols['sp'])[o]
    gf = np.concatenate(cols['gf'])[o]
    chi = np.concatenate(cols['chi'])[o]
    jl = np.concatenate(cols['jl'])[o]
    ju = np.concatenate(cols['ju'])[o]
    ll = np.concatenate(cols['ll'])[o]
    ul = np.concatenate(cols['ul'])[o]
    labs = sorted(set(ll) | set(ul))
    li = {v: k for k, v in enumerate(labs)}
    labbytes = '\n'.join(labs).encode()

    blob = b''.join([
        np.diff(np.rint(lam * 1e4).astype(np.int64),
                prepend=int(b0 * 1e4)).astype('<u4').tobytes(),
        sp.astype('u1').tobytes(),
        np.rint(gf * 1000).clip(-32768, 32767).astype('<i2').tobytes(),
        np.rint(chi * 1000).clip(0, 65535).astype('<u2').tobytes(),
        np.array([li[v] for v in ll], '<u2').tobytes(),
        np.array([li[v] for v in ul], '<u2').tobytes(),
        np.rint(jl * 2).clip(0, 255).astype('u1').tobytes(),
        np.rint(ju * 2).clip(0, 255).astype('u1').tobytes(),
        labbytes])
    n = len(lam)
    body = (b'LCAT' + np.array([n, len(labbytes)], '<u4').tobytes()
            + zlib.compress(blob, 9))
    (CAT / f'{int(b0)}.bin').write_bytes(body)
    cat_b += len(body)

    for t in tags:
        a = np.concatenate(dp[t])[o]
        b = np.concatenate(ds[t])[o]
        blob = (qdepth(a).tobytes()
                + np.rint(np.nan_to_num(b) * 255).clip(0, 255).astype('u1').tobytes())
        body = b'LDEP' + np.array([n], '<u4').tobytes() + zlib.compress(blob, 9)
        (DEP[t] / f'{int(b0)}.bin').write_bytes(body)
        dep_b[t] += len(body)
    index.append(int(b0))
    tot += n

(CAT / 'index.json').write_text(json.dumps(
    dict(block=BLOCK, blocks=index, species=allsp, cut=CUT, fmt='LCAT'),
    separators=(',', ':')))
print(f'\n{tot:,} lines in the union, {len(index)} blocks')
print(f'  shared catalogue {cat_b/1e6:.2f} MB  (web/data/lines/)')
for t in tags:
    print(f'  + {t:10s} depths {dep_b[t]/1e6:.2f} MB')
