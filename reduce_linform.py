"""Collapse a SYNTHE .linform to one formation depth per wavelength.

    python3 reduce_linform.py <chunkdir> <stem> <model.atm> [--keep]

.linform carries the monochromatic optical depth at every atmospheric layer
for every wavelength -- 909 bytes a wavelength, half a gigabyte a star over
the full band.  Almost all of that is only needed to answer one question:
where does this wavelength form?  So each record is reduced to the layer where
tau_lambda = 1, and both of the things you would want to read off there are
kept: the reference optical depth log tau_5000, which compares across stars,
and the temperature, which is what the line actually responds to.

The .linform is deleted afterwards unless --keep is given.
"""
import sys
from pathlib import Path

import numpy as np


def read_deck(model):
    """Layer count, T, and tau_5000 if the deck carries it."""
    txt = Path(model).read_text().split('\n')
    rows, n, reading, cols = [], None, False, []
    for line in txt:
        if line.startswith('READ DECK'):
            n = int(line.split()[2])
            cols = line.replace(',', ' ').split()[3:]
            reading = True
            continue
        if reading:
            p = line.split()
            if len(p) < 7:
                break
            rows.append([float(x) for x in p])
            if len(rows) == n:
                break
    w = max(len(r) for r in rows)
    a = np.array([r + [np.nan] * (w - len(r)) for r in rows])
    T = a[:, 1]
    t5 = a[:, cols.index('TAU5000')] if 'TAU5000' in cols and a.shape[1] > cols.index('TAU5000') else None
    return n, T, t5


def tau5000_from_deck(model):
    """For decks written before TAU5000 was a column: integrate it."""
    sys.path.insert(0, str(Path.home() / 'memos' / 'resline'))
    import resline_lib as R
    atm = R.Atmosphere(model)
    kc = R.continuum_opacity(atm, 5000.0)[0]
    return atm.tau(kc)


def main():
    d, stem, model = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    keep = '--keep' in sys.argv
    f = d / f'{stem}.linform'
    if not f.exists():
        print(f'  no linform in {d}')
        return
    n, T, t5 = read_deck(model)
    if t5 is None:
        t5 = tau5000_from_deck(model)
    per = 1 + -(-n // 10)                       # header + ceil(n/10) data lines

    txt = f.read_text().split('\n')
    while txt and not txt[-1].strip():
        txt.pop()
    nrec = len(txt) // per
    assert nrec * per == len(txt), f'{f}: {len(txt)} lines, {per} per record'

    lam = np.empty(nrec)
    tau = np.empty((nrec, n))
    for i in range(nrec):
        blk = txt[i * per:(i + 1) * per]
        lam[i] = float(blk[0].split()[0])
        v = [float(x) for l in blk[1:] for x in l.split()]
        tau[i] = v[:n]

    # where tau_lambda passes 1, interpolated in log tau against layer index
    lt = np.log10(np.maximum(tau, 1e-30))
    ge = tau >= 1.0
    j = np.clip(np.argmax(ge, axis=1), 1, n - 1)
    a0 = lt[np.arange(nrec), j - 1]
    a1 = lt[np.arange(nrec), j]
    frac = np.clip((0.0 - a0) / np.where(a1 - a0 == 0, 1e-30, a1 - a0), 0.0, 1.0)
    pos = np.where(ge.any(axis=1), (j - 1) + frac, n - 1.0)

    idx = np.arange(n)
    logt5 = np.interp(pos, idx, np.log10(np.maximum(t5, 1e-30)))
    Tform = np.interp(pos, idx, T)

    np.savez_compressed(d / 'form.npz', lam_vac=lam * 10.0,
                        logtau5000=logt5.astype('f4'), tform=Tform.astype('f4'))
    if not keep:
        f.unlink()
    print(f'  {d.name}: {nrec} wavelengths, log tau5000 '
          f'{logt5.min():.2f}..{logt5.max():.2f}, T {Tform.min():.0f}..{Tform.max():.0f} K')


main()
