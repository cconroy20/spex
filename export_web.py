"""Export the synthesis to a static web bundle.

Two tiers per series:
  ov/   overview, decimated, loaded up front so the default view is instant
  full/ native R = 300,000, fetched per series only when the user zooms in

Values are uint16 over a fixed range, which is 1.8e-5 in normalized flux --
far below anything visible, and half the size of float32.  The wavelength grid
is uniform in ln(lambda), so it is reconstructed in the browser from
(lam0_vac, dln, n) rather than stored.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.append(str(Path.home() / 'memos' / 'resline'))
import resline_lib as R
import binfmt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from starcfg import star, species, GROUP_ORDER                                             # noqa: E402

ST = star(sys.argv[1] if len(sys.argv) > 1 else 'sun')
OUT = HERE / 'web' / 'data' / ST['tag']
(OUT / 'ov').mkdir(parents=True, exist_ok=True)
(OUT / 'full').mkdir(parents=True, exist_ok=True)

NORM_MAX = 1.2                      # normalized flux quantization range
OV_BIN = 52                         # ~6000 points in the overview tier

# formation depth: where tau_lambda = 1, in both of the units worth reading it
# in.  Fixed quantization ranges so every star decodes the same way.
FTAU_LO, FTAU_HI = -8.0, 1.0        # log tau_5000
FTEMP_LO, FTEMP_HI = 0.0, 12000.0   # K

GROUPS = [(g, [n for _, n, gg, *_ in species() if gg == g])
          for g in GROUP_ORDER]


def q16(y, lo, hi):
    return np.clip((y - lo) / (hi - lo) * 65535.0, 0, 65535).astype('<u2')


def decimate(y, n=OV_BIN):
    """min, max and mean per bin.

    The envelope is what makes a zoomed-out spectrum honest: at 200 native
    points per screen pixel, plotting a mean washes the line forest away,
    while min/max shows the real depth range inside the pixel.  The mean is
    kept for the smoothed display modes.
    """
    m = (len(y) // n) * n
    b = y[:m].reshape(-1, n)
    return np.concatenate([b.min(axis=1), b.max(axis=1), b.mean(axis=1)])


d = np.load(HERE / 'cache' / f'{ST["tag"]}_species.npz')
s = np.load(HERE / 'cache' / f'{ST["tag"]}_synthe.npz')
lam_vac = d['lam_vac']
dln = float(np.median(np.diff(np.log(lam_vac))))
n = len(lam_vac)
air = R.vac_to_air(lam_vac)

series, meta_series = {}, []
for gname, members in GROUPS:
    for k in members:
        kk = k.replace(' ', '_')
        if kk not in d.files:
            continue
        series[k] = d[kk].astype(float)
        meta_series.append(dict(name=k, group=gname, kind='mol' if gname == 'molecules'
                                else 'atom',
                                absorption=float(1 - d[kk].mean())))

# panels (a) and (b), resampled onto the species grid (the assembled spectrum
# has a handful of one-point gaps at chunk joins)
flux = np.interp(air, s['lam_air'], s['flux'])
cont = np.interp(air, s['lam_air'], s['cont'])
norm = np.interp(air, s['lam_air'], s['norm'])
FLUX_MAX = float(np.nanmax(cont) * 1.02)

written = []
n_ov = len(norm) // OV_BIN


def write(name, y, lo, hi):
    binfmt.write_series(OUT / 'full' / f'{name}.bin', q16(y, lo, hi))
    # one chunk per overview channel, so a difference never crosses from the
    # min channel into the max channel
    binfmt.write_series(OUT / 'ov' / f'{name}.bin', q16(decimate(y), lo, hi),
                        chunk=n_ov)
    written.append(name)


for k, y in series.items():
    write(k.replace(' ', '_'), y, 0.0, NORM_MAX)
# _flux is not written: it is _norm x _cont, and reconstructing it in the
# browser from the two quantized arrays agrees with the real thing to two
# quantization steps -- below anything that can be drawn
write('_cont', cont, 0.0, FLUX_MAX)
write('_norm', norm, 0.0, NORM_MAX)

fp = HERE / 'cache' / f'{ST["tag"]}_form.npz'
have_form = fp.exists()
if have_form:
    fz = np.load(fp)
    fa = R.vac_to_air(fz['lam_vac'])
    write('_ftau', np.interp(air, fa, fz['logtau5000']), FTAU_LO, FTAU_HI)
    write('_ftemp', np.interp(air, fa, fz['tform']), FTEMP_LO, FTEMP_HI)
else:
    print('  (no formation depths yet)')

meta = dict(
    n=n, n_ov=len(norm) // OV_BIN, ov_bin=OV_BIN, ov_channels=3,
    lam0_vac=float(lam_vac[0]), dln=dln,
    lam_air_min=float(air[0]), lam_air_max=float(air[-1]),
    norm_max=NORM_MAX, flux_max=FLUX_MAX,
    fmt='SPC1', chunk=binfmt.CHUNK, derived=['_flux'],
    form=have_form,
    qrange={'_ftau': [FTAU_LO, FTAU_HI], '_ftemp': [FTEMP_LO, FTEMP_HI]},
    star=ST['tag'], name=ST['name'],
    model=dict(teff=ST['teff'], logg=ST['logg'], feh=ST['feh'],
               afe=ST.get('afe', 0.0), name='ATLAS12 / SYNTHE, Kurucz'),
    groups=[g for g, _ in GROUPS], series=meta_series,
    # what to show on first load: the species that actually do something in
    # THIS star, which is a different list for a G dwarf and a K giant
    default_on=[m['name'] for m in sorted(meta_series,
                                          key=lambda m: -m['absorption'])[:14]])
# drop anything left from an earlier export with a different series list --
# _flux stopped being written and its 621 kB files would otherwise linger
keep = {f'{n}.bin' for n in written}
for sub in ('full', 'ov'):
    for f in (OUT / sub).glob('*.bin'):
        if f.name not in keep:
            f.unlink()
            print(f'  removed stale {sub}/{f.name}')

(OUT / 'meta.json').write_text(json.dumps(meta, indent=1))
print(f'{len(written)} series, n={n} (overview {meta["n_ov"]} x3), '
      f'grid {air[0]:.1f}-{air[-1]:.1f} A air')
tot = sum(f.stat().st_size for f in OUT.rglob('*'))
print(f'bundle {tot/1e6:.1f} MB   (overview tier '
      f'{sum(f.stat().st_size for f in (OUT/"ov").glob("*"))/1e3:.0f} kB)')
