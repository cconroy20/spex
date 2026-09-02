"""Pedagogical figures: the solar spectrum from 3800 to 10000 A."""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

from plotstyle import apply_style, panel_label, BLUE, AMBER, PURPLE, TEAL, INK, GREY

sys.path.append(str(Path.home() / 'memos' / 'resline'))
import resline_lib as R

apply_style()
HERE = Path(__file__).resolve().parent
FIG = HERE / 'fig'
FIG.mkdir(exist_ok=True)
# 3600 A takes in the Balmer break at 3646 A; the synthesis runs to 3550
LO, HI = 3600.0, 10000.0

# Per-species spectra: one SYNTHE run per species, which has SYNTHE's own
# profile physics throughout -- in particular the Holtsmark quasi-static Stark
# profile for hydrogen, which a Voigt cannot reproduce.  species_flux_synthe
# is the 25-species 3550-10000 A set the committed figures were drawn from;
# sun_species is what run_star.sh + assemble_species.py write today, over
# 3550-25000 A, and the panels are cut to LO-HI either way.  specflux's own
# transfer (species_flux.npz) is the cross-check, not a figure source.
SPECIES_CACHE = next(
    (f for f in ('species_flux_synthe.npz', 'sun_species.npz')
     if (HERE / 'cache' / f).exists()), 'sun_species.npz')

s = np.load(HERE / 'cache' / 'sun_synthe.npz')
LAM, FLUX, CONT, NORM = s['lam_air'], s['flux'], s['cont'], s['norm']


def binned(x, y, nb, how='mean'):
    e = np.linspace(LO, HI, nb + 1)
    # Restrict to the plotted band first.  digitize + clip alone puts every
    # point outside it into the first or last bin, and the synthesis now runs
    # to 2.5 um: 48% of the grid would land in the last bin, which read as a
    # cliff at 10000 A rather than as data.
    m = (x >= LO) & (x <= HI)
    x, y = x[m], y[m]
    i = np.clip(np.digitize(x, e) - 1, 0, nb - 1)
    c = np.bincount(i, minlength=nb).astype(float)
    if how == 'mean':
        v = np.bincount(i, weights=y, minlength=nb) / np.maximum(c, 1)
    else:
        v = np.ones(nb)
        np.minimum.at(v, i, y)
    v[c == 0] = np.nan
    return 0.5 * (e[1:] + e[:-1]), v


# ======================================================================
# Figure 1 -- flux, continuum, normalized flux
# ======================================================================
def figure_overview():
    fig, ax = plt.subplots(2, 1, figsize=(10.5, 7.0), sharex=True)
    NB = 1550
    xb, fb = binned(LAM, FLUX / 1e7, NB)
    xc, cb = binned(LAM, CONT / 1e7, NB)
    ax[0].plot(xb, fb, lw=0.6, color=INK, label='emergent flux')
    ax[0].plot(xc, cb, lw=1.7, color=AMBER, label='continuum (no lines)')
    ax[0].set_ylabel(r'$F_\lambda$  [$10^{7}$ erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$]')
    ax[0].set_ylim(0, 1.35)
    ax[0].legend(loc='upper right')
    panel_label(ax[0], '(a)')

    xn, nb_ = binned(LAM, NORM, NB)
    ax[1].plot(xn, nb_, lw=0.6, color=INK)
    ax[1].axhline(1.0, lw=0.9, color=AMBER, zorder=0)
    ax[1].set_ylabel(r'$F_\lambda / F_\lambda^{\rm cont}$')
    ax[1].set_ylim(0, 1.06)
    ax[1].set_xlabel(r'wavelength [$\mathrm{\AA}$]')
    panel_label(ax[1], '(b)', 'lower right')
    for a in ax:
        a.set_xlim(LO, HI)
    fig.tight_layout(h_pad=0.5)
    fig.savefig(FIG / 'sun_overview.pdf', bbox_inches='tight')
    print('wrote fig/sun_overview.pdf')


# ======================================================================
# Figure 2 -- the spectrum species by species
# ======================================================================
# No 'all' row: panel (b) of the poster is already the full synthesis, and the
# two were the same array (rms 1.4e-4).  The freed row went to Mn I.
#
# Panels are grouped by family so the physics reads down the page, and ordered
# by contribution WITHIN each group, so the strongest member leads it.  Ca II
# sits with the neutrals rather than with the ions: H and K plus the infrared
# triplet make it a headline feature, and it would be buried otherwise.
SPECIES_GROUPS = [
    ['Fe I', 'Ca II', 'H I', 'Mg I', 'Ca I', 'Si I', 'Na I'],
    ['Cr I', 'Ni I', 'Ti I', 'Mn I'],
    ['Ti II', 'Fe II'],
    ['CH', 'CN', 'MgH'],
]


def species_order(d):
    """Flatten SPECIES_GROUPS, sorting each group by mean line absorption."""
    out = []
    for g in SPECIES_GROUPS:
        have = [k for k in g if k.replace(' ', '_') in d.files]
        out += sorted(have, key=lambda k: d[k.replace(' ', '_')].mean())
    return out


SHOW = [k for g in SPECIES_GROUPS for k in g]


def figure_species(nb=1550):
    d = np.load(HERE / 'cache' / SPECIES_CACHE)
    lam = R.vac_to_air(d['lam_vac'])
    have = ['all'] + species_order(d)
    fig, ax = plt.subplots(len(have), 1, figsize=(10.5, 1.05 * len(have) + 1.4),
                           sharex=True, sharey=True)
    for a, k in zip(ax, have):
        x, y = binned(lam, d[k.replace(' ', '_')], nb)
        a.fill_between(x, 1.0, np.minimum(y, 1.0), color=BLUE, lw=0)
        a.set_ylim(0.0, 1.06)
        a.set_yticks([0.0, 0.5, 1.0])
        a.axhline(1.0, lw=0.7, color=AMBER, zorder=0)
        lab = 'all lines' if k == 'all' else k
        a.annotate(lab, xy=(1, 0), xycoords='axes fraction',
                   textcoords='offset points', xytext=(-6, 5), ha='right',
                   va='bottom', fontsize=11.5,
                   bbox=dict(fc='white', ec='none', alpha=0.72, pad=1.2))
    ax[-1].set_xlabel(r'wavelength [$\mathrm{\AA}$]')
    ax[len(have) // 2].set_ylabel(r'$F_\lambda / F_\lambda^{\rm cont}$')
    ax[0].set_xlim(LO, HI)
    fig.tight_layout(h_pad=0.15)
    fig.savefig(FIG / 'sun_species.pdf', bbox_inches='tight')
    print(f'wrote fig/sun_species.pdf  ({len(have)} panels)')


# ======================================================================
# Figure 3 -- zooms, where the lines are actually resolved
# ======================================================================
WINDOWS = [(3925, 3975, 'Ca II H and K'), (4290, 4320, 'CH G band'),
           (5160, 5195, 'Mg b'), (5885, 5900, 'Na D'),
           (6550, 6575, r'H$\alpha$'), (8480, 8680, 'Ca II IR triplet')]


def figure_zoom():
    d = np.load(HERE / 'cache' / SPECIES_CACHE)
    lam = R.vac_to_air(d['lam_vac'])
    keys = [k for k in d.files if k not in ('lam_vac', 'all', 'other')]
    cols = [BLUE, AMBER, PURPLE]
    fig, axs = plt.subplots(3, 2, figsize=(11.0, 8.2))
    for a, (lo, hi, name) in zip(axs.ravel(), WINDOWS):
        m = (lam > lo) & (lam < hi)
        x = lam[m]
        try:
            wo, fo = R.read_iag(lo, hi)
            if len(wo) > 20:
                a.plot(wo, fo, lw=0.9, color=GREY, zorder=1, label='observed Sun')
        except Exception:
            pass
        rank = sorted(keys, key=lambda k: d[k][m].mean())[:3]
        for c, k in zip(cols, rank):
            if 1.0 - d[k][m].mean() < 0.004:
                continue
            a.plot(x, d[k][m], lw=1.2, color=c, zorder=2, label=k.replace('_', ' '))
        a.plot(x, d['all'][m], lw=0.9, color=INK, zorder=3, label='all lines')
        a.set_xlim(lo, hi)
        a.set_ylim(0, 1.62)
        a.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        h, l = a.get_legend_handles_labels()
        o = [l.index('all lines')] + [i for i, q in enumerate(l) if q != 'all lines']
        a.legend([h[i] for i in o], [l[i] for i in o], loc='upper left', ncol=3,
                 handlelength=1.3, columnspacing=1.0, borderaxespad=0.4)
        a.annotate(name, xy=(1, 1), xycoords='axes fraction',
                   textcoords='offset points', xytext=(-6, -6), ha='right',
                   va='top', fontsize=11.5)
    for a in axs[-1]:
        a.set_xlabel(r'wavelength [$\mathrm{\AA}$]')
    for a in axs[:, 0]:
        a.set_ylabel(r'$F_\lambda / F_\lambda^{\rm cont}$')
    fig.tight_layout(h_pad=0.8, w_pad=1.2)
    fig.savefig(FIG / 'sun_zoom.pdf', bbox_inches='tight')
    print('wrote fig/sun_zoom.pdf')


if __name__ == '__main__':
    figure_overview()
    if (HERE / 'cache' / SPECIES_CACHE).exists():
        figure_species()
        figure_zoom()


# ======================================================================
# Poster -- overview panels stacked on top of the species waterfall
# ======================================================================
MOLECULES = {'CH', 'CN', 'C2', 'MgH', 'SiH', 'NH', 'OH', 'CO', 'SiO'}


def smooth_to_R(lam, y, R):
    """Convolve with a Gaussian of FWHM = c/R.

    Both grids are uniform in ln(lambda), so a fixed-width kernel in index
    space is a fixed-width kernel in velocity.  The intrinsic sampling is
    R = 300,000, which adds in quadrature at the 0.06% level and is ignored.
    """
    dln = float(np.median(np.diff(np.log(lam))))
    return gaussian_filter1d(y, (1.0 / R) / dln / 2.35482, mode='nearest')


TITLE = (r'Solar Model:  $T_{\rm eff} = 5777$ K,  $\log g = 4.44$')


def figure_poster(width=30.811, row=1.5719, nb=60000, fs=2.4, lw=0.6,
                  smooth=None, title=TITLE, title_fs=11.0, fill=False):
    """One large vector figure sized for printing.

    Vector output, so there is no resolution to set -- what matters is the
    physical size and that the type scales with it.  The spectra are binned
    The defaults trim to exactly 30 x 40 in (3:4), a standard poster size, so
    it prints without scaling.  (`width` and `row` are solved for that trimmed
    size, so they move if the title size or the row count changes.)  The spectra are binned
    to `nb` points -- 0.10 A, about 2000 samples per inch at this width,
    which is well past what any printer resolves while keeping the file
    tractable.  Solar lines are ~0.15 A wide, so nothing is smoothed away.
    """
    d = np.load(HERE / 'cache' / SPECIES_CACHE)
    lam = R.vac_to_air(d['lam_vac'])
    have = species_order(d)
    nrow = 2 + len(have)

    # Flux and continuum are smoothed separately and then divided, which is
    # what an instrument does; for the species panels only the ratio is on
    # disk, but the continuum varies by <0.1% over one resolution element at
    # R = 10,000, so smoothing the ratio is the same thing.
    flux, cont, norm = FLUX, CONT, NORM
    spec = {k: d[k.replace(' ', '_')] for k in have}
    if smooth:
        flux = smooth_to_R(s['lam_vac'], FLUX, smooth)
        cont = smooth_to_R(s['lam_vac'], CONT, smooth)
        norm = flux / cont
        spec = {k: smooth_to_R(d['lam_vac'], v, smooth) for k, v in spec.items()}
    height = (4.5 + 3.5) * row + len(have) * row + 2.2

    with plt.rc_context({
            'font.size': 13 * fs, 'axes.labelsize': 14 * fs,
            'xtick.labelsize': 12 * fs, 'ytick.labelsize': 12 * fs,
            'legend.fontsize': 12 * fs, 'axes.linewidth': lw * fs,
            'xtick.major.size': 3.5 * fs, 'ytick.major.size': 3.5 * fs,
            'xtick.minor.size': 2.0 * fs, 'ytick.minor.size': 2.0 * fs,
            'xtick.major.width': lw * fs, 'ytick.major.width': lw * fs,
            'xtick.minor.width': lw * fs, 'ytick.minor.width': lw * fs,
            # embed TrueType rather than Type 3 -- print shops choke on Type 3
            'pdf.fonttype': 42, 'ps.fonttype': 42}):
        fig, ax = plt.subplots(
            nrow, 1, figsize=(width, height), sharex=True,
            gridspec_kw=dict(height_ratios=[4.5, 3.5] + [1.0] * len(have)))

        xb, fb = binned(LAM, flux / 1e7, nb)
        xc, cb = binned(LAM, cont / 1e7, nb)
        ax[0].plot(xb, fb, lw=lw * fs, color=INK, label='emergent flux')
        ax[0].plot(xc, cb, lw=lw * fs, color=AMBER, label='continuum (no lines)')
        ax[0].set_ylabel('$F_\\lambda$\n'
                         r'[$10^{7}$ erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$]',
                         labelpad=10 * fs, fontsize=11 * fs, linespacing=1.6)
        ax[0].set_ylim(0, 1.35)

        xn, yn = binned(LAM, norm, nb)
        ax[1].plot(xn, yn, lw=lw * fs, color=INK)
        ax[1].set_ylabel(r'$F_\lambda / F_\lambda^{\rm cont}$', labelpad=10 * fs)
        ax[1].set_ylim(0, 1.06)

        for a, k in zip(ax[2:], have):
            x, y = binned(lam, spec[k], nb)
            col = AMBER if k in MOLECULES else BLUE
            yy = np.minimum(y, 1.0)
            # drawn as a line, like panels (a) and (b): a shaded region reads
            # as a dark mass and makes the same data look unlike the
            # normalized panel above it.  The stroke weight is what keeps a
            # single-bin line from being a hairline on paper.
            if fill:
                a.fill_between(x, 1.0, yy, color=col, lw=0)
            a.plot(x, yy, lw=lw * fs, color=col, solid_joinstyle='round')
            a.set_ylim(0.0, 1.08)
            a.set_yticks([0.0, 0.5, 1.0])
            a.annotate('all lines' if k == 'all' else k, xy=(1, 0),
                       xycoords='axes fraction', textcoords='offset points',
                       xytext=(-7 * fs, 5 * fs), ha='right', va='bottom',
                       fontsize=13 * fs,
                       bbox=dict(fc='white', ec='none', alpha=0.75, pad=2.0))
        ax[2 + len(have) // 2].set_ylabel(r'$F_\lambda / F_\lambda^{\rm cont}$',
                                          labelpad=10 * fs)
        ax[-1].set_xlabel(r'wavelength [$\mathrm{\AA}$]', labelpad=8 * fs)
        ax[0].set_xlim(LO, HI)
        ax[0].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(500))
        ax[0].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(100))
        if smooth:
            ax[0].annotate(f'$R = {smooth:,}$'.replace(',', '{,}'), xy=(1, 1),
                           xycoords='axes fraction', textcoords='offset points',
                           xytext=(-8 * fs, -8 * fs), ha='right', va='top',
                           fontsize=13 * fs)
        fig.tight_layout(h_pad=0.25 * fs)
        if title:
            # placed above the axes area rather than inside it; bbox_inches
            # 'tight' expands the page to take it in
            fig.suptitle(title, y=1.0, va='bottom', fontsize=title_fs * fs)
        name = 'sun_poster.pdf' if not smooth else f'sun_poster_R{smooth//1000}k.pdf'
        fig.savefig(FIG / name, bbox_inches='tight')
    print(f'wrote fig/{name}  ({width:.0f} x {height:.1f} in, '
          f'{nrow} panels, {nb} bins'
          + (f', smoothed to R = {smooth}' if smooth else '') + ')')
