"""Per-species emergent spectra: radiative transfer with one species' lines only.

Panels 1-3 come from SYNTHE.  This module answers the different question
"what would the Sun look like if only Fe I (or CH, or Ca II) absorbed?",
which needs a formal solution per species, not opacity at one depth.

The continuum opacity is rescaled wavelength by wavelength so the line-free
flux reproduces SYNTHE's own continuum: resline_lib's continuum carries H-,
H I, Rayleigh and Thomson only and runs ~11% low against ATLAS12, which also
has metal bound-free.  Without that calibration every line depth would be
biased through kappa_line / kappa_cont.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import wofz

sys.path.append(str(Path.home() / 'memos' / 'resline'))
import resline_lib as R
import sunlib
import moldisp

HERE = Path(__file__).resolve().parent
PI = np.pi


# ----------------------------------------------------------------------
def calibrated_continuum(atm, lam_lo, lam_hi, n=200):
    """(lam_coarse, kappa_cont(lam, depth)) scaled to SYNTHE's continuum flux."""
    lg = np.geomspace(lam_lo, lam_hi, n)
    kc = np.array([R.continuum_opacity(atm, la)[0] for la in lg])
    s = np.load(HERE / 'cache' / 'sun_synthe.npz')
    target = np.interp(lg, s['lam_vac'], s['cont']) / PI * 1e8
    logc = np.empty(n)
    for i in range(n):
        lo, hi = -1.5, 1.5
        for _ in range(45):
            mid = 0.5 * (lo + hi)
            f = R.emergent_flux(atm, atm.tau(kc[i] * np.exp(mid)),
                                np.array([lg[i]]))[0]
            lo, hi = (mid, hi) if f > target[i] else (lo, mid)
        logc[i] = 0.5 * (lo + hi)
    return lg, kc * np.exp(logc)[:, None], logc


def continuum_at(lg, kc_cal, lam):
    x, xp = np.log(lam), np.log(lg)
    return np.exp(np.array([np.interp(x, xp, np.log(kc_cal[:, j]))
                            for j in range(kc_cal.shape[1])]))


# ----------------------------------------------------------------------
def build_lines(phot, potion):
    """Flat line table + species table with n/U at every depth."""
    atm = phot.atm
    parts, labels, nU = [], [], []

    def add(label, arr):
        labels.append(label)
        nU.append(arr)
        return len(labels) - 1

    g = sunlib.parse_gfall_full()
    gr, gs, gw, z, ic = sunlib.synthe_damping(g, potion)
    ok = (z >= 1) & (z <= len(sunlib.PT)) & (ic <= 2) & (g['aut'] == 0)
    code = np.round(z + ic / 100.0, 2)
    sidx = np.full(len(code), -1)
    for c in np.unique(code[ok]):
        zz = int(np.floor(c + 1e-6))
        el, ion = sunlib.PT[zz - 1], int(round((c - zz) * 100))
        try:
            U = np.array([phot.ad.Q(el, t, ion) for t in atm.T])
        except (KeyError, IndexError):
            continue
        # Saha, NOT the .mol dump.  SYNTHE's synthesis takes its populations
        # from run_xnfpelsyn / COMPUTE_ONE_POP, and for the low-ionization
        # metals those disagree with the NMOLEC densities the .mol file
        # reports -- K by 46x, Na by 15x, Al 2.2x, Mg 1.7x at tau = 1.
        # Verified on Na D: Saha reproduces SYNTHE's profile to 0.005 in
        # residual flux, the .mol densities are off by a factor of 5 in the
        # wing.  Molecules keep the .mol densities; there is no alternative
        # for them, and their constituents (H, C, N, O) agree either way.
        if el not in atm.abund:
            continue
        n = atm.abund[el] * atm.nH * np.array(
            [phot.ion_fraction_at(el, ion, j) for j in range(len(atm.T))])
        sidx[code == c] = add(sunlib.species_name(c), n / U / atm.nH)
    m = ok & (sidx >= 0)
    parts.append(dict(lam=g['wl_nm'][m] * 10.0, elo=g['elo'][m],
                      gf=10.0**(g['loggf'][m] + g['x1'][m] + g['x2'][m]),
                      mass=np.array([sunlib.MASS.get(sunlib.PT[i - 1], 2.5 * i)
                                     for i in z[m]], dtype=float),
                      gr=gr[m], gs=gs[m], gw=gw[m], sidx=sidx[m]))

    mm = sunlib.parse_mol_lines()
    disp = moldisp.parse_dispatch()
    broad = sunlib.read_mol_broad()
    Tq, Q = sunlib.read_mol_partition()
    icode, iso = mm['icode'].astype(int), mm['iso'].astype(int)
    ok = np.array([c in sunlib.MOL_NAME and sunlib.MOL_NAME[c] in Q for c in icode])
    ok &= np.array([(i, c) in disp or (i, None) in disp for i, c in zip(iso, icode)])
    icode, iso = icode[ok], iso[ok]
    k = {key: v[ok] for key, v in mm.items()}
    sh = np.array([moldisp.isotope_shift(disp, i, c) for i, c in zip(iso, icode)])
    sidx = np.full(len(icode), -1)
    for c in np.unique(icode):
        n = phot.nden_all.get(float(c))
        if n is None:
            continue
        sidx[icode == c] = add(
            sunlib.MOL_NAME[c],
            n / np.interp(atm.T, Tq, Q[sunlib.MOL_NAME[c]]) / atm.nH)
    gr = np.where(k['loggr'] == 0.0, 2.223e13 / k['wl_nm']**2,
                  10.0**(k['loggr'] * 0.01))
    isX = k['upperX'] > 0
    gs = np.where(isX, 3.0e-8, 3.0e-5)
    gw = np.where(isX, 1.0e-8, 1.0e-7)
    for c, (g0, dg) in broad.items():
        sel = icode == c
        if sel.any() and g0 > 0:
            gw[sel] = np.maximum(g0 - dg * k['jlo'][sel], 0.1 * g0)
    m = sidx >= 0
    parts.append(dict(lam=k['wl_nm'][m] * 10.0, elo=k['elo'][m],
                      gf=10.0**(k['loggf'][m] + sh[m, 0] + sh[m, 1]),
                      mass=sh[m, 2].astype(float),
                      gr=gr[m], gs=gs[m], gw=gw[m], sidx=sidx[m]))

    L = {key: np.concatenate([p[key] for p in parts]) for key in parts[0]}
    L['nu0'] = R.C / (L['lam'] * 1e-8)
    return L, np.array(nU), labels


# ----------------------------------------------------------------------
class Engine:
    def __init__(self, phot, L, nU, lam):
        self.p, self.L, self.nU, self.lam = phot, L, nU, lam
        atm = phot.atm
        self.nd = len(atm.T)
        self.T, self.ne, self.vturb = atm.T, atm.ne, atm.vturb
        na = phot.nden_all
        self.txnxn = ((na[1.00] + 0.42 * na[2.00] + 0.85 * na[101.00])
                      * (atm.T / 1e4)**0.3)

    def block(self, k):
        L, T = self.L, self.T[None, :]
        nu0 = L['nu0'][k][:, None]
        boltz = np.exp(-L['elo'][k][:, None] / (0.695034800 * T))
        stim = 1.0 - np.exp(-R.H * nu0 / (R.KB * T))
        vD = np.sqrt(2.0 * R.KB * T / (L['mass'][k][:, None] * R.AMU)
                     + self.vturb[None, :]**2)
        dnuD = nu0 * vD / R.C
        gamma = (L['gr'][k][:, None] + L['gs'][k][:, None] * self.ne[None, :]
                 + L['gw'][k][:, None] * self.txnxn[None, :])
        area = R.PIE2MC * L['gf'][k][:, None] * self.nU[L['sidx'][k]] * boltz * stim
        return area, gamma / (4.0 * PI * dnuD), dnuD, vD


def windows(eng, kc_lg, kc_cal, ix, dln, floor=3e-4, hw_max=8192,
            probe=(24, 34, 44, 52, 58, 64)):
    """Half-window per line (grid points) and a keep mask, from the largest
    line-centre / continuum opacity ratio over a set of probe depths."""
    L = eng.L
    n = len(L['lam'])
    eta = np.empty(n)
    amax = np.empty(n)
    vref = np.empty(n)
    pj = list(probe)
    step = 150000
    for s in range(0, n, step):
        k = np.arange(s, min(s + step, n))
        area, a, dnuD, vD = eng.block(k)
        peak = area / (np.sqrt(PI) * dnuD) * np.real(wofz(1j * a))
        kcl = continuum_at(kc_lg, kc_cal, L['lam'][k])[pj]      # (nprobe, nsel)
        eta[k] = (peak[:, pj] / kcl.T).max(axis=1)
        amax[k] = a[:, pj].max(axis=1)
        vref[k] = vD[:, pj[-2]]
    vg = np.sqrt(np.maximum(np.log(np.maximum(eta / floor, 1.0)), 0.0))
    vl = np.sqrt(np.maximum(amax * eta / (np.sqrt(PI) * floor), 0.0))
    vcut = np.maximum(np.maximum(vg, vl), 3.5)
    pts = vcut * vref / R.C / dln
    hw = np.clip(2**np.ceil(np.log2(np.maximum(pts, 4.0))), 4, hw_max).astype(np.int64)
    return hw, eta > floor


def run(phot, eng, kc_lg, kc_cal, ix, hw, keep, gid, ngroup, names,
        chunk=6000, budget=45000, nmu=8, verbose=True):
    """Normalized flux per group, chunked in wavelength.

    Groups PARTITION the lines (everything unlisted lands in 'other'), so the
    all-lines spectrum is the transfer through the summed group opacity and
    costs no extra Voigt evaluations.
    """
    import time
    atm, lam, nd = phot.atm, eng.lam, eng.nd
    nlam = len(lam)
    res = {g: np.ones(nlam) for g in list(names) + ['all']}
    t0 = time.time()
    for c0 in range(0, nlam, chunk):
        c1 = min(c0 + chunk, nlam)
        lc, nc = lam[c0:c1], c1 - c0
        kcc = continuum_at(kc_lg, kc_cal, lc)
        f_c = R.emergent_flux(atm, atm.tau(kcc), lc, nmu=nmu)
        acc = np.zeros((ngroup, nd, nc))
        sel = np.where(keep & (ix + hw >= c0) & (ix - hw < c1))[0]
        for g in range(ngroup):
            sg = sel[gid[sel] == g]
            if not len(sg):
                continue
            for h in np.unique(hw[sg]):
                s2 = sg[hw[sg] == h]
                nb = max(1, int(budget // (2 * h + 1)))
                for s in range(0, len(s2), nb):
                    k = s2[s:s + nb]
                    idx = ix[k][:, None] + np.arange(-h, h + 1)[None, :] - c0
                    good = (idx >= 0) & (idx < nc)
                    idxc = np.clip(idx, 0, nc - 1)
                    area, a, dnuD, _ = eng.block(k)
                    nu = R.C / (lc[idxc] * 1e-8)
                    v = ((nu[:, None, :] - eng.L['nu0'][k][:, None, None])
                         / dnuD[:, :, None])
                    phi = (np.real(wofz(v + 1j * a[:, :, None]))
                           / (np.sqrt(PI) * dnuD[:, :, None]))
                    con = (area[:, :, None] * phi) * good[:, None, :]
                    fi = idxc.ravel()
                    for j in range(nd):
                        acc[g, j] += np.bincount(fi, weights=con[:, j, :].ravel(),
                                                 minlength=nc)
        for g, nm in enumerate(names):
            if acc[g].any():
                res[nm][c0:c1] = R.emergent_flux(
                    atm, atm.tau(kcc + acc[g]), lc, nmu=nmu) / f_c
        tot = acc.sum(axis=0)
        if tot.any():
            res['all'][c0:c1] = R.emergent_flux(
                atm, atm.tau(kcc + tot), lc, nmu=nmu) / f_c
        if verbose and (c0 // chunk) % 5 == 0:
            print(f'  {lam[c0]:8.1f} A  {100*c1/nlam:5.1f}%  {len(sel):6d} lines '
                  f'{time.time()-t0:6.0f}s', flush=True)
    return res
