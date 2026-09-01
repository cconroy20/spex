"""Per-species line opacity in the solar photosphere, on the SYNTHE wavelength grid.

Panels 1-3 of the figure set (fluxed spectrum, continuum, normalized spectrum)
come from SYNTHE itself.  This module supplies panel 4: the line opacity of
each atomic species at a single reference depth, computed from the same
primary data SYNTHE uses --

  atmosphere    ~/kurucz/grids/THESUN/atm/ap00t5777g4.44at12.dat  (ATLAS12)
  line data     ~/kurucz/atlas12/data/gfallvac08oct17.dat         (Kurucz)
  partition fn  ~/kurucz/upgrade/raw_data/bc16_table8_vNov2022.dat (Barklem & Collet 2016)
  ion. energy   ~/kurucz/upgrade/raw_data/bc16_table4.dat
  isotopes      ~/kurucz/atlas12/data/isotopes.dat

The atmosphere, continuum opacity and Saha-Boltzmann machinery are reused
from ~/memos/resline/resline_lib.py, which was validated against the model's
own Rosseland mean and against ten solar equivalent widths.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.special import wofz

sys.path.insert(0, str(Path.home() / 'memos' / 'resline'))
sys.path.insert(0, str(Path.home() / 'kurucz' / 'atlas12' / 'tools'))

import resline_lib as R                                              # noqa: E402
from atomic_saha import AtomicData                                   # noqa: E402

HERE = Path(__file__).resolve().parent
GFALL_OPT = HERE / 'cache' / 'gfall_opt.dat'
ISOTOPES = Path.home() / 'kurucz' / 'atlas12' / 'data' / 'isotopes.dat'

PT = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
      'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
      'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',
      'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
      'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm',
      'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W',
      'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
      'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U']
ROMAN = ['I', 'II', 'III']


def isotope_masses():
    """Mass of the most abundant isotope of each element [amu], from ATLAS12.

    isotopes.dat holds ISOION(20, 265): for each element, ten isotope mass
    numbers followed by ten number fractions.  ATLAS12 takes the Doppler
    width from the most abundant isotope, so use the same choice here.
    """
    vals = []
    for line in ISOTOPES.read_text().split('\n'):
        if line.startswith('#') or not line.strip():
            continue
        vals.extend(float(x) for x in line.split())
    a = np.array(vals[:20 * 265]).reshape(265, 20)
    out = {}
    for z in range(1, len(PT) + 1):
        m, f = a[z - 1, :10], a[z - 1, 10:]
        if f.sum() <= 0:
            continue
        out[PT[z - 1]] = float(m[np.argmax(f)])
    return out


MASS = isotope_masses()


# ----------------------------------------------------------------------
# line list
# ----------------------------------------------------------------------
def parse_gfall(path=GFALL_OPT, cache=None):
    """Fixed-width gfall records -> numpy arrays.  Vacuum nm throughout."""
    cache = Path(cache) if cache else path.with_suffix('.npz')
    if cache.exists():
        d = np.load(cache)
        return {k: d[k] for k in d.files}

    def f(s, default=np.nan):
        try:
            return float(s)
        except ValueError:
            return default

    wl, gf, code, elo, gr, gs, gw = [], [], [], [], [], [], []
    with open(path, errors='replace') as fh:
        for line in fh:
            c = f(line[18:24])
            if c != c:
                continue
            # The two level blocks are not in energy order and a negative
            # energy marks a predicted level; the lower level is the one
            # with the smaller |E|.
            e1, e2 = abs(f(line[24:36], 0.0)), abs(f(line[52:64], 0.0))
            wl.append(f(line[0:11]))
            gf.append(f(line[11:18], -99.0))
            code.append(c)
            elo.append(min(e1, e2))
            gr.append(f(line[80:86], 0.0))
            gs.append(f(line[86:92], 0.0))
            gw.append(f(line[92:98], 0.0))
    out = dict(wl_nm=np.array(wl), loggf=np.array(gf), code=np.array(code),
               elo_cm=np.array(elo), loggr=np.array(gr),
               loggs=np.array(gs), loggw=np.array(gw))
    for k in out:
        out[k] = np.nan_to_num(out[k], nan=0.0)
    np.savez_compressed(cache, **out)
    return out


# ----------------------------------------------------------------------
# opacity at one depth
# ----------------------------------------------------------------------
class Photosphere:
    """The model atmosphere plus everything needed at one reference depth."""

    def __init__(self, lam_ref_A=5000.0, tau_ref=1.0, atm_path=None):
        self.atm = R.Atmosphere(atm_path) if atm_path else R.Atmosphere()
        self.ad = AtomicData()
        kc = R.continuum_opacity(self.atm, lam_ref_A)[0]
        self.tau_ref_scale = self.atm.tau(kc)
        self.j = int(np.argmin(abs(np.log10(self.tau_ref_scale) - np.log10(tau_ref))))
        self.T = float(self.atm.T[self.j])
        self.ne = float(self.atm.ne[self.j])
        self.nH = float(self.atm.nH[self.j])
        self.vturb = float(self.atm.vturb[self.j])
        self.mu_H = float(self.atm.mu_H)
        self._ionfrac = {}

    def ion_fraction(self, el, stage):
        key = (el, stage)
        if key not in self._ionfrac:
            T, ne = self.T, self.ne
            chi = self.ad._ie[el]
            U = [self.ad.Q(el, T, i) for i in (0, 1, 2)]
            r1 = (2.0 * R.SAHA_TR * T**1.5 * (U[1] / U[0])
                  * np.exp(-chi[0] / (R.KB_EV * T)) / ne)
            r2 = (2.0 * R.SAHA_TR * T**1.5 * (U[2] / U[1])
                  * np.exp(-chi[1] / (R.KB_EV * T)) / ne)
            den = 1.0 + r1 + r1 * r2
            for i, v in enumerate([1.0 / den, r1 / den, r1 * r2 / den]):
                self._ionfrac[(el, i)] = float(v)
        return self._ionfrac[key]

    def continuum(self, lam_A):
        """Continuum opacity per gram at the reference depth, vectorized in lam."""
        out = np.empty(len(lam_A))
        for i, la in enumerate(lam_A):
            out[i] = R.continuum_opacity(self.atm, la)[0][self.j]
        return out / self.mu_H


IONPOTS = Path.home() / 'kurucz' / 'atlas12' / 'data' / 'ionpots.dat'
MOLDUMP = HERE / 'run' / 'moltest' / 'ap00t5777g4.44at12.mol'


def read_ionpots(path=IONPOTS):
    """POTION(999): ionization potentials [cm^-1] on Kurucz's species index."""
    v = []
    for line in path.read_text().split('\n'):
        if line.startswith('#') or not line.strip():
            continue
        v.extend(float(x) for x in line.split())
    return np.array(v)


def potion_index(z, icharge):
    """Kurucz's POTION index for element z in charge state icharge (0 = neutral)."""
    return np.where(z <= 30, z * (z + 1) // 2 + icharge, z * 5 + 341 + icharge)


def read_moldump(path=MOLDUMP):
    """SYNTHE's own number densities [cm^-3] per species code, per depth."""
    txt = path.read_text().split('\n')
    out, i = {}, 0
    while i < len(txt):
        if txt[i].strip().startswith('J '):
            codes = [float(x) for x in txt[i].split()[1:]]
            rows = []
            i += 1
            while i < len(txt) and txt[i].strip() and not txt[i].strip().startswith('J '):
                p = txt[i].split()
                if len(p) != len(codes) + 1:
                    break
                rows.append([float(x) for x in p[1:]])
                i += 1
            a = np.array(rows)
            for k, c in enumerate(codes):
                out[round(c, 2)] = a[:, k]
        else:
            i += 1
    return out


def parse_gfall_full(path=GFALL_OPT):
    """gfall records with everything SYNTHE's reader uses."""
    cache = path.parent / 'gfall_full.npz'
    if cache.exists():
        d = np.load(cache)
        return {k: d[k] for k in d.files}

    def f(s, default=0.0):
        try:
            return float(s)
        except ValueError:
            return default

    cols = {k: [] for k in ('wl_nm', 'loggf', 'code', 'elo', 'eup',
                            'loggr', 'loggs', 'loggw', 'x1', 'x2', 'aut')}
    with open(path, errors='replace') as fh:
        for line in fh:
            c = f(line[18:24], -1.0)
            if c <= 0:
                continue
            flag = line[141:144] if len(line) > 144 else '   '
            if flag == 'COR':                     # SYNTHE skips these
                continue
            e1, e2 = abs(f(line[24:36])), abs(f(line[52:64]))
            cols['wl_nm'].append(f(line[0:11]))
            cols['loggf'].append(f(line[11:18], -99.0))
            cols['code'].append(c)
            cols['elo'].append(min(e1, e2))
            cols['eup'].append(max(e1, e2))
            cols['loggr'].append(f(line[80:86]))
            cols['loggs'].append(f(line[86:92]))
            cols['loggw'].append(f(line[92:98]))
            cols['x1'].append(f(line[109:115]))
            cols['x2'].append(f(line[118:124]))
            cols['aut'].append(1.0 if flag == 'AUT' else 0.0)
    out = {k: np.array(v) for k, v in cols.items()}
    np.savez_compressed(cache, **out)
    return out


def synthe_damping(g, potion):
    """gamma_rad, gamma_stark, gamma_vdW per line, with SYNTHE's defaults.

    Matches mod_mklinelist.f90 read_gfall: classical radiative damping when
    gr is absent, Cowley 1e-8 n*^5 for Stark, and the Unsold 4.5e-9 form for
    van der Waals, with Kurucz's iron-group override for Sc-Ni.
    """
    z = np.floor(g['code'] + 1e-6).astype(int)
    icharge = np.rint((g['code'] - z) * 100).astype(int)
    zeff = icharge + 1.0

    gr = np.where(np.abs(g['loggr']) < 1e-6,
                  2.223e13 / g['wl_nm']**2, 10.0**g['loggr'])

    idx = potion_index(z, icharge)
    idx = np.clip(idx, 1, len(potion))
    chi = potion[idx - 1]

    def effnsq(e, default):
        d = chi - e
        return np.where(d > 0.0, 109737.31 * zeff**2 / np.maximum(d, 1e-30), default)

    nsq_up = np.minimum(effnsq(g['eup'], 25.0), 1000.0)
    gs = np.where(np.abs(g['loggs']) < 1e-6,
                  1e-8 * nsq_up**2 * np.sqrt(nsq_up), 10.0**g['loggs'])

    nsq_lo = np.minimum(np.where(chi - g['elo'] > 0.0,
                                 109737.31 * zeff**2 / np.maximum(chi - g['elo'], 1e-30),
                                 25.0), 1000.0)
    rsq_up = 2.5 * (np.minimum(effnsq(g['eup'], 25.0), 1000.0) / zeff)**2
    rsq_lo = 2.5 * (nsq_lo / zeff)**2
    iron = (z > 20) & (z < 29)                       # Kurucz's Sc-Ni override
    rsq_up = np.where(iron, (45.0 - z) / zeff, rsq_up)
    rsq_lo = np.where(iron, 0.0, rsq_lo)
    rsq_up = np.where(rsq_up < rsq_lo, 2.0 * rsq_lo, rsq_up)
    gw = np.where(np.abs(g['loggw']) < 1e-6,
                  4.5e-9 * np.maximum(rsq_up - rsq_lo, 0.0)**0.4, 10.0**g['loggw'])
    return gr, gs, gw, z, icharge


def species_name(code):
    z = int(np.floor(code + 1e-6))
    ion = int(round((code - z) * 100))
    if z <= len(PT) and ion < 3:
        return f'{PT[z - 1]} {ROMAN[ion]}'
    return f'{code:.2f}'


class LineOpacity:
    """Line opacity per species at one depth, on a log-uniform wavelength grid."""

    def __init__(self, phot, lam_lo_A=3795.0, lam_hi_A=10010.0, resolu=300000.0,
                 mol_path=None):
        self.phot = phot
        self.dln = np.log1p(1.0 / resolu)
        n = int(np.log(lam_hi_A / lam_lo_A) / self.dln) + 1
        self.lam = lam_lo_A * np.exp(self.dln * np.arange(n))     # vacuum A
        self.n = n
        self.nden = read_moldump(mol_path) if mol_path else read_moldump()
        self.potion = read_ionpots()
        # continuum opacity per cm^3 at the reference depth, interpolated
        lg = np.geomspace(lam_lo_A, lam_hi_A, 400)
        kc = np.array([R.continuum_opacity(phot.atm, la)[0][phot.j] for la in lg])
        self.kcont = np.exp(np.interp(np.log(self.lam), np.log(lg), np.log(kc))) \
            * phot.atm.nH[phot.j]                                  # cm^-1

    # -- populations -----------------------------------------------------
    def n_over_U(self, code):
        """n(species)/U at the reference depth, from SYNTHE's EOS where it has it."""
        z = int(np.floor(code + 1e-6))
        ion = int(round((code - z) * 100))
        el = PT[z - 1]
        U = self.phot.ad.Q(el, self.phot.T, ion)
        key = round(code, 2)
        if key in self.nden:
            return self.nden[key][self.phot.j] / U
        # not in SYNTHE's molecular network (Z > 30 apart from Y, Zr, La):
        # fall back to Saha with the model's own electron density
        nel = self.phot.atm.abund.get(el)
        if nel is None:
            return None
        return (nel * self.phot.nH * self.phot.ion_fraction(el, ion)) / U

    # -- line quantities -------------------------------------------------
    def prepare(self, g):
        p = self.phot
        gr, gs, gw, z, icharge = synthe_damping(g, self.potion)
        ok = (z >= 1) & (z <= len(PT)) & (icharge <= 2) & (g['aut'] == 0)
        keep = {k: v[ok] for k, v in g.items()}
        gr, gs, gw = gr[ok], gs[ok], gw[ok]
        z, icharge = z[ok], icharge[ok]
        code = np.round(z + icharge / 100.0, 2)

        # SYNTHE's perturber sum for van der Waals
        nHI = self.nden[1.00][p.j]
        nHeI = self.nden[2.00][p.j]
        nH2 = self.nden[101.00][p.j]
        txnxn = (nHI + 0.42 * nHeI + 0.85 * nH2) * (p.T / 1e4)**0.3
        gamma = gr + gs * p.ne + gw * txnxn

        lam = keep['wl_nm'] * 10.0                       # vacuum A
        nu0 = R.C / (lam * 1e-8)
        mass = np.array([MASS.get(PT[i - 1], 2.5 * i) for i in z]) * R.AMU
        vD = np.sqrt(2.0 * R.KB * p.T / mass + p.vturb**2)
        dnuD = nu0 * vD / R.C
        stim = 1.0 - np.exp(-R.H * nu0 / (R.KB * p.T))
        gf = 10.0**(keep['loggf'] + keep['x1'] + keep['x2'])
        boltz = np.exp(-keep['elo'] / (0.695034800 * p.T))   # E in cm^-1

        nu_over_U = np.zeros(len(gf))
        for c in np.unique(code):
            m = code == c
            v = self.n_over_U(c)
            nu_over_U[m] = 0.0 if v is None else v

        # frequency-integrated line opacity [cm^-1 Hz], and the line-centre value
        area = R.PIE2MC * gf * nu_over_U * boltz * stim
        a = gamma / (4.0 * np.pi * dnuD)
        peak = area / (np.sqrt(np.pi) * dnuD) * np.real(wofz(1j * a))
        return dict(lam=lam, dnuD=dnuD, a=a, area=area, peak=peak,
                    code=code, vD=vD)

    # -- accumulation ----------------------------------------------------
    def accumulate(self, L, groups, floor=1e-4, hw_max=8192):
        """Sum Voigt profiles into one array per group.  groups: code -> label."""
        out = {lab: np.zeros(self.n) for lab in set(groups.values())}
        out['other'] = np.zeros(self.n)

        ix = np.rint(np.log(L['lam'] / self.lam[0]) / self.dln).astype(np.int64)
        inside = (ix >= 0) & (ix < self.n)
        kc_at = np.where(inside, self.kcont[np.clip(ix, 0, self.n - 1)], np.inf)
        eta0 = L['peak'] / kc_at

        # how far out to carry each line, in Doppler widths
        vg = np.sqrt(np.maximum(np.log(np.maximum(eta0 / floor, 1.0)), 0.0))
        vl = np.sqrt(np.maximum(L['a'] * eta0 / (np.sqrt(np.pi) * floor), 0.0))
        vcut = np.maximum(np.maximum(vg, vl), 3.0)
        pts = vcut * L['vD'] / R.C / self.dln                 # half-window, grid points
        hw = np.clip(2**np.ceil(np.log2(np.maximum(pts, 4.0))), 4, hw_max).astype(np.int64)

        alive = inside & (eta0 > floor)
        labels = np.array([groups.get(c, 'other') for c in L['code']])
        nused = 0
        for h in np.unique(hw[alive]):
            sel = np.where(alive & (hw == h))[0]
            block = max(1, int(4e6 // (2 * h + 1)))
            for s in range(0, len(sel), block):
                k = sel[s:s + block]
                idx = ix[k][:, None] + np.arange(-h, h + 1)[None, :]
                good = (idx >= 0) & (idx < self.n)
                idxc = np.clip(idx, 0, self.n - 1)
                lam = self.lam[idxc]
                nu = R.C / (lam * 1e-8)
                v = (nu - R.C / (L['lam'][k][:, None] * 1e-8)) / L['dnuD'][k][:, None]
                phi = np.real(wofz(v + 1j * L['a'][k][:, None])) \
                    / (np.sqrt(np.pi) * L['dnuD'][k][:, None])
                contrib = (L['area'][k][:, None] * phi) * good
                for lab in np.unique(labels[k]):
                    m = labels[k] == lab
                    out[lab] += np.bincount(idxc[m].ravel(),
                                            weights=contrib[m].ravel(),
                                            minlength=self.n)
                nused += len(k)
        return out, nused, int(alive.sum()), len(alive)


# ----------------------------------------------------------------------
# molecules
# ----------------------------------------------------------------------
MOL_OPT = HERE / 'cache' / 'mol_opt.dat'
BC16_T6 = Path.home() / 'kurucz' / 'upgrade' / 'raw_data' / 'bc16_table6.dat'
MOL_BROAD = Path.home() / 'kurucz' / 'atlas12' / 'data' / 'mol_broad.dat'

# Kurucz molecule code -> Barklem & Collet (2016) Table 6 label
MOL_NAME = {106: 'CH', 107: 'NH', 108: 'OH', 111: 'NaH', 112: 'MgH',
            113: 'AlH', 114: 'SiH', 120: 'CaH', 124: 'CrH', 126: 'FeH',
            606: 'C2', 607: 'CN', 608: 'CO', 812: 'MgO', 813: 'AlO',
            814: 'SiO', 822: 'TiO', 823: 'VO'}

MOLBROAD_TREF = 3000.0
_WH2, _PREF_DYN = 0.85, 1.01325e6


def read_mol_partition(path=BC16_T6):
    """BC16 Table 6: molecular partition functions Q(T) on their own T grid."""
    lines = path.read_text().split('\n')
    T = np.array([float(x) for x in lines[2].split()[2:]])
    q = {}
    for line in lines[3:]:
        p = line.split()
        if len(p) != len(T) + 1:
            continue
        q[p[0]] = np.array([float(x) for x in p[1:]])
    return T, q


def read_mol_broad(path=MOL_BROAD):
    """mol_broad.dat converted to SYNTHE's gamma_w convention (load_mol_broad)."""
    out = {}
    for line in path.read_text().split('\n'):
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        p = line.split()
        try:
            code, g2, n2, d2 = int(p[0]), float(p[2]), float(p[3]), float(p[6])
        except (ValueError, IndexError):
            continue
        conv = ((296.0 / MOLBROAD_TREF)**n2
                * (R.KB * MOLBROAD_TREF / _PREF_DYN) * 4.0 * np.pi * R.C
                / (_WH2 * (MOLBROAD_TREF / 1e4)**0.3))
        out[code] = (g2 * conv, d2 * conv)
    return out


def parse_mol_lines(path=MOL_OPT):
    """Kurucz molecular line records:
    (F10.4,F7.3,F5.1,F10.3,F5.1,F11.3,I4,A8,A8,I2,I4)."""
    cache = path.parent / 'mol_full.npz'
    if cache.exists():
        d = np.load(cache)
        return {k: d[k] for k in d.files}

    def f(s, default=0.0):
        try:
            return float(s)
        except ValueError:
            return default

    cols = {k: [] for k in ('wl_nm', 'loggf', 'elo', 'eup', 'jlo', 'icode',
                            'iso', 'loggr', 'upperX')}
    with open(path, errors='replace') as fh:
        for line in fh:
            e, ep = f(line[22:32]), f(line[37:48])
            xj, xjp = f(line[17:22]), f(line[32:37])
            # SYNTHE ignores the wavelength column for molecules and
            # recomputes the VACUUM wavelength from the level energies
            # (read_molec_ascii: wlvac = 1e7 / (|ep| - |e|)).  The column
            # itself is an air wavelength, 1.2 A off at 4300 A.
            lo_e, up_e = min(abs(e), abs(ep)), max(abs(e), abs(ep))
            dE = up_e - lo_e
            cols['wl_nm'].append(1.0e7 / dE if dE > 0 else f(line[0:10]))
            cols['loggf'].append(f(line[10:17], -99.0))
            cols['elo'].append(lo_e)
            cols['eup'].append(up_e)
            cols['jlo'].append(xj if abs(e) <= abs(ep) else xjp)
            cols['icode'].append(f(line[48:52]))
            cols['iso'].append(f(line[68:70]))
            cols['loggr'].append(f(line[70:74]))
            cols['upperX'].append(1.0 if line[60:61] == 'X' else 0.0)
    out = {k: np.array(v) for k, v in cols.items()}
    np.savez_compressed(cache, **out)
    return out


def _molecular_prepare(self, m, only=None):
    """Line quantities for molecular records, following SYNTHE's read_molec."""
    import moldisp
    p = self.phot
    disp = moldisp.parse_dispatch()
    broad = read_mol_broad()
    Tq, Q = read_mol_partition()

    icode = m['icode'].astype(int)
    iso = m['iso'].astype(int)
    ok = np.array([c in MOL_NAME and MOL_NAME[c] in Q and
                   (only is None or c in only) for c in icode])
    ok &= np.array([(i, c) in disp or (i, None) in disp
                    for i, c in zip(iso, icode)])
    k = {key: v[ok] for key, v in m.items()}
    icode, iso = icode[ok], iso[ok]

    shift = np.array([moldisp.isotope_shift(disp, i, c) for i, c in zip(iso, icode)])
    x1, x2, mass_amu = shift[:, 0], shift[:, 1], shift[:, 2]

    lam = k['wl_nm'] * 10.0
    nu0 = R.C / (lam * 1e-8)
    gr = np.where(k['loggr'] == 0.0, 2.223e13 / k['wl_nm']**2,
                  10.0**(k['loggr'] * 0.01))
    isX = k['upperX'] > 0
    gs = np.where(isX, 3.0e-8, 3.0e-5)
    gw = np.where(isX, 1.0e-8, 1.0e-7)
    for c, (g0, dg) in broad.items():
        sel = icode == c
        if sel.any() and g0 > 0:
            gw[sel] = np.maximum(g0 - dg * k['jlo'][sel], 0.1 * g0)

    nHI, nHeI, nH2 = self.nden[1.00][p.j], self.nden[2.00][p.j], self.nden[101.00][p.j]
    txnxn = (nHI + 0.42 * nHeI + 0.85 * nH2) * (p.T / 1e4)**0.3
    gamma = gr + gs * p.ne + gw * txnxn

    vD = np.sqrt(2.0 * R.KB * p.T / (mass_amu * R.AMU) + p.vturb**2)
    dnuD = nu0 * vD / R.C
    stim = 1.0 - np.exp(-R.H * nu0 / (R.KB * p.T))
    gf = 10.0**(k['loggf'] + x1 + x2)
    boltz = np.exp(-k['elo'] / (0.695034800 * p.T))

    n_over_U = np.zeros(len(gf))
    for c in np.unique(icode):
        n = self.nden.get(float(c))
        if n is None:
            continue
        n_over_U[icode == c] = n[p.j] / np.interp(p.T, Tq, Q[MOL_NAME[c]])

    area = R.PIE2MC * gf * n_over_U * boltz * stim
    a = gamma / (4.0 * np.pi * dnuD)
    peak = area / (np.sqrt(np.pi) * dnuD) * np.real(wofz(1j * a))
    code = np.array([1000 + c for c in icode], dtype=float)   # tag as molecular
    return dict(lam=lam, dnuD=dnuD, a=a, area=area, peak=peak, code=code, vD=vD)


LineOpacity.prepare_molecules = _molecular_prepare


def _nden_all(self):
    """SYNTHE's number densities for every species, at every depth."""
    if not hasattr(self, '_nda'):
        self._nda = read_moldump()
    return self._nda


def _ion_fraction_at(self, el, stage, j):
    """Saha ionization fraction at depth j (for species SYNTHE's EOS omits)."""
    T, ne = float(self.atm.T[j]), float(self.atm.ne[j])
    chi = self.ad._ie[el]
    U = []
    for i in (0, 1, 2):
        try:
            U.append(self.ad.Q(el, T, i))
        except KeyError:                 # H and He have no third stage
            U.append(0.0)
    r1 = 2.0 * R.SAHA_TR * T**1.5 * (U[1] / U[0]) * np.exp(-chi[0] / (R.KB_EV * T)) / ne
    r2 = (0.0 if U[2] == 0.0 or U[1] == 0.0 else
          2.0 * R.SAHA_TR * T**1.5 * (U[2] / U[1])
          * np.exp(-chi[1] / (R.KB_EV * T)) / ne)
    den = 1.0 + r1 + r1 * r2
    return [1.0 / den, r1 / den, r1 * r2 / den][stage]


Photosphere.nden_all = property(_nden_all)
Photosphere.ion_fraction_at = _ion_fraction_at
