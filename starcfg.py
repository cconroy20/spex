"""The stars the explorer ships, and where each one's inputs live.

Adding a star is a row here plus `./run_star.sh <tag> <model.atm>`.
Parameters are the ones the model was computed at; the placeholders in the
web app carry the same numbers.
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent
GRIDS = Path.home() / 'kurucz' / 'grids'

STARS = {
    'sun': dict(
        name='Sun', sp='G2 V', teff=5777, logg=4.44, feh=0.0, afe=0.0,
        model=GRIDS / 'THESUN' / 'atm' / 'ap00t5777g4.44at12.dat'),
    'procyon': dict(
        name='Procyon', sp='F5 IV-V', teff=6530, logg=3.96, feh=0.0, afe=0.0,
        model=HERE / 'run' / 'atm' / 'procyon' / 'procyon.atm'),
    'hd122563': dict(
        name='HD 122563', sp='G8 III', teff=4587, logg=1.61, feh=-2.64, afe=0.4,
        model=HERE / 'run' / 'atm' / 'hd122563_a' / 'hd122563.atm'),
    'arcturus': dict(
        name='Arcturus', sp='K1.5 III', teff=4286, logg=1.66, feh=-0.52, afe=0.3,
        model=HERE / 'run' / 'atm' / 'arcturus_a' / 'arcturus.atm'),
    'barnard': dict(
        name="Barnard's Star", sp='M4 V', teff=3220, logg=5.05, feh=-0.40,
        afe=0.0, model=HERE / 'run' / 'atm' / 'barnard' / 'barnard.atm'),
}


def star(tag):
    if tag not in STARS:
        raise SystemExit(f'unknown star {tag!r}; have {", ".join(STARS)}')
    s = dict(STARS[tag], tag=tag)
    s['stem'] = s['model'].stem                  # SYNTHE names outputs from this
    s['run'] = HERE / 'run' / tag
    s['mol'] = s['run'] / 'eos' / f'{s["stem"]}.mol'
    return s


# ----------------------------------------------------------------------
# species
#
# (run tag, display name, group, kind, source files).  The run tag is also
# the directory each SYNTHE run writes into, so it has to stay stable.
# `gfall` files are the per-species splits under cache/species/ made by
# split_gfall.py; `mol` files live in the ATLAS12 data/mol directory.
# ----------------------------------------------------------------------
LIGHT = 'light elements'
IRON = 'iron peak'
IONS = 'ions'
NCAP = 'neutron capture'
MOLEC = 'molecules'

SPECIES = [
    ('a1.00', 'H I', LIGHT), ('a3.00', 'Li I', LIGHT), ('a6.00', 'C I', LIGHT),
    ('a7.00', 'N I', LIGHT), ('a8.00', 'O I', LIGHT), ('a11.00', 'Na I', LIGHT),
    ('a12.00', 'Mg I', LIGHT), ('a13.00', 'Al I', LIGHT),
    ('a14.00', 'Si I', LIGHT), ('a16.00', 'S I', LIGHT),
    ('a19.00', 'K I', LIGHT), ('a20.00', 'Ca I', LIGHT),

    ('a21.00', 'Sc I', IRON), ('a22.00', 'Ti I', IRON), ('a23.00', 'V I', IRON),
    ('a24.00', 'Cr I', IRON), ('a25.00', 'Mn I', IRON),
    ('a26.00', 'Fe I', IRON), ('a27.00', 'Co I', IRON),
    ('a28.00', 'Ni I', IRON), ('a29.00', 'Cu I', IRON),
    ('a30.00', 'Zn I', IRON),

    ('a12.01', 'Mg II', IONS), ('a14.01', 'Si II', IONS),
    ('a20.01', 'Ca II', IONS), ('a21.01', 'Sc II', IONS),
    ('a22.01', 'Ti II', IONS), ('a23.01', 'V II', IONS),
    ('a24.01', 'Cr II', IONS), ('a25.01', 'Mn II', IONS),
    ('a26.01', 'Fe II', IONS),

    ('a38.01', 'Sr II', NCAP), ('a39.01', 'Y II', NCAP),
    ('a40.01', 'Zr II', NCAP), ('a56.01', 'Ba II', NCAP),
    ('a57.01', 'La II', NCAP), ('a60.01', 'Nd II', NCAP),
    ('a63.01', 'Eu II', NCAP),

    ('mCH', 'CH', MOLEC, ['chjorg.dat']),
    ('mCN', 'CN', MOLEC, ['cnaxbrooke.dat', 'cnbxbrooke.dat', 'cnxx12brooke.dat']),
    ('mC2', 'C2', MOLEC, ['c2ax.dat', 'c2ba.dat', 'c2dabrooke.dat', 'c2ea.dat']),
    ('mMgH', 'MgH', MOLEC, ['mgh24_owens22.dat', 'mgh25_owens22.dat',
                            'mgh26_owens22.dat', 'mghax.dat', 'mghbx.dat']),
    ('mSiH', 'SiH', MOLEC, ['sihaxsightly.asc', 'sihxxsightly.asc']),
    ('mCaH', 'CaH', MOLEC, ['cah_owens22.dat']),
    ('mCrH', 'CrH', MOLEC, ['crhax.dat']),
    ('mFeH', 'FeH', MOLEC, ['fehfx.dat']),
    ('mTiO', 'TiO', MOLEC, ['tiototo2024.bin']),
    ('mVO', 'VO', MOLEC, ['voax.asc', 'vobx.asc', 'vocx.asc']),
    ('mCO', 'CO', MOLEC, ['coax.dat', 'coxx.dat']),
    ('mOH', 'OH', MOLEC, ['ohax.dat', 'ohxx.dat']),
    ('mH2O', 'H2O', MOLEC, ['h2opokazatel.bin'], 'h2o'),
]

GROUP_ORDER = [LIGHT, IRON, IONS, NCAP, MOLEC]


def species():
    """(tag, name, group, kind, files) with the atomic files filled in."""
    out = []
    for row in SPECIES:
        tag, name, group = row[:3]
        if tag.startswith('a'):
            out.append((tag, name, group, 'gfall', [f'gf_{tag[1:]}.dat']))
        else:
            # lines.list has separate reader types; h2o is not read as `mol`
            out.append((tag, name, group, row[4] if len(row) > 4 else 'mol', row[3]))
    return out


def write_species_list(path):
    """The same table, in the tab-separated form run_star.sh reads."""
    lines = []
    for tag, name, group, kind, files in species():
        lines.append('\t'.join([tag, kind] + files))
    Path(path).write_text('\n'.join(lines) + '\n')
    return len(lines)
