"""Parse SYNTHE's molec_dispatch isotope table straight out of the Fortran."""
import re
from pathlib import Path

SRC = Path.home() / 'kurucz' / 'atlas12' / 'src' / 'mod_mklinelist.f90'


def parse_dispatch(src=SRC):
    """{(iso, icode_or_None): (x1, x2)} exactly as molec_dispatch assigns it."""
    txt = src.read_text().split('\n')
    lo = next(i for i, l in enumerate(txt) if 'SUBROUTINE molec_dispatch' in l)
    hi = next(i for i, l in enumerate(txt[lo:], lo) if 'END SUBROUTINE molec_dispatch' in l)
    body = txt[lo:hi]

    out, iso, icode, depth = {}, None, None, 0
    for raw in body:
        s = raw.strip()
        if s.startswith('SELECT CASE (iso)'):
            depth = 1
            continue
        if s.startswith('SELECT CASE (icode)'):
            depth = 2
            continue
        if s.startswith('END SELECT'):
            depth = max(0, depth - 1)
            icode = None
            continue
        m = re.match(r'CASE \((\d+)\)', s)
        if m:
            if depth <= 1:
                iso, icode = int(m.group(1)), None
            else:
                icode = int(m.group(1))
        elif s.startswith('CASE DEFAULT'):
            icode = None
        elif re.match(r'IF \(icode \.EQ\. (\d+)\)', s):
            icode = int(re.match(r'IF \(icode \.EQ\. (\d+)\)', s).group(1))
        elif s.startswith('ELSE'):
            icode = None
        mx = re.search(r'iso1=\s*(\d+);\s*iso2=\s*(\d+);\s*x1=\s*(-?[\d.]+);\s*x2=\s*(-?[\d.]+)', s)
        if mx and iso is not None:
            out[(iso, icode)] = (float(mx.group(3)), float(mx.group(4)),
                                 int(mx.group(1)) + int(mx.group(2)))
    return out


def isotope_shift(table, iso, icode):
    """x1 + x2 for this (iso, icode), falling back to the CASE DEFAULT branch."""
    if (iso, icode) in table:
        return table[(iso, icode)]
    return table.get((iso, None))


if __name__ == '__main__':
    t = parse_dispatch()
    print(f'{len(t)} (iso, icode) entries parsed')
    for iso, icode, nm in [(12, 106, 'C12-H  (CH)'), (13, 106, 'C13-H  (CH)'),
                           (24, 112, 'Mg24-H (MgH)'), (25, 112, 'Mg25-H'),
                           (26, 112, 'Mg26-H'), (12, 607, 'C12-N14 (CN)'),
                           (13, 607, 'C13-N (CN)'), (15, 607, 'C-N15 (CN)'),
                           (12, 606, 'C12-C12 (C2)'), (16, 108, 'O16-H (OH)'),
                           (28, 114, 'Si28-H (SiH)')]:
        print(f'  {nm:14s} iso={iso:3d} icode={icode:4d} -> x1,x2 = {isotope_shift(t, iso, icode)}')
