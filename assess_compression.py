"""What the shipped bundle would cost per star under each candidate scheme.

    python3 assess_compression.py [tag ...]

Schemes measured:
  raw            what is shipped today
  gzip           plain deflate, whole file
  d+split        first difference, byte planes separated, deflate
  d+split/chunk  the same but per 16,384-point chunk, so HTTP Range still works
"""
import json
import sys
import zlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CHUNK = 16384


def pack(a):
    d = np.diff(a.astype(np.int32), prepend=a[0]).astype('<i2').view('<u2')
    return (d >> 8).astype('u1').tobytes() + (d & 255).astype('u1').tobytes()


def sizes(a, chunked):
    if not chunked:
        return len(zlib.compress(pack(a), 9))
    return sum(len(zlib.compress(pack(a[s:s + CHUNK]), 9))
               for s in range(0, len(a), CHUNK))


tags = sys.argv[1:] or ['sun']
tot = {}
for tag in tags:
    D = HERE / 'web' / 'data' / tag
    if not (D / 'meta.json').exists():
        print(f'{tag}: no bundle'); continue
    M = json.loads((D / 'meta.json').read_text())
    r = g = c = k = 0
    for f in (D / 'full').glob('*.bin'):
        a = np.frombuffer(f.read_bytes(), '<u2')
        r += len(a) * 2
        g += len(zlib.compress(a.tobytes(), 9))
        c += sizes(a, False)
        k += sizes(a, True)
    ovr = ovc = 0
    n = M['n_ov']
    for f in (D / 'ov').glob('*.bin'):
        a = np.frombuffer(f.read_bytes(), '<u2'); ovr += len(a) * 2
        ovc += sum(sizes(a[i * n:(i + 1) * n], False) for i in range(3))
    lr = lz = 0
    for f in (D / 'lines').glob('[0-9]*.json'):
        b = f.read_bytes(); lr += len(b); lz += len(zlib.compress(b, 9))
    tot[tag] = dict(full_raw=r, full_gzip=g, full_ds=c, full_chunk=k,
                    ov_raw=ovr, ov_ds=ovc, lines_raw=lr, lines_gzip=lz)
    print(f'{tag:10s} full {r/1e6:6.2f} -> gzip {g/1e6:5.2f} / d+split {c/1e6:5.2f}'
          f' / chunked {k/1e6:5.2f} MB   ({r/k:.1f}x)')
    print(f'{"":10s} ov   {ovr/1e6:6.2f} -> {ovc/1e6:5.2f} MB      '
          f'lines {lr/1e6:5.2f} -> {lz/1e6:5.2f} MB (json gzip)')

if len(tot) > 1:
    R = sum(v['full_raw'] + v['ov_raw'] + v['lines_raw'] for v in tot.values())
    C = sum(v['full_chunk'] + v['ov_ds'] + v['lines_gzip'] for v in tot.values())
    print(f'\nall {len(tot)} stars: {R/1e6:.1f} MB now -> {C/1e6:.1f} MB '
          f'({R/C:.1f}x), before sharing the line index')
