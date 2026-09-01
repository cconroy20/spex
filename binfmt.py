"""The on-disk format for the compressed series.

A series is a uint16 array on the log-lambda grid.  Neighbouring points differ
by very little -- the grid is critically sampled at R = 300,000 -- so a first
difference makes the high byte nearly constant, and separating the byte planes
lets deflate find that.  Compressing each chunk independently instead of the
whole file costs about 5% and keeps HTTP Range working, which is the whole
reason the native tier is usable at all.

    magic   'SPC1'                       4 bytes
    n       uint32                       points in the series
    chunk   uint32                       points per chunk
    nchunk  uint32
    offset  uint32 x (nchunk + 1)        byte offsets, from the start of file
    payload deflate(hi_bytes || lo_bytes) per chunk

The difference is taken mod 2**16 and reconstructed by a running sum with the
same wraparound, so nothing is lost: this is exact, not lossy.
"""
import struct
import zlib

import numpy as np

MAGIC = b'SPC1'
CHUNK = 16384


def encode_chunk(a):
    d = np.diff(a.astype(np.uint16), prepend=np.uint16(0)).astype('<u2')
    return zlib.compress((d >> 8).astype('u1').tobytes()
                         + (d & 255).astype('u1').tobytes(), 9)


def decode_chunk(blob):
    raw = np.frombuffer(zlib.decompress(blob), dtype='u1')
    h = len(raw) // 2
    d = (raw[:h].astype('<u2') << 8) | raw[h:].astype('<u2')
    return np.cumsum(d, dtype='<u8').astype('<u2')      # wraps, as intended


def write_series(path, a, chunk=CHUNK):
    a = np.ascontiguousarray(a, dtype='<u2')
    parts = [encode_chunk(a[s:s + chunk]) for s in range(0, len(a), chunk)]
    head = 16 + 4 * (len(parts) + 1)
    off, p = [head], head
    for b in parts:
        p += len(b)
        off.append(p)
    with open(path, 'wb') as f:
        f.write(MAGIC + struct.pack('<III', len(a), chunk, len(parts)))
        f.write(np.array(off, '<u4').tobytes())
        for b in parts:
            f.write(b)
    return p


def read_series(path):
    """Reader for verification -- the browser does this incrementally."""
    raw = open(path, 'rb').read()
    assert raw[:4] == MAGIC, path
    n, chunk, nchunk = struct.unpack('<III', raw[4:16])
    off = np.frombuffer(raw[16:16 + 4 * (nchunk + 1)], '<u4')
    out = np.empty(n, '<u2')
    for c in range(nchunk):
        v = decode_chunk(raw[off[c]:off[c + 1]])
        out[c * chunk:c * chunk + len(v)] = v
    return out
