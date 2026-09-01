'use strict';
/* The Solar Spectrum, species by species.
   Static site: no server logic.  Two data tiers per series -- a decimated
   overview loaded up front (min/max/mean per bin) and the native R = 300,000
   arrays fetched in 32 kB chunks by HTTP range request only when zoomed in. */

const DPR = Math.min(window.devicePixelRatio || 1, 2);
const CHUNK = 16384;          // native points per fetched chunk
const FULL_MAX = 80000;       // window size below which we use native data
const ATOM = '#004C8C', MOL = '#B26B00', INK = '#1a1a1a', GREY = '#a9a49b',
      AMBER = '#B26B00';
let M = null, OV = {}, CH = {}, WHOLE = {};
const state = { i0: 0, i1: 1, on: new Set(), R: 300000, ymin: 0, mode: 'individual' };
const $ = s => document.querySelector(s);
const R_MIN = 1000, R_NATIVE = 300000;
const posToR = p => Math.round(10 ** (Math.log10(R_MIN)
  + (p / 1000) * (Math.log10(R_NATIVE) - Math.log10(R_MIN))));
const rToPos = R => Math.round(1000 * (Math.log10(R || R_NATIVE) - Math.log10(R_MIN))
  / (Math.log10(R_NATIVE) - Math.log10(R_MIN)));
/* display names carry a space ('Fe I'); files use an underscore */
const fileOf = n => n.replace(/ /g, '_');

/* Which star is on show.  Everything the app fetches hangs off DB, so adding
   a star is a directory of the same shape plus a row in STARS.  The three
   without a `dir` are placeholders: their parameters are nominal literature
   values, to be replaced by whatever the models are actually computed at. */
const STARS = [
  { name: 'Sun', sp: 'G2 V', teff: 5777, logg: 4.44, feh: 0.0, afe: 0.0, dir: 'data/sun' },
  { name: 'Procyon', sp: 'F5 IV-V', teff: 6530, logg: 3.96, feh: 0.0, afe: 0.0, dir: 'data/procyon' },
  { name: 'HD 122563', sp: 'G8 III', teff: 4587, logg: 1.61, feh: -2.64, afe: 0.4, dir: 'data/hd122563' },
  { name: 'Arcturus', sp: 'K1.5 III', teff: 4286, logg: 1.66, feh: -0.52, afe: 0.3, dir: 'data/arcturus' },
  { name: "Barnard's Star", sp: 'M4 V', teff: 3220, logg: 5.05, feh: -0.40, afe: 0.0, dir: 'data/barnard' },
];
let DB = STARS[0].dir;
function fail(err) {
  console.error(err);
  let el = document.getElementById('err');
  if (!el) {
    el = document.createElement('div');
    el.id = 'err';
    document.body.appendChild(el);
  }
  el.textContent = 'Error: ' + (err && err.message ? err.message : err);
  el.style.display = 'block';
}
window.addEventListener('error', e => fail(e.error || e.message));
window.addEventListener('unhandledrejection', e => fail(e.reason));

/* ---------- wavelength scale ---------- */
const vac2air = lv => {
  if (lv <= 2000) return lv;
  const s2 = (1e4 / lv) ** 2;
  return lv / (1 + 8.34254e-5 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2));
};
const air2vac = la => { let v = la; for (let k = 0; k < 4; k++) v = la * (v / vac2air(v)); return v; };
const idxAir = i => vac2air(M.lam0_vac * Math.exp(i * M.dln));
const airIdx = a => Math.log(air2vac(a) / M.lam0_vac) / M.dln;

/* ---------- data ---------- */
const dec = (u, name) => {
  const hi = (name === '_flux' || name === '_cont') ? M.flux_max : M.norm_max;
  return hi / 65535;
};
/* Series arrive as SPC1: the array first-differenced, its byte planes split,
   and each 16,384-point chunk deflated on its own so a Range request still
   lands on a chunk boundary.  See binfmt.py for the container. */
const HDR = {}, INFLIGHT = {};

function inflate(buf) {
  const st = new Blob([buf]).stream().pipeThrough(new DecompressionStream('deflate'));
  return new Response(st).arrayBuffer().then(b => new Uint8Array(b));
}
function undelta(raw) {
  const h = raw.length >> 1, out = new Uint16Array(h);
  let acc = 0;
  for (let i = 0; i < h; i++) {
    acc = (acc + ((raw[i] << 8) | raw[h + i])) & 0xFFFF;
    out[i] = acc;
  }
  return out;
}
function parseHdr(b) {
  if (b.length < 16 || b[0] !== 83 || b[1] !== 80 || b[2] !== 67 || b[3] !== 49) return null;
  const dv = new DataView(b.buffer, b.byteOffset, b.byteLength);
  const n = dv.getUint32(4, true), chunk = dv.getUint32(8, true), nchunk = dv.getUint32(12, true);
  if (b.length < 16 + 4 * (nchunk + 1)) return null;
  const off = new Uint32Array(nchunk + 1);
  for (let i = 0; i <= nchunk; i++) off[i] = dv.getUint32(16 + 4 * i, true);
  return { n, chunk, nchunk, off };
}
async function decodeAll(bytes, h) {
  const out = new Uint16Array(h.n);
  for (let c = 0; c < h.nchunk; c++) {
    const v = undelta(await inflate(bytes.subarray(h.off[c], h.off[c + 1])));
    out.set(v.subarray(0, Math.min(v.length, h.n - c * h.chunk)), c * h.chunk);
  }
  return out;
}

/* _flux is not shipped: it is _norm x _cont, and rebuilding it from the two
   quantized arrays lands within two quantization steps of the real thing. */
const DERIVED = { _flux: ['_norm', '_cont'] };
function fluxFrom(no, co) {
  const n = Math.min(no.length, co.length), out = new Uint16Array(n);
  const k = M.norm_max / 65535;
  for (let i = 0; i < n; i++) {
    const v = Math.round(no[i] * co[i] * k);
    out[i] = v > 65535 ? 65535 : v;
  }
  return out;
}

async function fetchSeries(dir, name) {
  const r = await fetch(`${DB}/${dir}/${fileOf(name)}.bin`);
  if (!r.ok) return null;
  const b = new Uint8Array(await r.arrayBuffer());
  const h = parseHdr(b);
  if (!h) return null;
  HDR[dir + '/' + name] = h;
  return decodeAll(b, h);
}

async function loadOv(name) {
  if (OV[name]) return OV[name];
  const key = 'O' + name;
  if (INFLIGHT[key]) return INFLIGHT[key];
  INFLIGHT[key] = (async () => {
    let a;
    if (DERIVED[name]) {
      const [p, q] = DERIVED[name];
      const [A, B] = await Promise.all([loadOv(p), loadOv(q)]);
      if (!A || !B) return null;
      a = { min: fluxFrom(A.min, B.min), max: fluxFrom(A.max, B.max),
            mean: fluxFrom(A.mean, B.mean) };
      OV[name] = a;
      return a;
    }
    const v = await fetchSeries('ov', name);
    if (!v) return null;
    const n = M.n_ov;
    OV[name] = { min: v.subarray(0, n), max: v.subarray(n, 2 * n),
                 mean: v.subarray(2 * n, 3 * n) };
    return OV[name];
  })().catch(() => null).finally(() => { delete INFLIGHT[key]; });
  return INFLIGHT[key];
}

async function getWhole(name) {
  if (WHOLE[name]) return WHOLE[name];
  const key = 'W' + name;
  if (INFLIGHT[key]) return INFLIGHT[key];
  INFLIGHT[key] = (async () => {
    if (DERIVED[name]) {
      const [p, q] = DERIVED[name];
      await Promise.all([getWhole(p), getWhole(q)]);
      if (WHOLE[p] && WHOLE[q]) WHOLE[name] = fluxFrom(WHOLE[p], WHOLE[q]);
      return WHOLE[name];
    }
    const v = await fetchSeries('full', name);
    if (v) WHOLE[name] = v;
    return v;
  })().catch(() => null).finally(() => { invalidateChunkCache(); delete INFLIGHT[key]; });
  return INFLIGHT[key];
}
async function getFull(name, i0, i1) { return WHOLE[name] ? WHOLE[name] : getWhole(name); }

let _cName = '', _cIdx = -1, _cBuf = null;
function fullAt(name, i) {
  const whole = WHOLE[name];
  return whole ? whole[i] : -1;
}
function invalidateChunkCache() { _cIdx = -1; _cName = ''; _cBuf = null; }

/* ---------- per-pixel reduction ---------- */
function blank(W) {
  const a = new Float32Array(W).fill(NaN);
  return { lo: a, hi: a, mid: a };
}
function envelope(name, i0, i1, W, useFull) {
  // refresh() paints before awaiting, so a series may not be here yet; draw
  // nothing for it rather than throwing and killing the whole frame
  if (!useFull && !OV[name]) return blank(W);
  const lo = new Float32Array(W), hi = new Float32Array(W), mid = new Float32Array(W);
  const k = dec(null, name);
  for (let x = 0; x < W; x++) {
    const a = i0 + (i1 - i0) * x / W, b = i0 + (i1 - i0) * (x + 1) / W;
    let mn = Infinity, mx = -Infinity, sum = 0, cnt = 0;
    if (useFull) {
      for (let i = Math.floor(a); i < Math.max(Math.ceil(b), Math.floor(a) + 1); i++) {
        const v = fullAt(name, i); if (v < 0) continue;
        if (v < mn) mn = v; if (v > mx) mx = v; sum += v; cnt++;
      }
    } else {
      const o = OV[name];
      const j0 = Math.floor(a / M.ov_bin), j1 = Math.max(Math.ceil(b / M.ov_bin), j0 + 1);
      for (let j = j0; j < j1 && j < M.n_ov; j++) {
        if (o.min[j] < mn) mn = o.min[j]; if (o.max[j] > mx) mx = o.max[j];
        sum += o.mean[j]; cnt++;
      }
    }
    lo[x] = mn === Infinity ? NaN : mn * k;
    hi[x] = mx === -Infinity ? NaN : mx * k;
    mid[x] = cnt ? sum * k / cnt : NaN;
  }
  return { lo, hi, mid };
}
/* Effective sigma of what actually reaches the eye: the instrumental Gaussian
   of FWHM = c/R, combined in quadrature with the box of one screen pixel.
   Without the pixel term the kernel could be far narrower than the pixel, and
   the curve was point-sampled -- one bin in four at low zoom -- which aliased
   badly and made the selector look broken.  Returns sigma in tier points. */
function kernelSigma(i0, i1, W, step, R) {
  const sigR = R ? (1 / R) / M.dln / 2.35482 / step : 0;
  const px = (i1 - i0) / W / step;
  return Math.sqrt(sigR * sigR + px * px / 12);
}


/* ---------- min/max band ----------
   A spectrum panel shows the min-max range within each screen sample, drawn as
   vertical strokes: one object, not a band outline (which would read as two
   lines) and not a mean (which is invariant to any smoothing finer than a
   pixel, so it made the resolution control look dead when zoomed out). */
/* Running-sum box blur: O(n) whatever the width, which is what lets the
   resolution slider cost the same at R = 2,000 as at R = 100,000.  Measured at
   ~0.3 ms for 80,000 points, so sixteen species is ~5 ms a frame.  A hoisted-
   clamp, buffer-reusing version was tried and was 2x SLOWER at the moderate
   sigma that matters most, so this stays. */
function boxBlur(src, half) {
  const n = src.length, out = new Float32Array(n), w = 2 * half + 1;
  let acc = 0;
  for (let i = -half; i <= half; i++) acc += src[Math.min(n - 1, Math.max(0, i))];
  for (let i = 0; i < n; i++) {
    out[i] = acc / w;
    acc += src[Math.min(n - 1, i + half + 1)] - src[Math.min(n - 1, Math.max(0, i - half))];
  }
  return out;
}
function directGauss(src, sig) {
  const half = Math.max(1, Math.ceil(3 * sig)), n = src.length;
  const w = []; let s = 0;
  for (let t = -half; t <= half; t++) { const q = Math.exp(-0.5 * (t / sig) ** 2); w.push(q); s += q; }
  const out = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    let a = 0;
    for (let t = -half; t <= half; t++) a += src[Math.min(n - 1, Math.max(0, i + t))] * w[t + half];
    out[i] = a / s;
  }
  return out;
}
/* Three box passes approximate a Gaussian, but equal integer widths quantise
   sigma in steps of about one point -- which is why dragging the resolution
   used to jump.  Mixing two widths (Kutskir) makes the effective sigma
   near-continuous; below ~1 point a direct kernel is cheap and exact. */
function boxesForGauss(sig, n) {
  const wIdeal = Math.sqrt((12 * sig * sig / n) + 1);
  let wl = Math.floor(wIdeal); if (wl % 2 === 0) wl--;
  const wu = wl + 2;
  const m = Math.round((12 * sig * sig - n * wl * wl - 4 * n * wl - 3 * n) / (-4 * wl - 4));
  return Array.from({ length: n }, (_, i) => (i < m ? wl : wu));
}
function gaussApprox(src, sig) {
  if (sig < 0.05) return src;
  if (sig < 1.2) return directGauss(src, sig);
  let cur = src;
  for (const w of boxesForGauss(sig, 3)) cur = boxBlur(cur, (w - 1) / 2);
  return cur;
}

/* Convolve, then min/max.  Always.
   Smoothing needs the real samples, so when a resolution is asked for we pull
   the native array and convolve that -- no shrink factors, no tier-dependent
   approximation, and nothing that can produce a value the data cannot reach.
   The convolved array depends only on (species, R), so it is cached and
   panning or zooming costs nothing.
   With R off, the overview's stored min/max IS the exact min/max of the raw
   samples over whole bins, so zoomed-out views need no native data at all. */
const SM = {};
function smoothedNative(name, R) {
  const w = WHOLE[name];
  if (!w) return null;
  const hit = SM[name];
  if (hit && hit.R === R) return hit.arr;
  const k = dec(null, name);
  const src = new Float32Array(w.length);
  for (let i = 0; i < w.length; i++) src[i] = w[i] * k;
  const arr = gaussApprox(src, (1 / R) / M.dln / 2.35482);
  SM[name] = { R, arr: arr === src ? src : Float32Array.from(arr) };
  return SM[name].arr;
}

let CB = { key: null, arr: null };
function combinedNative(names, R) {
  const key = R + '|' + names.join(',');
  if (CB.key === key) return CB.arr;
  const n = M.n, out = new Float32Array(n).fill(1);
  for (const nm of names) {
    const w = WHOLE[nm];
    if (!w) return null;                       // still loading
    const k = dec(null, nm);
    for (let i = 0; i < n; i++) out[i] *= w[i] * k;
  }
  const sm = gaussApprox(out, (1 / R) / M.dln / 2.35482);
  CB = { key, arr: sm === out ? out : Float32Array.from(sm) };
  return CB.arr;
}

function minMax(get, n, i0, i1, NS) {
  const lo = new Float32Array(NS), hi = new Float32Array(NS), mid = new Float32Array(NS);
  for (let x = 0; x < NS; x++) {
    const a = i0 + (i1 - i0) * x / NS, b = i0 + (i1 - i0) * (x + 1) / NS;
    let mn = Infinity, mx = -Infinity, s = 0, c = 0;
    const j1 = Math.min(n, Math.max(Math.ceil(b), Math.floor(a) + 1));
    for (let j = Math.max(0, Math.floor(a)); j < j1; j++) {
      const v = get(j);
      if (v < mn) mn = v; if (v > mx) mx = v; s += v; c++;
    }
    lo[x] = c ? mn : NaN; hi[x] = c ? mx : NaN; mid[x] = c ? s / c : NaN;
  }
  return { lo, hi, mid };
}

function band(name, i0, i1, NS, useFull, R) {
  const arr = smoothedNative(name, R);
  if (!arr) return envelope(name, i0, i1, NS, useFull);   // native still loading
  return minMax(j => arr[j], arr.length, i0, i1, NS);
}

/* Right-aligned text with subscript runs: [{t:'F'}, {t:'\u03bb', sub:true}, ...].
   Canvas has no markup, so measure the whole run first, then place each piece. */
function runsWidth(g, runs, fs) {
  let w = 0;
  for (const r of runs) {
    g.font = `${(r.sub ? fs * 0.72 : fs)}px Charter, Georgia, serif`;
    w += g.measureText(r.t).width;
  }
  return w;
}
function subText(g, runs, xRight, y, fs) {
  const font = s => `${s}px Charter, Georgia, serif`;
  const small = fs * 0.72, drop = fs * 0.2;
  let w = 0;
  for (const r of runs) { g.font = font(r.sub ? small : fs); w += g.measureText(r.t).width; }
  const prev = g.textAlign;
  g.textAlign = 'left';
  let x = xRight - w;
  for (const r of runs) {
    g.font = font(r.sub ? small : fs);
    g.fillText(r.t, x, r.sub ? y + drop : y);
    x += g.measureText(r.t).width;
  }
  g.textAlign = prev;
}

function tickLabel(g, txt, x, y, lo, hi) {
  const w = g.measureText(txt).width / 2;
  g.textAlign = (x - w < lo) ? 'left' : (x + w > hi) ? 'right' : 'center';
  g.fillText(txt, g.textAlign === 'left' ? lo : g.textAlign === 'right' ? hi : x, y);
  g.textAlign = 'center';
}

function yTicks(lo, hi) {
  const step = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1].find(s => (hi - lo) / s <= 4) || 1;
  const out = [];
  for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + 1e-9; v += step) out.push(v);
  return { ticks: out, dp: step < 0.1 ? 2 : 1 };
}

function axisTicks(a0, a1) {
  const span = a1 - a0;
  const step = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000].find(s => span / s < 11) || 1000;
  const out = [];
  for (let v = Math.ceil(a0 / step) * step; v <= a1; v += step) out.push(v);
  return out;
}

/* ---------- drawing ---------- */
const PANEL_H = 132;              // every panel the same height
function panels() {
  const p = [{ key: '_fluxpanel', h: PANEL_H, label: 'flux' },
             { key: '_norm', h: PANEL_H, label: 'normalized' }];
  if (state.mode === 'combined') {
    const sel = M.series.filter(s => state.on.has(s.name)).map(s => s.name);
    if (sel.length) p.push({ key: '_combined', h: PANEL_H, combined: sel });
  } else {
    for (const s of M.series) if (state.on.has(s.name)) p.push({ key: s.name, h: PANEL_H, sp: s });
  }
  return p;
}
function fmt(x) { return x.toLocaleString('en-US'); }

function draw() {
  const wrap = $('#plotwrap'), cv = $('#plot');
  const W = Math.max(320, wrap.clientWidth - 44);
  const P = panels();
  const padL = 100, padR = 14, padT = 6, padB = 4, gap = 17;
  const H = P.reduce((a, q) => a + q.h + gap, 0) + padT + padB;
  /* WebKit caps canvas area near 16.8 M px^2 and hands back a dead context
     past it -- with every species on we would ask for 19.2 M.  Back the
     supersampling off until it fits rather than losing the canvas. */
  let ss = DPR;
  while (ss > 1 && (W + padL + padR) * ss * H * ss > 12e6) ss -= 0.25;
  cv.width = (W + padL + padR) * ss; cv.height = H * ss;
  cv.style.height = H + 'px';
  const g = cv.getContext('2d');
  g.setTransform(ss, 0, 0, ss, 0, 0);
  g.clearRect(0, 0, W + padL + padR, H);
  /* The min/max band is a zigzag that reverses direction at every sample.
     Canvas's default miter join extends such a reversal by up to miterLimit
     (10) times the line width, which showed up as spikes shooting above the
     continuum -- pure geometry, not data: the convolution never exceeds
     0.99999. Round joins have no such overshoot. */
  g.lineJoin = 'round';
  g.lineCap = 'round';
  const css = getComputedStyle(document.body);
  const ink = css.getPropertyValue('--ink').trim() || INK;
  const line = css.getPropertyValue('--line').trim() || '#ddd';
  const atom = css.getPropertyValue('--atom').trim() || ATOM;
  const mol = css.getPropertyValue('--mol').trim() || MOL;

  /* Reduce one sample per DEVICE pixel, not per CSS pixel: the backing store
     is already DPR times wider, so sampling at CSS resolution threw away half
     the detail on a retina display. */
  const NS = Math.max(1, Math.round(W * ss));
  const px = i => padL + (i + 0.5) / ss;    // x0 === padL for every panel
  const useFull = (state.i1 - state.i0) <= FULL_MAX;
  let y = padT;
  const FS = 15;
  g.font = `${FS}px Charter, Georgia, serif`;
  for (const q of P) {
    const x0 = padL, y0 = y, h = q.h;
    g.strokeStyle = line; g.lineWidth = 1;
    g.strokeRect(x0 + .5, y0 + .5, W, h);
    const isFlux = q.key === '_fluxpanel';
    const vmax = isFlux ? M.flux_max / 1e7 : 1.08;
    const vmin = isFlux ? 0 : state.ymin;
    const Y = v => y0 + h - ((v - vmin) / (vmax - vmin)) * h;

    // y ticks
    g.fillStyle = css.getPropertyValue('--muted').trim() || '#777';
    g.textAlign = 'right'; g.textBaseline = 'middle';
    const tk = isFlux ? { ticks: [0, 0.5, 1.0], dp: 1 } : yTicks(vmin, 1.0);
    for (const t of tk.ticks) {
      const yy = Y(t); if (yy < y0 - 1 || yy > y0 + h + 1) continue;
      g.fillText(t.toFixed(tk.dp), x0 - 8, yy);
      g.strokeStyle = line; g.beginPath(); g.moveTo(x0, yy); g.lineTo(x0 + 4, yy); g.stroke();
    }

    // rotated y-axis label: F_lambda on the flux panel, "norm." on the rest
    g.save();
    g.fillStyle = css.getPropertyValue('--muted').trim() || '#777';
    g.translate(x0 - 56, y0 + h / 2);
    g.rotate(-Math.PI / 2);
    g.textBaseline = 'middle';
    const ylab = isFlux ? [{ t: 'F' }, { t: '\u03bb', sub: true }] : [{ t: 'norm.' }];
    subText(g, ylab, runsWidth(g, ylab, FS) / 2, 0, FS);
    g.restore();

    const series = isFlux ? ['_flux', '_cont'] : [q.key];
    const comb = q.combined || null;
    // clip to the PANEL BOX: the polyline deliberately runs a couple of points
    // past each edge so the trace enters and leaves cleanly, and those would
    // otherwise hang outside the axes.  Not a clip at y = 1.0 -- that halved
    // any trace sitting exactly on the continuum.
    g.save();
    g.beginPath(); g.rect(x0, y0, W, h); g.clip();
    for (const nm of series) {
      const col = nm === '_cont' ? AMBER
        : (q.sp ? (q.sp.kind === 'mol' ? mol : atom) : ink);
      const scale = (nm === '_flux' || nm === '_cont') ? 1 / 1e7 : 1;
      g.strokeStyle = col; g.lineWidth = nm === '_cont' ? 1.4 : 1;
      // the combined trace is already convolved, so it must not be smoothed
      // again on the polyline path
      const arrC = comb ? combinedNative(comb, state.R) : null;
      if (comb && !arrC) continue;                 // still loading
      if (nm === '_cont') {                       // the continuum is smooth
        const e = envelope(nm, state.i0, state.i1, NS, useFull);
        g.beginPath();
        let st = false;
        for (let i = 0; i < NS; i++) {
          if (isNaN(e.mid[i])) { st = false; continue; }
          const yy = Y(e.mid[i] * scale);
          st ? g.lineTo(px(i), yy) : g.moveTo(px(i), yy); st = true;
        }
        g.stroke();
      } else if (useFull && (state.i1 - state.i0) < 4 * NS) {
        /* Draw the data points themselves as a polyline whenever there are
           only a few per screen sample.  Reducing to min/max in that regime
           adds a vertical excursion at every sample, which reads as a ragged,
           undersampled line even though the grid carries 9-16 points across a
           line FWHM.  At 4x the sample count this is ~10,000 segments: free. */
        const kk = dec(null, nm);
        const a = Math.max(0, Math.floor(state.i0) - 1);
        const b = Math.min(M.n, Math.ceil(state.i1) + 2);
        let vals = new Float32Array(b - a);
        if (arrC) {
          for (let j = a; j < b; j++) vals[j - a] = arrC[j];
        } else {
          let last = 1;
          for (let j = a; j < b; j++) {
            let v = fullAt(nm, j);
            if (v < 0) v = last; else last = v;
            vals[j - a] = v * kk;
          }
          vals = gaussApprox(vals, (1 / state.R) / M.dln / 2.35482);
        }
        g.beginPath();
        for (let j = a; j < b; j++) {
          const x = padL + (j - state.i0) / (state.i1 - state.i0) * W;
          const yy = Y(vals[j - a] * scale);
          j === a ? g.moveTo(x, yy) : g.lineTo(x, yy);
        }
        g.stroke();
      } else {
        /* min-max within each screen sample, as a CONNECTED zigzag: isolated
           vertical strokes left gaps wherever the value moved more than the
           stroke between samples, which read as chunky. */
        const e = arrC
          ? minMax(j => arrC[j], arrC.length, state.i0, state.i1, NS)
          : band(nm, state.i0, state.i1, NS, useFull, state.R);
        g.beginPath();
        let st = false;
        for (let i = 0; i < NS; i++) {
          if (isNaN(e.lo[i])) { st = false; continue; }
          const yTop = Y(e.hi[i] * scale), yBot = Y(e.lo[i] * scale);
          st ? g.lineTo(px(i), yTop) : g.moveTo(px(i), yTop);
          g.lineTo(px(i), yBot);
          st = true;
        }
        g.stroke();
      }
    }
    g.restore();
    // label
    g.fillStyle = ink; g.textBaseline = 'bottom';
    const runs = isFlux ? null
      : q.key === '_norm' ? [{ t: 'all lines' }]
      : comb ? [{ t: `${comb.length} species multiplied` }]
      : [{ t: q.key }];
    if (runs) subText(g, runs, x0 + W - 8, y0 + h - 5, FS + 1);
    g.font = `${FS}px Charter, Georgia, serif`;
    y += h + gap;
  }
  const a0 = idxAir(state.i0), a1 = idxAir(state.i1);
  syncWave(a0, a1);
  const round = v => (v >= 1000 ? (Math.round(v / 100) * 100).toLocaleString('en-US')
                                : String(Math.round(v)));
  $('#yminlabel').textContent = `${state.ymin.toFixed(2)} – 1.0`;
  const waiting = P.some(q => (q.sp && !WHOLE[q.key])
    || (q.combined && q.combined.some(n => !WHOLE[n])));
  $('#reslabel').textContent = `R = ${round(state.R)}` + (waiting ? ' · loading…' : '');
  cv._geom = { padL, W, P, padT, gap, LW: W + padL + padR, LH: H };
  drawRuler(padL, W, padR);
  redrawTip();
}

function drawRuler(padL, W, padR) {
  const cv = $('#axis'), H = 64;
  cv.width = (W + padL + padR) * DPR; cv.height = H * DPR;
  cv.style.height = H + 'px';
  const g = cv.getContext('2d');
  g.setTransform(DPR, 0, 0, DPR, 0, 0);
  const css = getComputedStyle(document.body);
  const ink = css.getPropertyValue('--ink').trim() || INK;
  const line = css.getPropertyValue('--line').trim() || '#ddd';
  g.clearRect(0, 0, W + padL + padR, H);
  g.font = '18px Charter, Georgia, serif';
  g.fillStyle = css.getPropertyValue('--muted').trim() || '#777';
  g.textAlign = 'center'; g.textBaseline = 'top';
  const a0 = idxAir(state.i0), a1 = idxAir(state.i1);
  g.strokeStyle = line; g.lineWidth = 1;
  for (const v of axisTicks(a0, a1)) {
    const x = padL + (airIdx(v) - state.i0) / (state.i1 - state.i0) * W;
    g.beginPath(); g.moveTo(x, 0); g.lineTo(x, 7); g.stroke();
    tickLabel(g, fmt(v), x, 10, 2, W + padL + padR - 2);
  }
  g.textAlign = 'center';                 // same face, size and colour as the ticks
  g.fillText('wavelength  (air, \u00c5)', padL + W / 2, 40);
}

function drawBrush() {
  const cv = $('#brush'), W = cv.clientWidth, H = 54;
  cv.width = W * DPR; cv.height = H * DPR;
  const g = cv.getContext('2d'); g.setTransform(DPR, 0, 0, DPR, 0, 0);
  const css = getComputedStyle(document.body);
  g.clearRect(0, 0, W, H);
  const NB = Math.max(1, Math.round(W * DPR));
  const e = envelope('_norm', 0, M.n, NB, false);
  g.strokeStyle = css.getPropertyValue('--ink').trim() || INK; g.lineWidth = 1 / DPR;
  g.beginPath();
  for (let i = 0; i < NB; i++) {
    if (isNaN(e.lo[i])) continue;
    const x = (i + 0.5) / DPR;
    g.moveTo(x, H - e.hi[i] / 1.08 * H); g.lineTo(x, H - e.lo[i] / 1.08 * H);
  }
  g.stroke();
  const sel = brushSel || [state.i0, state.i1];
  let xa = sel[0] / M.n * W, xb = sel[1] / M.n * W;
  if (xb - xa < 4) { const m = (xa + xb) / 2; xa = m - 2; xb = m + 2; }  // stay grabbable
  g.fillStyle = 'rgba(120,120,120,.30)';
  g.fillRect(0, 0, xa, H); g.fillRect(xb, 0, W - xb, H);
  g.strokeStyle = css.getPropertyValue('--accent').trim() || '#8E2F6E';
  g.lineWidth = 1.6; g.strokeRect(xa, 1, xb - xa, H - 2);
}

/* ---------- interaction ---------- */
async function ensure() {
  const useFull = (state.i1 - state.i0) <= FULL_MAX;
  const need = ['_flux', '_cont', '_norm'];
  for (const s of M.series) if (state.on.has(s.name)) need.push(s.name);
  await Promise.all(need.map(n => loadOv(n)));
  await Promise.all(need.map(n => getWhole(n)));   // every R convolves real samples
}
let running = false, dirty = false, raf = null, hashT = null;
function scheduleDraw() {
  if (raf) return;
  raf = requestAnimationFrame(() => {
    raf = null;
    try { draw(); drawBrush(); } catch (err) { fail(err); }
  });
}
/* Coalesce, never drop.  The old guard returned early while a fetch was in
   flight, so a view change during the wait never triggered another ensure() --
   the draw that followed could be missing data with nothing to re-request it,
   which is what made zooming back out feel stuck. */
async function refresh() {
  dirty = true;
  if (running) return;
  running = true;
  while (dirty) {
    dirty = false;
    scheduleDraw();                       // paint with whatever is already here
    try { await ensure(); } catch (err) { /* keep the loop alive */ }
    scheduleDraw();
  }
  running = false;
}
function deferHash() {
  clearTimeout(hashT);
  hashT = setTimeout(writeHash, 250);
}
/* The URL carries the view, so a link can point at one line. */
function writeHash() {
  const a = idxAir(state.i0).toFixed(2), b = idxAir(state.i1).toFixed(2);
  const on = M.series.filter(s => state.on.has(s.name)).map(s => s.name).join(',');
  const p = [`w=${a}-${b}`];
  if (DB !== STARS[0].dir) p.push('star=' + encodeURIComponent(STARS.find(x => x.dir === DB).name));
  if (on) p.push('s=' + encodeURIComponent(on));
  p.push('R=' + state.R);
  if (state.mode === 'combined') p.push('m=c');
  if (state.ymin) p.push('y=' + state.ymin.toFixed(2));
  // Safari caps replaceState at 100 calls per 10 s and throws past it; the
  // URL is a convenience, so never let it take the render down with it
  try { history.replaceState(null, '', '#' + p.join('&')); } catch (e) { /* ignore */ }
}
function readHash() {
  const h = new URLSearchParams(location.hash.slice(1));
  if (h.has('s')) {
    state.on.clear();
    decodeURIComponent(h.get('s')).split(',').filter(Boolean).forEach(n => state.on.add(n));
  }
  if (h.has('R')) { state.R = +h.get('R'); const el = $('#res'); if (el) el.value = rToPos(state.R); }
  if (h.get('m') === 'c') state.mode = 'combined';
  if (h.has('y')) { state.ymin = +h.get('y'); const el = $('#ymin'); if (el) el.value = Math.round(state.ymin * 100); }
  if (h.has('w')) {
    const m = h.get('w').split('-').map(Number);
    if (m.length === 2 && isFinite(m[0]) && isFinite(m[1]) && m[1] > m[0]) {
      return [airIdx(m[0]), airIdx(m[1])];
    }
  }
  return null;
}

function setView(i0, i1) {
  // 40 native points is ~0.67 A at 5000 A: a single line fills the panel and
  // the individual samples become visible as vertices
  const span = Math.max(40, Math.min(M.n, i1 - i0));
  i0 = Math.max(0, Math.min(M.n - span, i0));
  state.i0 = i0; state.i1 = i0 + span;
  deferHash();          // WebKit throws past 100 replaceState calls / 10 s
  refresh();
}

function syncWave(a0, a1) {
  const el = $('#wrange'); if (!el) return;
  const w = el.clientWidth;
  const a = state.i0 / M.n, b = state.i1 / M.n;
  $('#hlo').style.left = (a * w) + 'px';
  $('#hhi').style.left = (b * w) + 'px';
  const f = $('#wfill');
  f.style.left = (a * w) + 'px';
  f.style.width = Math.max(1, (b - a) * w) + 'px';
  // don't overwrite a field while it is being typed into
  const lo = $('#wlo'), hi = $('#whi');
  if (document.activeElement !== lo) lo.value = a0.toFixed(a0 > 1000 ? 1 : 2);
  if (document.activeElement !== hi) hi.value = a1.toFixed(a1 > 1000 ? 1 : 2);
}

function applyWaveInputs() {
  const lo = parseFloat($('#wlo').value), hi = parseFloat($('#whi').value);
  if (!isFinite(lo) || !isFinite(hi)) { draw(); return; }     // put the old values back
  let a = airIdx(Math.min(lo, hi)), b = airIdx(Math.max(lo, hi));
  a = Math.max(0, Math.min(M.n - 40, a));
  b = Math.max(a + 40, Math.min(M.n, b));
  setView(a, b);
}

function wireWave() {
  const el = $('#wrange');
  let mode = null, grab = 0;
  const posOf = ev => {
    const r = el.getBoundingClientRect();
    return Math.min(Math.max((ev.clientX - r.left) / r.width, 0), 1) * M.n;
  };
  const start = m => ev => {
    mode = m;
    grab = posOf(ev) - state.i0;
    el.setPointerCapture(ev.pointerId);
    ev.preventDefault(); ev.stopPropagation();
  };
  $('#hlo').addEventListener('pointerdown', start('lo'));
  $('#hhi').addEventListener('pointerdown', start('hi'));
  $('#wfill').addEventListener('pointerdown', start('pan'));
  el.addEventListener('pointermove', ev => {
    if (!mode) return;
    const p = posOf(ev);
    if (mode === 'lo') setView(Math.min(p, state.i1 - 40), state.i1);
    else if (mode === 'hi') setView(state.i0, Math.max(p, state.i0 + 40));
    else { const span = state.i1 - state.i0; setView(p - grab, p - grab + span); }
  });
  const end = () => { mode = null; };
  el.addEventListener('pointerup', end);
  el.addEventListener('pointercancel', end);

  for (const id of ['#wlo', '#whi']) {
    const f = $(id);
    f.addEventListener('focus', () => f.select());
    f.addEventListener('change', applyWaveInputs);          // fires on Enter and on blur
    f.addEventListener('keydown', ev => {
      if (ev.key === 'Enter') { applyWaveInputs(); f.blur(); }
      if (ev.key === 'Escape') { f.blur(); draw(); }
    });
  }
}

let brushSel = null;              // live selection while dragging out a new range
function wireBrush() {
  const cv = $('#brush');
  let mode = null, anchor = 0, grab = 0;
  const EDGE = 6;                 // px within which an edge is grabbable
  const wid = () => cv.getBoundingClientRect().width;
  const xToIdx = ev => {
    const r = cv.getBoundingClientRect();
    return Math.min(Math.max((ev.clientX - r.left) / r.width, 0), 1) * M.n;
  };
  // hit box, widened so both edges stay separately grabbable at deep zoom
  const boxPx = () => {
    let xa = state.i0 / M.n * wid(), xb = state.i1 / M.n * wid();
    if (xb - xa < 16) { const m = (xa + xb) / 2; xa = m - 8; xb = m + 8; }
    return [xa, xb];
  };
  const hit = ev => {
    const x = ev.clientX - cv.getBoundingClientRect().left;
    const [xa, xb] = boxPx();
    if (Math.abs(x - xa) <= EDGE) return 'lo';
    if (Math.abs(x - xb) <= EDGE) return 'hi';
    return (x > xa && x < xb) ? 'pan' : 'new';
  };

  cv.addEventListener('pointermove', ev => {
    if (!mode) {                                   // hover: show what is grabbable
      const h = hit(ev);
      cv.style.cursor = (h === 'lo' || h === 'hi') ? 'ew-resize'
                      : h === 'pan' ? 'grab' : 'crosshair';
      return;
    }
    const p = xToIdx(ev);
    if (mode === 'lo') setView(Math.min(p, state.i1 - 40), state.i1);
    else if (mode === 'hi') setView(state.i0, Math.max(p, state.i0 + 40));
    else if (mode === 'pan') {
      const span = state.i1 - state.i0;
      setView(p - grab, p - grab + span);
    } else {                                       // dragging out a fresh range
      if (Math.abs(p - anchor) / M.n * wid() < 4) return;
      brushSel = [Math.min(anchor, p), Math.max(anchor, p)];
      drawBrush();
    }
  });

  cv.addEventListener('pointerdown', ev => {
    mode = hit(ev);
    anchor = xToIdx(ev);
    grab = anchor - state.i0;
    cv.setPointerCapture(ev.pointerId);
    ev.preventDefault();
    if (mode === 'pan') cv.style.cursor = 'grabbing';
  });

  const finish = ev => {
    if (!mode) return;
    if (mode === 'new') {
      if (brushSel) setView(brushSel[0], brushSel[1]);
      else { const span = state.i1 - state.i0, p = xToIdx(ev); setView(p - span / 2, p + span / 2); }
    }
    brushSel = null; mode = null; cv.style.cursor = '';
    drawBrush();
  };
  cv.addEventListener('pointerup', finish);
  cv.addEventListener('pointercancel', () => { mode = null; brushSel = null; drawBrush(); });
}

/* ---------------------------------------------------------------------
   line identification

   Every line SYNTHE could place at 0.2% central depth or deeper, from the
   same gfall records the synthesis read.  Split into 500 A blocks so a
   hover pulls only the block it is over.  Blocks arrive asynchronously;
   the tooltip redraws itself when one lands.
   ------------------------------------------------------------------ */
/* The canvas is stretched to the wrapper by CSS while its logical width is
   padL + W + padR, so a client point is not a canvas point: everything that
   reads the mouse has to divide the difference out first. */
function local(ev) {
  const cv = $('#plot'), gm = cv._geom;
  if (!gm) return null;
  const r = cv.getBoundingClientRect();
  const sx = r.width / gm.LW, sy = r.height / gm.LH;
  return { x: (ev.clientX - r.left) / sx - gm.padL,
           y: (ev.clientY - r.top) / sy,
           sx, sy, r, gm };
}

let LIDX = null;
const LCAT = {}, LDEP = {}, LREQ = {};

/* The catalogue -- wavelength, species, log gf, chi, level labels -- is a
   property of the atom, so it is shared by every star and survives a switch.
   Only the two depth bytes per line are per star. */
function loadLineIndex() {
  fetch('data/lines/index.json').then(r => r.ok ? r.json() : null)
    .then(j => { LIDX = j; redrawTip(); }).catch(() => {});
}

function slice(u8, off, len) {
  return u8.buffer.slice(u8.byteOffset + off, u8.byteOffset + off + len);
}

async function fetchCat(b) {
  const r = await fetch(`data/lines/${b}.bin`);
  if (!r.ok) { LCAT[b] = null; return; }
  const raw = new Uint8Array(await r.arrayBuffer());
  const dv = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
  const n = dv.getUint32(4, true), lablen = dv.getUint32(8, true);
  const z = await inflate(raw.subarray(12));
  let o = 0;
  const dlam = new Uint32Array(slice(z, o, 4 * n)); o += 4 * n;
  const sp = z.slice(o, o + n); o += n;
  const gf = new Int16Array(slice(z, o, 2 * n)); o += 2 * n;
  const chi = new Uint16Array(slice(z, o, 2 * n)); o += 2 * n;
  const ll = new Uint16Array(slice(z, o, 2 * n)); o += 2 * n;
  const ul = new Uint16Array(slice(z, o, 2 * n)); o += 2 * n;
  const jl = z.slice(o, o + n); o += n;
  const ju = z.slice(o, o + n); o += n;
  const lab = new TextDecoder().decode(z.subarray(o, o + lablen)).split('\n');
  const lam = new Float64Array(n);
  let acc = b * 1e4;
  for (let i = 0; i < n; i++) { acc += dlam[i]; lam[i] = acc / 1e4; }
  LCAT[b] = { n, lam, sp, gf, chi, ll, ul, jl, ju, lab };
  redrawTip();
}

async function fetchDep(b, key) {
  const r = await fetch(`${DB}/lines/${b}.bin`);
  if (!r.ok) { LDEP[key] = null; return; }
  const raw = new Uint8Array(await r.arrayBuffer());
  const n = new DataView(raw.buffer, raw.byteOffset, raw.byteLength).getUint32(4, true);
  const z = await inflate(raw.subarray(8));
  LDEP[key] = { dp: z.slice(0, n), ds: z.slice(n, 2 * n) };
  redrawTip();
}

function getBlock(b) {
  const dk = DB + '/' + b;
  if (b in LCAT && dk in LDEP) {
    return (LCAT[b] && LDEP[dk]) ? { c: LCAT[b], d: LDEP[dk] } : null;
  }
  if (!LREQ[b]) { LREQ[b] = 1; fetchCat(b).catch(() => { LCAT[b] = null; }); }
  if (!LREQ[dk]) { LREQ[dk] = 1; fetchDep(b, dk).catch(() => { LDEP[dk] = null; }); }
  return undefined;                       // not here yet
}

/* lines within tol of a, restricted to the species in `allow` (null = any).

   Ranking is on the predicted central depth, not the measured one: two lines
   inside one blend share a measured depth and would tie, while the predicted
   depth is per line and separates them.  It is stored on a log scale, since
   ranking cares about ratios.  Distance from the cursor is folded in with a
   Gaussian of half the tolerance, so pointing at a line picks that line
   rather than whatever deeper line sits a few pixels away. */
function findLines(a, tol, allow) {
  if (!LIDX) return null;
  const step = LIDX.block, out = [];
  const L = Math.log10(LIDX.cut);
  let pending = false;
  const b0 = Math.floor((a - tol) / step) * step;
  const b1 = Math.floor((a + tol) / step) * step;
  for (let b = b0; b <= b1; b += step) {
    if (LIDX.blocks.indexOf(b) < 0) continue;
    const k = getBlock(b);
    if (k === undefined) { pending = true; continue; }
    if (!k) continue;
    const c = k.c, d = k.d;
    let lo = 0, hi = c.n;
    while (lo < hi) { const m = (lo + hi) >> 1; c.lam[m] < a - tol ? lo = m + 1 : hi = m; }
    for (let i = lo; i < c.n && c.lam[i] <= a + tol; i++) {
      const q = d.dp[i];
      if (!q) continue;                   // below this star's cut
      const name = LIDX.species[c.sp[i]];
      if (allow && !allow.has(name)) continue;
      out.push({ name, lam: c.lam[i], gf: c.gf[i] / 1000, chi: c.chi[i] / 1000,
                 dp: Math.pow(10, L * (1 - q / 255)), ds: d.ds[i] / 255,
                 jl: c.jl[i] / 2, ju: c.ju[i] / 2,
                 lo: c.lab[c.ll[i]], up: c.lab[c.ul[i]] });
    }
  }
  for (const h of out) {
    const u = (h.lam - a) / (0.5 * tol);
    h.score = h.dp * Math.exp(-0.5 * u * u);
  }
  out.sort((p, q) => q.score - p.score);
  return { hits: out, pending };
}

const esc = t => String(t).replace(/&/g, '&amp;').replace(/</g, '&lt;')
  .replace(/>/g, '&gt;');
const sgn = v => (v >= 0 ? '+' : '\u2212') + Math.abs(v).toFixed(3);
function levels(h) {
  if (!h.up) return '';
  // molecules carry a band designation with no lower label: one rotational
  // line out of thousands is not what anyone wants named
  if (!h.lo) return `<div class="lev">${esc(h.up)}</div>`;
  const j = v => ` <span class="j">J=${v % 1 ? v.toFixed(1) : v}</span>`;
  return `<div class="lev">${esc(h.lo)}${j(h.jl)} &#8594; ${esc(h.up)}${j(h.ju)}</div>`;
}

let lastHover = null;
function redrawTip() { if (lastHover) showTip(lastHover); }

function showTip(ev) {
  const L = local(ev), tip = $('#tip'), mark = $('#lmark');
  if (!L) return;
  const { x, y, sx, sy, r, gm } = L;
  if (x < 0 || x > gm.W || y < 0 || y > gm.LH) {
    lastHover = null; tip.style.display = 'none'; mark.style.display = 'none';
    return;
  }
  lastHover = { clientX: ev.clientX, clientY: ev.clientY };

  // which panel, and so which species may be named
  let yy = gm.padT, hit = null, ytop = yy;
  for (const q of gm.P) {
    if (y >= yy && y <= yy + q.h) { hit = q; ytop = yy; break; }
    yy += q.h + gm.gap;
  }
  const allow = !hit ? null
    : hit.sp ? new Set([hit.sp.name])
    : hit.combined ? new Set(hit.combined)
    : null;

  const i = state.i0 + x / gm.W * (state.i1 - state.i0);
  const lam = idxAir(i);
  const tol = Math.abs(idxAir(i + 3 * (state.i1 - state.i0) / gm.W) - lam);
  const res = findLines(lam, tol, allow);

  let txt = `<b>${lam.toFixed(3)} &#8491;</b>`;
  if (hit && hit.sp) txt += ` <span class="j">${esc(hit.sp.name)}</span>`;
  mark.style.display = 'none';
  if (res && res.hits.length) {
    const h = res.hits[0];
    const d = h.ds > 0 ? h.ds : h.dp;   // measured where we have it
    txt = `<div class="hd"><b>${esc(h.name)}</b> ${h.lam.toFixed(3)} &#8491;</div>`
        + `<div class="par">log <i>gf</i> ${sgn(h.gf)} &middot; `
        + `&chi; ${h.chi.toFixed(3)} eV &middot; depth ${d.toFixed(2)}</div>`
        + levels(h);
    const rest = res.hits.slice(1, 4);
    if (rest.length) {
      txt += '<div class="also">also here: ' + rest.map(q =>
        `${esc(q.name)} ${q.lam.toFixed(3)}`).join(' &middot; ') + '</div>';
    }
    // a rule at the identified line, so it is obvious which feature was named
    if (hit) {
      const mx = (gm.padL + (airIdx(h.lam) - state.i0) / (state.i1 - state.i0) * gm.W) * sx;
      const wrap = $('#plotwrap').getBoundingClientRect();
      mark.style.display = 'block';
      mark.style.left = (r.left - wrap.left + mx) + 'px';
      mark.style.top = (r.top - wrap.top + ytop * sy) + 'px';
      mark.style.height = (hit.h * sy) + 'px';
    }
  } else if (res && res.pending) {
    txt += '<div class="par">identifying\u2026</div>';
  }
  tip.innerHTML = txt;
  tip.style.display = 'block';
  const tw = tip.offsetWidth || 220;
  tip.style.left = Math.min(ev.clientX + 14, window.innerWidth - tw - 10) + 'px';
  tip.style.top = Math.min(ev.clientY + 16,
                           window.innerHeight - tip.offsetHeight - 8) + 'px';
}

function wirePlot() {
  const cv = $('#plot');
  // Plain wheel zooms.  This does mean the page will not scroll while the
  // cursor is over the plot -- scroll from the rail or the margins instead.
  cv.addEventListener('wheel', ev => {
    ev.preventDefault();
    const L = local(ev);
    if (!L) return;
    const f = Math.min(Math.max(L.x / L.gm.W, 0), 1);
    const c = state.i0 + f * (state.i1 - state.i0);
    // deltaMode 1 is lines and 2 is pages; without normalising, a line-mode
    // wheel reports ~3 per notch and zooming out took ~1,900 events
    const unit = ev.deltaMode === 1 ? 16 : ev.deltaMode === 2 ? 400 : 1;
    const k = Math.exp(Math.max(-0.6, Math.min(0.6, ev.deltaY * unit * 0.006)));
    const span = (state.i1 - state.i0) * k;
    setView(c - f * span, c - f * span + span);
  }, { passive: false });
  cv.addEventListener('dblclick', () => setView(0, M.n));
  let drag = null;
  cv.addEventListener('mousedown', ev => { drag = { x: ev.clientX, i0: state.i0, i1: state.i1 }; });
  window.addEventListener('mouseup', () => { drag = null; });
  window.addEventListener('mousemove', ev => {
    const gm = cv._geom; if (!gm) return;
    if (drag) {
      const sx = cv.getBoundingClientRect().width / gm.LW;
      const d = (ev.clientX - drag.x) / sx / gm.W * (drag.i1 - drag.i0);
      lastHover = { clientX: ev.clientX, clientY: ev.clientY };
      setView(drag.i0 - d, drag.i1 - d);   // draw() re-runs the identification
      return;
    }
    showTip(ev);
  });
  cv.addEventListener('mouseleave', () => {
    lastHover = null;
    $('#tip').style.display = 'none';
    $('#lmark').style.display = 'none';
  });
}

function buildChips() {
  const box = $('#species');
  box.innerHTML = '';
  for (const gname of M.groups) {                 // grouped in order, unlabelled
    for (const s of M.series.filter(q => q.group === gname)) {
      const el = document.createElement('label');
      el.className = 'chip' + (s.kind === 'mol' ? ' mol' : '') + (state.on.has(s.name) ? ' on' : '');
      el.innerHTML = `<span class="sw"></span>${s.name}`;
      el.addEventListener('click', () => {
        state.on.has(s.name) ? state.on.delete(s.name) : state.on.add(s.name);
        el.classList.toggle('on');
        deferHash(); refresh();
      });
      box.appendChild(el);
    }
  }
}

const num = v => v.toFixed(2).replace('-', '\u2212');   // typographic minus
const rat = v => (v > 0 ? '+' : '') + num(v);          // abundance ratios signed

function buildStars() {
  const box = $('#stars');
  box.innerHTML = '';
  STARS.forEach((st, i) => {
    const b = document.createElement('button');
    const on = st.dir === DB;
    b.className = 'star' + (on ? ' on' : '') + (st.dir ? '' : ' soon');
    b.setAttribute('role', 'radio');
    b.setAttribute('aria-checked', on ? 'true' : 'false');
    // the Sun's numbers come from the synthesis, not from this table
    const m = on && M ? M.model : st;
    b.innerHTML = `${st.dir ? '' : '<span class="tag">soon</span>'}`
      + `<span class="nm">${st.name}</span><span class="sp">${st.sp}</span>`
      // two lines: all four numbers on one wraps at this card width
      + `<div class="pr">T<sub>eff</sub> ${m.teff} K &middot; log <i>g</i> `
      + `${num(m.logg)}</div>`
      + `<div class="pr2">[Fe/H] ${rat(m.feh)} &middot; [&alpha;/Fe] `
      + `${rat(m.afe || 0)}</div>`;
    b.addEventListener('click', () => selectStar(i));
    box.appendChild(b);
  });
}

/* Switching star swaps the whole data root, so every cache keyed on a series
   name has to go with it -- otherwise the new star draws the old one's arrays.
   The view and the species selection carry over: the wavelength grid and the
   species set are the same for every star. */
async function selectStar(i) {
  const st = STARS[i];
  if (!st.dir) {
    $('#starnote').textContent =
      `${st.name} has not been computed yet.`;
    return;
  }
  $('#starnote').textContent = '';
  if (st.dir === DB) return;
  DB = st.dir;
  OV = {}; CH = {}; WHOLE = {};
  for (const k of Object.keys(INFLIGHT)) delete INFLIGHT[k];
  for (const k of Object.keys(SM)) delete SM[k];
  CB = { key: null, arr: null };
  invalidateChunkCache();
  // LCAT is star-independent and deliberately kept across the switch
  for (const k of Object.keys(HDR)) delete HDR[k];
  const f = state.i0 / M.n, g = state.i1 / M.n;      // keep the window
  M = await (await fetch(`${DB}/meta.json`)).json();
  state.i0 = f * M.n; state.i1 = g * M.n;
  buildStars(); buildChips(); loadLineIndex();
  deferHash();
  await refresh();
}

/* the rail sticks below the controls bar, whose height depends on wrapping */
function syncRail() {
  const h = $('#controls').getBoundingClientRect().height;
  document.documentElement.style.setProperty('--ctrlh', h + 'px');
}

async function main() {
  // the star has to be settled before anything is fetched
  const h0 = new URLSearchParams(location.hash.slice(1));
  if (h0.has('star')) {
    const st = STARS.find(x => x.name === decodeURIComponent(h0.get('star')));
    if (st && st.dir) DB = st.dir;
  }
  M = await (await fetch(`${DB}/meta.json`)).json();
  M.default_on.forEach(n => state.on.add(n));
  state.i0 = 0; state.i1 = M.n;
  buildStars();
  const hv = readHash();
  buildChips();
  wireBrush(); wirePlot(); loadLineIndex();
  if (hv) { state.i0 = Math.max(0, hv[0]); state.i1 = Math.min(M.n, hv[1]); }
  $('#reset').addEventListener('click', () => setView(0, M.n));
  wireWave();
  $('#ymin').addEventListener('input', e => {
    state.ymin = +e.target.value / 100;      // upper limit stays pinned at 1
    scheduleDraw(); deferHash();
  });
  $('#res').addEventListener('input', e => {
    state.R = posToR(+e.target.value);
    refresh(); deferHash();       // may still need the native arrays

  });
  const setMode = m => {
    state.mode = m;
    $('#mindiv').classList.toggle('on', m === 'individual');
    $('#mcomb').classList.toggle('on', m === 'combined');
    deferHash(); refresh();
  };
  $('#mindiv').addEventListener('click', () => setMode('individual'));
  $('#mcomb').addEventListener('click', () => setMode('combined'));
  if (state.mode === 'combined') setMode('combined');

  $('#toggleall').addEventListener('click', () => {
    M.series.forEach(s => state.on.add(s.name)); buildChips(); refresh();
  });
  $('#togglenone').addEventListener('click', () => {
    state.on.clear(); buildChips(); refresh();
  });
  window.addEventListener('resize', () => { syncRail(); draw(); drawBrush(); });
  syncRail();
  await refresh();
}
main();
