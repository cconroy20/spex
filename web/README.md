# The Spectral Explorer

An interactive atlas of model stellar spectra, computed line by line and split
into the contribution of each atom and molecule — what a star would look like if
only Fe I, or Ca II, or CH absorbed. The Sun (3550–10000 Å) is the first model
loaded; the layout is built to take others.

**Live: https://USER.github.io/REPO/**

- Toggle any of 25 species on or off.
- Drag the strip under the header to zoom; scroll or drag on the plot to pan.
- Degrade to any resolving power from R = 300,000 down to R = 3,000.
- The URL tracks the view, so a link points at one line:
  `#w=5885-5900&s=Fe%20I,Na%20I`

## What is being plotted

A 1D LTE synthesis with Kurucz's ATLAS12 and SYNTHE on a solar model atmosphere
(Teff 5777 K, log g 4.44, [Fe/H] 0, v_turb 2 km/s), at R = 300,000, using the
`gfall` atomic line list and sixteen molecular line lists.

Each species panel is a **separate run of the radiative transfer** with only
that species' lines in the line list — the panels are computed, not estimated by
scaling or differencing. SYNTHE has no species filter, but it reads whatever
file `lines.list` names, so `gfall` was split by species code and each species
given a line list of its own.

## Caveats

- Strict LTE with S = B. The Ca II H and K cores go to zero and Hα is too
  shallow; both form in the chromosphere.
- Species panels do not sum to the total, because line blending is not additive.
- The H I panel steps at the Balmer limit (3646 Å): SYNTHE partitions hydrogen
  opacity there, lines redward and bound-free continuum blueward. Neither half
  is physical alone; their sum is continuous.

## How the site works

No server logic — it is static files, one directory per star under `data/`,
plus a shared `data/lines/` catalogue. Each series has two tiers:

| tier | contents |
|---|---|
| `data/<star>/ov/` | decimated ~5,974 bins, min/max/mean per bin, loaded up front |
| `data/<star>/full/` | native 310,691 points, fetched per series on zoom |

The wavelength grid is uniform in ln λ, so it is reconstructed in the browser
from `(lam0_vac, dln, n)` rather than stored, and converted to air wavelengths
with the IAU (Morton 2000) formula.

## Compression

Series are stored in the `SPC1` container written by `binfmt.py`. The grid is
critically sampled at R = 300,000, so neighbouring points differ by very
little: taking a first difference makes the high byte of each `uint16` nearly
constant, separating the byte planes puts those constants together, and
deflate finds them. The difference is taken mod 2¹⁶ and undone by a running
sum with the same wraparound, so this is exact, not lossy.

    magic   'SPC1'                      4 bytes
    n       uint32                      points in the series
    chunk   uint32                      points per chunk (16,384)
    nchunk  uint32
    offset  uint32 x (nchunk + 1)       byte offsets from the start of file
    payload deflate(hi_bytes || lo_bytes) per chunk

Each chunk is deflated on its own, which costs about 5% against compressing
the whole file and keeps HTTP Range landing on a chunk boundary. The browser
inflates with `DecompressionStream('deflate')` — Safari 16.4+, Chrome 80+,
Firefox 113+; there is no uncompressed fallback.

Measured on the shipped bundles: 13.1x for the Sun, 15.8x Procyon, 10.5x
Arcturus, 32.8x HD 122563. The ratio rises with the number of species,
because most species are flat in any one star and flat costs nothing.

`_flux` is not shipped at all. It is `_norm` x `_cont`, and rebuilding it in
the browser from the two quantized arrays agrees with the real thing to two
quantization steps — below anything that can be drawn. See `DERIVED` in
`app.js`.

`.nojekyll` is required: several data files begin with an underscore, and
GitHub Pages' Jekyll step would otherwise drop them.

`.nojekyll` is required: several data files begin with an underscore, and
GitHub Pages' Jekyll step would otherwise drop them.

## Adding a star

The header carries one card per star. Everything the app fetches hangs off
`DB`, the data root, so a second star is a directory of the same shape as
`data/` plus a row in the `STARS` table in `app.js` giving it a `dir`. The
three cards without one are placeholders, and their parameters are nominal
literature values — replace them with whatever the models are actually
computed at. The Sun's card reads its numbers from `meta.json`, not from the
table.

## Line identification

Hovering a panel names the line under the cursor. The index behind it is built
by `build_lines.py` and `pack_lines.py` from the same Kurucz records the
synthesis read — no separate catalogue, so a line the tooltip names is a line
that is actually in the spectrum.

There are 6.0 million lines in the band, far too many to ship, so each is
given a predicted central depth against the continuum, computed with the
single-depth machinery in `sunlib.py`, and only those above 0.2% survive.

The catalogue half — wavelength, species, log gf, χ, level labels — is a
property of the atom and identical for every star, so it is built once over
the union of every star's visible lines and shipped at `data/lines/`; each
star adds only two depth bytes per line at `data/<star>/lines/`. The union
turns out to be barely larger than the richest single star (229,759 against
Arcturus's 223,802), because the strong lines of an F dwarf and a metal-poor
giant are close to a subset of a cool giant's. That makes sharing nearly
free: 1.5 MB once, plus 0.08–0.32 MB per star, against 20.7 MB of per-star
JSON before.

The depth is evaluated at **two** reference depths, τ = 1 and τ = 0.1, and the
larger is kept. Ranking on τ = 1 alone (T = 6518 K) dissociates the molecules
and threw away CN lines that SYNTHE puts at 10% depth — CN completeness went
from 54% to 99.6% when the second depth was added. Against the synthesized
per-species spectra, every species now has ≥ 98.5% of its minima deeper than
2% matched to an indexed line within 0.03 Å.

Two depths are carried per line. The **predicted** depth ranks candidates:
lines inside one blend share a measured depth and would tie, while the
predicted depth is per line and separates them. The **measured** depth, read
off that species' synthesized spectrum, is what the tooltip reports, so the
number agrees with the trace on screen. Distance from the cursor is folded in
with a Gaussian of half the pick tolerance, so pointing at a line picks that
line rather than a deeper one a few pixels away.

Molecules are named by band, not by rotational line: `A–X (1,1)`, decoded from
the state letter and vibrational quantum numbers in the level codes. Hydrogen
records carry the literal text `AVERAGE ENERGIES` where the level labels go,
so those are named from the energies instead — `Balmer α (n = 2 → 3)`.

One trap: `#plot` is stretched to the wrapper by CSS while its logical width is
`padL + W + padR`, about 6% wider. Every handler that reads the mouse has to
divide that scale out first (`local()`), or the wavelength under the cursor is
wrong by up to 6% of the span at the right-hand edge.

## Rendering

Below ~4 data points per screen sample the panels draw the data points
themselves as a polyline. Reducing to min/max in that regime puts a vertical
excursion at every sample, which reads as an undersampled line even though the
grid carries 9-16 points across a line FWHM (measured: Fe I 6252 has 10,
Fe I 5232 has 16, against a Nyquist requirement of 2). Above that threshold,
each panel shows the **min-max range within one screen sample**, drawn as
vertical strokes: one object, rather than a band outline (which reads as two
lines) or a mean (which is invariant to any smoothing finer than a pixel, so
the resolution control appeared to do nothing when zoomed out). Samples are
taken per *device* pixel, not per CSS pixel, so a retina display gets the
factor of two it has paid for.

Choosing a resolution smooths before that reduction. There is only ONE
synthesis behind the slider -- nothing is interpolated between precomputed
resolutions; every R is a live convolution of the native R = 300,000 data with
a Gaussian of FWHM = c/R. Three box passes of mixed width stand in for the
Gaussian, which is O(n) whatever the width: measured in node on the real Fe I
array, 0.25 ms for 80,000 points and flat from sigma = 1.3 to 63.7 points, so
sixteen species cost ~4 ms a frame. The approximation costs accuracy -- up to
8e-3 in normalized flux against a true Gaussian at sigma = 6 -- which is
invisible on screen but would matter if these curves were ever measured;
a fourth pass would take it to ~2e-3.

Do not time this in headless Chrome with `--virtual-time-budget`: the clock is
frozen, `performance.now()` never advances, and every measurement comes back
0.00 ms (or worse, a plausible-looking wrong number). Use node. The
readout reports what is actually on screen: the limit is the data tier, not
the pixel, because a min-max band does show sub-pixel structure. The overview
tier bins 52 native points, so it caps at R ~ 5,800; below that the app must
zoom in far enough to pull native data.

## Publishing

    git init && git add -A && git commit -m "Interactive solar spectrum atlas"
    git branch -M main
    git remote add origin git@github.com:USER/REPO.git
    git push -u origin main

Then in the repository settings, Pages → Deploy from a branch → `main` / root.

## Credits

Synthesis: ATLAS12 and SYNTHE, R. L. Kurucz; `gfall` line list; Barklem &
Collet (2016) partition functions. Built from the `sunspec` pipeline.
