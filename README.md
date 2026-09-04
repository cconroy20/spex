# The Spectral Explorer

An interactive atlas of model stellar spectra, computed line by line and split
into the contribution of each atom and molecule — what a star would look like if
only Fe I, or Ca II, or CH absorbed. Five stars from an F dwarf to an M dwarf,
each over 3549–24993 Å.

Live: https://cconroy20.github.io/spex/

- Toggle any of 51 species on or off.
- Drag the strip under the header to zoom; scroll or drag on the plot to pan.
- Degrade to any resolving power from R = 300,000 down to R = 1,000.
- Overlay telluric transmission, a blackbody, formation depth, or photon noise
  at a chosen signal-to-noise.
- The URL tracks the view, so a link points at one line:
  `#w=5885.00-5900.00&s=41`

## What is being plotted

A 1D LTE synthesis with Kurucz's ATLAS12 and SYNTHE, at R = 300,000, using the
`gfall` atomic line list and sixteen molecular line lists. The five models:

| star | Teff | log g | [Fe/H] | [α/Fe] |
|---|---|---|---|---|
| Sun | 5777 | 4.44 | 0.00 | 0.0 |
| Procyon | 6530 | 3.96 | 0.00 | 0.0 |
| HD 122563 | 4587 | 1.61 | −2.64 | 0.4 |
| Arcturus | 4286 | 1.66 | −0.52 | 0.3 |
| Barnard's Star | 3220 | 5.05 | −0.40 | 0.0 |

Each species panel is a separate run of the radiative transfer with only that
species' lines in the line list — the panels are computed, not estimated by
scaling or differencing. SYNTHE has no species filter, but it reads whatever
file `lines.list` names, so `gfall` was split by species code and each species
given a line list of its own.

## The sequences

Two cards are not stars but grids along one axis, computed with the same
settings as each other (`mlt=1.25`, `vturb=2.0`, solar abundances) so they are
comparable. Everything but the named axis is held fixed, which is the whole
point of them. Neither has species panels, because the per-species runs have
not been done for the grids.

| card | axis | range | models | held |
|---|---|---|---|---|
| Teff Sequence | Teff | 2520–12650 K, 0.025 dex | 27 | log g 5.0 |
| log g Sequence | log g | 2.0–5.5, 0.5 dex | 8 x 3 | Teff 3500 / 4500 / 5500 K |

Both cards read the same way: the rail names the parameter being held and
shows the value it is pinned at, then the slider for the one that moves. The
log g card has three temperatures to choose between and the Teff card has one
gravity, which is still shown as a card so the two look like the same
instrument. The three log g sequences are three data roots sharing one card,
so switching between them is a star switch in everything but name.

Grid spacing was chosen by dropping every other model and interpolating it
from its neighbours. The two axes need very different densities:

| axis | bracket tested | rms |
|---|---|---|
| Teff | 0.05 dex | 0.0084 (0.0131 below 3570 K) |
| log g | 0.5 dex | 0.0025 at 3500 K, 0.0007 at 5500 K |

Gravity interpolates an order of magnitude better than temperature and needs a
grid twenty times coarser, because changing it rescales the pressure structure
without switching species on and off. It is not a small effect — rms 0.157
between log g 5.5 and 2.0 at 3500 K, against 0.32 for 3570 to 5650 K in
temperature — it is simply a smooth one.

## The temperature sequence

The sixth card is not a star but a grid of models along one axis. Teff runs
from 2520 to 12650 K at **fixed log g = 5.0 and solar abundances**, so the
slider changes exactly one thing, which is the whole point of it. There are 27
computed models, even in log T at a spacing of 0.025 dex, and the browser
interpolates between the two that bracket the slider.

The cool end is where the sequence earns its keep and where interpolation
might have been expected to fail: at 2520 K lines remove 86% of the flux,
against 21% for the Sun. It does not fail. The same held-out test gives rms
0.0131 over 2520-3570 K against 0.0084 over the rest, so the molecular forest
costs about half again as much error and nothing worse.

The two quantities are interpolated differently, because they behave
differently, and both choices were measured against models held out of the
grid rather than assumed:

| quantity | scheme | rms at 0.05 dex |
|---|---|---|
| normalized flux | linear in log T | 0.0084 |
| | log-linear in log T | 0.0097 |
| continuum | log-linear in log T | 0.0110 |
| | linear in log T | 0.0444 |

Lines saturate, so residual flux is not the exponential in temperature that a
single line's opacity is, and interpolating it in the log makes things worse.
The continuum spans 7000x across the sequence and has to be interpolated in
the log, which is why `_cont` ships as log10 already: lerping the stored value
IS the log interpolation, and the browser exponentiates once at the end.

Every model ships, not every second one, so the brackets the app actually
interpolates across are half the width of the ones in that table. The 0.0084
is an upper bound on what you see. The same held-out test at 0.075 dex gives
0.0171, so the error grows roughly with the square of the bracket and 0.025
dex should land near 0.002, which is a quarter of a pixel on a panel.

Two models between 7110 and 8460 K failed to converge in ATLAS12, so that one
bracket is three times wider than the rest. The app says so in the rail when
the slider is inside it.

Above 12650 K the grid stops. A 17870 K model exists and everything between it
and 12650 K segfaulted, leaving a 0.15 dex hole; `export_tgrid.py` drops any
node whose nearest neighbour is further than 0.10 dex rather than let the
slider interpolate across a gap that wide and draw a spectrum belonging to no
star. Below 2520 K the cost curve, not the physics, is what stops it: the
models converge but each iteration slows by more than an order of magnitude,
and a 2380 K model was abandoned after two hours at nine iterations of
thirty.

The sequence ships at R = 50,000 rather than 300,000, which is 3.8 MB for all
21 models against roughly 23 MB at native resolution, and the resolution
slider caps itself accordingly rather than offering a resolution the data does
not have. There are no species panels, because the per-species runs have not
been done for the grid.

## Caveats

- Strict LTE with S = B. The Ca II H and K cores go to zero and Hα is too
  shallow; both form in the chromosphere.
- Species panels do not sum to the total, because line blending is not additive.
- The H I panel steps at the Balmer limit (3646 Å): SYNTHE partitions hydrogen
  opacity there, lines redward and bound-free continuum blueward. Neither half
  is physical alone; their sum is continuous.

## Layout

Two lines of bar and a strip, then the spectrum: the first panel starts 197 px
down the page, against 580 before. The bar holds the title, one line of
framing and Reset, with the five stars as tabs beneath it — only the star on
view spells out its four parameters, since five cards of four numbers is a lot
of ink for something clicked once a session. Under that sit the navigator and
the two controls worth having to hand while reading, the wavelength range and
the resolution.

Everything else lives in the left rail under a heading: Species first, because
it is what you touch constantly, then Display, then Simulated observation.
Nothing is behind a popover or a disclosure — the rail is taller than the
window, so the settings are reached by scrolling it while the bar and the
navigator stay put.

The species list is grouped into its five families, and `buildChips` reads them
from `meta.json` rather than naming them itself, so the headings follow the
species table. Fifty-one identical pills was a list you read; five named
families is a list you scan.

## How the site works

No server logic — it is static files, one directory per star under `data/`,
plus a shared `data/lines/` catalogue and a shared `data/telluric*.bin`. Each
star carries 55 series (51 species, the continuum, the full synthesis, and two
formation-depth arrays) in two tiers:

| tier | contents |
|---|---|
| `data/<star>/ov/` | decimated to 11,261 bins of 52 points, min/max/mean per bin |
| `data/<star>/full/` | native 585,579 points |

The overview tier paints the first frame and carries a zoomed-out view on its
own, since its stored min/max is the exact min/max of the raw samples over
whole bins. It cannot represent a resolution finer than 52 native points,
which caps it at R ≈ 5,800.

Everything past that reads the native tier, and `ensure()` fetches it whole for
every selected series rather than by range: each setting of the resolution
control is a live convolution of the real samples, so the array is needed in
full whatever the zoom. First load is 2.5 MB for the Sun and 4.3 MB for
Barnard's Star, on top of the overview tier.

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
the whole file and leaves the chunk boundaries addressable. The browser
inflates with `DecompressionStream('deflate')` — Safari 16.4+, Chrome 80+,
Firefox 113+; there is no uncompressed fallback.

Measured on the shipped bundles:

| star | native tier | overview | line depths | ratio |
|---|---|---|---|---|
| Sun | 4.20 MB | 564 kB | 0.34 MB | 15.3x |
| Procyon | 3.61 MB | 513 kB | 0.22 MB | 17.9x |
| Arcturus | 4.95 MB | 621 kB | 0.53 MB | 13.0x |
| HD 122563 | 1.63 MB | 274 kB | 0.14 MB | 39.5x |
| Barnard's Star | 6.95 MB | 648 kB | 0.71 MB | 9.3x |

plus 3.56 MB of shared catalogue, for 30 MB and 828 files in all. The ratio
tracks how many species are flat in that star: HD 122563 is metal-poor and
compresses by 39.5x, Barnard's Star has a spectrum full of molecules and
compresses by 9.3x.

`_flux` is not shipped at all. It is `_norm` x `_cont`, and rebuilding it in
the browser from the two quantized arrays agrees with the real thing to two
quantization steps — below anything that can be drawn. See `DERIVED` in
`app.js`.

`.nojekyll` is required: several data files begin with an underscore, and
GitHub Pages' Jekyll step would otherwise drop them.

## Adding a star

The bar carries one tab per star. Everything the app fetches hangs off `DB`,
the data root, so a second star is a directory of the same shape as the five
under `data/` plus a row in the `STARS` table in `app.js` giving it a `dir`.
Build the directory with `run_star.sh`, `assemble.py`, `assemble_species.py`,
`assemble_form.py`, `export_web.py` and `pack_lines.py`; see the pipeline
README one level up. The tabs read their parameters from the `STARS` table,
except the star on show, which reads them from its `meta.json` and spells them
out inline.

## Simulated observation

A continuum signal-to-noise per pixel and a fixed number of pixels per
resolution element, drawn as a noisy curve on the flux and all-lines panels --
the two that correspond to something a telescope could record.  A species
panel is a decomposition, not an observable, so it stays clean.

A detector pixel is not the model grid.  The model is one sample per
resolution element at R = 300,000, so a spectrograph at resolution R with p
pixels per element steps `300000/(R*p)` model points.  That is only a
coarsening while `R*p < 300,000`; above it the pixel would be finer than the
model, neighbouring pixels would share model points, and the "noise" would
come out correlated and read as smooth wiggles.  `pixelGridOK` refuses that
case and the control says so rather than drawing something untrue.  One model
point is sampled per pixel rather than averaging the span -- cheaper, and the
difference is invisible once noise is on top.

The noise is photon noise, so the quoted S/N is per pixel IN THE CONTINUUM and
falls as sqrt(F) into a line: at S/N 100, a core at 10% of the continuum has
S/N 31.  A flat S/N across the line would make deep cores look far better
measured than they are.

The two noisy panels grow their upper limit as 4/(S/N) once that beats the 8%
they already leave above the continuum, so a noise excursion is never clipped:
with a thousand pixels across a panel the largest is routinely three and a
half sigma, and clipping it would hide exactly the scatter the control exists
to show.  Nothing else moves -- the formation and species panels keep their
own ranges.

Within one realization the noise is deterministic in the pixel's absolute
index, so panning and zooming carry it with the spectrum instead of
reshuffling it every frame.  The seed itself is random on every load and on
every press of `resample noise`, and is deliberately NOT carried in the URL:
a shared link reproduces the configuration, not the particular draw.

## Line identification

Hovering a panel names the line under the cursor. The index behind it is built
by `build_lines.py` and `pack_lines.py` from the same Kurucz records the
synthesis read — no separate catalogue, so a line the tooltip names is a line
that is actually in the spectrum.

There are millions of lines in the band, far too many to ship, so each is
given a predicted central depth against the continuum, computed with the
single-depth machinery in `sunlib.py`, and only those above 0.2% in at least
one of the five stars survive. That union runs to 488,665 lines, against
360,782 for Barnard's Star, the richest single star.

The catalogue half — wavelength, species, log gf, χ, level labels — is a
property of the atom and identical for every star, so it is built once over the
union and shipped at `data/lines/` in 44 blocks of 500 Å; each star adds only
two depth bytes per line at `data/<star>/lines/`. That makes sharing nearly
free: 3.56 MB once, plus 0.14–0.71 MB per star.

Coverage is 49 of the 51 species. TiO and H₂O are the two left out, and
deliberately: both ship as multi-gigabyte binaries with no text form, and
naming one rotational line out of 97 million is not the point.

Sr II, Ba II, Nd II and Eu II were missing until 2026-09-02, and the reason is
worth keeping. `LineOpacity.n_over_U` returned None for Z = 38, 56, 60 and 63,
so their predicted depth came out exactly zero for every line and
`pack_lines.py`, which culls on that quantity, dropped the species whole — even
though Sr II 4077.70 Å sits at 97.7% depth in its own panel. The cause was in
`atmlib.Atmosphere`: the deck lists all 99 elements, but the closing loop of
`_abundances` walked a 30-element table that stopped at Zn and discarded the
rest. `ELEMENTS` now runs to Z = 92. The bulk sums that build `mu_H` and
`ntot_over_nH` still cover only the elements with a tabulated mass, so `rho`
and every optical depth are unchanged to the last bit; the heavier elements
are carried for line opacity alone.

The depth is evaluated at two reference depths, τ = 1 and τ = 0.1, and the
larger is kept. Ranking on τ = 1 alone (T = 6518 K) dissociates the molecules
and threw away CN lines that SYNTHE puts at 10% depth — CN completeness went
from 54% to 99.6% when the second depth was added.

Two depths are carried per line. The measured depth, read off that species'
synthesized spectrum, ranks the candidates, and the predicted depth breaks the
ties: lines inside one blend share a measured depth, while the predicted depth
is per line and separates them. Predicted depth alone let a line the synthesis
puts at zero outrank a real one — SiH beat the CO bandhead at 2.3 μm that way.
Distance from the cursor is folded in with a Gaussian of half the pick
tolerance, so pointing at a line picks that line rather than a deeper one a few
pixels away.

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
each panel shows the min-max range within one screen sample, drawn as
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
the pixel, because a min-max band does show sub-pixel structure.

## Local testing

    python3 serve.py 8731

`serve.py` sends no-cache headers, without which a browser holds a stale
`app.js` against new markup and the page collapses in a way that looks like a
bug in the code being edited. It also honours Range, so the container's chunk
boundaries can still be exercised.

## Publishing

    ./deploy.sh

The bundle is ~30 MB of binaries rewritten wholesale on every export, so the
script rebuilds `gh-pages` from scratch as a single orphan commit and
force-pushes it: the branch never accumulates history, and `web/data/` stays
out of `main` entirely. It also stamps the current commit onto the `app.js`
and `style.css` URLs, because Pages serves HTML and assets with the same
ten-minute max-age and caches them independently — without the stamp a browser
can hold new markup against an old script, and the new controls appear with
nothing wired to them. Set Pages to deploy from the `gh-pages` branch at root.

## Credits

Synthesis: ATLAS12 and SYNTHE, R. L. Kurucz; `gfall` line list; Barklem &
Collet (2016) partition functions. Telluric transmission: ESO SkyCalc (Noll et
al. 2012; Jones et al. 2013).
