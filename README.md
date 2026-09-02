# The solar spectrum, 3600-10000 A

Pedagogical figure set, plus the pipeline behind the interactive atlas in
`web/`. The synthesis itself now runs 3549-24993 A for five stars; the figures
here are cut to 3600-10000 A and are the Sun alone.

Four quantities, all for the same ATLAS12 solar model (Teff 5777, log g 4.44,
v_turb 2 km/s; `starcfg.py` names the deck, and all five stars were regenerated
with the current ATLAS12 so every deck carries 80 layers and a TAU5000 column):

1. the fluxed spectrum,
2. the continuum-normalized spectrum,
3. the line-free continuum,
4. the spectrum species by species -- what the Sun would look like if only
   Fe I, or Ca II, or CH absorbed.

## Figures

| file | contents |
|---|---|
| `fig/sun_overview.pdf` | (a) F_lambda with the continuum over it, (b) F/F_cont |
| `fig/sun_species.pdf`  | 17 stacked panels: all lines, then one species each |
| `fig/sun_zoom.pdf`     | six named windows at full resolution, with the observed Sun |
| `fig/sun_poster.pdf`   | print version: the two overview panels stacked on top of the species waterfall, **30 x 40 in exactly (3:4)**, native R = 300,000 |
| `fig/sun_poster_R10k.pdf` | the same sheet smoothed to R = 10,000 |

`figure_poster()` takes `width` (inches), `row` (height of one species panel),
`nb` (bins), `fs` (type scale), `lw` (stroke weight), `title` and `title_fs`,
so the sheet
can be re-cut to any size in a few seconds. The species panels are drawn as lines,
like panels (a) and (b); a shaded region reads as a dark mass rather than as a
spectrum. `fill=True` restores the filled version. There is no `all lines` row
on the poster: panel (b) is already the full synthesis, and the two were the
same array to an rms of 1.4e-4 (`assemble_species.py` builds `all` by
interpolating panel (b) onto the species grid). That row went to Mn I instead.
Panels are grouped by family -- neutrals and Ca II, iron-peak neutrals, ions,
molecules -- and ordered by contribution within each group, so the strongest
member leads it. Every drawn line, including the continuum, the spines and the
tick marks, carries the same weight. The defaults (width 30.811, row
1.5719) are solved so that the page trims to exactly 30.00 x 40.00 in -- a
standard poster size, so it prints without scaling. The canvas is slightly
larger than the trimmed page, which is why those numbers are not round: the
title sits above the axes and `bbox_inches='tight'` pulls the edges in.
Height is not a free parameter at fixed `row` -- it follows from the row count. Each species profile is stroked as well as
filled: at 0.10 A bins a single-bin line is a hairline on paper, and the stroke
gives every feature a visible minimum width without diluting its depth.

`smooth=R` convolves everything with a Gaussian of FWHM = c/R first. Both grids
are uniform in ln(lambda), so a fixed-width kernel in index space is a fixed
width in velocity. Flux and continuum are smoothed separately and then divided,
which is what an instrument does; the species panels hold only the ratio, but
the continuum varies by <0.1% across one resolution element at R = 10,000, so
smoothing the ratio is the same thing. `nb` drops to 30,000 for the smoothed
version -- three samples per resolution element is all it can carry.
Checks: equivalent widths are conserved to 0.1-0.3% (Na D 1.6192 -> 1.6152 A,
Mg b 8.1726 -> 8.1571 A), and isolated lines broaden as intrinsic width and
lambda/R in quadrature. Deepest points go Ca II K 0.000 -> 0.029, Mg b 0.073 ->
0.266, Na D 0.124 -> 0.412, Halpha 0.327 -> 0.435. It is vector, so there is no resolution to set -- what matters is
the physical size and that the type scales with it. The spectra are binned to
60,000 points (0.10 A, ~1500 samples per inch), past what any printer resolves,
which keeps the file at 12 MB rather than 60. Fonts are embedded as subsetted
TrueType, not Type 3. Molecular species are drawn in amber, atomic in blue.

## How each panel is computed

**Panels 1-3 are SYNTHE.** `run_star.sh` drives `$ATLAS12/bin/synthe.exe`
over 355-2500 nm at R = 300,000 in 86 chunks of 25 nm (chunking is for memory:
the full line list runs to ~10 million lines per chunk in the blue, 3.2 GB).
The blue end reaches past the Balmer break at 3646 A, where the continuum steps
by a factor of 1.211; the figures are cut at 3600 A, leaving margin below it.
Full line list -- gfall, Kurucz's predicted lines, and every molecule in
`data/lines.list`. `assemble.py` stitches the chunks, converts to air
wavelengths and to F_lambda = 4 pi H_nu c / lambda^2, and validates.

**Panel 4 is one SYNTHE run per species** (the second half of `run_star.sh`;
`run_species_synthe.sh` is the earlier solar-only driver that produced the
25-species cache the committed figures were drawn from). SYNTHE's
line list cannot be filtered by species through its interface, but it reads
whatever file `lines.list` names, so splitting gfall by species code into
`cache/species/` and pointing a private `lines.list` at each one isolates a
species exactly. Molecules were already one file per species in `data/mol/`.
Each run covers the whole band in one pass -- a single species is only
~10^4-10^5 lines, so it needs no chunking and takes 8 s.  Extending the band
therefore means extra chunks for panels 1-3 and a re-run of the species; the
per-species gfall files are split by species code over the whole line list, not
by wavelength, so nothing needs re-extracting. `assemble_species.py` collects
them into `cache/<tag>_species.npz`.

`makefigs.py` prefers `cache/species_flux_synthe.npz`, the 25-species
3550-10000 A set the committed PDFs were drawn from, and falls back to
`cache/<tag>_species.npz` from the current pipeline. Nothing in the repository
regenerates the former.

This matters most for **hydrogen**. gfall holds only 82 H I records in the band
(the Balmer series plus a Paschen pile-up near the limit) with no Stark or van
der Waals constants, but SYNTHE routes species code 1.00 to a complex-profile
branch and gives it the Holtsmark quasi-static Stark profile. A Voigt gets the
saturated cores about right (Halpha 0.344 against 0.327) and misses two thirds
of the wing area: H I supplies 0.62% of all optical line absorption with the
proper profile, 0.19% with a Voigt (0.19% vs 0.62% measured over 3800-10000 A;
over the full 3600-10000 A band, where the higher Balmer members crowd toward
the limit, H I rises to 1.18% and is second only to Fe I).

(`hilines.bin`, despite the name, is *high-ionization* lines -- Fe VI, Ni VI --
not H I. Hydrogen comes from gfall.)

### The independent python transfer

`specflux.py` does the same job with its own formal solution on the same model
atmosphere, and is kept as a cross-check rather than as the source of the
figures. It agrees with SYNTHE to 1-3% on every species except the molecules
(6-7%, from a line-strength floor) and hydrogen, where its Voigt profile is
3.5x too weak. Building it is what surfaced two traps worth remembering:

* **Molecular wavelengths.** SYNTHE ignores the wavelength column of the
  molecular line files and recomputes the vacuum wavelength from the level
  energies; the column itself is an air wavelength, 1.2 A off at 4300 A.
  Using the column put every CH line in the wrong place (G-band rms 0.29;
  it is 0.04 with the energies).
* **Populations.** SYNTHE's `.mol` diagnostic dump and the populations its
  synthesis actually uses disagree for the low-ionization metals -- K by 46x,
  Na 15x, Al 2.2x, Mg 1.7x at tau = 1, while every other element agrees to a
  few percent. Saha is what the synthesis uses: it reproduces SYNTHE's Na D
  profile to 0.005 in residual flux, the `.mol` densities are 5x off in the
  wing.

## Validation

The python transfer against SYNTHE on the same (gfall + molecules) line list:

| window | mean depression, SYNTHE | python | rms |
|---|---|---|---|
| 5161-5179 A (Mg b)   | 0.2244 | 0.2238 | 0.033 |
| 4296-4314 A (G band) | 0.5303 | 0.5063 | 0.040 |

Species by species, over the whole band, python / SYNTHE mean depression:
Fe I 0.99, Ca II 0.97, Mg I 1.03, Na I 1.05, Ca I 0.97, Ti I 0.99, Cr I 0.98,
Ni I 0.99, Fe II 0.99, CH 0.93, CN 0.94, MgH 1.16 -- and H I 0.29, the
Voigt-versus-Holtsmark gap.

Panels 1-3 against the Sun: over the figure band, 3600-10000 A, the synthesis
holds 65.6% of sigma Teff^4 and line blanketing removes 14.0% of the flux; over
the full 3549-25043 A synthesis those become 93.1% and 11.1%. The continuum at
5000 A is 1.069e7 erg/s/cm^2/A. Ratios to the IAG solar flux atlas are 0.97
(Mg b), 1.04 (Na D), 1.00 (Ca II IR triplet).

## Caveats

* Strict LTE with S = B. Ca II H and K cores go to zero (the Sun has ~5%
  residual plus an emission reversal) and Halpha is shallow (0.30 against 0.18
  observed) -- both form in the chromosphere.
* The H I panel cliffs at the Balmer limit: line absorption is exactly 1.0000
  up to 3644 A, then steps to 0.826 and sits on a merged pedestal. That is
  SYNTHE partitioning hydrogen opacity at the series limit with a hard
  boundary -- lines redward (high-n members merged by the Inglis-Teller
  tapering, `inglis = 1600/n_e^(2/15)`), bound-free continuum blueward, where
  the continuum jumps by 1.2103. Neither half is physical alone; the sum is
  continuous. The real artifact is that Stark wings of the near-limit lines
  should bleed a few A shortward of 3646 and are set to zero instead. The
  pedestal depth is not the cause of the blue over-blanketing: against the
  Kitt Peak flux atlas the synthesis runs ~0.07 low in mean residual flux
  equally on both sides of the limit (3560-3640 A: 0.44-0.51 vs 0.51-0.58;
  3650-3720 A: 0.48-0.50 vs 0.56-0.57).
* Panel (b), the full synthesis, carries Kurucz's predicted line list and the
  exotic molecules (TiO, H2O, VO, AlO, CaOH) that have no panel of their own.
  The species rows account for about 87% of it, so they do not sum to it -- nor
  would they exactly in any case, since blending is not additive.

## Interactive version

`web/` is a self-contained static site of the same synthesis, live at
https://cconroy20.github.io/spex/ : five stars, 51 species each over
3549-24993 A, drag-to-zoom, a resolution selector from R = 300,000 down to
1,000, telluric and blackbody overlays, formation depths, simulated photon
noise, line identification on hover, and a URL that tracks the view so a link
can point at one line. 30 MB, 828 files. See `web/README.md`. Rebuild its data
with `export_web.py`, `pack_lines.py` and `export_telluric.py`, and publish it
with `deploy.sh`.

## Files

    starcfg.py         the five stars and the 51 species, in one table
    split_gfall.py     gfall -> one file per species code, in cache/species/
    run_star.sh        every SYNTHE run one star needs: EOS, chunks, species
    assemble.py        stitch + validate -> cache/<tag>_synthe.npz
    assemble_species.py  per-species runs -> cache/<tag>_species.npz
    reduce_linform.py  .linform -> one formation depth per wavelength
    assemble_form.py   stitch those -> cache/<tag>_form.npz
    sunlib.py          model, continuum, line lists, SYNTHE's damping defaults
    moldisp.py         parses SYNTHE's molec_dispatch isotope table
    specflux.py        per-species formal solution, kept as a cross-check
    run_species.py     driver -> cache/species_flux.npz
    accumulate.py      optional: line OPACITY per species at tau_5000 = 1
    run_diag.sh        line-list ablation (gfall / +predicted / +molecules)
    makefigs.py        the figures
    plotstyle.py       house style

    binfmt.py          the SPC1 container the browser reads
    export_web.py      one star's spectra -> web/data/<tag>/
    export_telluric.py ESO SkyCalc transmission -> web/data/telluric*.bin
    build_lines.py     predicted central depths -> cache/<tag>_lineindex.npz
    pack_lines.py      shared catalogue + per-star depths -> web/data/lines/
    assess_lines.py    how much the catalogue grows across stars
    assess_compression.py  container variants, measured
    deploy.sh          publish web/ to gh-pages as a single commit
