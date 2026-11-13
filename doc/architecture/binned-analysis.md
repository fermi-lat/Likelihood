# Likelihood: binned analysis

Part of the [Likelihood architecture docs](README.md). This file covers
`BinnedLikelihood`, counts maps, source maps, PSF integration, energy
dispersion and weighted likelihood. Binned exposure maps and the mean
PSF's livetime weighting are in
[exposure-and-livetime.md](exposure-and-livetime.md).

## The statistic

`BinnedLikelihood : LogLike` (`Likelihood/BinnedLikelihood.h`, 806
lines; `src/BinnedLikelihood.cxx`) evaluates

    log L = Σ_{k in [kmin,kmax)} Σ_{filled pixels j} n_kj · log m_kj  −  Npred_total

`value()`:
1. `computeModelMap_internal(weighted=true)` fills `m_model` for the
   **filled pixels only** (pixels with data) and returns the total
   predicted counts. Fixed sources contribute through one precomputed
   array (`m_fixedModelCounts`) unless the fixed set has changed.
2. Loops over energy planes `k` and the filled pixels of that plane
   (`BinnedCountsCache::firstPixels()`/`filledPixels()`), adding
   `data · log(model)`. With a weights map, `data` is the **weighted**
   counts.
3. Subtracts Npred, adds priors, records the best fit.

`set_klims(kmin, kmax)` restricts the energy planes used.
`set_ebounds`/`unset_ebounds` throw for the binned case.

### Class map

| Class | Role |
|-------|------|
| `CountsMapBase : evtbin::DataProduct` | Counts cube read from a `gtbin` file or built in memory. `CountsMap` (WCS, `CCUBE`/`CMAP`) and `CountsMapHealpix` (all-sky or partial `HEALPIX`). Energies come from `EBOUNDS`. Also hosts the weights-map arithmetic (`makeBkgEffMap`, `makeAlphaMap`, `makeWtsMap`). |
| `BinnedCountsCache` | Data-side cache: counts, the list of filled pixels per energy plane, counts spectra (plain and weighted), log energy ratios, and the `WeightMap`. |
| `SourceMap` | One source's model cube on the counts-map grid **without its spectrum**: exposure × PSF-convolved spatial model, at the energy-bin **edges** (`num_ebins + 1` planes, plus 2·`edisp_bins` with energy dispersion). Stored dense or as a `SparseVector<float>`. Also holds the per-plane `npreds` (differential counts at the edges) and spectral weights. |
| `SourceMapCache` | Name → `SourceMap*`. Creates, loads from file, saves, and recomputes source maps. Applies energy-dispersion corrections (`Drm_Cache`). |
| `Drm`, `Drm_Cache` | Detector response matrix at the counts-map centre, and per-source cached true/measured counts spectra with the correction factors ξ. |
| `WeightMap` | The likelihood weights resampled onto the counts-map grid. |
| `BinnedLikeConfig`, `PsfIntegConfig` | All switches (see below). The constructor reads the `BinnedConfig.cxx` environment variables **every time** and lets them override its arguments. |
| `ModelMap` | Model counts map for `gtmodel`. |
| `CountsSpectra` | Observed and model counts spectra (`gtlike specfile`). |

Counts in an energy bin come from integrating `spectrum × npreds` between
the bin edges, with log–log trapezoids (`Source::pixelCounts`; linear with
`USE_LINEAR_QUADRATURE`). This is why source maps are evaluated at the
edges, and why the binned exposure map should also be at the edges
(`gtexpcube2 bincalc=EDGE`, the default).

## Source maps

### Computing a source map

`PSFUtils::makeModelMap` dispatches on source type. Which `MeanPsf` is
used:
- **Point sources** get a `MeanPsf` at their own position
  (`PSFUtils::build_psf`), unless `USE_SINGLE_PSF` is set.
- **Diffuse sources** use the `Observation`'s `MeanPsf`, computed at the
  counts-map centre (built by `AppHelpers`).

**Point sources** (`makePointSourceMap_wcs`/`_healpix`):
- For each pixel and energy, the PSF value is estimated by
  `psfValueEstimate` according to `PsfIntegConfig::integ_type`:
  `adaptive` (default; adaptive quadrature over the pixel with
  `psfEstimatorFtol`/`psfEstimatorPeakTh`), `annular`, or `pixel_center`.
- **PSF corrections** (`psfcorr`, default yes): the map is rescaled so
  that each energy plane integrates to the analytic PSF containment within
  the map radius. This is done only when the source is inside the map by
  at least 4 pixels (`withinBounds(..., 4)`).
- Multiplied by the point-source exposure from the binned exposure map
  (unless `apply_exposure="false"`).
- With `ptsrc=no` (`computePointSources=false`), point-source maps are
  not precomputed into the srcmaps file. They are computed on demand when
  first needed.

**Diffuse sources** (`makeDiffuseMap_wcs`/`_healpix`/`_native`):
- **WCS counts map:**
  - The template is projected onto a finer grid: the resampling factor
    comes from `rfactor`/`resample` and the template's own pixel size
    (`computeResampFactor`), and `minbinsz` limits rebinning of fine
    templates.
  - For maps made by `gtbin` (`CountsMap::conformingMap()`), the grid is
    padded by the 99.9% PSF radius (clamped to 1–10°) to capture PSF
    leakage from outside the ROI.
  - The padded grid is multiplied by exposure and **convolved with the
    PSF by FFT** (`Convolve::convolve2d`, FFTW), then summed back onto
    the counts-map pixels.
- **HEALPix counts map:** if the region needs more than 45°
  (`neededMapSize`), the template is convolved on the sphere with
  spherical harmonics (`ConvolveHealpix`, healpix_cxx `alm`). Otherwise it
  is convolved in a native WCS projection and resampled to HEALPix.
- **Analytic `SpatialFunction`s** (`RadialDisk`, `RadialGaussian`) use
  `spatialResponse` with a precomputed PSF-convolved profile instead of
  an FFT.
- `apply_psf="false"` skips convolution; `convol=no` skips it for every
  source.

**Composite sources** sum the maps of their components.

### Source maps file

`gtsrcmaps` (and `BinnedLikelihood::saveSourceMaps`) write:

| HDU | Content |
|-----|---------|
| Primary | The counts map (copied from `cmap`) with DSS keywords (including `IRF_VERSION`) and the WCS. |
| `EBOUNDS`, `GTI` | From the counts map. |
| `<source name>` (one per source) | WCS: 3-D image with `num_ebins+1` (+2·`edisp_bins`) planes. HEALPix: table, all-sky, partial, or sparse (`FileUtils::SrcMapType` = `WCS`, `WCS_Sparse`, `HPX_AllSky`, `HPX_Sparse`, `HPX_Partial`). `SourceMapCache::copyWcsToExtensions` copies the WCS keywords to each image. |
| `__weights__` | The weights map on the counts-map grid (weighted likelihood). |
| `__FIXED__`, `__FIXED_PAR__` | Optional single summed map of all fixed sources, and the parameters it was made with (`saveTotalFixedSourceMap`). |

**Reuse.** When `BinnedLikelihood` is created with `srcMapsFile` and
`load_existing_srcmaps` is true (the default), `SourceMapCache::getSourceMap`
reads a source's map from the HDU **with the same name** if one exists.
Otherwise it computes it. Nothing checks that the stored map matches the
current position, spatial model, IRFs, PSF settings or `edisp_bins`.
`AppHelpers::makeBinnedLikelihood` passes the **`cmap` file** as the
srcmaps file, so `gtlike cmap=srcmaps.fits` reuses everything `gtsrcmaps`
wrote, and `gtsrcmaps cmap=<old srcmaps>` would also reuse old maps for
matching names.

### Fixed sources

With `use_single_fixed_map` (the default in `BinnedLikeConfig`; `gtsrcmaps`
turns it off), sources whose spectra have no free parameters are summed
once into `m_fixedModelCounts` and their fixed Npreds. Their `SourceMap`
models are then cleared from memory (`clear_model`, unless `save_model`)
and reloaded from the file if the source is freed again
(`reloadIfCleared`). `fixedModelUpdated()` detects changes to the fixed
set and triggers `buildFixedModelWts`.

### Memory

Source maps are the dominant memory cost: sources × pixels × (energy
planes + 1) floats. The knobs:

| Knob | Where | Effect |
|------|-------|--------|
| `computePointSources` (`ptsrc`) | config / `gtsrcmaps` | Don't precompute point-source maps. |
| `use_single_fixed_map` | config | Sum fixed sources once and drop their maps. |
| `delete_local_fixed` | config, pyLikelihood `BinnedConfig` | Also free locally computed (not file-backed) fixed-source maps. |
| `no_cached_weightmaps` | config, pyLikelihood `BinnedConfig` (added 2025, LK-127 work) | Delete the original (unprojected) weights map after projection. |
| `save_all_srcmaps` / `SAVE_ALL_SRCMAPS` | config / env | Keep all models in memory (the opposite: faster, more memory). |
| Sparse maps | automatic for sparse HEALPix / WCS | `SparseVector<float>` storage. |

The regression suite's `test_memory*` scripts measure the effect of
these; `notes.md` there has the history.

## Energy dispersion

Controlled by one integer, `edisp_val` (`BinnedLikeConfig::edisp_val()`),
set from the `edisp_bins` tool parameter, the pyLikelihood config, or the
`USE_BL_EDISP` environment variable (which overrides everything):

| `edisp_val` | Meaning |
|-------------|---------|
| 0 | No energy dispersion. |
| N > 0 | **Per-source** dispersion. Source maps, the binned exposure map and the `MeanPsf` are computed on an energy grid extended by N log-spaced bins on **each side** (`FitUtils::expand_energies`). Predicted counts in measured-energy bin k are summed over true-energy bins with the DRM (`FitUtils::npred_edisp`). |
| N < 0 | **Old-style** correction with \|N\| DRM bins. The total true counts spectrum of each source is convolved with the DRM (`Drm::convolve`), and the ratios ξ_k = measured/true rescale the model. No extra source-map planes. |

- Sources with `apply_edisp="false"` (models that already include
  dispersion) get `edisp_val = 0` (`BinnedLikelihood::edisp_val(name)`).
- `Drm` is built at the counts-map reference direction from the livetime
  cube and IRF `IEdisp`. It turns phi dependence off while it integrates
  and then restores it.
- `gtexpcube2`, `gtsrcmaps` and `gtlike`/`gtmodel` must use the **same**
  `edisp_bins` and counts-map energy binning. Binned-exposure lookups need
  an exact energy-plane match (to 1e-5 fractional), so a mismatch throws
  "The energy … is not available" when the source maps are built (see
  [exposure-and-livetime.md](exposure-and-livetime.md#binned-exposure-map-gtexpcube2--binnedexposure)).
  `AppHelpers::checkExposureMap` checks only spatial coverage of the
  counts-map edge, at the first energy.
- `expand_energies` treats `edisp_bins ≤ 0` as a no-op; negative values
  from the tools are clamped to 0 for the energy expansion
  (`AppHelpers`, `gtexpcube2`).
- `Drm_Cache::update` clamps negative measured counts to 0 and prints a
  one-time message ("Measured counts < 0"; the code says "Not sure why
  this can happen"). `kref` guards were added for GCC 15 in September
  2026 (`0fb89cf2`).

Unbinned analysis has no energy dispersion; see
[unbinned-analysis.md](unbinned-analysis.md#gotchas).

## Weighted likelihood

The weighted likelihood down-weights pixels where systematic
uncertainties in the background model dominate. The theory is in
`doc/weighted_like/weighted_like.tex`.

Chain of tools:

```
gtmodel ─► model cube ─► gteffbkg (cmap=model, efact) ─► effective background B_eff
             [several components: gtalphabkg inputs=list ─► α]
                         gtwtsmap effbkgfile=B_eff [alphafile=α] epsilon=ε ─► wmap
             gtsrcmaps wmap=… / gtlike wmap=… / pyLikelihood BinnedAnalysis(wmap=…)
```

- `gteffbkg`: `CountsMap::makeBkgEffMap(MeanPsf, efact)` convolves the
  model cube with the PSF and normalizes each plane by the PSF peak value
  at the bin's geometric-centre energy, giving the background counts
  inside an effective PSF solid angle. Each plane then accumulates the
  planes **above** it whose mean energy is within a factor `efact` (so
  B_eff at E counts background from E to about `efact`·E).
- `gtalphabkg`: for each pixel, over the input B_eff maps (one per data
  component) with minimum b_min, α = (1 + ε² b_min) / (1 + ε² b_min Σ_i
  (b_min/B_i)²) (`getAlphaWts`). Where b_min ≤ 0 the sum is set to 0.
- `gtwtsmap`: w = α / (1 + ε² B_eff) per pixel (`CountsMapBase::getWts`).
  Written with `LatWeightsMapTemplate`/`LatHealpixWeightsTemplate`.
- Using it: `BinnedLikelihood` gets the weights as a `ProjMap`. It is read
  through `WcsMapLibrary` with **interpolation off and extrapolation on**,
  resampled onto the counts grid as `WeightMap`, and the source maps and
  data are weighted. Weighted source maps are written to `__weights__`.
  An existing `__weights__` HDU is reused unless `overwriteWeights`.
- Npred: `value()` uses the weighted Npred; `NpredValue(name, weighted)`
  can return either.

## Gotchas

- **Stale source maps:** see *Reuse* above. This is the most common
  cause of "my change had no effect".
- **Pixels with model ≤ 0 are skipped** in the data term (`value()`
  `continue`s), so data in such pixels does not penalize the model.
- **Only filled pixels** are evaluated in the data term. Pixels with zero
  counts contribute only through Npred, which is correct for Poisson
  statistics. Changing the counts with `setCountsMap` triggers a rebuild
  of the filled-pixel list and the fixed-model weights.
- **HEALPix counts maps require an all-sky binned exposure map**
  (`checkExposureMap_healpix` throws otherwise).
- **Non-`gtbin` WCS counts maps** (`conformingMap()` false) get no padding
  for PSF leakage when diffuse maps are convolved.
- **`emapbnds`** (`gtsrcmaps`): with the default `yes`, a diffuse source
  that needs exposure outside the binned exposure map throws ("Request
  for exposure at a sky position that is outside of the map boundaries").
  `gtsrcmaps` explains the choice: use an all-sky exposure map, or
  `emapbnds=no` to drop the contribution from outside.
- **Progress dots divide by `npts/20`** (`PSFUtils`) and `n/20`
  (`CountsMapBase::getWts`, `getAlphaVector`). For maps with fewer than
  20 values (×energies) this is a modulo by zero; `test_Likelihood`
  builds with `TRAP_FPE` (from reading the code; not checked at runtime).
- **Two deprecated constructors** of `BinnedLikelihood` print a warning to
  `stderr` and hard-code `edisp_val = 0` and the PSF integration settings.
- **Environment variables override configuration.** A `BinnedLikeConfig`
  built with `edisp_val=0` will still have dispersion on if `USE_BL_EDISP`
  is set.
