# Likelihood: unbinned analysis

Part of the [Likelihood architecture docs](README.md). This file covers
`LogLike`, events, diffuse responses and `Npred`. For exposure (the
`Npred` input) see [exposure-and-livetime.md](exposure-and-livetime.md).

## The statistic

`LogLike` (`Likelihood/LogLike.h`, `src/LogLike.cxx`) evaluates the
extended unbinned log-likelihood

    log L = Σ_events log( Σ_sources ε_j · M_i(event_j) )  −  Σ_sources Npred_i

where `M_i(event)` is `Source::fluxDensity(event)` (spectrum × response
at the event's measured energy and direction, for the event's type) and
ε_j is the event's efficiency factor (livetime-fraction correction from
the IRF's `IEfficiencyFactor`, evaluated at event time).

`value()`:
1. If `set_ebounds` was called, events outside [emin, emax] are skipped.
   If the new range excludes all the ROI data, it returns 0 (requested by
   J. Ballet; `UnbinnedAnalysis.setEnergyRange` relies on it).
2. Data term: `logSourceModel(event)` for each event, using the
   per-event `ResponseCache`.
3. Model term: subtracts the cached `m_npredValues` (default) or, with
   `USE_OLD_LOGLIKE`, recomputes `Npred` for every source on every call.
4. Terms go through a Kahan accumulator; priors
   (`Parameter::log_prior_value`, added with `addPrior`) are added when
   `include_priors` is true.
5. `saveBestFit` remembers the best parameters seen so far.
   `restoreBestFit()` puts them back. A NaN total throws.

Derivatives (`getFreeDerivs`) follow the same structure with
`fluxDensityDeriv` and `NpredDeriv`.

## Events

| Class | Role |
|-------|------|
| `EventContainer` | Owns `std::vector<Event>`. `getEvents(file)` reads FT1 through `tip`: `RA`, `DEC`, `ENERGY`, `TIME`, `ZENITH_ANGLE`, `CONVERSION_TYPE`, `EVENT_CLASS` (bit array for Pass 8, integer for P7/none), `EVENT_TYPE`. Rejects events outside `RoiCuts` (`accept`). Reads any diffuse-response columns. |
| `Event` | One photon: direction, energy, time, spacecraft z/x axes at that time (interpolated from `ScData`), cos(zenith), event type, efficiency. Holds per-diffuse-source responses (`m_respDiffuseSrcs`) and the running `modelSum`. |
| `ResponseCache` / `CachedResponse` (`EventSourceCache.h`) | Per (event, source) cached responses, so that only the spectrum is re-evaluated when a source's spectral parameters change. |

Event type: `EVENT_TYPE` is masked with the `EVENT_TYPE` bit-mask cut from
the file's DSS keywords (default mask 3 = FRONT/BACK) and converted to an
index with `log2`. Without an `EVENT_TYPE` column, `CONVERSION_TYPE` is
used.

## Diffuse responses

For a diffuse source the response at an event is the PSF-convolved,
exposure-free integral of the spatial model around the event direction.
It does not depend on the spectral parameters, so it is computed once per
(event, source).

- **Precomputed (`gtdiffrsp`):** adds one `1E` column per diffuse source
  to the FT1 file **in place**. Column names are generic (`DIFRSP0`,
  `DIFRSP1`, …). Header keywords `NDIFRSP` and `DIFRSPn` map them to the
  real names, which have the form `<irfs>__<source name>`, lowercased
  (`Event::diffuseSrcName`, `DiffRespNames`). Older files use the real
  name as the column name; `gtdiffrsp convert=yes` rewrites the header.
- **On the fly:** `LogLike::computeEventResponses()` (called by `gtlike`
  after reading the model) computes the missing responses with Gaussian
  quadrature (`Event::computeResponseGQ`). Integration extends to
  4× the 2σ (95.4%) PSF containment radius, or further if needed to reach
  72.5° off-axis.
  `MAP_BASED_DIFFRSP` switches to map-pixel integration in some cases.
- Which sources are "new" is decided from the **first event only**
  (`getNewDiffuseSrcs` on `m_events[0]`).

A precomputed column is used only when its key matches the current IRF
name (with any `::` event-type qualifier stripped) and the source name.
Otherwise the response is recomputed, which is correct but slow.

`gtsrcprob` uses the same machinery: for each event and each source it
writes `fluxDensity_i / Σ fluxDensity` as a new column in a copy of the
FT1 file.

## Npred and exposure

`Npred` (`src/Npred.cxx`) is an `optimizers::Function` wrapper around
`Source::Npred()` / `NpredDeriv()`, so that Npred can be treated like any
other function of the parameters.

- **Point sources:** `Source::Npred` integrates spectrum × `m_exposure`
  over `m_energies` with log–log trapezoids. `m_energies` is
  `RoiCuts::energies()`, 100 log-spaced points spanning the ROI energy
  cut. `m_exposure` is computed once per source, when it is created or
  moved; see [exposure-and-livetime.md](exposure-and-livetime.md#point-source-exposure).
- **Diffuse sources:** exposure comes from the **unbinned exposure map**
  (`gtexpmap`, `ExposureMap::integrateSpatialDist`), or from map-pixel
  quadrature (`map_based_integral="true"` / `MAP_BASED_NPRED`).

`set_ebounds(emin, emax)` recomputes every source's exposure on a new
100-point grid from emin to emax. `unset_ebounds()` only clears the flag
and recomputes Npred: the exposure arrays stay on the narrowed grid
(from reading the code; pyLikelihood never calls `unset_ebounds`, it
calls `set_ebounds` with the full range instead).

## Inputs a run needs

| Input | Why | Checks done by `gtlike` |
|-------|-----|------------------------|
| FT1 (one file or `@list`) | Events; DSS keywords define the ROI (`RoiCuts`), energy range, GTIs, IRF selection. | `AppHelpers::checkCuts` between files (event-class cuts skipped unless `irfs=DSS`). |
| FT2 | `ScData`: pointing and livetime fraction for each event. | — |
| Livetime cube | Point-source exposure. Without it, exposure is summed over `ScData` intervals (slow). | `checkTimeCuts` vs FT1 GTIs (warning only). |
| Unbinned exposure map | Diffuse-source Npred. | `checkExpMapCuts`: DSS keywords must match exactly except energy; the map's energy range must cover the data's. |
| XML model | Sources. | — |

## Gotchas

- **Energy dispersion is not available in unbinned analysis.** Any
  attempt throws ("Attempt to use energy dispersion handling in unbinned
  analysis", or "energy dispersion handling disabled" in
  `PointSource::Aeff`). `gtlike`'s `edisp` parameter is never read.
- **Non-positive model density** at an event returns −1e30 from
  `logSourceModel` (or 0 if there is a single source) instead of failing.
  `LOGLIKE_CATCH_NEG_PROB` turns this into an exception. Fits that wander
  into negative-normalization regions therefore see a huge penalty rather
  than an error.
- **Diffuse-response names are lowercased**, so two diffuse sources whose
  names differ only in case share one response.
- **`gtdiffrsp` modifies its input FT1 file.** Work on a copy if the
  original must stay pristine.
- **`EVENT_TYPE` → index via `log2(type & mask)`** assumes exactly one bit
  of the mask is set per event. If `type & mask` is 0, the result is
  `log(0)` cast to an integer (undefined); such events should already have
  been removed by `gtselect` (from reading the code).
- **`ScData` window:** FT2 rows are read for the ROI time range ±60 s.
  The code assumes FT2 intervals are at most 30 s and that FT2 files are
  given in time order; out-of-order files throw. See
  [exposure-and-livetime.md](exposure-and-livetime.md#scdata).
- **Response cache memory** grows as events × sources with free spectra.
  The regression suite's `test_memory*` scripts track it.
- `likelihood::readEventData` has a `BINNED` branch that reads events
  into the `EventContainer` (commented "Why?"). It is a no-op, because
  the event-file list is only filled in `UNBINNED` mode.
