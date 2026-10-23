# Likelihood: exposure and livetime

Part of the [Likelihood architecture docs](README.md). This file covers
the spacecraft data, ROI, livetime cubes, exposure maps, the mean PSF, and
how IRFs are loaded and used.

## Overview

```
FT2 ──► ScData (unbinned: per-event pointing and livetime fraction)
FT2 ──► gtltcube: LikeExposure (map_tools::Exposure) ──► ltcube.fits
                                         │ read by
                                         ▼
                                    ExposureCube ──► PointSource exposure (unbinned)
                                         │      ──► ExposureMap   (gtexpmap: unbinned diffuse)
                                         │      ──► BinnedExposure / BinnedHealpixExposure (gtexpcube2)
                                         │      ──► MeanPsf (binned PSF), Drm (energy dispersion)
ResponseFunctions (irfs) ────────────────┘
```

The livetime cube stores, for every sky direction (HEALPix, 1° by
default), the livetime accumulated in each cos θ bin (θ = off-axis angle),
and optionally in φ bins. Exposure for any direction and energy is then
Σ_bins livetime(cos θ[, φ]) × A_eff(E, θ[, φ]), with no further pass over
FT2.

## ScData

`ScData` (`src/ScData.cxx`) reads FT2 `START`, `STOP`, `LIVETIME`,
`RA_SCZ`, `DEC_SCZ`, `RA_SCX`, `DEC_SCX` with a cfitsio row filter for the
requested time range **±60 s**. It does not use `astro`'s pointing
history.
- `zAxis(t)`/`xAxis(t)` interpolate the spacecraft axes linearly between
  interval starts.
- `livetimefrac(t)` = `LIVETIME/(STOP−START)` of the interval containing
  t, or 0 if t is in a gap.
- Used for each event in unbinned analysis (`Event` stores the axes and
  the efficiency factor), and for point-source exposure when there is no
  livetime cube.

Gotchas:
- The ±60 s pad assumes FT2 intervals of at most 30 s (`maxIntervalSize`).
  Standard 30 s and 1 s FT2 files are fine.
- Multiple FT2 files must be given **in time order**. Otherwise it throws
  "start times … are not monotonically increasing". The check starts at
  the third row read (`m_start.size() > 1`).
- `zAxis`/`xAxis` need at least two intervals; with one row,
  `size() - 2` underflows (from reading the code).
- Times outside the loaded range throw from `time_index`.

## RoiCuts

`RoiCuts` holds the ROI cone, energy range, time ranges and GTIs. It
reads them from the **DSS keywords** of FT1 or counts-map files via
`dataSubselector::Cuts` (see
`src/dataSubselector/doc/architecture/dss-keywords.md`). `energies()` is a
100-point log grid over the energy cut (`makeEnergyVector`), used for
unbinned exposure and Npred. `writeGtiExtension`/`writeDssTimeKeywords`
put the time selection into `gtltcube` outputs.

## Livetime cubes

### Writing: `gtltcube` / `LikeExposure`

`LikeExposure : map_tools::Exposure` (see
`src/map_tools/doc/architecture/` and `src/healpix/doc/architecture/` for
the `HealpixArray<CosineBinner>` storage and its **process-wide static
binning**).

`gtltcube` (`src/makeExposureCube/makeExposureCube.cxx`):
1. Takes the GTIs from the FT1 file's DSS keywords (`RoiCuts::readCuts`).
   With `evfile=""` it uses `tmin`/`tmax` instead.
2. Reads FT2 rows in the GTI range ±60 s.
3. `LikeExposure::load` walks the rows. For each interval it computes
   the fraction inside the GTIs (`acceptInterval`) and adds
   `LIVETIME × fraction` to every sky pixel whose zenith angle is within
   [`zmin`, `zmax`] (`map_tools::Exposure::fill`). It also fills a second
   cube weighted by the livetime fraction (`LIVETIME/(STOP−START)`), which
   becomes `WEIGHTED_EXPOSURE`.
4. With `phibins > 0`, it calls `CosineBinner::setPhiBins` and
   `fill_zenith`, which also bins in instrument φ using `RA_SCX`/`DEC_SCX`.
5. Writes `EXPOSURE`, `WEIGHTED_EXPOSURE`, `CTHETABOUNDS` and `GTI`, from
   `data/LivetimeCubeTemplate`. Header keys: `THETABIN`, `NBRBINS`,
   `PHIBINS`, `TSTART`/`TSTOP`, and the zenith cut.

`dcostheta` (0.025) and `binsz` (1°) set the binning. The `file_version`
parameter is not read.

### Summing: `gtltsum`

`gtltsum` (`src/gtaddlivetime/`) adds the `COSBINS` of two cubes or a list
of cubes, row by row, for both `EXPOSURE` and `WEIGHTED_EXPOSURE`. It
merges GTIs and takes the earliest `TSTART` and latest `TSTOP`. It
requires identical `CTHETABOUNDS`, row counts, `COSBINS` lengths and
`ZENMAX`. It does **not** check that the inputs cover disjoint times, so
summing overlapping cubes double-counts livetime silently (from reading
the code).

### Reading: `ExposureCube`

`ExposureCube::readExposureCube(file)` loads `EXPOSURE` and (if present)
`WEIGHTED_EXPOSURE` into `map_tools::Exposure` objects. It detects φ
binning from the `PHIBINS` keyword.

`value(dir, aeff, energy)` is the exposure used everywhere:

    exposure = f1(E) · Σ L(cos θ) A(cos θ)  +  f2(E) · Σ L_w(cos θ) A(cos θ)

where `f1`, `f2` come from the IRF's `IEfficiencyFactor::getLivetimeFactors`
(the livetime-fraction correction), evaluated at the **midpoint MET**
`(tstart + tstop)/2` of the cube. If the IRFs have no efficiency factor,
f1 = 1 and f2 = 0. A negative result throws.

### Phi dependence

- If the cube has φ bins, `value` integrates A_eff(θ, φ) over them
  (`AeffWrapper` converts φ to degrees for the IRF).
- If it does not, **constructing an `ExposureCube::Aeff` turns phi
  dependence off on every loaded IRF** (`setPhiDependence(false)`), and it
  stays off for the rest of the process.
- `gtexpcube2 ignorephi=yes` does the same deliberately.
- `PointSource::Aeff` saves the IRF's phi state, sets its own, and
  restores it for each call.
- `Drm` turns phi dependence off while it computes and restores it
  afterwards.

## Point-source exposure

`PointSource::computeExposureWithHyperCube` evaluates, at each of the
source's energies, `ExposureCube::value(srcDir, PointSource::Aeff, E)`.
`PointSource::Aeff` returns A_eff × (PSF fraction inside the ROI
acceptance cone), summed over the loaded event types. So for unbinned
analysis, "exposure" includes the PSF containment of the ROI.

Without a livetime cube, `PointSource::computeExposure` sums over `ScData`
intervals instead, which is much slower.

Gotcha (from reading the code): in `PointSource::Aeff::value`, if
**any** event type's A_eff is below 0.1 cm² at that angle, the function
returns 0 for the whole sum. The comment calls this a kluge because the
PSF is ill-defined there. With PSF0–3 or FRONT/BACK types this can drop
other types' contributions at large off-axis angles.

## Unbinned exposure map: `gtexpmap` / `ExposureMap`

`ExposureMap::computeMap` builds a CAR grid of `nlong × nlat` pixels
centred on the ROI, spanning `srcrad`. `srcrad` must be at least the ROI
radius (otherwise it throws), and a message suggests ROI radius + 10° so
that sources whose PSF leaks in are covered. The map is computed at `nenergies` energies over the ROI
range. Each pixel's value is a point-source exposure at that direction
(livetime cube or `ScData`), so it includes the ROI's PSF containment.
`submap`/`nlongmin…` compute a piece of the map (for parallel runs by
hand).

Output: primary 3-D image plus an `ENERGIES` table. DSS keywords come
from the FT1 ROI cuts plus `IRF_VERSION`. The GTI extension is taken from
the **livetime cube** when one is given (requested by the catalog group,
who use hand-made cubes). `DiffuseSource` uses it through
`integrateSpatialDist` for Npred.

`gtexpmap` prints two warnings on every run: maps are for **unbinned**
analysis only, and half-degree pixels are a nominal choice that must be
refined for small extended sources.

## Binned exposure map: `gtexpcube2` / `BinnedExposure*`

`BinnedExposureBase` with two implementations:
- `BinnedExposure` (WCS; primary image + `ENERGIES` table).
- `BinnedHealpixExposure` (`HPXEXPOSURES` table + `ENERGIES`).

For every pixel centre and energy, it sums over event types
`ExposureCube::value(dir, BinnedExposureBase::Aeff, E)`. This `Aeff` is
A_eff only (no PSF), restricted to `thmin ≤ θ ≤ thmax`. The result is in
cm² s.

Geometry and energies:
- Taken from the counts map, or from `nxpix`/`nypix`/`binsz`/`xref`/
  `yref`/`proj`/`coordsys` and `emin`/`emax`/`enumbins`/`ebinfile` when
  given. By default an all-sky map is recommended so that source maps can
  include sources outside the ROI.
- `bincalc=EDGE` (default) puts planes at the energy-bin edges (what
  `SourceMap` needs); `CENTER` puts them at log centres.
- With `edisp_bins = N`, N extra planes are added on each side.

Lookup (`operator()`, `get_exposures_for_dir`):
- The energy must match a plane to 1e-5 fractional
  (`BinnedExposureBase::findNearest`), otherwise it throws "The energy …
  is not available". So the source maps' energies (counts-map edges,
  expanded by `edisp_bins`) must be exactly the exposure map's planes:
  same counts-map binning, `bincalc=EDGE`, same `edisp_bins`.
- The sky position is **not interpolated**: the pixel index is
  `static_cast<int>(pixel − 1)` of the 1-based fractional pixel
  coordinate. That is a truncation, while the map was computed at
  integer pixel centres, so lookups can come from the neighbouring pixel
  up to one pixel away (from reading the code; small for smooth exposure
  and fine pixels).
- Outside the map: 0, or a throw when `setBoundaryFlag(true)`
  (`gtsrcmaps emapbnds=yes`, and `checkExposureMap`).

## MeanPsf

`MeanPsf(ra, dec, energies, observation)` is the exposure-weighted PSF at
one direction:

    PSF_mean(δ; E) = Σ_types ∫ L(cos θ) A(E,θ) PSF(δ; E,θ) / exposure(E)

tabulated on a fixed δ grid. It provides `value(δ, E)`, `integral`,
`containmentRadius`, `peakValue` and `exposure(E)`. `gtpsf` writes it out
(`PSF` table with `Energy`, `Exposure` and the PSF values; `THETA` table).
The binned analysis uses one at the counts-map centre (built by
`AppHelpers`) for diffuse sources, and one per point source. See
[binned-analysis.md](binned-analysis.md#computing-a-source-map).

## ResponseFunctions and IRF selection

`ResponseFunctions` maps event-type index → `irfInterface::Irfs*`. The
irfs package docs describe the factory and naming:
`src/irfs/doc/architecture/loading-and-names.md`.

`AppHelpers::createResponseFuncs` decides which IRFs to load:
- **`irfs=CALDB`:** read DSS keywords (`dataSubselector::Cuts`) from the
  exposure map or FT1 (unbinned), or from `bexpmap`, falling back to
  `cmap`/`srcmaps` (binned). If there is no `EVENT_TYPE` cut, add one from
  the `evtype` parameter. Then `Cuts::CALDB_implied_irfs()` gives the name
  (e.g. `P8R3_SOURCE_V3`) and the selected event types.
- **Explicit name:** `::` qualifiers are rejected. The event class comes
  from the name (`Cuts::setIrfs`) and the event types from `evtype`
  (`evtype=3` prints "Using evtype=3 (i.e., FRONT/BACK irfs)").
- `ResponseFunctions::load(name, "", evtTypes)` creates one `Irfs` per
  event type in the selection. It throws listing all valid names if the
  name is unknown.
- `ResponseFunctions::efficiencyFactor()` returns the **first** event
  type's efficiency factor, and that one is used for the whole livetime
  cube.

The `irfs=DC2` special case in `AppHelpers::responseFuncs` (pick a DC2
class from `EVENT_CLASS` cuts) is a leftover from 2006-era data
challenges.

## Gotchas

- **Phi dependence is global, mutable IRF state** (see above). In a
  pyLikelihood session, loading a φ-less livetime cube switches φ
  dependence off for all later exposure calculations with those IRFs.
- **Livetime correction uses the cube's mid-time.** For long data sets
  the efficiency factors vary with time, but only `(tstart+tstop)/2` is
  used.
- **`LikeExposure` keeps references** to the `timeCuts` and `gtis`
  vectors passed to its constructor. `gtltcube` passes locals of
  `createDataCube`. They are used only in `load()`, which runs inside that
  function, so this is safe today, but calling `load()` later would read
  freed memory.
- **`gtltcube` with no matching FT2 rows** warns ("All livetimes will be
  identically zero") and writes an all-zero cube with exit status 0.
- **Exposure-map consistency checks are partial.** `checkTimeCuts`
  between FT1/counts map and livetime cube only warns (the
  "Aggregate absolute differences in GTI start times" message).
  `checkExposureMap` checks spatial coverage of the counts-map edge at
  one energy.
- **`gtexpcube2 irfs=CALDB` fails** (confirmed at runtime); see
  [tools.md](tools.md#gtexpcube2).
