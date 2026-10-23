# Likelihood: tools

Part of the [Likelihood architecture docs](README.md). This file lists
every gt tool this package builds, the shared setup in `AppHelpers`, and
exit-status behaviour. User-facing help text is in `fermitools-fhelp`.
Parameter-file mechanics are in `src/ape/doc/architecture/parameter-files.md`,
and the tool framework (startup, standard hidden parameters `chatter`,
`clobber`, `debug`, `gui`, `mode`) is in `src/st_app/doc/architecture/`.

All 20 tools are `st_app::StApp` subclasses in `src/<dir>/`, each linked
`PRIVATE` to `Likelihood` (plus `st_graph` for `gtlike`, `st_stream` for
`gtexpmap`). Par files are in `pfiles/`.

## Tool table

| Tool | Source | Purpose | Main inputs → outputs | Key parameters | Drives |
|------|--------|---------|-----------------------|----------------|--------|
| `gtltcube` | `makeExposureCube/` | Livetime cube | FT1 (for GTIs) + FT2 → ltcube | `dcostheta`, `binsz`, `phibins`, `zmin`/`zmax`, `tmin`/`tmax` (only if `evfile=""`) | `LikeExposure`, `RoiCuts`. Does **not** use `AppHelpers`. |
| `gtltsum` | `gtaddlivetime/` | Sum livetime cubes | 2 cubes or `@list` → ltcube | `table` (`Exposure`), `table2` (`WEIGHTED_EXPOSURE`) | `tip`, `dataSubselector::Cuts::mergeGtis` |
| `gtexpmap` | `expMap/` | Unbinned exposure map | FT1, FT2, ltcube → CAR cube | `srcrad`, `nlong`, `nlat`, `nenergies`, `submap` | `ExposureMap`, `PointSource` exposure |
| `gtexpcube2` | `gtexpcube2/` | Binned exposure map | ltcube + counts map (or explicit geometry) → WCS or HEALPix cube | `evtype` (**must be given**), `irfs` (see below), `edisp_bins`, `bincalc`, `ignorephi`, `thmin`/`thmax`, geometry and energy pars, `hpx_*` | `BinnedExposure`, `BinnedHealpixExposure` |
| `gtdiffrsp` | `diffuseResponses/` | Diffuse responses into FT1 (**in place**) | FT1, FT2, model → same FT1 + `DIFRSPn` columns | `evclsmin`/`evclass` (legacy), `convert` | `Event::computeResponseGQ`, `DiffRespNames`, `xmlBase` parser |
| `gtsrcmaps` | `gtsrcmaps/` | Source maps | counts map, model, ltcube, bexpmap [, wmap] → srcmaps | `ptsrc`, `psfcorr`, `convol`, `resample`/`rfactor`, `minbinsz`, `emapbnds`, `edisp_bins`, `copyall` | `AppHelpers::makeBinnedLikelihood`, `BinnedLikelihood::saveSourceMaps` |
| `gtlike` | `likelihood/` | Fit | UNBINNED: FT1, FT2, ltcube, expmap; BINNED: srcmaps (`cmap`), bexpmap, ltcube [, wmap] → `results.dat`, `counts_spectra.fits`, fitted XML | `statistic`, `optimizer`, `ftol`/`toltype`, `tsmin`, `refit`, `sfile`, `plot`, `edisp_bins` (binned) | `LogLike`/`BinnedLikelihood`, `optimizers::OptimizerFactory`, `CountsSpectra` |
| `gtmodel` | `gtmodelmap/` | Model counts map | srcmaps (`srcmaps`), model, ltcube, bexpmap → CMAP or CCUBE | `outtype`, `edisp_bins`, `psfcorr` | `BinnedLikelihood`, `ModelMap` |
| `gttsmap` | `TsMap/` | TS map | as `gtlike`, plus a grid → TS image | `statistic`, grid pars, `optimizer` | `LogLike`/`BinnedLikelihood`; adds a `PowerLaw` test source at each pixel and refits |
| `gttscube` | `TsCube/TsCube.cxx` | TS cube (WCS) with per-energy likelihood scans | srcmaps, bexpmap, ltcube, model → TS cube | `nnorm`, `nsigma`, `covscale_bb`, `covscale`, `stlevel`, `maxiter`, `remakesrc`, `target` | `FitScanner`, Newton fitter in `FitUtils` |
| `gthealcube` | `TsCube/HealCube.cxx` | Same on a HEALPix grid | as `gttscube` | `hpx_order`, `hpx_ordering_scheme`, `hpx_region` | `FitScanner` |
| `gtfindsrc` | `gtfindsrc/` | Point-source position | FT1, FT2, ltcube, expmap, model → trial points (text) | `target`, `ra`/`dec` or `l`/`b`, `atol`, `posacc`, `reopt` | `LogLike`, `optimizers::Amoeba` (Nelder–Mead); unbinned only |
| `gtsrcprob` | `gtsrcprob/` | Per-event source probabilities | FT1 (with diffuse responses), FT2, model → copy of FT1 with one column per source | `srclist` | `SourceModel`, `Source::fluxDensity` |
| `gtpsf` | `meanPsf/` | Mean PSF at a position | ltcube → `PSF` and `THETA` tables | `ra`, `dec`, `emin`/`emax`/`nenergies`, `thetamax`/`ntheta`; `irfs` default **`P8R2_SOURCE_V6`**, `evtype` default 3 | `MeanPsf` |
| `gtdrm` | `gtdrm/` | Detector response matrix | srcmaps/cmap, ltcube, bexpmap → `DRM` table | `edisp_bins` (the other binned pars are marked "Not used") | `BinnedLikelihood::drm()`, `FileUtils::write_drm_to_table` |
| `gtbkg` | `backfile/` | Background PHA for XSPEC | PHA1 (from `gtbin`), model, FT2, ltcube [, expmap] → PHA copy with model counts | `target` (excluded source) | `LogLike` Npred per channel |
| `gteffbkg` | `gteffbkg/` | Effective background | model or counts cube, ltcube, bexpmap → B_eff cube | `efact`, `edisp_bins` | `CountsMap::makeBkgEffMap`, `MeanPsf` |
| `gtalphabkg` | `gtalphabkg/` | α map for several components | text list of B_eff maps → α cube | `epsilon` | `CountsMapBase::makeAlphaMap` |
| `gtwtsmap` | `gtwtsmap/` | Likelihood weights | B_eff [, α] → weights map | `epsilon`, `alphafile` | `CountsMapBase::makeWtsMap` |
| `gtebl` | `gtebl/` | EBL optical depth table | — → ASCII (E, τ) | `model_id` (0–19, `IRB::EblModel`; not range-checked), `redshift` | `IRB::EblAtten` |

Not built: `src/exposureMap/` (old exposure tool) and `src/gtcntsmap/`
(counts maps are now made by `gtbin`).

Recent tool changes:
- `gtexpcube2`'s `evtype` became a required prompt (December 2025).
- `CREATOR` keywords with the tool version were added to `gtsrcmaps`,
  `gtmodel` and `gtalphabkg` outputs (October 2026).
- `gtmodel` now updates existing header keywords instead of duplicating
  them (`lk-123`, PR #132).
- `gtlike plot=yes` draws through embedded matplotlib
  (`matplotlibcpp.h`) since 2024; a Linux segfault was fixed in October
  2026.
- `gtexpmap` gained the pixel-size warning (March 2026).
- GTI and DATE keyword ordering in `gtexpmap` outputs was fixed (July
  2026, PR #145).

## AppHelpers

`AppHelpers` (`Likelihood/AppHelpers.h`, `src/AppHelpers.cxx`) is the
shared setup for most tools. `gtltcube`, `gtltsum` and `gtebl` don't
construct one; `gtwtsmap` and `gtalphabkg` use only its static
`readCountsMap`.

`AppHelpers(pars, analysisType)` with `"UNBINNED"`, `"BINNED"` or
`"NONE"`:
1. Builds a `FunctionFactory` with all prototypes
   (`addFunctionPrototypes`; see
   [models-and-sources.md](models-and-sources.md#function-registry)).
2. `createResponseFuncs(analysisType)`: works out the IRF name and event
   types (see
   [exposure-and-livetime.md](exposure-and-livetime.md#responsefunctions-and-irf-selection)).
   With `irfs=CALDB` the DSS keywords come from:
   - `UNBINNED`: `expmap`, else `evfile`.
   - `BINNED`: `bexpmap`, else `cmap`, falling back to `cmap`/`srcmaps`
     if the file has no `EVENT_CLASS` cut.
   - `NONE`: no file, which is the root of the `gtexpcube2` bug below.
3. Creates `RoiCuts`, `ScData`, `ExposureCube` (reads `expcube` if the
   parameter exists and is not `none`), `ExposureMap`, `EventContainer`.
4. `BINNED`: reads `bexpmap` (WCS or HEALPix, via `readBinnedExposure`)
   and `phased_expmap` (via `WcsMapLibrary`). Builds the `Observation`.
   If an `expcube` parameter exists, builds a `MeanPsf` at the counts-map
   reference direction, on the counts-map energies expanded by
   `edisp_bins`.
5. Missing parameters are tolerated: each `pars["x"]` is wrapped in
   `catch (hoops::Hexception&)`, which is how one constructor serves 15
   tools with different par files.

Helpers used by tools and pyLikelihood:
- `setRoi(file, ext, strict)`, `readScData()`, `readExposureMap()`.
- `checkCuts` / `checkTimeCuts` / `checkExpMapCuts` / `checkExposureMap`:
  consistency checks between files' DSS keywords and geometry. GTI
  differences only **warn** by default (`gtiWarningOnly=true`).
- `readCountsMap(file)` and `readBinnedExposure(file)`: pick the WCS or
  HEALPix class from the file.
- `makeBinnedLikelihood(pars, helper, cmap_par)`: requires a livetime
  cube, reads the counts map, checks exposure-map coverage, builds a
  `BinnedLikeConfig` from `ptsrc`, `psfcorr`, `convol`, `resample`,
  `rfactor`, `minbinsz`, `edisp_bins` (each with a default if the
  parameter is absent), reads `wmap`, and returns a `BinnedLikelihood`
  whose srcmaps file is the `cmap` file.
- `param<T>(pars, name, default, suppress_warn)`: read a parameter or
  return a default.
- `checkOutputFile(clobber, file)`: if the file exists and `clobber=no`,
  prints a message and calls `std::exit(1)`.

`AppHelpers`'s destructor deletes everything it created, including the
`phased_expmap` pointer it got from `WcsMapLibrary`. The library still
holds that pointer, so a later `wcsmap()` call for the same file would
return freed memory (from reading the code; tools exit right after).

## Exit status

The `st_app` framework turns an uncaught exception into exit status 1
(`src/st_app/doc/architecture/README.md`). Deviations in this package:

| Tool | Situation | Exit status |
|------|-----------|-------------|
| `gtlike` | The optimizer throws during `find_min` or `getUncertainty`. The message is printed, then results, XML and spectra are written as usual. | **0** (from reading the code) |
| `gtsrcprob` | `outfile` equals `evfile`, or exists with `clobber=no`. The second message wrongly says "clobber=yes". | **0** (`std::exit(0)`) |
| Any tool using `checkOutputFile`; `gtltcube` | Output exists and `clobber=no`. | 1 (`std::exit(1)`, bypassing `st_app`) |
| `gtltcube` | No FT2 rows match the FT1 GTIs. Writes an all-zero cube with a warning. | 0 |
| `gttsmap` | A fit at one grid point throws an `optimizers::Exception`: that pixel gets TS 0 and the scan continues. If the **null** fit throws (any exception), the reference log-likelihood is set to 0, so every TS value is then 2·logL, which is meaningless, with no message. | 0 (from reading the code) |
| `gtebl` | The constructor fails (`setVersion`). | 1 |

## Tool-specific notes

### gtexpcube2

**`irfs=CALDB` does not work** (confirmed at runtime with fermitools 2.5.3;
the code path is unchanged at this commit):

```
$ gtexpcube2 infile=ltcube.fits cmap=srcmap.fits outfile=bexp.fits irfs=CALDB evtype=3
Caught St13runtime_error at the top level: File not found:
$ gtexpcube2 ... cmap=none irfs=CALDB ...
Caught St13runtime_error at the top level: File not found: none
```

With an explicit name (`irfs=P8R2_SOURCE_V6` on the same files) it
works. Cause: `ExpCube::run` creates `AppHelpers(&m_pars, "NONE")` when a
counts map **is** given and `"BINNED"` when it is not. Commit `c4fd32dd`
(2017, "Allow for user defined healpix exposure maps") inverted the
original logic, which was `cmap != "none"` → `"BINNED"`. With `"NONE"`
there is no file to read DSS keywords from. With `"BINNED"` and
`cmap=none`, it tries to open "none". The par file still defaults `irfs`
to `CALDB`.

Other notes:
- `cmap=none` (WCS) or `cmap=hpx` (HEALPix) means "no counts map; use the
  geometry parameters". If a counts map is given **and** geometry
  parameters are set (not `INDEF`), the geometry parameters win, and the
  counts map supplies only what is missing.
- It checks GTIs against the livetime cube (`checkTimeCuts`, warning
  only) and copies the livetime cube's GTI and header keywords.
- A local `bool useEbounds(true)` in `run()` is unused; the member
  `m_useEbounds` is what `bincalc` sets.
- `gtexpcube2.par` has no `gui` parameter.

### gtsrcmaps

- Reads the model with `requireExposure=false` and `loadMaps=false`.
  Turns `use_single_fixed_map` off, so every source gets its own HDU.
- Sets `RoiCuts` to a 20° cone around the counts-map reference point (WCS)
  or the whole sky (HEALPix). This only affects the "far from ROI"
  warnings.
- `copyall=yes` copies the whole input `cmap` file (including any
  existing source-map HDUs) to `outfile` first. With `copyall=no` (the
  default) only the counts map is written. The `if (!copyall) { bool
  copyall = true; }` block is dead code (it shadows the variable).
- Because the input `cmap` is also passed as the srcmaps file to
  `BinnedLikelihood`, running `gtsrcmaps` on an existing srcmaps file
  **reuses** the maps already in it for matching source names; see
  [binned-analysis.md](binned-analysis.md#source-maps-file).
- Writes DSS keywords (adding `IRF_VERSION`, and bit-mask cuts when
  `irfs` is not `CALDB`) and the GTI.

### gtlike

- `statistic` is a par-file switch: `BINNED` enables `cmap`, `bexpmap`,
  `wmap`, `psfcorr`, `phased_expmap`, `edisp_bins`; `UNBINNED` enables
  `evfile`, `evtable`, `scfile`, `sctable`, `expmap`.
- `edisp` is in the par file but never read; unbinned energy dispersion
  is not supported. Binned dispersion uses `edisp_bins`.
- `wmap` is read twice: once in `run()` into an unused member, once in
  `makeBinnedLikelihood`.
- TS values in `results.dat`: for each free source, the source is
  removed and the null model renormalized (`tsmin=no`, default) or refit
  (`tsmin=yes`).
- `chatter` is passed to the optimizer minus 1.
- `refit=yes` loops, re-reading parameter values from the XML each time.
- Plotting needs a working Python with matplotlib at run time (embedded
  through `libpython`). The `st_graph` includes are leftovers.

### gtdiffrsp

- Opens the FT1 file with `editTable` and **appends columns in place**.
  Files without `NDIFRSP` need `convert=yes` (rewrites the header to the
  `DIFRSPn` scheme), otherwise it throws.
- `evclsmin`/`evclass` are pre-Pass 8 class selections.

### gtsrcprob

- Copies FT1 (`FitsUtil::fcopy`) and appends one `1E` column per source
  in `srclist` (default: all sources). Values are
  `fluxDensity_i / Σ_all fluxDensity`.
- Diffuse responses must already be in the FT1 file (`gtdiffrsp`) for the
  same IRF name; otherwise they are computed on the fly.
- Commit `5cac2f8c` (September 2026) fixed a large slowdown: point
  sources had been computing exposure even though `gtsrcprob` passes
  `requireExposure=false`.

### gtpsf

- The `irfs` default is the specific, superseded name `P8R2_SOURCE_V6`,
  not `CALDB`.
