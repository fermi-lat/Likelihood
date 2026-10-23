# Likelihood: models and sources

Part of the [Likelihood architecture docs](README.md). This file covers
the sky model: `SourceModel`, the `Source` hierarchy, spectral and spatial
functions, the XML source-model format, and function registration.

## The model container: `SourceModel`

`SourceModel` (`Likelihood/SourceModel.h`) derives from
`optimizers::Statistic`. It holds `std::map<std::string, Source*>
m_sources` and a reference to the `Observation`. `LogLike` (unbinned) and
`BinnedLikelihood` (binned) derive from it.

How parameters reach the optimizer:
- `syncParams()` rebuilds `m_parameter` by walking `m_sources` **in map
  order (alphabetical by source name)** and appending each source's
  **spectrum** parameters in function order.
- `setParams_()` / `setFreeParamValues_()` walk the same order to push
  values back into each `Source`'s spectrum.
- **Only spectral parameters are fit.** Spatial parameters (RA/Dec of a
  `SkyDirFunction`, the width of a `RadialGaussian`, …) never enter the
  parameter vector. Positions are fit outside the likelihood by
  `gtfindsrc` (Nelder–Mead around `LogLike`) or by pyLikelihood scripts
  that move sources and refit.
- `getParam(paramName, funcName, srcName)` and `setParam(...)` address a
  single parameter by triple; `funcName` is `"Spectrum"`, `"Position"`
  or `"SpatialDist"` (the keys of `Source::getSrcFuncs()`).

Other `SourceModel` responsibilities:
- `addSource` / `deleteSource` / `getSource` / `steal_source` /
  `give_source` (move a source between two models with its `SourceMap`).
- `mergeSources` / `splitCompositeSource`: build or break up a
  `CompositeSource`.
- `readXml` / `reReadXml` (parameters only) / `writeXml` /
  `write_fluxXml`.
- `createCountsMap`: model counts map (overridden by `BinnedLikelihood`).

`set_ebounds`/`unset_ebounds` throw in `SourceModel` and
`BinnedLikelihood`; only `LogLike` implements them (see
[unbinned-analysis.md](unbinned-analysis.md)).

## Sources

```
Source (abstract; spectrum + map of named functions; exposure arrays)
├── PointSource      "Point"      functions: Spectrum, Position (SkyDirFunction)
├── DiffuseSource    "Diffuse"    functions: Spectrum, SpatialDist
└── CompositeSource  "Composite"  one shared Spectrum scaling a nested SourceModel
```

| Class | Key points |
|-------|-----------|
| `Source` | Holds `m_spectrum` (an `optimizers::Function*`, cloned in), `m_functions` (name → function), the energy grid `m_energies` and the exposure `m_exposure` used for `Npred`. Flags: `use_edisp()` (XML `apply_edisp`), `use_psf()` (`apply_psf`), `use_exposure()` (`apply_exposure`). `fluxDensity(Event)` is the response-folded model density at an event (unbinned). `Npred()` integrates spectrum × exposure over `m_energies` with log–log trapezoids (`Source::pixelCounts`, linear if `USE_LINEAR_QUADRATURE`). `flux()`/`energyFlux()` integrate the spectrum only. |
| `PointSource` | Position in a `SkyDirFunction`. `computeExposure` integrates effective area × the PSF fraction inside the ROI (acceptance cone) over the livetime cube (`computeExposureWithHyperCube`) or, without a cube, over `ScData` intervals. `setDir(ra, dec, updateExposure)` moves it. See [exposure-and-livetime.md](exposure-and-livetime.md). |
| `DiffuseSource` | Spatial model is any `optimizers::Function` taking a `SkyDirArg`. Unbinned use needs either the exposure map (`ExposureMap::integrateSpatialDist`) or map-based integration. `mapBaseObject()` gives the `MapBase` side of map-based models. |
| `CompositeSource` | Has its own `SourceModel` of sub-sources and a single spectrum that scales them all (`m_functions["Spectrum"]`). Binned only in practice: source maps are summed from the components (`PSFUtils::makeCompositeMap`). Read from XML as `<source type="CompositeSource">` with a nested `<source_library>`. |

`SourceFactory` (`SourceFactory.h`) parses XML into prototype `Source`
objects keyed by name. `create(name)` returns a clone. `SourceModel::readXml`
uses a `SourceFactory` and then `addSource`s each clone.

## Function registry

`AppHelpers::addFunctionPrototypes(FunctionFactory*)` is **the** list of
model names. Tools call it through `AppHelpers`; pyLikelihood calls it
through its `SourceFactory.funcFactory()` `%extend`.
`optimizers::FunctionFactory`'s constructor already registers
`PowerLaw`, `BrokenPowerLaw`, `Gaussian`, `AbsEdge` and `ConstantValue`
(used as the isotropic spatial model). Likelihood adds:

| Kind | Registered names |
|------|-----------------|
| Spectral | `PowerLaw2`, `BrokenPowerLaw2`, `BrokenPowerLaw3`, `MultipleBrokenPowerLaw`, `PiecewisePowerLaw`, `SmoothBrokenPowerLaw`, `SmoothDoubleBrokenPowerLaw`, `LogParabola`, `LogGaussian`, `LogNormal`, `LogNormalLog`, `BandFunction`, `GaussianError`, `ExpCutoff`, `ExpCutoffSEDPeak`, `BPLExpCutoff` (class `BrokenPowerLawExpCutoff`), `PLSuperExpCutoff`, `PLSuperExpCutoff2`, `PLSuperExpCutoff3`, `PLSuperExpCutoff4`, `FileFunction`, `DMFitFunction` |
| Decorators | `EblAtten::<base>`, `EnergyBand::<base>` and `ScaleFactor::<base>` for most of the spectral models above (not `MultipleBrokenPowerLaw` or `PiecewisePowerLaw`), plus `ScaleFactor::Gaussian` |
| Spatial | `SkyDirFunction`, `SpatialMap`, `MapCubeFunction` (class `MapCubeFunction2`), `RadialProfile`, `RadialDisk`, `RadialGaussian` |

Notes on specific functions:
- **Decorators** wrap a prototype of the base function and add their own
  parameters:
  - `EblAtten` multiplies by exp(−τ(E, z)) from `IRB::EblAtten`
    (`celestialSources/eblAtten`). Its parameters are `tau_norm`,
    `redshift` and `ebl_model` (the `IRB::EblModel` enum index, 0–19;
    always fixed).
  - `EnergyBand` zeroes the function outside [`Emin`, `Emax`] (both
    always fixed).
  - `ScaleFactor` multiplies by a free `ScaleFactor` parameter.
- **`FileFunction`** reads an ASCII table of (E, dN/dE) from the XML
  `file=` attribute and interpolates in log–log. `ScaleFactor::FileFunction`
  passes the file to the inner function.
- **`DMFitFunction`** (dark-matter annihilation spectra) is C++ glue
  (`src/dmfit/DMFitFunction.cxx`) over three f2c-translated Fortran files
  (`dmfit_comm.c`, `dmfit_func.c`, `dmfit_load.c`, built as the
  `dmfit_obj` OBJECT library). `dmfit_load` reads a 72,576-value table
  into a **Fortran COMMON block** (`hasim_`), so there is one table per
  process. Two `DMFitFunction`s with different `file=` attributes
  overwrite each other's table. The table ships as
  `data/gammamc_dif.dat` (installed to `refdata/fermi/Likelihood`). If
  the XML has no `file=`, the default path is
  `$(LIKELIHOODDATAPATH)/gammamc_dif.dat`, and nothing sets
  `LIKELIHOODDATAPATH` any more, so a `file=` attribute is effectively
  required. The `.f` originals, `test.f` and `test.py` are kept for
  reference and are not built.
- **Any function** can take a `scaling_file=` attribute on `<spectrum>`.
  This attaches a `FileFunction` via `optimizers::Function::setScalingFunction`.

Spatial models:

| Name | Class | Notes |
|------|-------|-------|
| `SkyDirFunction` | `SkyDirFunction` | Point-source position. **Required** for `type="PointSource"`, otherwise `SourceFactory` throws. |
| `SpatialMap` | `SpatialMap` (`Function` + `MapBase`) | 2-D template (WCS or HEALPix) from `file=`, via `WcsMapLibrary` (cached per process). `map_based_integral="true"` integrates exposure over map pixels instead of the exposure-map grid. |
| `MapCubeFunction` | `MapCubeFunction2` | 3-D template (energy planes) from `file=`, interpolated in energy. The Galactic diffuse model. |
| `RadialProfile` | `RadialProfile` | ASCII radial profile from `file=`. |
| `RadialDisk`, `RadialGaussian` | `SpatialFunction` subclasses | Analytic extended sources. `RadialGaussian` uses GSL `gsl_sf_bessel_I0_scaled`. They get their parameters through `SpatialFunction::setParams(xml_node*)`, which takes a non-const RapidXML node (hence the `const_cast` in `SourceFactory`). |
| `ConstantValue` | from `optimizers` | Isotropic background. |

`SpatialFunction` (the base for analytic extended sources) also provides
`diffuseResponse` (unbinned) and `spatialResponse` (binned, PSF-convolved)
through `ResponseFunctor`s.

## XML source-model format

Reading: `SourceFactory::readXml` (`src/SourceFactory.cxx`).
- The file path goes through `facilities::Util::expandEnvVar`. The whole
  file is read into a buffer owned by the factory and parsed with
  **RapidXML** (vendored in `xmlBase`; see `src/xmlBase/doc/architecture/`).
  The xerces → RapidXML port landed in February–April 2026.
- The root element must be `<source_library>`. If it has a
  `function_library="...xml"` attribute, that file is loaded into the
  `FunctionFactory` first (`optimizers::FunctionFactory::readXml`).
- Each `<source type=… name=…>` needs exactly one `<spectrum>` and (except
  for `CompositeSource`) exactly one `<spatialModel>`.
- `<parameter>` elements are applied by name with
  `optimizers::Parameter::extractDomData` (value, scale, min, max, free,
  error).
- Attributes:

| Element | Attribute | Effect |
|---------|-----------|--------|
| `spectrum` | `type` | Registered function name. |
| `spectrum` | `file` | Table for `FileFunction`, `ScaleFactor::FileFunction`, `DMFitFunction`. |
| `spectrum` | `scaling_file` | Multiplicative `FileFunction`. |
| `spectrum` | `apply_edisp` | `true`/`false`; default true. False means "this model already includes energy dispersion" (e.g. the Galactic diffuse model). Any other value throws. |
| `spatialModel` | `type` | Registered spatial name. |
| `spatialModel` | `file` | Template for `SpatialMap`, `MapCubeFunction`, `RadialProfile` (env vars expanded). |
| `spatialModel` | `map_based_integral` | `true` → map-pixel exposure quadrature. |
| `spatialModel` | `apply_psf`, `apply_exposure` | `true`/`false`; default true. False skips PSF convolution / exposure when building source maps (binned). |

Options passed by callers:
- `requireExposure` — when false, `PointSource`s are built without
  computing exposure (`setDir(..., updateExposure=false)`) and
  `DiffuseSource`s skip it too. Binned tools and `gtsrcprob`, `gtmodel`,
  `gttscube`, `gthealcube` pass false. The fix in commit `5cac2f8c`
  (September 2026) made point sources honour this flag again; before it,
  `gtsrcprob` computed unneeded exposures (a large slowdown).
- `addPointSources` — when false, point sources are skipped (used for
  diffuse-only setups).
- `loadMaps` — when false, `SpatialMap`/`MapCubeFunction` defer reading
  the FITS pixels until first use (`MapBase::readFitsFile(..., loadMap)`).
  `gtlike` and `gtsrcmaps` pass false.

Writing:
- `SourceModelBuilder` (`writeXml`) writes the same format, including
  the `file`, `apply_*` attributes and parameter errors. Serialization of
  the `error` attribute was added in September 2026 (`cad8d93a`).
- `FluxBuilder` (`write_fluxXml`) writes a **flux-package** XML (the
  `gtobssim` format; see `src/flux/doc/architecture/source-xml.md`) so a
  fitted model can be simulated.
- Both derive from `XmlBuilder` and use RapidXML.

`SourceModel::reReadXml` re-reads only parameter values for sources that
already exist. `gtlike refit=yes` uses it.

## Combining likelihoods

| Class | Use |
|-------|-----|
| `SummedLikelihood` | Sum of `LogLike`/`BinnedLikelihood` components that have **identical source models** (joint fit of several data sets, e.g. PSF event types). The parameter vector is that of the first component and is pushed to all of them. pyLikelihood's `SummedLikelihood.py` wraps it. `FitScanner` has a summed-model wrapper. |
| `Composite2` | Several likelihoods with explicitly tied parameters (`TiedParameter`). |
| `CompositeLikelihood` | Stacking: several `LogLike`s whose target sources share the same spectral parameters (written for the dwarf-spheroidal dark-matter analysis). |
| `OptEM` | Expectation–maximization fitter over a `LogLike`. Used only from Python. |
| `OneSourceFunc` | Extended likelihood for one source (helper `Statistic`). |

## Gotchas

**Silent drops and name mismatches** (all confirmed at runtime with
fermitools 2.5.3 using pyLikelihood's `SourceFactory`; the logic is the same
at this commit)
- A `<source>` with an unrecognised `type` (e.g. `Pointsource`) is
  **skipped silently**: `makeSources` only acts on the three exact
  strings. The error appears later, when something asks for the source by
  name. A point source is also skipped silently when `addPointSources` is
  false.
- **`MultipleBrokenPowerLaw` cannot be read from XML.** The class's
  `genericName` is `MultipleBPL`, but it is registered as
  `MultipleBrokenPowerLaw`, and `SourceFactory::setSpectrum` adds the
  `Index*`/`Break*` parameters only when `type == "MultipleBPL"`.
  - `type="MultipleBrokenPowerLaw"`: the prototype has no parameters, so
    the first `<parameter>` fails with "A Parameter named Normalization is
    not a Parameter of Function".
  - `type="MultipleBPL"`: "Cannot create Function named 'MultipleBPL'".
  - Writing uses `genericName`, so it cannot round-trip either.
  - This has been broken since at least 2016 (`d1388981`).
  `test_MultipleBrokenPowerLaw` tests the function directly, not through
  XML.
- `PiecewisePowerLaw` from XML always gets the default outer indices
  (−2, −3) from `addParamsToPiecewisePL`; the `dNdE*` values come from the
  later generic parameter pass. Its behaviour also depends on
  `PIECEWISEPOWERLAW_DECOUPLED`.
- For point sources, `RA`/`DEC` are read with
  `getAttributeValueAsDouble(…, "value")`: a `scale` attribute on these
  parameters is ignored.

**Ownership**
- `Source::setSpectrum` clones the function. `SourceFactory::setSpectrum`
  deletes its temporary afterwards. Holding on to a `Function*` you passed
  in does **not** let you change the source.
- `SourceModel::addSource(src, fromClone=true)` clones `src`. With
  `fromClone=false` it takes ownership.
- `Observation` owns nothing. Destroying the `AppHelpers` (or the Python
  objects) that own the pieces leaves a `SourceModel` with dangling
  references.

**Spatial templates**
- `WcsMapLibrary` caches templates by `filename::extension` for the life
  of the process. Two sources using the same file share one map object.
  Changes to the file on disk are not seen.
- `checkRoiDist` only **warns** when a point source lies more than
  ROI radius + 10° from the ROI centre.

**DMFit**
- One table per process (Fortran COMMON block); see above.
- The installed `libLikelihood.so` (fermitools 2.5.3) leaves the f2c I/O
  symbols `f_open`/`f_clos` undefined and has no `NEEDED` entry for
  `libf2c`. `DMFitFunction::readFunction` works from pyLikelihood
  (confirmed at runtime), so the symbols are resolved there. This was not
  checked for the C++ tools.

## Extension points

See the README's [Extension points](README.md#extension-points) for
adding spectral and spatial functions. For a new **source type**:
1. Subclass `Source`. Implement `fluxDensity`, its derivatives, `flux*`,
   `energyFlux*` and `clone`. Set `m_srcType`.
2. Add a branch in `SourceFactory::makeSources` and in
   `SourceModelBuilder`.
3. Teach `PSFUtils::makeModelMap` and `SourceMap` how to build its binned
   map, and `LogLike`/`Event` how to compute its unbinned response.
