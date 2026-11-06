# Proposal: per-pixel threading for `gtltcube`

Status: plan agreed 2026-11-06, not implemented. Nothing here has been
benchmarked yet. The work order, tests and benchmarks are in
[Implementation plan](#implementation-plan).

## Summary

`gtltcube` spends almost all of its time in `map_tools::Exposure::fill`.
Each call adds one FT2 interval to every HEALPix pixel of the livetime
cube. The work for each pixel is independent of every other pixel, so we
can split the pixels across threads with OpenMP. The result is
**bit-identical** to the serial code for any thread count, because each
bin still gets the same additions in the same row order.

We propose two stages:

1. **Stage 1 (small change, in `map_tools`):** parallelize the pixel loop
   inside each `Exposure::fill*` call. This is about 30 lines in
   `map_tools/src/Exposure.cxx` plus a CMake link line. No API changes.
2. **Stage 2 (larger change, in `map_tools` and `Likelihood`):** have
   `LikeExposure::load` collect blocks of accepted FT2 rows and fill both
   cubes (`EXPOSURE` and `WEIGHTED_EXPOSURE`) in one parallel pass over
   pixels. This cuts thread start-up overhead and memory traffic by
   roughly the block size.

Ship Stage 1 first. It is also the baseline for deciding whether Stage 2
is worth doing.

The alternative is splitting by time (a partial cube per time chunk,
summed at the end, as `gtltsum` and `GtBurst/gtapps_mp/gtltcube_mp.py`
do). That is also correct, but it needs one copy of both cubes per
thread, and the float32 bins round differently, so the output is not
bit-identical. See [Why not split by time](#why-not-split-by-time).

## How the code works today

`gtltcube` (`src/makeExposureCube/makeExposureCube.cxx`) builds a
`LikeExposure`, which holds two `map_tools::Exposure` objects: itself
(`EXPOSURE`) and `m_weightedExposure` (`WEIGHTED_EXPOSURE`). Each one is
a `HealpixArray<CosineBinner>`. A `CosineBinner` is a
`std::vector<float>` of cos θ bins, plus φ bins if `phibins > 0`.

`LikeExposure::load` (`src/LikeExposure.cxx`) walks the FT2 rows. For
each row that overlaps the GTIs it makes **two** calls:

```cpp
fill(dirz, zenith, deltat*fraction);                              // EXPOSURE
m_weightedExposure->fill(dirz, zenith, deltat*fraction*weight);   // WEIGHTED_EXPOSURE
```

or the `fill_zenith(dirz, dirx, zenith, ...)` versions when
`phibins > 0`.

Each `fill*` (`map_tools/src/Exposure.cxx`) runs `std::for_each` over
`m_dir_cache`, a `std::vector<std::pair<CosineBinner*, Simple3Vector>>`
with one entry per pixel. `Exposure::Filler::operator()`, for one pixel:

1. If a zenith cut is active, it skips the pixel when
   `pixeldir·zenith` is outside (`m_zcut`, `m_zmaxcut`).
2. It computes cos θ: `pixeldir·dirz`. With φ binning it instead rotates
   `pixeldir` into instrument coordinates and takes `z()` and `phi()`.
3. It calls `CosineBinner::fill(costheta[, phi], deltat)`, which does
   `bin += deltat`. The φ version adds `deltat` to two bins: the cos θ
   bin and the matching φ bin.
4. It adds `deltat` to the Filler's `m_total` (accepted) or `m_lost`
   (cut). These only feed `Exposure::total()` and `lost()`, which are
   printed by `map_tools`' `exposure_cube` but **not written to the
   ltcube file**.

At the default `binsz=1` there are 49,152 pixels (nside 64); at
`binsz=0.5` there are 196,608 (nside 128). Each `CosineBinner` holds
`nbins × (1 + phibins)` floats: 40 with no φ bins, 440 at `phibins=10`.
The pointing transform is built once per row, in the `Filler`
constructor. All per-pixel work only reads the row's vectors and writes
that pixel's own `CosineBinner`.

`map_tools`' own `Exposure::load` / `processEntry` (used by
`exposure_cube`) also bins through `fill_zenith`, so it picks up Stage 1
too.

## Stage 1: parallel pixel loop inside `Exposure::fill*`

Keep `Filler` exactly as it is. Give each thread its own copy of the
row's `Filler`, split `m_dir_cache` across threads with a static
schedule, and add up the per-thread `total`/`lost` at the end.

```cpp
// map_tools/src/Exposure.cxx

#include <exception>

namespace {
   // Apply one row's Filler to every pixel, splitting pixels across
   // OpenMP threads. Each pixel's CosineBinner is touched by exactly one
   // thread, so no locking is needed.
   template <class Cache, class F>
   void parallel_fill(Cache & cache, const F & proto, int nthreads,
                      double & total, double & lost) {
      const long npix(static_cast<long>(cache.size()));
      double my_total(0), my_lost(0);
      std::exception_ptr error;
      // With nthreads == 1 (the default) the region runs serially on the
      // calling thread, whatever OMP_NUM_THREADS says.
#pragma omp parallel num_threads(nthreads) if(nthreads > 1) \
                     reduction(+:my_total, my_lost)
      {
         F filler(proto);
#pragma omp for schedule(static)
         for (long i = 0; i < npix; i++) {
            try {
               filler(cache[i]);
            } catch (...) {
               // An exception must not leave the parallel region (that
               // calls std::terminate). Keep the first one and rethrow
               // it on the calling thread.
#pragma omp critical(exposure_fill_error)
               if (!error) error = std::current_exception();
            }
         }
         my_total += filler.total();
         my_lost += filler.lost();
      }
      if (error) std::rethrow_exception(error);
      total = my_total;
      lost = my_lost;
   }
}

void Exposure::fill(const astro::SkyDir& dirz, const astro::SkyDir& zenith, double deltat)
{
    double total, lost;
    parallel_fill(m_dir_cache, Filler(deltat, dirz, zenith, m_zcut, m_zmaxcut),
                  numThreads(), total, lost);
    addtotal(total);
    m_lost += lost;
}
```

The thread count belongs to each `Exposure` object, not to the process,
so that no other tool changes behaviour:

```cpp
// map_tools/map_tools/Exposure.h
public:
    /// Threads used by fill*() and fill_block(). The default is 1
    /// (serial). 0 means the OpenMP default (OMP_NUM_THREADS, otherwise
    /// every core).
    void setNumThreads(int n) { m_nthreads = n; }
    /// The resolved count: 0 becomes omp_get_max_threads().
    int numThreads() const;
private:
    int m_nthreads; ///< set to 1 in every constructor
```

Both constructors initialise `m_nthreads(1)`, including the one that
reads a file, and the implicit copy constructor copies it.
`LikeExposure` adds a `setNumThreads(int)` that sets both itself and
`m_weightedExposure`. `gtltcube` calls it with its `nthreads` parameter
before `load()`. Nothing calls `omp_set_num_threads`, so other OpenMP
code in the same process is not affected either.

`Filler` is a private nested class of `Exposure`. A free function
template can still use it, because `F` is deduced and the private name
`Exposure::Filler` is never written outside the class. If a compiler
objects, make `parallel_fill` a private static member template instead.
The body stays the same. `fill(dirz, deltat)` and
`fill_zenith(...)` change the same way. `fill(dirz, deltat)` keeps its
current behaviour of calling `addtotal(deltat)` and ignoring the Filler's
totals.

The `Filler` copy is the per-thread state. It holds the row's
`HepRotation` and vectors, which are read-only, and its own
`m_total`/`m_lost` accumulators.

### Build change

The top-level `CMakeLists.txt` already does
`find_package(OpenMP REQUIRED COMPONENTS CXX)`, and `Likelihood` already
links `OpenMP::OpenMP_CXX`. `map_tools` does not, so add:

```cmake
# map_tools/CMakeLists.txt
target_link_libraries(map_tools PUBLIC tip healpix astro hoops st_app st_stream irfLoader
                      OpenMP::OpenMP_CXX)
```

This affects every user of `map_tools::Exposure`: `Likelihood`
(`gtltcube`, `ExposureCube`), `map_tools`' own `exposure_cube` and
`exposure_map`, `SolarSystemTools` (`ExposureSun`, `LikeExposureSun`),
`skymaps` and the pyLikelihood SWIG wrapper. None of them change
behaviour, apart from `total()`/`lost()` rounding (see
[Determinism](#determinism)).

### Expected cost and limits

- **Fork/join for every row, twice.** Every accepted FT2 row now opens
  two parallel regions, one per cube. Each region does 49,152 pixel
  updates, which should be much more work than the few microseconds an
  OpenMP region costs. This needs measuring, especially for larger
  `binsz`, where there are fewer pixels.
- **Memory bandwidth.** Each row streams through every pixel's
  `CosineBinner`:

  | Case | Floats per pixel | Per cube | Both cubes |
  |------|------------------|----------|------------|
  | `binsz=1`, no φ | 40 (160 B) | 7.9 MB | 16 MB |
  | `binsz=0.5`, no φ | 40 (160 B) | 31 MB | 63 MB |
  | `binsz=1`, `phibins=10` | 440 (1.76 KB) | 86 MB | 173 MB |

  With φ bins each fill also writes two floats per pixel, not one.
  Whether the cubes fit in cache depends on the machine. The
  development laptop (i7-11800H) has 1.25 MB of L2 per core and 24 MB of
  shared L3, so both default cubes fit in L3 there, and Stage 1 may
  scale well in the default case. The 0.5° and φ cases don't fit
  anywhere, and a server part with less L3 per core may not fit even
  the default case. Once the cubes don't fit, Stage 1 may stop scaling
  when it saturates memory bandwidth. This is the main reason for
  Stage 2, and it is why the benchmarks need more than one machine.
- **Duplicate geometry.** The two cubes recompute the same cos θ, φ and
  zenith test for every pixel. Stage 2 removes this.

## Stage 2: blocked rows, both cubes in one pass

Turn the loops inside out: for a **block** of FT2 rows, give each thread
a range of pixels, and have it apply every row in the block, in order,
to its pixels in both cubes. While a thread works through a block, each
pixel's bins stay in its cache.

### New `map_tools::Exposure` API

```cpp
// map_tools/map_tools/Exposure.h
public:
    /// One accepted FT2 interval, ready to be binned.
    struct Interval {
        astro::SkyDir dirz, dirx, zenith;
        double deltat;          ///< added to *this
        double weighted_deltat; ///< added to the weighted cube, if any
    };

    /// Add every interval in rows, in order, to *this and (if non-null)
    /// to weighted. weighted must have been built with the same pixel
    /// size, cos-theta binning and zenith cuts as *this.
    void fill_block(const std::vector<Interval>& rows,
                    bool use_phi,
                    Exposure* weighted = 0);
```

```cpp
// map_tools/src/Exposure.cxx
void Exposure::fill_block(const std::vector<Interval>& rows, bool use_phi,
                          Exposure* weighted)
{
    // Serial pre-pass: do the astro/CLHEP work once per row, on this
    // thread, and keep the results in plain arrays.
    struct RowGeom {
        Simple3Vector dirz, zenith;
        CLHEP::HepRotation rot;
        double deltat, wdeltat;
    };
    std::vector<RowGeom> geom;
    geom.reserve(rows.size());
    for (const Interval& r : rows) {
        RowGeom g;
        g.dirz = Simple3Vector(r.dirz());
        g.zenith = Simple3Vector(r.zenith());
        if (use_phi) {
            g.rot = astro::PointingTransform(r.dirz, r.dirx).localToCelestial().inverse();
        }
        g.deltat = r.deltat;
        g.wdeltat = r.weighted_deltat;
        geom.push_back(g);
    }

    const bool cut(!(m_zcut == -1 && m_zmaxcut == 1));
    const long npix(static_cast<long>(m_dir_cache.size()));
    double total(0), lost(0);

    const int nthreads(numThreads());
#pragma omp parallel for schedule(static) num_threads(nthreads) \
                         if(nthreads > 1) reduction(+:total, lost)
    for (long i = 0; i < npix; i++) {
        CosineBinner& pix(*m_dir_cache[i].first);
        CosineBinner* wpix(weighted ? weighted->m_dir_cache[i].first : 0);
        const Simple3Vector& pixeldir(m_dir_cache[i].second);
        for (const RowGeom& g : geom) {        // row order preserved per pixel
            if (cut) {
                double z(pixeldir.dot(g.zenith));
                if (!(z > m_zcut && z < m_zmaxcut)) { lost += g.deltat; continue; }
            }
            if (use_phi) {
                CLHEP::Hep3Vector inst(pixeldir.transform(g.rot));
                double costheta(inst.z()), phi(inst.phi());
                pix.fill(costheta, phi, g.deltat);
                if (wpix) wpix->fill(costheta, phi, g.wdeltat);
            } else {
                double costheta(pixeldir.dot(g.dirz));
                pix.fill(costheta, g.deltat);
                if (wpix) wpix->fill(costheta, g.wdeltat);
            }
            total += g.deltat;
        }
    }
    addtotal(total);
    m_lost += lost;
    // The weighted cube's total()/lost() are not updated here. Nothing
    // reads them for the weighted cube today; if that matters, track
    // them separately.
}
```

The pixel at index `i` in `weighted->m_dir_cache` must be the same
pixel as index `i` in `this->m_dir_cache`. That holds when both objects
come from the same constructor arguments, which is how `LikeExposure`
builds them. `fill_block` should check that the two cache sizes match
and throw if they don't.

### `LikeExposure::load` changes

Replace the two `fill*` calls per row with appending to a block, and
flush the block when it is full and after the loop:

```cpp
// Likelihood/src/LikeExposure.cxx (inside load)
std::vector<map_tools::Exposure::Interval> block;
block.reserve(s_blockSize);
const bool use_phi(healpix::CosineBinner::nphibins() != 0);
...
if (acceptInterval(start, stop, m_timeCuts, m_gtis, fraction)) {
   ...
   map_tools::Exposure::Interval r;
   r.dirz = astro::SkyDir(ra, dec);
   r.dirx = astro::SkyDir(rax, decx);
   r.zenith = astro::SkyDir(ra_zenith, dec_zenith);
   r.deltat = deltat*fraction;
   r.weighted_deltat = deltat*fraction*weight;   // same expression as today
   block.push_back(r);
   if (block.size() == s_blockSize) {
      fill_block(block, use_phi, m_weightedExposure);
      block.clear();
   }
   m_numIntervals++;
}
...
if (!block.empty()) {
   fill_block(block, use_phi, m_weightedExposure);
}
```

`s_blockSize` is a compile-time constant, not a user parameter. The
pixel loop reads the whole block of rows once per pixel, so the block
should stay in L2. A `RowGeom` is about 140 bytes: two 3-vectors, a
3×3 rotation and two doubles. So 4,096 rows is about 560 KB, almost
half the L2 on the development laptop (1.25 MB), and it shares that L2
with the pixels' bins. Even 256 rows cuts the traffic to the cubes by
256×. The benchmark sweeps 128, 256, 512, 1024 and 4096 rows and picks
the knee of the curve (see [Implementation plan](#implementation-plan)).

The per-pixel body of `fill_block` needs the same try/catch and rethrow
as `parallel_fill` (see [Exceptions](#exceptions)).

The progress dots still print from the row loop on the main thread, so
they no longer line up exactly with the binning work. That is
cosmetic.

## Thread safety

What the parallel loops touch, and why each one is safe:

| Item | Access in the parallel region | Why it is safe |
|------|-------------------------------|----------------|
| `CosineBinner` per pixel | read/write | Each pixel index belongs to exactly one thread. |
| `CosineBinner::s_nbins`, `s_cosmin`, `s_sqrt_weight`, `s_phibins` (process-wide statics) | read | Set before any filling: `setPhiBins` in `makeExposureCube::createDataCube`, `setBinning` in the `Exposure` constructor. Nothing changes them during `load`. |
| `m_dir_cache[i].second` (pixel direction) | read | Built once in `create_cache`. |
| Row geometry (`Filler` copy or `RowGeom`) | read | Built on the calling thread before the region. |
| `m_zcut`, `m_zmaxcut` | read | Constant after construction. |
| `total` / `lost` | write | OpenMP `reduction`. |
| `astro::SkyDir`, `astro::PointingTransform` | not touched | All `astro` calls stay on the calling thread, so we don't depend on `astro` being thread-safe. |

**False sharing.** Each `CosineBinner`'s data is its own heap block.
With `schedule(static)`, each thread gets one contiguous range of
pixels, so two threads can only write to the same cache line at the
edges of their ranges, if neighbouring blocks were allocated next to
each other. That effect is negligible.

**Existing hazard to fix first.** `LikeExposure` stores `m_timeCuts` and
`m_gtis` as **const references** to locals in
`makeExposureCube::createDataCube`. These references already dangle:
`m_exposure` outlives `createDataCube`, because `ExposureCube::run`
calls `m_exposure->writeFile()` after `createDataCube` returns. Today
this is harmless only because nothing after `load()` reads them;
`writeFile` uses the copied `m_tmin`/`m_tmax`. Any later change that
reads them afterwards, or loads on another thread, would read freed
memory. Changing both members to stored copies is step 0 of the plan.
(Also listed in the Gotchas of
[../architecture/exposure-and-livetime.md](../architecture/exposure-and-livetime.md).)

### Exceptions

An exception that escapes an OpenMP parallel region calls
`std::terminate`. Today, an exception thrown while filling reaches
`gtltcube`'s top-level catch and is reported. Both stages must catch
inside the region, keep the first exception, and rethrow it on the
calling thread after the region (see `parallel_fill` above). On the
normal path this costs nothing.

At least one throw can happen in practice.
`CosineBinner::fill(costheta, phi, value)` (`healpix/src/CosineBinner.cxx`)
clamps the cos θ index but not the φ index. `phi_index(phi)` returns
`s_phibins` when `folded_phi(phi)` is exactly 1, which happens when
`fmod(phi, π/2)` is exactly π/4, for example when a pixel's instrument
x and y components are exactly equal. The bin index
`(phibins+1)*nbins + ci` is then one past the end, and `at()` throws
`std::out_of_range`. It is rare, but it is reachable. Clamping the φ
index the way `ci` is clamped would not change any output the current
code produces, because those runs abort today. So it can be fixed
separately, in `healpix`. Record it in `healpix`'s known issues.

## Determinism

- **COSBINS (what goes in the file): bit-identical** to the serial code
  for any thread count and any block size. Each bin receives the same
  sequence of `float += double` updates, in FT2 row order, with values
  computed by the same expressions. Threads only change *which core*
  does a pixel, never the order of additions *within* a pixel.
- **`Exposure::total()` / `lost()`: not bit-identical.** These are sums
  over pixels, and the reduction changes their order. They are only
  printed by `map_tools`' `exposure_cube` and are not in the ltcube
  file.

- **Compiler contraction (Stage 2).** Bit-identity also assumes that the
  new code does the same floating-point operations as `Filler`. Stage 1
  reuses `Filler` unchanged. Stage 2 rewrites the dot product and the
  rotation inside `fill_block`. If the compiler fuses a multiply and an
  add (FMA) in one place but not the other, the last bits differ.
  - x86-64 without `-march` flags that enable FMA cannot emit FMA, so
    the two paths cannot differ there. Our RelWithDebInfo build uses
    `-O2 -g` with no `-march`.
  - arm64 (Apple silicon) has FMA in its base instruction set, and
    clang contracts within an expression by default.

  To keep the two paths the same, `fill_block` should call the same
  `Simple3Vector::dot` and `transform` members that `Filler` uses,
  rather than writing the arithmetic out again. Don't build
  `Exposure.cxx` with `-ffp-contract=off` to force agreement: that would
  also change the serial results on arm64 relative to the current arm64
  release.

The goal is bit-identical output on every platform, but only x86-64
Linux can be tested on the development machine. Someone with a Mac must
run the same comparison on arm64 before release.

So the regression test is strict: run `gtltcube` with `nthreads=1` and
`nthreads=N`, and compare both `EXPOSURE` and `WEIGHTED_EXPOSURE`
`COSBINS` columns for exact equality, with `phibins=0` and `phibins>0`,
and with and without a `zmax` cut. Both stages must also reproduce the
gold-standard cubes made by the unmodified code exactly. The tests and
reference outputs are in [Implementation plan](#implementation-plan).

## Controlling the thread count

OpenMP uses all available cores by default and honours
`OMP_NUM_THREADS`. Decided (2026-11-06):

- **`gtltcube` gets a hidden `nthreads` parameter, default 1.**
  `ExposureCube::createDataCube`
  (`src/makeExposureCube/makeExposureCube.cxx`) passes it to
  `LikeExposure::setNumThreads` before `load()`. `0` means "OpenMP
  default", which honours `OMP_NUM_THREADS` and otherwise uses every
  core. gt tools often run on shared batch nodes, where grabbing every
  core is unfriendly, so users opt in. The diffrsp and srcmaps proposals
  use the same convention. Nobody calls `omp_set_num_threads`; the count
  is per `Exposure` object (see Stage 1). The parameter goes in
  `pfiles/gtltcube.par` and in
  `fermitools-fhelp/fhelp_files/gtltcube.txt`.
- **Process-level wrappers.** `GtBurst/gtapps_mp/gtltcube_mp.py` runs
  one `gtltcube` per core. Because `nthreads` defaults to 1, its child
  processes stay single-threaded, and it needs no change as long as it
  doesn't pass `nthreads`.
- **Other `map_tools::Exposure` users stay single-threaded.** This
  covers `exposure_cube`, `SolarSystemTools`, `skymaps` and
  pyLikelihood. Each `Exposure` defaults to one thread, so these users
  run exactly as before, even when `OMP_NUM_THREADS` is set, until their
  own code calls `setNumThreads`. The goal for now is no impact on any
  tool other than `gtltcube`. Later, once we've looked at how they are
  used, they might default to every core or to some reasonable count
  above 1.
- **`SolarSystemTools`' `ExposureSun`** has separate fill loops with the
  same structure. It is out of scope here and becomes a follow-up once
  Stage 1 has numbers.

## Why not split by time

Filling a separate partial cube per time chunk and adding them up is
correct. The cube is a linear sum over FT2 rows, and `acceptInterval`'s
GTI-overlap fractions split a boundary-straddling row exactly. But:

- It needs one extra pair of cubes per thread: about 16 MB at nside 64
  with no φ, and much more with φ bins or finer `binsz`.
- It changes the order of float32 additions, so the output differs from
  the serial code in the last bits. Every regression comparison then
  needs a tolerance.
- It still needs a final reduction pass over all pixels.

Splitting by pixel avoids all three.

## Implementation plan

The original open questions are answered here and in
[Controlling the thread count](#controlling-the-thread-count). The
guiding rule: **before any code changes, run the unmodified tool to get
both a timing baseline and a gold-standard output.** Every later step
is checked against both. Output is compared with the gold standard;
timing is compared with the baseline at every step, not just at the
end (see [Monitoring at every step](#monitoring-at-every-step)). We may
later decide that some output differences are acceptable, but the
decision needs the unmodified baseline to compare against.

### Environment and inputs

- **Scripts run in whatever conda environment is active.** They find
  everything through `$CONDA_PREFIX` and the tools on `PATH`, and never
  hard-code an environment. To get results from another version, switch
  environments and run the same scripts again.
  - **Development build:** the `fermidev` environment. After
    `conda activate fermidev`, run `source recipe/activate.sh` from the
    top of the ScienceTools checkout. That sets `CALDB`, `PFILES`,
    `FERMI_INST_DIR` and the other variables the tools need. The
    `recipe/` directory is at the top level of every ScienceTools
    checkout.
  - **Current public release:** Fermitools 2.5.3 (the `FT-2.5.3`
    environment on the development machine).
- **Large input:**
  `$CONDA_PREFIX/share/fermitools/data/Likelihood/outref/ft2_10years_reduced.fits`.
  It is installed by the `fermitools-test-scripts-data` package, so
  anyone can use it. Its SHA-256 is in that package's `MANIFEST.sha256`.
- **Small input:** `data/oneday_scData_0000.fits` in this package. It is
  simulated, has 2,882 rows and covers MET 0–86430.
- **New data files** go in the same `outref` directory, where
  `fermitools-test-scripts-data` will eventually install them. This
  covers the gold-standard cubes, their manifests and any pre-cut FT2
  slice. Each one is listed in
  [`../../src/test/likelihood_test_suite/outref-additions.md`](../../src/test/likelihood_test_suite/outref-additions.md).
  Until the files are uploaded to the package's web host and the
  package is rebuilt, they exist only locally.
- **Timing logs** are specific to each machine and are not package data.
  They go in `src/test/ltcube_bench/results/` in this repository, one
  directory per machine.
- **Python:** the harness uses `numpy` and `astropy`. Adding them to the
  development environment and to the package requirements for
  development is a separate change.
- **Not used as a reference:** `ltcube_12years_zmax105.fits` in the same
  directory. Nobody knows the parameters it was made with. We make our
  own cubes from `ft2_10years_reduced.fits`, with the parameters
  recorded.
- **Platform:** only x86-64 Linux can be tested on the development
  machine (8-core i7-11800H under WSL2). See
  [Determinism](#determinism) for arm64.

### Configuration matrix

| Parameter | Values |
|-----------|--------|
| `binsz` | 1, 0.5 |
| `phibins` | 0, 10 |
| `zmax` | 180 (no cut), 105 |
| `dcostheta`, `zmin` | 0.025, 0 (defaults) |
| `nthreads` (changed code only) | 1, 2, 4, 8 |
| Stage 2 block size | 128, 256, 512, 1024, 4096 |

Since `gtltcube` has no event file here, `tmin`/`tmax` come from the FT2
file's own time range unless a time slice is given. Every run writes a
manifest next to its output that records the full parameter list.

**Check the cost before running the full matrix.** A rough look at the
one-day file with the current build suggested about 0.5 ms per FT2 row
at `binsz=1` without φ, and several milliseconds per row with φ or at
0.5°. If `ft2_10years_reduced.fits` still has 30-second rows (about 10
million), one serial run would take hours, and the φ and 0.5° cases much
longer. So the first harness step records the file's row count and
times a short slice of each configuration. From that we decide whether
the gold standard covers the whole file or a fixed time slice. The slice
boundaries would then be recorded in the manifest.

### Harness pieces (no library changes)

These live in `src/test/ltcube_bench/`, except the driver.

1. **`run_ltcube_matrix.py`: full-tool runs.**
   - Runs `gtltcube` for each selected configuration in the active
     environment.
   - Each run is repeated 3 times. It records the median and spread of
     wall time, and peak RSS from `os.wait4`'s rusage.
   - Results go to a JSON file that records:
     - the environment (`$CONDA_DEFAULT_ENV`, `$CONDA_PREFIX`) and the
       `CREATOR` keyword of the output;
     - the CPU model;
     - the input file and its SHA-256;
     - the configuration and the timings.
   - Timed runs go one at a time on an otherwise idle machine. A
     `--no-timing --jobs N` mode can produce gold-standard files in
     parallel when only the output matters.
2. **`bench_ltcube`: a C++ driver** (`src/test/bench_ltcube.cxx`, a
   CMake target that is built but not installed).
   - It mirrors `ExposureCube::createDataCube`: the same FT2 filter, and
     `setPhiBins` before the `LikeExposure` is built.
   - It times the FT2 read, `LikeExposure::load` and `writeFile`
     separately, so that the parallel part can be measured on its own.
   - It runs one configuration per process, because `CosineBinner`'s
     binning is process-wide static state.
   - It only builds in a development tree. The `load()` baseline
     therefore comes from the unmodified master build, and the 2.5.3
     release gets full-tool numbers only.
3. **`compare_ltcubes.py`: exact comparison of two ltcube files.**
   - Compares the `COSBINS`, `RA` and `DEC` columns of `EXPOSURE` and
     `WEIGHTED_EXPOSURE`.
   - Compares the binning keywords (`NSIDE`, `NBRBINS`, `COSMIN`,
     `THETABIN`, `PHIBINS`, `TSTART`, `TSTOP`, `ZENMAX`/`ZENMIN`), the
     `CTHETABOUNDS` table and the `GTI` extension.
   - Ignores `DATE`, `CREATOR`, `FILENAME`, `CHECKSUM` and `DATASUM`.
   - If the files are not identical, it reports for each HDU the number
     of differing bins and the maximum absolute and relative difference.
     That gives us what we need to decide whether a variation is
     acceptable.

### Tests

- **CppUnit, in `test_Likelihood` (`src/test/test.cxx`).** A new
  `test_LikeExposure_threads` builds a `LikeExposure` from
  `oneday_scData_0000.fits` at `binsz=4` (nside 16, 3,072 pixels, so it
  is fast).
  - It loads with 1 thread and with N threads, and compares both cubes
    in memory for exact equality.
  - It covers `phibins` 0 and 10, each with and without a zenith cut.
  - It never compares against a stored file, so it works the same on
    every platform.
  - It checks that a newly built `Exposure` reports `numThreads() == 1`,
    with `OMP_NUM_THREADS` set to something larger, so the
    single-thread default for other tools is tested.
  - It must restore the `CosineBinner` statics (`setPhiBins(0)`, and the
    binning) afterwards, because they are process-wide and other tests
    in the same executable depend on them.
- **`likelihood_test_suite` entries.**
  - Two new scripts: `scripts/test_ltcube.py` (`nthreads=1`) and
    `scripts/test_ltcube_threaded.py` (`nthreads=4`). Both run
    `gtltcube` on the same time slice of `ft2_10years_reduced.fits`, at
    the default `binsz=1`, `phibins=0`, `zmax=105`. Both use the same
    key, made from the gold standard. The threaded test uses a fixed 4
    threads rather than `nthreads=0`, so its nominal time doesn't depend
    on how many cores the machine has. Most development machines have at
    least 4.
  - **Slice size.** The slice must be large enough that the speed-up
    shows clearly, but no single test may take more than 1–2 minutes.
    Each run's time is roughly a fixed part (start-up, reading and
    filtering the FT2 file, writing the cube) plus a per-row part (the
    `load()` loop). Only the per-row part gets faster with threads. So
    we size the slice from the cost check in harness step 1:
    - Pick the number of rows that makes the serial (`nthreads=1`) run
      about 60–90 s.
    - The rough estimate of 0.5 ms per row at `binsz=1` without φ gives
      about 120,000–180,000 rows. With 30-second rows that is roughly
      6–8 weeks of data.
    - If the fixed part is a few seconds, a 4-thread run that is 3×
      faster is 20–30 s, so the speed-up is unmistakable.
    - **Risk: the FT2 read may not shrink with the slice.** `gtltcube`
      reads the FT2 file through a CFITSIO row filter, which may scan
      the whole ten-year file whatever the slice. If the cost check
      shows that this read alone takes a large share of the 1–2 minute
      budget, the speed-up gets hidden behind it. The fix is a
      pre-cut FT2 file containing only the slice, added to
      `fermitools-test-scripts-data` next to the ten-year file, and
      recorded in its `MANIFEST.sha256`.
  - **How the speed-up shows up.** The test output can't include
    timings, because it is compared against the key. The suite already
    reports each test's time against its nominal time, so the speed-up
    is read from the results table:
    - Both tests start with `NOMINAL_TIME` set to the baseline serial
      time.
    - After Stage 1, `test_ltcube` must stay within the usual 10% band.
      That shows the threading adds no overhead at one thread.
    - `test_ltcube_threaded` should come in well under its nominal
      time, which the suite shows in amber.
    - Once the speed-up is confirmed, set the threaded test's
      `NOMINAL_TIME` to the measured value, so a later loss of the
      speed-up shows up red.
    - The suite's time calibration is a single-thread CPU measurement.
      It doesn't account for core count or memory bandwidth, so expect
      the threaded nominal to transfer less well between machines.
  - The suite compares numbers with a 1% tolerance, which is too loose
    for this. So each test prints a SHA-256 of each HDU's `COSBINS`,
    prefixed `sha256:` so the token never parses as a number and must
    match exactly. Each test also prints per-HDU sums, so a failure
    shows how large the difference is.
  - The suite measures memory for the Python process only, and
    `gtltcube` runs as a child process. Add an opt-in module flag (for
    example `MEASURE_CHILDREN = True`) that makes the suite call
    `memory_usage(..., include_children=True)` for that test.
  - `NOMINAL_TIME` and `NOMINAL_MEMORY` come from the gold-standard
    runs.

### Monitoring at every step

Timing is watched at every step of the work order, not just at the
milestones. A change that shouldn't affect speed (step 0, for example)
can still slow the code down, and it is far easier to find at the step
that caused it.

- **Timing set.** All eight `binsz` × `phibins` × `zmax` configurations
  run on one fixed time slice, through both `bench_ltcube` (FT2 read,
  `load()`, write) and `run_ltcube_matrix.py` (full tool). The slice is
  chosen in harness step 1 so that the whole set finishes in a
  reasonable time, roughly under an hour serially. After Stage 1 the
  set also runs at 1, 2, 4 and 8 threads.
- **Full-file runs** happen only at the milestones: the baseline, after
  Stage 1 and after Stage 2.
- **Every step:**
  - rerun the timing set and the likelihood test suite;
  - compare the output against the gold standard;
  - compare the timings against the baseline, and append them to a
    results log kept next to the gold standard: one JSON file per step,
    plus a summary table.
- **Pass criteria.**
  - A step that isn't meant to change speed must stay within the
    run-to-run spread measured at baseline. Use the spread from the 3
    repeats rather than a fixed percentage, since WSL2 timings are
    noisy.
  - With threading added, the `nthreads=1` timings must not be slower
    than baseline by more than that spread.
  - A slowdown outside the spread is investigated before the next step
    starts.

### Work order

1. **Harness.** Write the three harness pieces and the cost check
   described above. Nothing in the libraries changes. The cost check
   measures three things on the unmodified build:
   - the FT2 file's row count;
   - the fixed costs (start-up and the filtered FT2 read for a short
     slice of the ten-year file);
   - the per-row `load()` cost for each configuration.

   From these it picks the timing-set slice and the test-suite slice,
   and says whether the test suite needs a pre-cut FT2 file.
2. **Gold standard and baseline.**
   - On unmodified master in `fermidev`: run the matrix
     (there is no `nthreads` yet) to get the gold-standard cubes and the
     full-tool and `load()` timings.
   - In `FT-2.5.3`: run the same matrix, full-tool only.
   - Compare release against master. We expect them to be identical. If
     they are not, record the differences. Master stays the reference
     for the changed code.
3. **Step 0.** Change `LikeExposure::m_timeCuts` and `m_gtis` to stored
   copies. Compare against the gold standard: the result must be
   identical.
4. **Test-suite entries.** Add `test_ltcube` with its key from the gold
   standard. The threaded variant is added in the next step, once
   `nthreads` exists.
5. **Stage 1.**
   - Add `parallel_fill` with its exception handling, the `map_tools`
     OpenMP link, `Exposure::setNumThreads` (default 1), and the
     `nthreads` parameter (pfile, fhelp).
   - Add the CppUnit test and `test_ltcube_threaded`.
   - Compare against the gold standard at 1, 2, 4 and 8 threads, then
     benchmark.
6. **Decide on Stage 2.** Suggested rule: go ahead if, at 8 threads,
   Stage 1's `load()` speedup is below 4× in any configuration, and
   `load()` is still most of the full tool's time. If FT2 reading or
   FITS writing dominates instead, Stage 2 won't help, and that serial
   part is the next target.
7. **Stage 2.**
   - Add `fill_block` and the `LikeExposure::load` changes.
   - Sweep the block sizes and fix the chosen one as `s_blockSize`.
   - Compare against the gold standard, then benchmark.
8. **Documentation.** Update `../architecture/exposure-and-livetime.md`,
   the known-issues index, `gtltcube.txt` (fhelp) and the release
   notes. Record the `healpix` φ-index bug from
   [Exceptions](#exceptions).

### Still open

- **A representative batch or cluster node** for benchmarking. Not
  available yet. Results from the development laptop alone may be
  misleading: WSL2 timings are noisy, and its 24 MB L3 holds both
  default cubes.
- **arm64 bit-identity check.** It needs a Mac.
- **Which gold-standard cubes to upload.** A full-file gold standard for
  all eight configurations would be about 2 GB; the 0.5° φ cubes are
  about 690 MB each. The test-suite keys are hashes in `key_data/` and
  need no FITS file. If step 2 shows that 2.5.3 and master give
  identical output, anyone can regenerate the gold standard by running
  the harness in a 2.5.3 environment. So the upload might only need the
  pre-cut FT2 slice, if one is needed, plus the small manifests. Decide
  when the list is complete.
- **Thread defaults for other tools.** They stay at one thread for now.
  Review them once we've looked at how each is used.
