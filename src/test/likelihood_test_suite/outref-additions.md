# Files to add to `fermitools-test-scripts-data`

This lists the data files we create for Likelihood tests and benchmarks
that belong in the `fermitools-test-scripts-data` package. Until they
are uploaded, they exist only locally in
`$CONDA_PREFIX/share/fermitools/data/Likelihood/outref/`, which is where
the package installs them.

The package's large files are hosted on a web server and downloaded when
the conda package is built (`likelihood_url` and `likelihood_files` in
the package's `recipe/build.sh`). To publish these files:

1. Upload the files to the web host.
2. Add each name to `likelihood_files` in `recipe/build.sh` and to
   `payload.list`, and add its SHA-256 to `MANIFEST.sha256`.
3. Rebuild and release the package.

Add a row for each file as soon as it is created. "Status" is `planned`,
`local` (created, not uploaded) or `uploaded`.

| File | Status | Made by | Inputs and parameters | Size | SHA-256 | Used by |
|------|--------|---------|-----------------------|------|---------|---------|
| *(none yet)* | | | | | | |

## Expected additions

From [ltcube pixel threading](../../../doc/proposals/ltcube-pixel-threading.md#implementation-plan):

- **Gold-standard ltcubes** from the unmodified `gtltcube`, made from
  `ft2_10years_reduced.fits`, one for each `binsz` × `phibins` × `zmax`
  configuration. Each comes with a manifest giving the exact parameters,
  the time slice, the tool version and the input's SHA-256. Which of
  these get uploaded is still open; see "Still open" in the proposal.
- **A pre-cut FT2 slice** of `ft2_10years_reduced.fits`, only if the
  harness cost check shows that reading the whole ten-year file hides
  the speed-up in the test-suite tests.
