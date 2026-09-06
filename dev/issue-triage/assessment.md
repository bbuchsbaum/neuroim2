# GitHub issue triage, 2026-09-06

Reviewed all 14 open issues against live GitHub and default branch
`77b1ddb41a30d251e7751dee84bc81eb93efa07f`. Do not confuse existing unmerged
PR fixes with default-branch behavior.

| Issue | Disposition and evidence |
| --- | --- |
| #1 | Closed as obsolete. The old MNI_SPACE_1MM dataset is no longer shipped. A fresh-process synthetic world-coordinate -> grid -> purrr spherical ROI workflow succeeds. Plotting is now available. |
| #2 | Closed as resolved. Sparse matrix conversion matches an independent dense reference, including zeros outside the mask. Existing regression file passes. Original report omitted its input object; this verifies the currently supported construction/conversion path. |
| #5 | Closed as resolved. Three-timepoint data with a million-timepoint space is rejected. Existing spatial, time, and mask-cardinality validity tests pass. |
| #6, #7 | Closed as obsolete 0.8.0 release checklists, superseded by current 0.19.0 development. This is not certification of unchecked historic CRAN/publicity steps. |
| #19 | Reproduced quantile clipping with a synthetic skull-stripped background: default limits 0..100 clipped bright tissue at 200. Patch uses data limits by default (0..200), retaining explicit robust and numeric limits. No original template was attached; no claim of image-specific visual validation. |
| #20 | Valid documentation request. Details now recommend a shared fixed intensity bandwidth for batch workflows. Clarifies that bilateral_filter_4d pools within each 4-D call, whereas bilateral_filter_vec estimates per volume. No automatic messages or pooled-estimator policy added. |
| #21 | Valid enhancement. Patch adds independent soft-alpha knee/cap/floor and auto-gamma policy controls, exports the curve helper, and records resolved parameters in both return forms. Actual raster alpha bytes are tested. Defaults preserve the existing soft curve. Neurosurf integration is outside this repository; minor argument-renaming and proportional-mode redesign are deferred. |
| #29 | Closed as resolved. Existing file-backed sequence tests pass 18 expectations, including missing backing files during validation and a dense-materialization guard. |
| #30 | Closed as implemented in 77b1ddb. Existing index-only searchlight contract tests pass 403 expectations. |
| #31 | Reproduced the square-matrix convention; patch documents it and adds explicit matrix orientation. Auto retains compatibility; explicit time-by-voxel input round-trips square and rectangular matrices. |
| #32 | Keep open. PR #33 already supplies the reference-index fix and its hosted pkgdown job passed. Default-branch pkgdown still fails; the PR is unmerged. |
| #34 | Keep open. PR #35 already supplies upload support; live Codecov API still reports activated=false, 55.26% coverage, and no latest_commit. Activation and exact default-SHA processing remain unproven. |
| #36 | Keep open. PR #37 already fixes sparse NIfTI writing; five hosted R CMD check platforms passed, but it is unmerged. No duplicate writer patch included here. |

## Reproduction and validation

From the package root:

```sh
LC_ALL=C LANG=C Rscript --vanilla dev/issue-triage/probe.R 77b1ddb
LC_ALL=C LANG=C Rscript --vanilla dev/issue-triage/probe.R
LC_ALL=C LANG=C Rscript --vanilla -e 'devtools::test(reporter="summary", stop_on_failure=TRUE)'
```

R 4.5.1, Matrix 1.7-3, macOS arm64. Full suite: zero failures; 30 new
assertions; three existing skips (empty test and two opt-in performance tests).
Ten warnings: dependency build-version notices, two parallel-worker CPU-budget
warnings, a k-means Quick-TRANSfer warning, and a multiple-loaded-DLL warning.
No tests were weakened. Synthetic rendering was also inspected; the pixel-level
alpha oracle, rather than visual appearance alone, establishes the curve values.

Build/check and hosted verification are separate. The local package check omits
vignette rebuilding, test reruns (full suite run separately), and the PDF manual.
The existing PRs also have docs-check/AppVeyor failures; their successful R checks
are not a claim that all gates are green. No PR was merged.

Local tarball check: zero errors, two warnings, zero notes. Both warnings are
missing built vignette artifacts caused by `--no-build-vignettes`; installation,
namespace, compiled code, Rd consistency, and examples passed. The 30 new
assertions also pass against the package installed from this tarball.

```sh
R CMD build . --no-build-vignettes --no-manual
R CMD check --no-manual --no-vignettes --no-tests neuroim2_0.19.0.tar.gz
```

Tarball SHA256: `1207a820e909988bd62afc0aad7f0f92222d0f189a291418420da5cf7de3f18e`.
Local logs: `/tmp/neuroim2-triage/` (temporary, not a hosted evidence source).
