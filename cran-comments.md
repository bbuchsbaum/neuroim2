## Release candidate 0.19.1 (not submitted)

This candidate updates CRAN 0.13.0 with the changes recorded in NEWS through
0.19.1. NEWS includes migration guidance for result-changing smoothing,
reorientation, local-maxima, mapped scaling, sequence and sparse-result
orientation, and plotting behavior. Dense-to-sparse conversions also preserve
singleton voxel/time dimensions and keep values at their original voxels for
unsorted or repeated numeric masks. Maintainer and license are unchanged.

The three help topics reported in the 2026-10-03 CRAN check snapshot now have
usage sections generated from their roxygen sources: image, as.raster and
as-ClusteredNeuroVol-DenseNeuroVol.

The nine vignettes share bundled local fonts to keep installed documentation
below CRAN's general 5 MB limit. All content and styling are retained; the
preflight checks packaged and installed resources for offline resolution.

## Validation status

Evidence is recorded in dev/cran-release-0.19.1.md and the draft PR's exact-head
workflow artifacts. Replace this status with final observed results before
any submission; the obsolete 0.8.5 check claim has been removed.

Full local R CMD check has not run: the initial preflight was blocked by the
repository's 20 GiB disk-space guard. After space recovered, the singleton
fix was verified with an isolated local install and focused tests. The release-preflight
workflow builds the source tarball, vignettes and PDF manual with release R,
then checks that exact tarball using Linux R-devel. The regular CI matrix also
checks macOS, Windows and Linux.

The downstream workflow compares current CRAN reverse dependencies against
the CRAN baseline and candidate and runs pinned fmridataset integration tests.
Winbuilder, R-hub and macbuilder have not been run. The coordinating parent
owns the separately authorized CRAN submission after final validation and
independent review. Use the final handoff's submission comments and exact
checked tarball; do not reuse results from a superseded candidate.
