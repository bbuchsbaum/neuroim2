## Release candidate 0.19.1 (not submitted)

This candidate updates CRAN 0.13.0 with the changes recorded in NEWS through
0.19.1. NEWS includes migration guidance for result-changing smoothing,
reorientation, local-maxima, mapped scaling, sequence and sparse-result
orientation, and plotting behavior. Maintainer and license are unchanged.

The three help topics reported in the 2026-07-16 CRAN check snapshot now have
usage sections generated from their roxygen sources: image, as.raster and
as-ClusteredNeuroVol-DenseNeuroVol.

## Validation status

Evidence is recorded in dev/cran-release-0.19.1.md and the draft PR's exact-head
workflow artifacts. Replace this status with final observed results before
any submission; the obsolete 0.8.5 check claim has been removed.

Full local R CMD check has not run: the connected Mac has less than the 20 GiB
free space required by its repository instructions. The release-preflight
workflow builds vignettes and the PDF manual and runs R CMD check --as-cran on
the retained source tarball using Linux R-devel. The regular CI matrix also
checks macOS, Windows and Linux.

The downstream workflow compares current CRAN reverse dependencies against
the CRAN baseline and candidate and runs pinned fmridataset integration tests.
Winbuilder, R-hub, macbuilder and CRAN submission have not been requested.
