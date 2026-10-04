## Release candidate 0.19.1 (resubmission preparation)

The previous submission 357569 was archived after pretest on 2026-10-03.
Although Windows and Debian R CMD check both ended Status: OK, the incoming
summary flagged Windows overall runtime of 37 minutes above ten minutes.
Vignette rebuilding accounted for 31 minutes on Windows and 19 on Debian.

Profiling identified repeated extraction of the complete S4 data array in
dense vectors() iteration. The iterator now extracts it once; voxel order and
values are unchanged. All vignette examples, full datasets and tests remain
enabled. Release checks now enforce a ten-minute complete-check budget and
explicit review of every NOTE, including official win-builder R-devel and
R-hub v2 clang/UBSan checks of the exact release-built tarball.

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

During initial preparation, full local checking was blocked by the repository's
20 GiB disk-space guard. During runtime remediation, a local linker crash was
isolated to the toolchain and worked around with a per-command Apple clang and
classic-linker configuration; Bradley's global configuration was unchanged.
The release-preflight
workflow builds the source tarball, vignettes and PDF manual with release R,
then checks that exact tarball using Linux R-devel. The regular CI matrix also
checks macOS, Windows and Linux.

The downstream workflow compares current CRAN reverse dependencies against
the CRAN baseline and candidate and runs pinned fmridataset integration tests.
Final win-builder and R-hub results must be attached to the exact new artifact;
the first submission's passing status and excessive runtime are superseded.
Macbuilder has not been run. The coordinating parent
owns the separately authorized CRAN submission after final validation and
independent review. Use the final handoff's submission comments and exact
checked tarball; do not reuse results from a superseded candidate.
