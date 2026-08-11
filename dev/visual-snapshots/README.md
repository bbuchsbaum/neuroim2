# Visual snapshot harness

`test-plot-snapshots.R` compares rendered plots against the golden SVGs in
`_snaps_plot-snapshots/` using vdiffr.

It lives here rather than in `tests/testthat/` because these comparisons are
byte-exact on SVG text, which encodes the font metrics of the machine that
produced them. They cannot pass across the R-CMD-check matrix (Windows, macOS
and Linux), and testthat deletes any `_snaps/` subdirectory whose test file did
not run, so they cannot simply be skipped either.

The structural equivalents in `tests/testthat/test-plot-structure.R` cover the
same plotting calls in a platform-independent way and run on every check.

To use this harness when reviewing an intentional visual change, copy both
entries back into `tests/testthat/` (the snapshot directory must be named
`_snaps/plot-snapshots`), run the suite on the reference machine, review with
`testthat::snapshot_review()`, then move them back here.
