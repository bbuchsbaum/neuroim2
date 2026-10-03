# CRAN release candidate 0.19.1

Review candidate only; do not submit, merge, tag, release or deploy as part of
this preparation.

## Scope and version

Base: `f1f7f00264b0566beb73d17654ff1e874e66bf3c` on remote master.
CRAN's package page reports 0.13.0 (published 2026-04-16); its retrieved check
snapshot is dated 2026-07-16 and is not evidence for this candidate.

0.19.1 advances the existing master development version 0.19.0.9000 and keeps
0.20.0 available for the separate lazy-iteration feature branch. It does not
imply a small change from CRAN 0.13.0: NEWS retains the intermediate history
and adds an upgrade guide. Bradley should confirm this provisional version
before submission.

The local `feat/plot-hillclimb` commit
`b7b874c4aee9b656ae8bf333f3c0a13fc7274f9b` also contains `vec_blocks()` and
accessor expansion. Only the reproduced mapped-scaling and square-sequence
correctness fixes, their generated help and regression fixtures are selectively
backported here. No new iterator is exported. The original branch is unchanged.
Website-theme PR #45 is separate and has not been merged.

## Reproduce checks and submission artifact

Use a clean checkout of the reviewed candidate SHA with current R-devel,
R package dependencies, Pandoc, a working LaTeX installation and at least 20 GiB
free space on the connected Mac. From the checkout root:

Use roxygen2 7.3.3, matching DESCRIPTION, for regeneration; newer roxygen2
releases change unrelated generated output and are not part of this release.

```sh
export RCPP_PARALLEL_NUM_THREADS=2 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
export MAKEFLAGS=-j2 _R_CHECK_LIMIT_CORES_=true
export RELEASE_OUT="$(mktemp -d)"
Rscript -e 'roxygen2::roxygenize(".", roclets=c("rd", "namespace"), load_code=roxygen2::load_installed)'
git diff --exit-code -- man NAMESPACE
Rscript tools/release/baseline.R
Rscript tools/release/preflight.R
export DOWNSTREAM_OUT="$(mktemp -d)"
Rscript tools/release/downstream.R
```

Install this exact checkout into an isolated R library before the documentation
and before/after checks (for example `R CMD INSTALL --library="$R_LIBS_USER" .`
with `R_LIBS_USER` pointing to a newly created temporary directory). The workflow
does this via `setup-r-dependencies`'s `local::.` entry. Do not use a different
installed neuroim2 version for those checks.

`preflight.R` uses `pkgbuild::build(vignettes=TRUE, manual=TRUE)` and then
`rcmdcheck::rcmdcheck(tarball, args="--as-cran")`. The uploadable artifact is
`neuroim2_0.19.1.tar.gz` in RELEASE_OUT; `SOURCE_SHA`, `SHA256SUMS`, the complete
check directory, session information and parsed results identify what was
checked. The workflow retains these as `release-source-<head SHA>` for 30 days.
Rebuild and rerun checks if that artifact expires or the source changes.

The focused before/after script requires both independent regression fixtures
to fail on an isolated build of the base SHA and pass on the candidate. The
full suite exercises dense, sparse, BigNeuroVec, mapped, file-backed and
sequence backends, FLOAT/SHORT/UBYTE decoding, affine/NIfTI geometry,
smoothing, resampling and plot structure. Visual reference snapshots retain
the existing opt-in environment guard; a green check does not claim those ran.

## Downstream coverage and limits

`downstream.R` saves the current CRAN package index and downloaded sources,
discovers direct Depends/Imports/LinkingTo/Suggests reverse dependencies and
checks them with isolated baseline and candidate libraries. It requires bidser
to appear, retains existing findings, and fails on newly observed errors or
warnings. These checks omit the downstream PDF manuals.

It also runs the `nifti-array-source` and `feature-space` tests from fmridataset
`ef141af6623b699168fa4329f304992c8b782ed0`. This is focused lab integration,
not a full lab-wide check. Other lab packages, recursive reverse dependencies,
reference-machine visual snapshots, sanitizers and external CRAN test services
remain unverified. If a downstream API break is found, resolve it or obtain
maintainer coordination before submission; no communications are sent here.

## Evidence status

Pending exact-head workflow completion. Full local check is blocked by the
disk-space guard (18 GiB free at the start, required minimum 20 GiB). Preliminary
small mixed-source reproductions confirmed the two defects, but only the
workflow's isolated builds establish the before/after package evidence.

Release decision: no-go for submission until workflow results and any findings
are reviewed and the version is confirmed. Draft PR review can proceed.

Policy references checked during preparation:

* <https://cran.r-project.org/web/packages/policies.html>
* <https://cran.r-project.org/web/packages/submission_checklist.html>
* <https://cran.r-project.org/package=neuroim2>
* <https://cran.r-project.org/web/checks/check_results_neuroim2.html>
