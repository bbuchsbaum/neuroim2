# CRAN release candidate 0.19.1

The preparation task produces a checked candidate. The coordinating parent
thread owns any separately authorized CRAN upload; this task does not submit,
merge, tag, release or deploy.

## Scope and version

Initial base: `f1f7f00264b0566beb73d17654ff1e874e66bf3c` on remote master.
PR #46 and the separate website-theme PR #45 were subsequently merged by
another actor. The follow-up candidate in PR #47 starts from current master
`f11f6469b32e9143e4a605c8aa1aa5e47ac9a80a`; this task performed no merge.
CRAN's package page reports 0.13.0 (published 2026-04-16). Its current check
snapshot is dated 2026-10-03 and is not evidence for this candidate.

0.19.1 advances the existing master development version 0.19.0.9000 and keeps
0.20.0 available for the separate lazy-iteration feature branch. It does not
imply a small change from CRAN 0.13.0: NEWS retains the intermediate history
and adds an upgrade guide. The coordinating parent owns the separately
authorized submission of 0.19.1 after final validation and independent review.

The local `feat/plot-hillclimb` commit
`b7b874c4aee9b656ae8bf333f3c0a13fc7274f9b` also contains `vec_blocks()` and
accessor expansion. Only the reproduced mapped-scaling and square-sequence
correctness fixes, their generated help and regression fixtures are selectively
backported here. No new iterator is exported. The original branch is unchanged.
Website-theme PR #45 was kept separate during release preparation; its changes
are now in the base because it was merged independently.

Additional isolated reproductions exposed the same square-matrix ambiguity
in sparse downsampling, scaling, arithmetic, concatenation and ROI conversion.
Those result producers now state their known matrix orientation explicitly;
the public constructor still warns when callers supply ambiguous square data.

Independent review also reproduced a baseline dimension-dropping bug in both
dense-to-sparse conversion methods. Matrix subsetting now preserves dimensions
for a single selected voxel or a single time point, with numeric indices and
LogicalNeuroVol masks. This changes two subsetting expressions and retains the
existing matrix-orientation contract. Follow-up independent review found that
unsorted numeric masks could assign values to the wrong voxels. Numeric
conversion now derives both data order and support from the same logical
selection. Unsorted and repeated indices, negative exclusions, zeros and empty
selections are tested against independent expected voxel values and equivalent
logical masks. Invalid indices are rejected instead of silently disappearing.

## Reproduce checks and submission artifact

Use a clean checkout of the reviewed candidate SHA with current release R
for the source build and current R-devel for the tarball check,
R package dependencies, Pandoc, a working LaTeX installation and at least 20 GiB
free space on the connected Mac. From the checkout root:

Use roxygen2 7.3.3, matching DESCRIPTION, for regeneration; newer roxygen2
releases change unrelated generated output and are not part of this release.
The LaTeX installation must include Courier, Helvetica and Times metrics
and MakeIndex (`tlmgr install courier helvetic times makeindex` for TinyTeX).
The workflow checks the
manual early so missing TeX dependencies fail before the long vignette build.

```sh
export RCPP_PARALLEL_NUM_THREADS=2 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
export MAKEFLAGS=-j2 _R_CHECK_LIMIT_CORES_=true
export RELEASE_OUT="$(mktemp -d)"
Rscript -e 'roxygen2::roxygenize(".", roclets=c("rd", "namespace"), load_code=roxygen2::load_installed)'
git diff --exit-code -- man NAMESPACE
R CMD Rd2pdf --no-preview --force --output="$RELEASE_OUT/neuroim2-manual.pdf" .
Rscript tools/release/baseline.R
Rscript tools/release/build-source.R
# Switch to a separate R-devel installation and dependency library for this step.
Rscript tools/release/preflight.R
# Use current release R for downstream comparisons.
export DOWNSTREAM_OUT="$(mktemp -d)"
Rscript tools/release/downstream.R
```

Install this exact checkout into an isolated R library before the documentation
and before/after checks (for example `R CMD INSTALL --library="$R_LIBS_USER" .`
with `R_LIBS_USER` pointing to a newly created temporary directory). The workflow
does this via `setup-r-dependencies`'s `local::.` entry. Do not use a different
installed neuroim2 version for those checks.

`build-source.R` uses `pkgbuild::build(vignettes=TRUE, manual=TRUE)` with release
R. `preflight.R` checks the recorded source SHA and tarball checksum, then runs
`rcmdcheck::rcmdcheck(tarball, args="--as-cran")` with R-devel. Separate CI jobs
keep their R libraries isolated. The uploadable artifact is
`neuroim2_0.19.1.tar.gz` in RELEASE_OUT; `SOURCE_SHA`, `SHA256SUMS`, the complete
check directory, both build/check session records and parsed results identify what was
checked. The workflow retains these as `release-source-<head SHA>` for 30 days.
Rebuild and rerun checks if that artifact expires or the source changes.

The nine HTML vignettes share their unchanged local font files rather than
embedding the same fonts nine times. Plots, scripts and layout CSS remain
embedded. `vignettes/.install_extras` packages the shared font stylesheet,
all seven fonts and their license. A scratch two-vignette render verified
identical scripts, images and non-font CSS and distinct embedded figures.
The preflight also validates both the built tarball's and installed package's
documentation: every resource resolves offline, shared assets match source
bytes, and the complete documentation directory is below 5 MB. It records
these checks and measured sizes in `vignette-assets.csv`. This is a package
size fix independent of the separately merged theme work; article content is
unchanged.

The focused before/after script requires all seven independent regression cases
to fail on an isolated build of the base SHA and pass on the candidate. The
singleton-conversion probe separately exercises 32 numeric/logical-mask cases:
14 singleton cases fail on master and 18 controls pass; all 32 must pass on
the candidate. Its logs are retained alongside the other regression evidence.
The exact public-writer reproduction from issue #40 also runs against the
candidate with its SHORT assertion corrected to require dense/mapped agreement;
FLOAT remains an identity control. The test suite additionally checks the
zero-slope convention. These tests concern uncompressed native-endian files
and do not establish behavior for every compressed or endian-conversion path.
The full suite exercises dense, sparse, BigNeuroVec, mapped, file-backed and
sequence backends, FLOAT/SHORT/UBYTE decoding, affine/NIfTI geometry,
smoothing, resampling and plot structure. Visual reference snapshots retain
the existing opt-in environment guard; a green check does not claim those ran.

## Downstream coverage and limits

`downstream.R` saves the current CRAN package index and downloaded sources,
discovers direct Depends/Imports/LinkingTo/Suggests reverse dependencies and
checks them with isolated baseline and candidate libraries. It requires bidser
to appear, retains existing findings, and fails on newly observed errors or
warnings. These checks omit the downstream PDF manuals.
The comparison removes elapsed times from diagnostic headers only; complete
raw diagnostics remain in the artifacts. Rechecking a current CRAN package
can produce an existing incoming-feasibility warning in both versions.

It also runs the `nifti-array-source` and `feature-space` tests from fmridataset
`ef141af6623b699168fa4329f304992c8b782ed0`. This is focused lab integration,
not a full lab-wide check. Other lab packages, recursive reverse dependencies,
reference-machine visual snapshots, sanitizers and external CRAN test services
remain unverified. If a downstream API break is found, resolve it or obtain
maintainer coordination before submission; no communications are sent here.

## Evidence status

The authoritative completion record is the exact-head results and artifact
links in [draft PR #47](https://github.com/bbuchsbaum/neuroim2/pull/47). Match
its head SHA to each workflow and the artifact's `SOURCE_SHA`; results from a
superseded commit do not establish the final candidate's status.

The initial full local check was blocked by the disk-space guard (18 GiB free,
required minimum 20 GiB); full release checks run on GitHub. Local validation
includes a 279-topic Rd audit and a successful 257-page PDF manual. After free
space recovered above the guard, the singleton fix was installed into an
isolated local library for focused regression and independent runtime probes.
Small mixed-source reproductions provided
initial evidence; the workflow's isolated source installs establish the seven
before/after package comparisons. Coverage, test counts and any check notes
belong in the PR evidence record after the final workflows finish.

Release decision: no-go for submission until the final candidate's workflow
results and independent review findings are reconciled. Earlier passing
artifacts are superseded whenever runtime code changes. Coverage against the
internal 90% target is reported separately; it is not a CRAN policy gate.

Policy references checked during preparation:

* <https://cran.r-project.org/web/packages/policies.html>
* <https://cran.r-project.org/web/packages/submission_checklist.html>
* <https://cran.r-project.org/package=neuroim2>
* <https://cran.r-project.org/web/checks/check_results_neuroim2.html>
