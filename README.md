# neuroim2 <img src="https://raw.githubusercontent.com/bbuchsbaum/neuroim2/master/docs/hex-neuroim2.png" align="right" width="120" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/neuroim2)](https://cran.r-project.org/package=neuroim2)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/neuroim2)](https://cran.r-project.org/package=neuroim2)
[![R-CMD-check](https://github.com/bbuchsbaum/neuroim2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/neuroim2/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/bbuchsbaum/neuroim2/actions/workflows/pkgdown.yaml/badge.svg)](https://bbuchsbaum.github.io/neuroim2/)
[![codecov](https://codecov.io/gh/bbuchsbaum/neuroim2/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bbuchsbaum/neuroim2)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
<!-- badges: end -->

Data structures and I/O for volumetric brain imaging with a focus on
fMRI. This is the successor to
[`neuroim`](https://github.com/bbuchsbaum/neuroim) with improved S4
classes, sparse/dense 3D–4D representations, and fast
resampling/filtering.

**Website:** <https://bbuchsbaum.github.io/neuroim2/>  
**CRAN:** <https://cran.r-project.org/package=neuroim2>  
**Cheatsheet:** [neuroim2_cheatsheet.md](neuroim2_cheatsheet.md)

## Installation

### CRAN

``` r
install.packages("neuroim2")
```

### R-universe (daily builds)

``` r
install.packages("neuroim2",
  repos = c("https://bbuchsbaum.r-universe.dev", "https://cloud.r-project.org"))
```

### Development version (GitHub)

``` r
# install.packages("remotes")
remotes::install_github("bbuchsbaum/neuroim2")
```

## Usage

Read a NIFTI image and perform simple operations:

``` r
library(neuroim2)
```

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'neuroim2'

    ## The following object is masked from 'package:base':
    ## 
    ##     scale

``` r
fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
vol <- read_vol(fname)

v1 <- vol[1,1,1]
vol2 <- vol + vol
all(vol == (vol2 - vol))
```

    ## [1] TRUE

Create a 4D image from volumes:

``` r
vec <- vec_from_vols(list(vol, vol, vol2))
series1 <- vec[1,1,1,]
length(series1)
```

    ## [1] 3

## Vignettes

See examples of use of `neuroim2` in the
[vignettes](https://bbuchsbaum.github.io/neuroim2/articles/index.html).


## Albers theme
This package uses the albersdown theme. Vignettes are styled with `vignettes/albers.css` and a local `vignettes/albers.js`; the palette family is provided via `params$family` (default 'red'). The pkgdown site uses `template: { package: albersdown }`.

<!-- albersdown:theme-note:start -->
## Albers theme
This package uses the albersdown theme. Existing vignette theme hooks are replaced so `albers.css` and local `albers.js` render consistently on CRAN and GitHub Pages. The palette family is provided via `params$family` (default 'red'). The pkgdown site uses `template: { package: albersdown }`.
<!-- albersdown:theme-note:end -->
