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

neuroim2 reads, manipulates and writes volumetric brain images in R. Its
objects behave like ordinary arrays – you can index them, do arithmetic
on them, compare them – but they never forget where they are in space:
what comes back from an operation is still an image, with its affine
attached.

It covers the ground between a file on disk and a statistical model:
NIfTI and AFNI I/O, 3D volumes and 4D time-series in dense, sparse and
memory-mapped forms, coordinate systems, regions of interest and
searchlights, resampling, spatial filtering, and slice plots. It is the
successor to [`neuroim`](https://github.com/bbuchsbaum/neuroim).

**Website:** <https://bbuchsbaum.github.io/neuroim2/>  ·  **CRAN:**
<https://cran.r-project.org/package=neuroim2>

## Installation

``` r
install.packages("neuroim2")
```

Daily builds are on R-universe, and the development version is on
GitHub:

``` r
install.packages("neuroim2",
  repos = c("https://bbuchsbaum.r-universe.dev", "https://cloud.r-project.org"))

# install.packages("remotes")
remotes::install_github("bbuchsbaum/neuroim2")
```

## Five minutes with neuroim2

Everything below runs as written; the two images ship with the package.

### Read an image and look at it

``` r
library(neuroim2)

anat <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
plot_montage(anat, zlevels = seq(-30, 60, by = 18), unit = "mm", ncol = 6, style = "dark")
```

<img src="man/figures/README-anatomy-1.png" alt="Six axial slices of the MNI152 template, from the cerebellum to the vertex." width="100%" />

### An array that knows where it is

A volume indexes and computes like an array, and carries its geometry
with it:

``` r
dim(anat)
#> [1] 48 57 48
spacing(anat)                      # voxel size, mm
#> [1] 4.02083 4.01754 4.02083
affine_to_axcodes(trans(anat))     # which way the grid axes run
#> [1] "R" "A" "S"

class(anat > 5000)[1]              # a comparison is still an image
#> [1] "LogicalNeuroVol"
```

That geometry is what turns a millimetre coordinate into a voxel, and a
voxel back into millimetres. Voxel indices are 1-based, as everywhere in
R.

``` r
round(coord_to_grid(anat, c(0, 0, 0)), 1)      # the voxel at the MNI origin
#> [1] 24.9 33.9 20.4
round(index_to_coord(anat, which.max(anat)))   # where the brightest voxel is, mm
#>      [,1] [,2] [,3]
#> [1,]  -44   -3  -30
```

### From a time-series to a map

A `NeuroVec` is a 4D image: a time-series at every voxel. With real data
you would call `read_vec("bold.nii.gz")`; here `simulate_fmri()` makes a
spatially smooth one on the template grid so the example is
self-contained.

``` r
brain <- as.mask(anat > 3000)
bold  <- simulate_fmri(brain, n_time = 80, spatial_fwhm = 14, seed = 1)
dim(bold)
#> [1] 48 57 48 80
```

A seed-based correlation map is four steps: put a sphere at a millimetre
coordinate, average the time-series inside it, correlate that with every
voxel in the brain, and put the result back on the grid.

``` r
centre  <- round(coord_to_grid(anat, c(-40, -24, 52)))
seed    <- spherical_roi(brain, centre, radius = 8)        # radius in mm
seed_ts <- rowMeans(values(series_roi(bold, seed)))        # scans x voxels -> scans

vox   <- which(brain)
r     <- cor(series(bold, vox), seed_ts)                   # one r per brain voxel
r_map <- NeuroVol(as.vector(r), space(anat), indices = vox)

round(index_to_coord(r_map, which.max(r_map)))   # the peak, back in mm
#>      [,1] [,2] [,3]
#> [1,]  -40  -24   51
```

``` r
plot_overlay(anat, r_map, zlevels = c(40, 46, 52, 58, 64), unit = "mm", ncol = 5,
             ov_thresh = 0.5, cbar_title = "r", style = "dark")
```

<img src="man/figures/README-overlay-1.png" alt="Five axial slices of the template with the thresholded correlation map overlaid around the seed in the left hemisphere." width="100%" />

The data are simulated, so the map shows only what was put in – voxels
near the seed share its signal because the simulation is spatially
smooth. The steps are the ones a real analysis uses.

### Write it out

``` r
out <- tempfile(fileext = ".nii.gz")
write_vol(r_map, out)

all.equal(trans(read_vol(out)), trans(r_map))    # geometry survives the round trip
#> [1] TRUE
```

## What is in the package

- **Reading and writing** – NIfTI (`.nii`, `.nii.gz`) and AFNI
  (`.HEAD`/`.BRIK`), with header inspection before you load:
  `read_vol()`, `read_vec()`, `read_header()`, `write_vol()`,
  `write_vec()`
- **Containers** – `NeuroVol` (3D) and `NeuroVec` (4D), each in dense,
  sparse, memory-mapped, file-backed and clustered (parcel) forms that
  share one interface
- **Coordinates** – conversion between linear index, grid index and
  world millimetres, orientation codes, reorientation and deobliquing
- **Regions and searchlights** – spherical, cuboid and parcel ROIs;
  exhaustive, random and clustered searchlights, in parallel if you want
- **Resampling** – onto another image’s grid, to a new voxel size, or to
  a new orientation: `resample_to()`, `downsample()`, `reorient()`
- **Filtering** – Gaussian, bilateral and guided filters, Laplacian
  sharpening, connected components: `gaussian_blur()`,
  `bilateral_filter()`, `conn_comp()`
- **Plotting** – slice montages, orthogonal views, statistical overlays
  and registration checks, built on ggplot2: `plot_montage()`,
  `plot_ortho()`, `plot_overlay()`

## Coming from another package

neuroim2 exchanges images with RNifti and oro.nifti objects directly, so
it can sit in the middle of a pipeline that uses RNiftyReg, neurobase or
fslr without writing files in between:

``` r
img <- RNifti::readNifti(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))

vol  <- as_neuro(img)     # RNifti or oro.nifti image -> NeuroVol / NeuroVec
img2 <- as_nifti(vol)     # and back, with the affine

all.equal(trans(vol), trans(anat))
#> [1] TRUE
```

If you know nibabel or nilearn, most of what you reach for has a direct
counterpart:

- `nib.load(f)` → `read_vol(f)`, `read_vec(f)`, or `read_image(f)` to
  let the header decide
- `img.header` → `read_header(f)`, `header(x)`
- `img.get_fdata()` → `as.array(x)`
- `img.affine` → `trans(x)`
- `img.header.get_zooms()` → `spacing(x)`
- `nib.aff2axcodes(img.affine)` → `affine_to_axcodes(trans(x))`
- `nib.as_closest_canonical(img)` → `as_canonical(x)`
- `apply_affine(img.affine, ijk)` → `grid_to_coord(x, ijk + 1)`
- `resample_to_img(src, target)` → `resample_to(src, target)`
- `nib.save(img, f)` → `write_vol(x, f)`, `write_vec(x, f)`

The one habit to change is in that `+ 1`: grid indices are **1-based**
here. Voxel `(1, 1, 1)` sits at the origin of the affine, where nibabel
has `(0, 0, 0)`. World coordinates are RAS+ millimetres in both.

## Learn more

The
[articles](https://bbuchsbaum.github.io/neuroim2/articles/index.html) on
the website are written to be read in order, then by task.

**Learn** – about 45 minutes end to end

- [A tour](https://bbuchsbaum.github.io/neuroim2/articles/neuroim2.html)
  – the package in one sitting
- [Spaces and
  coordinates](https://bbuchsbaum.github.io/neuroim2/articles/spaces-and-coordinates.html)
  – the affine, orientation codes, and voxel-to-millimetre conversion
- [Volumes and
  vectors](https://bbuchsbaum.github.io/neuroim2/articles/volumes-and-vectors.html)
  – building, slicing and combining the containers
- [Reading and
  writing](https://bbuchsbaum.github.io/neuroim2/articles/reading-and-writing.html)
  – file formats, headers, and working with other packages

**Do** – read the one you need

- [Regions and
  searchlights](https://bbuchsbaum.github.io/neuroim2/articles/regions-and-searchlights.html)
- [Resampling and
  orientation](https://bbuchsbaum.github.io/neuroim2/articles/resampling-and-orientation.html)
- [Smoothing and
  filtering](https://bbuchsbaum.github.io/neuroim2/articles/smoothing-and-filtering.html)
- [Visualization](https://bbuchsbaum.github.io/neuroim2/articles/visualization.html)

**Scale** – when the data stops fitting in memory

- [Large
  data](https://bbuchsbaum.github.io/neuroim2/articles/large-data.html)

Every function is listed in the
[reference](https://bbuchsbaum.github.io/neuroim2/reference/index.html).
Bugs and questions go to the [issue
tracker](https://github.com/bbuchsbaum/neuroim2/issues).

<!-- albersdown:theme-note:start -->

## Albers theme

This package uses the albersdown theme. Existing vignette theme hooks
are replaced so `albers.css` and local `albers.js` render consistently
on CRAN and GitHub Pages. The defaults are configured via
`params$family` and `params$preset` (family = `red`, preset = `homage`).
The pkgdown site uses `template: { package: albersdown }`.
<!-- albersdown:theme-note:end -->
