# Overview: Getting Started with neuroim2

## Introduction

`neuroim2` gives you a small set of data structures for 3D and 4D
neuroimaging data, plus the spatial tools you need to move between file
I/O, coordinate systems, regions of interest, and resampling. The
package is broad, so this overview is intentionally narrow: it shows the
first objects and workflows to learn, then points you to the focused
vignettes that carry the rest.

## Quick start

Start by reading one image and inspecting its spatial metadata.

``` r
img <- read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))

dim(img)
#> [1] 64 64 25
spacing(img)
#> [1] 3.5 3.5 3.7
origin(img)
#> [1]  112.00 -108.50  -46.25
```

The most important thing to notice is that a `NeuroVol` is not just an
array. It also carries a `NeuroSpace`, which tracks voxel spacing,
origin, and affine transforms.

## What should you read next?

The recommended path through the package is:

1.  [`vignette("ChoosingBackends", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/ChoosingBackends.md)
    for dense, sparse, mapped, file-backed, and hyper-vector backends.
2.  [`vignette("coordinate-systems", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/coordinate-systems.md)
    for voxel, grid, and world-coordinate conversions.
3.  [`vignette("VolumesAndVectors", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/VolumesAndVectors.md)
    for the core manipulation story.
4.  [`vignette("Resampling", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/Resampling.md)
    for
    [`resample()`](https://bbuchsbaum.github.io/neuroim2/reference/resample-methods.md),
    [`downsample()`](https://bbuchsbaum.github.io/neuroim2/reference/downsample-methods.md),
    [`reorient()`](https://bbuchsbaum.github.io/neuroim2/reference/reorient-methods.md),
    and
    [`deoblique()`](https://bbuchsbaum.github.io/neuroim2/reference/deoblique.md).
5.  [`vignette("AnalysisWorkflows", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/AnalysisWorkflows.md)
    for ROIs, searchlights, and map-reduce style analyses.

If you only read one follow-on article after this overview, make it
[`vignette("VolumesAndVectors", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/VolumesAndVectors.md).

## The core objects

Most work in `neuroim2` starts with three ideas:

- `NeuroVol` for 3D images such as anatomical volumes, masks, and single
  summary maps.
- `NeuroVec` for 4D data such as fMRI time-series or stacked volumes.
- `ROI` objects for region-based extraction and local analyses.

Here is the smallest possible example of each.

``` r
mask <- img > 0
sum(mask)
#> [1] 29532

vec <- read_vec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
dim(vec)
#> [1] 64 64 25  4

roi <- spherical_roi(space(vec), c(45, 45, 20), radius = 4)
length(roi)
#> [1] 7
```

That is the core mental model for the package:

- read or construct a spatial object
- operate in image or vector form as needed
- define ROIs or neighborhoods
- extract, transform, or summarize

## A small end-to-end workflow

The next common step is to move from a 4D image to a region-level
summary.

``` r
roi_ts <- series_roi(vec, roi)
roi_mat <- values(roi_ts)
mean_ts <- rowMeans(roi_mat)

stopifnot(
  nrow(roi_mat) == dim(vec)[4],
  ncol(roi_mat) == length(roi),
  all(is.finite(mean_ts))
)

head(mean_ts)
#> [1] 0 0 0 0
```

This is a deliberately small example, but it shows the typical
`neuroim2` workflow:

1.  Load a spatial object.
2.  Define a spatial support such as an ROI.
3.  Extract values with the correct geometry preserved.
4.  Compute summaries at the level you care about.

For broader ROI and searchlight patterns, move directly to
[`vignette("AnalysisWorkflows", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/AnalysisWorkflows.md).

## Spatial operations come next

Once you are comfortable reading data and extracting values, the next
important layer is spatial transformation.

``` r
img_down <- downsample(img, spacing = c(2, 2, 2))

dim(img)
#> [1] 64 64 25
dim(img_down)
#> [1] 112 112  46
spacing(img_down)
#> [1] 2.00000 2.00000 2.01087
```

For the full story, including orientation handling and affine-aware
transforms, use:

- [`vignette("coordinate-systems", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/coordinate-systems.md)
- [`vignette("Resampling", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/Resampling.md)

## When should you change backends?

You do not need a special backend to start. Use the default dense path
first, then switch when the workload demands it.

- Use dense objects when the data fits comfortably in memory.
- Use sparse objects when most voxels are absent and should be treated
  as missing support, not stored zeros.
- Use file-backed or mapped objects when the array is too large to
  materialize eagerly.

``` r
big_vec <- read_vec(
  system.file("extdata", "global_mask_v4.nii", package = "neuroim2"),
  mode = "filebacked"
)

series(big_vec, 45, 45, 20)
```

The details and tradeoffs belong in
[`vignette("ChoosingBackends", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/ChoosingBackends.md).

## Where to go next

### Core path

- [`vignette("ChoosingBackends", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/ChoosingBackends.md)
- [`vignette("coordinate-systems", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/coordinate-systems.md)
- [`vignette("VolumesAndVectors", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/VolumesAndVectors.md)
- [`vignette("Resampling", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/Resampling.md)
- [`vignette("AnalysisWorkflows", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/AnalysisWorkflows.md)

### Advanced and specialized articles

- [`vignette("ImageVolumes", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/ImageVolumes.md)
- [`vignette("NeuroVector", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/NeuroVector.md)
- [`vignette("regionOfInterest", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/regionOfInterest.md)
- [`vignette("clustered-neurovec", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/clustered-neurovec.md)
- [`vignette("pipelines", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/pipelines.md)
- [`vignette("slice-visualization", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/slice-visualization.md)
- [`vignette("Cookbook", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/Cookbook.md)

### Reference and help

``` r
help(package = "neuroim2")
help.search("roi", package = "neuroim2")
```

## Summary

The package becomes much easier to navigate if you treat this overview
as a map, not a manual. Learn `NeuroVol`, `NeuroVec`, and ROI extraction
here, then move into the focused workflow vignettes for backend choice,
spatial transforms, and analysis patterns.
