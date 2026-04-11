# Working with Volumes and Vectors

`neuroim2` has two basic data containers you reach for constantly:
`NeuroVol` for one 3D image and `NeuroVec` for a 4D stack of volumes.
This article shows how to move between them, inspect geometry, and
extract the small pieces you actually analyze.

## Start with a single volume

Use
[`read_vol()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol.md)
when the file is a 3D image or when you want one volume at a time.

``` r
vol_file <- system.file("extdata", "global_mask2.nii.gz", package = "neuroim2")
vol <- read_vol(vol_file)

dim(vol)
#> [1] 64 64 25
spacing(vol)
#> [1] 3.5 3.5 3.7
origin(vol)
#> [1]  112.00 -108.50  -46.25
```

``` r
stopifnot(length(dim(vol)) == 3L)
stopifnot(all(spacing(vol) > 0))
```

The returned object behaves like an array, but it also carries a
`NeuroSpace` so you can keep voxel values tied to real spatial
coordinates.

``` r
space(vol)
#> <NeuroSpace> [3D] 
#> ── Geometry ──────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Spacing       : 3.5 x 3.5 x 3.7 mm
#>   Origin        : 112, -108.5, -46.25
#>   Orientation   : LAS
#>   Voxels        : 102,400
vol[1, 1, 1]
#> [1] 0
```

## Read a 4D time series

Use
[`read_vec()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md)
when the file contains a time series or another stack of aligned
volumes.

``` r
vec_file <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
vec <- read_vec(vec_file)

dim(vec)
#> [1] 64 64 25  4
vec
#> <DenseNeuroVec> [3.1 Mb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25 (4 timepoints)
#>   Spacing       : 3.5 x 3.5 x 3.7
#>   Origin        : 112, -108, -46.2
#>   Orientation   : LAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Mean +/- SD   : 0.288 +/- 0.453 (t=1)
#>   Label         : /home/runner/work/_temp/Library/neuroim2/extdata/global_mask_v4.nii
```

``` r
stopifnot(length(dim(vec)) == 4L)
stopifnot(dim(vec)[4] > 1L)
```

A `NeuroVec` still has one shared 3D spatial frame. The fourth dimension
is the axis you usually split, subset, and summarize.

## Pull out the volumes you need

[`sub_vector()`](https://bbuchsbaum.github.io/neuroim2/reference/sub_vector-methods.md)
is the direct way to keep only selected timepoints.

``` r
first_two <- sub_vector(vec, 1:2)

dim(first_two)
#> [1] 64 64 25  2
```

``` r
stopifnot(dim(first_two)[4] == 2L)
```

## Convert to a matrix when you need table-like work

[`as.matrix()`](https://bbuchsbaum.github.io/neuroim2/reference/as.matrix.md)
is useful when you want one row per voxel and one column per timepoint.

``` r
mat <- as.matrix(vec)

dim(mat)
#> [1] 102400      4
mat[1:4, 1:2]
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0
#> [3,]    0    0
#> [4,]    0    0
```

``` r
stopifnot(nrow(mat) == prod(dim(vec)[1:3]))
stopifnot(ncol(mat) == dim(vec)[4])
```

That representation is convenient for temporal summaries, scaling, or
passing the data into general R modeling code.

## Extract one voxel or one ROI time series

For a single voxel, use
[`series()`](https://bbuchsbaum.github.io/neuroim2/reference/series-methods.md).

``` r
voxel_ts <- series(vec, 12, 12, 12)
voxel_ts
#> [1] 0 0 0 0
```

For a spatial region, create an ROI and use
[`series_roi()`](https://bbuchsbaum.github.io/neuroim2/reference/series_roi.md).

``` r
roi <- spherical_roi(drop(first_two[[1]]), c(12, 12, 12), radius = 6)
roi_ts <- series_roi(vec, roi)

dim(values(roi_ts))
#> [1]  4 19
```

``` r
stopifnot(length(roi) > 0L)
stopifnot(nrow(values(roi_ts)) == dim(vec)[4])
```

## Switch to sparse storage when the support matters

If only a subset of voxels should be considered present, read through a
mask and keep the result sparse.

``` r
mask_vol <- read_vol(vec_file) > 0
sparse_vec <- read_vec(vec_file, mask = mask_vol)

class(sparse_vec)
#> [1] "SparseNeuroVec"
#> attr(,"package")
#> [1] "neuroim2"
dim(sparse_vec)
#> [1] 64 64 25  4
```

``` r
stopifnot(inherits(sparse_vec, "SparseNeuroVec"))
stopifnot(identical(dim(sparse_vec), dim(vec)))
```

This is the right move when absent voxels should be treated as missing
support, not as literal zeros.

## Where to go next

- [`vignette("ChoosingBackends")`](https://bbuchsbaum.github.io/neuroim2/articles/ChoosingBackends.md)
  for storage tradeoffs
- [`vignette("coordinate-systems")`](https://bbuchsbaum.github.io/neuroim2/articles/coordinate-systems.md)
  for affine and orientation details
- [`vignette("Resampling")`](https://bbuchsbaum.github.io/neuroim2/articles/Resampling.md)
  for grid changes and reorientation
- [`vignette("AnalysisWorkflows")`](https://bbuchsbaum.github.io/neuroim2/articles/AnalysisWorkflows.md)
  for ROIs, searchlights, and split-map-reduce patterns
