# Convenience shape generators for `resampled_searchlight()`

Helpers that return ready-to-use `shape_fun` callbacks for
[`resampled_searchlight()`](https://bbuchsbaum.github.io/neuroim2/reference/resampled_searchlight.md),
covering a few sensible non-spherical defaults.

## Usage

``` r
ellipsoid_shape(scales = c(1, 1, 1), jitter = 0)

cube_shape()

blobby_shape(drop = 0.3, edge_fraction = 0.7)
```

## Arguments

- scales:

  Length-3 positive numeric vector scaling the x/y/z axes relative to a
  sphere (for `ellipsoid_shape`). Values \>1 stretch; \<1 compress.

- jitter:

  Non-negative numeric; standard deviation of multiplicative Gaussian
  noise applied to `scales` each draw (ellipsoid).

- drop:

  Numeric in \[0,1\]; probability of dropping a voxel (blobby).

- edge_fraction:

  Numeric in (0,1\]; fraction of farthest voxels (by Euclidean distance
  from the center, in voxel units) considered "edge" and eligible for
  random dropping (blobby).

## Value

A function suitable for the `shape_fun` argument of
[`resampled_searchlight()`](https://bbuchsbaum.github.io/neuroim2/reference/resampled_searchlight.md).

## Details

Each returned function has signature
`function(mask, center, radius, iter, nonzero)` and should return an \\n
\times 3\\ integer coordinate matrix. The coordinates are later
converted to a `ROIVolWindow` internally.

## Examples

``` r
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Ellipsoid stretched along z with modest per-iteration jitter
sl_ellip <- resampled_searchlight(mask, radius = 6,
                                   shape_fun = ellipsoid_shape(scales = c(1, 1, 1.4),
                                                              jitter = 0.1))

# Simple axis-aligned cube (Chebyshev ball)
sl_cube <- resampled_searchlight(mask, radius = 5, shape_fun = "cube")

# Blobby sphere with 40% dropout on boundary voxels
sl_blob <- resampled_searchlight(mask, radius = 6,
                                 shape_fun = blobby_shape(drop = 0.4, edge_fraction = 0.6))
```
