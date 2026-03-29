# Create a resampled searchlight iterator

This function generates a resampled searchlight iterator by sampling
regions from within a brain mask. By default it builds spherical
searchlights, but users can provide a custom `shape_fun` to return
ellipsoids, cubes, or arbitrary irregular searchlight shapes. Centers
are drawn with replacement, so the same voxel (and its neighborhood) may
appear multiple times. Each searchlight can also draw its radius from a
user-specified set of radii.

## Usage

``` r
resampled_searchlight(
  mask,
  radius = 8,
  iter = 100,
  shape_fun = NULL,
  nonzero = TRUE
)

bootstrap_searchlight(mask, radius = 8, iter = 100)
```

## Arguments

- mask:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  object representing the brain mask.

- radius:

  A numeric scalar or vector specifying candidate radii (in voxel units)
  for the searchlight sphere. If a vector is supplied, a radius is
  sampled uniformly (with replacement) for each searchlight. All radii
  must be positive. Default is 8.

- iter:

  An integer specifying the total number of searchlights to sample (with
  replacement). Default is 100.

- shape_fun:

  Either `NULL` (default spherical kernel), a character keyword
  (`"sphere"`, `"ellipsoid"`, `"cube"`, `"blobby"`), or a custom
  function. Custom functions are called as
  `shape_fun(mask, center, radius, iter, nonzero)` and must return
  either a
  [`ROIVolWindow`](https://bbuchsbaum.github.io/neuroim2/reference/ROIVolWindow-class.md)
  or an `n x 3` integer matrix of voxel coordinates. This enables
  anisotropic or irregular searchlights.

- nonzero:

  Logical; if `TRUE` (default), the generated searchlight is intersected
  with the non-zero voxels of `mask`. Applies to both the default sphere
  and any `shape_fun` that returns coordinates.

## Value

A `deferred_list` object containing
[`ROIVolWindow`](https://bbuchsbaum.github.io/neuroim2/reference/ROIVolWindow-class.md)
objects, each representing a sampled searchlight region drawn from
within the mask.

## Details

Searchlight centers are sampled with replacement, so the same center
(and its surrounding voxels) can be selected multiple times. When
multiple radii are provided, each searchlight independently samples one
radius from the supplied values. Supplying `shape_fun` lets you draw
non-spherical searchlights (e.g., ellipsoids, cubes, blobby
deformations, or task-specific kernels). Built-in shortcuts are
available via `shape_fun = "ellipsoid"`, `"cube"`, and `"blobby"`;
`"sphere"` or `NULL` uses the default spherical kernel.

## Examples

``` r
# Load an example brain mask
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Generate a resampled searchlight iterator with radii drawn from {4,6,8}
searchlights <- resampled_searchlight(mask, radius = c(4, 6, 8))

# Use a custom shape: random ellipsoid scaled along each axis
ellipsoid_fun <- function(mask, center, radius, iter, nonzero) {
  scales <- runif(3, 0.5, 1.5)        # axis-wise stretch/compress
  vox <- spherical_roi(mask, center, radius, nonzero = FALSE)@coords
  ctr_mat <- matrix(center, nrow(vox), 3, byrow = TRUE)
  keep <- rowSums(((vox - ctr_mat) * scales)^2) <= radius^2
  vox[keep, , drop = FALSE]
}
ellip_searchlights <- resampled_searchlight(mask, radius = c(4, 6),
                                            iter = 50, shape_fun = ellipsoid_fun)

# Or use built-in named shapes
ellip_builtin <- resampled_searchlight(mask, radius = 6, shape_fun = "ellipsoid")
cube_builtin  <- resampled_searchlight(mask, radius = 6, shape_fun = "cube")

```
