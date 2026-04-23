# Estimate an image clip level

Computes an AFNI-inspired clip threshold for separating foreground from
low-intensity background. For 4D images, the threshold is computed from
a representative 3D volume, typically the voxelwise median across time.

## Usage

``` r
clip_level(x, mfrac = 0.5, gradual = FALSE, representative = "median")

# S4 method for class 'FileBackedNeuroVec'
clip_level(x, mfrac = 0.5, gradual = FALSE, representative = "median")

# S4 method for class 'MappedNeuroVec'
clip_level(x, mfrac = 0.5, gradual = FALSE, representative = "median")

# S4 method for class 'DenseNeuroVec'
clip_level(x, mfrac = 0.5, gradual = FALSE, representative = "median")

# S4 method for class 'AbstractSparseNeuroVec'
clip_level(x, mfrac = 0.5, gradual = FALSE, representative = "median")

# S4 method for class 'NeuroVol'
clip_level(x, mfrac = 0.5, gradual = FALSE, representative = "median")
```

## Arguments

- x:

  A neuroimaging object.

- mfrac:

  Fraction used in the median-update step. Values outside `(0, 0.99)`
  fall back to `0.5`.

- gradual:

  If `FALSE`, return a scalar clip level. If `TRUE`, return a 3D clip
  map interpolated across image octants.

- representative:

  For multi-volume inputs, the 3D summary image used for thresholding.
  Supported values are `"median"` and `"mean_abs"`.

## Value

If `gradual = FALSE`, a numeric scalar clip level. If `gradual = TRUE`,
a `DenseNeuroVol` with voxelwise clip levels.

## Examples

``` r
sp <- NeuroSpace(c(8, 8, 8), spacing = c(1, 1, 1))
vol <- NeuroVol(array(abs(rnorm(512)), c(8, 8, 8)), sp)
clip_level(vol)
#> [1] 0.5240544
```
