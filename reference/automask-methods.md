# Compute a brain-like mask from image intensities

Builds a spatial mask from image content rather than requiring an
external mask. The implementation is AFNI-inspired: it computes a clip
level, applies thresholding to a representative 3D image, retains the
largest connected component, and optionally applies peel/unpeel cleanup.

## Usage

``` r
automask(
  x,
  mfrac = 0.5,
  gradual = TRUE,
  representative = "mean_abs",
  peels = 1L,
  peel_threshold = 17L,
  connect = c("26-connect", "18-connect", "6-connect")
)

# S4 method for class 'FileBackedNeuroVec'
automask(
  x,
  mfrac = 0.5,
  gradual = TRUE,
  representative = "mean_abs",
  peels = 1L,
  peel_threshold = 17L,
  connect = c("26-connect", "18-connect", "6-connect")
)

# S4 method for class 'MappedNeuroVec'
automask(
  x,
  mfrac = 0.5,
  gradual = TRUE,
  representative = "mean_abs",
  peels = 1L,
  peel_threshold = 17L,
  connect = c("26-connect", "18-connect", "6-connect")
)

# S4 method for class 'DenseNeuroVec'
automask(
  x,
  mfrac = 0.5,
  gradual = TRUE,
  representative = "mean_abs",
  peels = 1L,
  peel_threshold = 17L,
  connect = c("26-connect", "18-connect", "6-connect")
)

# S4 method for class 'AbstractSparseNeuroVec'
automask(
  x,
  mfrac = 0.5,
  gradual = TRUE,
  representative = "mean_abs",
  peels = 1L,
  peel_threshold = 17L,
  connect = c("26-connect", "18-connect", "6-connect")
)

# S4 method for class 'NeuroVol'
automask(
  x,
  mfrac = 0.5,
  gradual = TRUE,
  representative = "mean_abs",
  peels = 1L,
  peel_threshold = 17L,
  connect = c("26-connect", "18-connect", "6-connect")
)
```

## Arguments

- x:

  A neuroimaging object.

- mfrac:

  Fraction used in the clip-level update step.

- gradual:

  If `TRUE`, use a spatially varying clip map across octants. If
  `FALSE`, use a single global clip threshold.

- representative:

  For multi-volume inputs, the 3D summary image used for thresholding.
  Supported values are `"mean_abs"` and `"median"`.

- peels:

  Number of AFNI-style peel/unpeel iterations. Use `0` to disable.

- peel_threshold:

  Minimum number of 18-neighbors required for a voxel to survive
  peeling.

- connect:

  Connectivity used when retaining the largest component.

## Value

A `LogicalNeuroVol`.

## Examples

``` r
sp <- NeuroSpace(c(8, 8, 8), spacing = c(1, 1, 1))
arr <- array(abs(rnorm(512)), c(8, 8, 8))
arr[3:6, 3:6, 3:6] <- arr[3:6, 3:6, 3:6] + 10
vol <- NeuroVol(arr, sp)
mask <- automask(vol, gradual = FALSE, peels = 0)
```
