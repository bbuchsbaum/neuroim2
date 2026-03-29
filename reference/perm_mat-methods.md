# Extract permutation matrix associated with an image

A permutation matrix defines how the native voxel coordinates can be
transformed to standard (LPI) orientation.

## Usage

``` r
perm_mat(x, ...)

# S4 method for class 'AxisSet2D'
perm_mat(x, ...)

# S4 method for class 'AxisSet3D'
perm_mat(x, ...)

# S4 method for class 'NeuroSpace'
perm_mat(x, ...)
```

## Arguments

- x:

  A NeuroSpace object

- ...:

  Additional arguments (not used)

## Value

A numeric N x N `matrix` representing the permutation transform, where N
is the dimensionality of the image.

A matrix representing the axis directions

A matrix representing the axis directions

A numeric N x N `matrix` representing the permutation transform, where N
is the dimensionality of the image.

## Details

a permutation matrix can be used to convert between cardinal image
orientations. For example, if an image is stored in "RPI"
(Right-Posterior-Inferior) format, a coordinate in this space can be
converted to LPI (Left-Posterior-Inferior) by multiplying a coordinate
vector by the permutation matrix.

## Examples

``` r
fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
vol <- read_vol(fname)
pmat <- perm_mat(space(vol))

vox <- c(12,12,8)
pvox <- vox %*% perm_mat(space(vol))

stopifnot(all(pvox == c(-12,12,8)))
```
