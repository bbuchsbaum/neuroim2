# Temporal Mean of a NeuroVec

Computes the voxel-wise mean across the 4th dimension (time), returning
a 3D
[`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md)
or
[`SparseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVol-class.md).

## Usage

``` r
# S4 method for class 'DenseNeuroVec'
mean(x, ..., na.rm = FALSE)

# S4 method for class 'SparseNeuroVec'
mean(x, ..., na.rm = FALSE)

# S4 method for class 'NeuroVec'
mean(x, ..., na.rm = FALSE)
```

## Arguments

- x:

  A
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
  object.

- ...:

  Additional arguments (currently unused).

- na.rm:

  Logical; if `TRUE`, missing values are omitted before computing the
  temporal mean at each voxel.

## Value

A
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
containing the temporal mean at each voxel.

## Examples

``` r
bspace <- NeuroSpace(c(10, 10, 10, 20), c(1, 1, 1))
dat <- array(rnorm(10 * 10 * 10 * 20), c(10, 10, 10, 20))
vec <- DenseNeuroVec(dat, bspace)
mean_vol <- mean(vec)
dim(mean_vol)  # 10 10 10
#> [1] 10 10 10
```
