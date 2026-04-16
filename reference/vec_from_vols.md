# Create NeuroVec from list of NeuroVol objects

Factory function to create a NeuroVec object from a list of NeuroVol
objects. This is a convenience wrapper around the NeuroVec constructor
that combines multiple 3D volumes into a single 4D NeuroVec.

## Usage

``` r
vec_from_vols(vols, mask = NULL)
```

## Arguments

- vols:

  A list of
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  objects. All volumes must have identical spatial dimensions.

- mask:

  An optional logical array or
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
  object defining the subset of voxels to include. If provided, a
  SparseNeuroVec will be created.

## Value

A
[`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
object (either DenseNeuroVec or SparseNeuroVec depending on whether a
mask is provided).

## See also

[`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md),
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol.md)

## Examples

``` r
# Create a simple NeuroVec from list of volumes
spc <- NeuroSpace(c(10, 10, 10))
vol1 <- NeuroVol(rnorm(10*10*10), spc)
vol2 <- NeuroVol(rnorm(10*10*10), spc)
vec <- vec_from_vols(list(vol1, vol2))
print(dim(vec))  # Should be c(10, 10, 10, 2)
#> [1] 10 10 10  2
```
