# Apply a spatial mask to an image

Zeroes voxels outside a 3D spatial mask while preserving the geometry of
the input object. Unlike
[`mask`](https://bbuchsbaum.github.io/neuroim2/reference/mask-methods.md),
which returns an object's spatial domain, `apply_mask()` modifies image
values.

## Usage

``` r
apply_mask(x, mask)

# S4 method for class 'FileBackedNeuroVec'
apply_mask(x, mask)

# S4 method for class 'MappedNeuroVec'
apply_mask(x, mask)

# S4 method for class 'NeuroHyperVec'
apply_mask(x, mask)

# S4 method for class 'DenseNeuroVec'
apply_mask(x, mask)

# S4 method for class 'AbstractSparseNeuroVec'
apply_mask(x, mask)

# S4 method for class 'NeuroVol'
apply_mask(x, mask)

# S4 method for class 'LogicalNeuroVol'
apply_mask(x, mask)

# S4 method for class 'SparseNeuroVol'
apply_mask(x, mask)
```

## Arguments

- x:

  A neuroimaging object.

- mask:

  A 3D mask supplied as a `LogicalNeuroVol`, numeric `NeuroVol`
  (thresholded at `> 0`), logical array/vector, integer voxel indices,
  or `NULL` for no masking.

## Value

A masked neuroimaging object. Dense inputs remain dense; sparse inputs
return a sparse object with the intersected mask.

## Examples

``` r
sp <- NeuroSpace(c(4, 4, 4), spacing = c(1, 1, 1))
vol <- NeuroVol(array(rnorm(64), c(4, 4, 4)), sp)
msk <- LogicalNeuroVol(array(runif(64) > 0.5, c(4, 4, 4)), sp)
masked <- apply_mask(vol, msk)
```
