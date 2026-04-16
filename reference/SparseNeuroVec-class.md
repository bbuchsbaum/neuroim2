# SparseNeuroVec Class

A class representing a sparse four-dimensional brain image, optimized
for efficient storage and access of large, sparse neuroimaging data.

Constructs a SparseNeuroVec object for efficient representation and
manipulation of sparse neuroimaging data with many zero or missing
values.

## Usage

``` r
SparseNeuroVec(data, space, mask, label = "", volume_labels = character())
```

## Arguments

- data:

  A matrix or a 4-D array containing the neuroimaging data. The
  dimensions of the data should be consistent with the dimensions of the
  provided NeuroSpace object and mask.

- space:

  A
  [NeuroSpace](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace.md)
  object representing the dimensions and voxel spacing of the
  neuroimaging data.

- mask:

  A 3D array, 1D vector of type logical, or an instance of type
  [LogicalNeuroVol](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md),
  which specifies the locations of the non-zero values in the data.

- label:

  Optional character string providing a label for the vector

- volume_labels:

  Optional character vector of length `dim(space)[4]` giving per-volume
  labels.

## Value

A SparseNeuroVec object, containing the sparse neuroimaging data, mask,
and associated NeuroSpace information.

## Details

SparseNeuroVec objects store data in a compressed format, where only
non-zero values are retained. This approach significantly reduces memory
usage for sparse brain images. The class leverages the mask and mapping
from its parent class
[`AbstractSparseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/AbstractSparseNeuroVec-class.md)
to efficiently manage the spatial structure of the data.

## Slots

- `data`:

  A `matrix` where each column represents a non-zero vector spanning the
  fourth dimension (e.g., time series for each voxel). Rows correspond
  to voxels in the sparse domain defined by the mask.

## Inheritance

`SparseNeuroVec` inherits from:

- [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md):
  Base class for 4D brain images

- [`AbstractSparseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/AbstractSparseNeuroVec-class.md):
  Provides sparse representation framework

- `ArrayLike4D`: Interface for 4D array-like operations

## See also

[`AbstractSparseNeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/AbstractSparseNeuroVec-class.md)
for the parent sparse representation class.
[`NeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
for the base 4D brain image class.

## Examples

``` r
# Create a sparse 4D brain image
mask <- LogicalNeuroVol(array(runif(64*64*32) > 0.7, c(64,64,32)), NeuroSpace(c(64,64,32)))
data <- matrix(rnorm(sum(mask) * 100), nrow=sum(mask), ncol=100)
sparse_vec <- SparseNeuroVec(data=data, mask=mask, space=NeuroSpace(dim=c(64,64,32,100)))

# Access a subset of the data
subset <- sparse_vec[,,, 1:10]


bspace <- NeuroSpace(c(10,10,10,100), c(1,1,1))
mask <- array(rnorm(10*10*10) > .5, c(10,10,10))
mat <- matrix(rnorm(sum(mask)), 100, sum(mask))
svec <- SparseNeuroVec(mat, bspace, mask)
length(indices(svec)) == sum(mask)
#> [1] TRUE
```
