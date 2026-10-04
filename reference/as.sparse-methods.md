# Convert DenseNeuroVec to sparse representation using mask

This method converts a DenseNeuroVec object to a sparse representation
using a given LogicalNeuroVol mask.

This method converts a DenseNeuroVec object to a sparse representation
using a given numeric mask.

## Usage

``` r
# S4 method for class 'DenseNeuroVec,LogicalNeuroVol'
as.sparse(x, mask)

# S4 method for class 'DenseNeuroVec,numeric'
as.sparse(x, mask)

# S4 method for class 'DenseNeuroVol,LogicalNeuroVol'
as.sparse(x, mask)

# S4 method for class 'DenseNeuroVol,numeric'
as.sparse(x, mask)

# S4 method for class 'ROIVol,ANY'
as.sparse(x)
```

## Arguments

- x:

  A DenseNeuroVec object to convert to a sparse representation.

- mask:

  A numeric vector representing the mask to apply during conversion.

## Value

A SparseNeuroVec object resulting from the conversion.

A SparseNeuroVec object resulting from the conversion.

## Details

For a `DenseNeuroVec`, numeric masks use R's vector-indexing rules:
positive indices retain voxels, negative indices exclude voxels, and
zeros are ignored. Repeated indices retain a voxel once. Retained voxels
are stored in ascending spatial-index order, regardless of the order
supplied in `mask`. An empty selection produces an all-zero sparse image
with the original dimensions. Missing or non-finite indices and positive
indices outside the spatial extent are rejected.
