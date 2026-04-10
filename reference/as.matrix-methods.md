# convert a `NeuroVec` to a matrix

convert a `NeuroVec` to a matrix

Converts a
[`NeuroVecSeq`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSeq-class.md)
to a dense voxel-by-time matrix.

convert a
[`ROIVec`](https://bbuchsbaum.github.io/neuroim2/reference/ROIVec-class.md)
to a matrix

Convert to Matrix

## Usage

``` r
# S4 method for class 'ClusteredNeuroVec'
as.matrix(x, by = c("cluster", "voxel"))

# S4 method for class 'MappedNeuroVec'
as.matrix(x)

# S4 method for class 'NeuroVec'
as.matrix(x)

# S4 method for class 'DenseNeuroVec'
as.matrix(x)

# S4 method for class 'NeuroVecSeq'
as.matrix(x, ...)

# S4 method for class 'ROIVec'
as.matrix(x)

# S4 method for class 'AbstractSparseNeuroVec'
as.matrix(x, ...)
```

## Arguments

- x:

  The object to convert to a matrix

- by:

  For ClusteredNeuroVec: controls the conversion target. Defaults to
  "cluster" to return a TxK matrix of cluster time-series. "voxel" is
  reserved for future use.

- ...:

  Additional arguments

## Value

A matrix representation of the object
