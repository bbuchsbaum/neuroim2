# Coerce SparseNeuroVol to DenseNeuroVol

Convert a sparse volumetric image to a dense representation with the
same spatial geometry. Non-zero values stored in the sparse vector are
placed at their corresponding linear indices in the dense array; all
other voxels are 0.

Identity method: returns a `DenseNeuroVol` (or subclass such as
`LogicalNeuroVol`) unchanged.

This function provides a method to coerce an object of class `ROIVol` to
a `DenseNeuroVol` using the `as.dense` method.

## Usage

``` r
# S4 method for class 'ClusteredNeuroVol'
as.dense(x)

# S4 method for class 'SparseNeuroVol'
as.dense(x)

# S4 method for class 'DenseNeuroVol'
as.dense(x)

# S4 method for class 'ROIVol'
as.dense(x)

# S4 method for class 'SparseNeuroVec'
as.dense(x)
```

## Arguments

- x:

  An object of class `ROIVol` to be coerced to a `DenseNeuroVol`.

## Value

A
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
object representing the dense version of the clustered volume.

A DenseNeuroVol with identical spatial dimensions and values expanded
from the sparse representation.

A `DenseNeuroVol` object obtained by coercing the `ROIVol` object.
