# Arithmetic Operations

Methods for performing arithmetic operations on neuroimaging objects

Perform an arithmetic operation between two DenseNeuroVec objects. The
input DenseNeuroVec objects must have the same dimensions and NeuroSpace
objects. The method computes the elementwise arithmetic operation and
returns a new DenseNeuroVec object.

Perform an arithmetic operation between a SparseNeuroVol object and a
NeuroVol object. The input SparseNeuroVol and NeuroVol objects must have
the same dimensions. The method performs the arithmetic operation on the
non-zero values of the SparseNeuroVol and the corresponding values of
the NeuroVol. The result is returned as a new DenseNeuroVol object.

Perform an arithmetic operation between a NeuroVol object and a
SparseNeuroVol object. The input NeuroVol and SparseNeuroVol objects
must have the same dimensions. The method performs the arithmetic
operation on the values of the NeuroVol and the non-zero values of the
SparseNeuroVol. The result is returned as a new DenseNeuroVol object.

Perform an arithmetic operation between two NeuroVec objects. The input
NeuroVec objects must have the same dimensions. The method performs the
arithmetic operation on the elements of the NeuroVec objects. The result
is returned as a new DenseNeuroVec object.

This function performs arithmetic operations on a NeuroVec object and a
NeuroVol object.

This function performs arithmetic operations on a NeuroVol object and a
NeuroVec object.

## Usage

``` r
# S4 method for class 'SparseNeuroVol,SparseNeuroVol'
Arith(e1, e2)

# S4 method for class 'ROIVol,ROIVol'
Arith(e1, e2)

# S4 method for class 'DenseNeuroVol,DenseNeuroVol'
Arith(e1, e2)

# S4 method for class 'DenseNeuroVec,DenseNeuroVec'
Arith(e1, e2)

# S4 method for class 'SparseNeuroVol,NeuroVol'
Arith(e1, e2)

# S4 method for class 'NeuroVol,SparseNeuroVol'
Arith(e1, e2)

# S4 method for class 'SparseNeuroVec,SparseNeuroVec'
Arith(e1, e2)

# S4 method for class 'NeuroVec,NeuroVec'
Arith(e1, e2)

# S4 method for class 'NeuroVec,NeuroVol'
Arith(e1, e2)

# S4 method for class 'NeuroVol,NeuroVec'
Arith(e1, e2)

# S4 method for class 'DenseNeuroVol,numeric'
Arith(e1, e2)

# S4 method for class 'numeric,DenseNeuroVol'
Arith(e1, e2)

# S4 method for class 'SparseNeuroVol,numeric'
Arith(e1, e2)

# S4 method for class 'numeric,SparseNeuroVol'
Arith(e1, e2)

# S4 method for class 'ClusteredNeuroVol,ClusteredNeuroVol'
Arith(e1, e2)

# S4 method for class 'ClusteredNeuroVol,numeric'
Arith(e1, e2)

# S4 method for class 'numeric,ClusteredNeuroVol'
Arith(e1, e2)

# S4 method for class 'ClusteredNeuroVol,NeuroVol'
Arith(e1, e2)

# S4 method for class 'NeuroVol,ClusteredNeuroVol'
Arith(e1, e2)
```

## Arguments

- e1, e2:

  Neuroimaging operands or numeric values.

## Value

A SparseNeuroVol object representing the result of the arithmetic
operation.

For the `ROIVol` method, an `ROIVol` with the same coordinates as `e1`,
containing the element-wise result.

A DenseNeuroVec object representing the result of the arithmetic
operation.

A DenseNeuroVol object representing the result of the arithmetic
operation.

A DenseNeuroVol object representing the result of the arithmetic
operation.

A DenseNeuroVec object representing the result of the arithmetic
operation.

A DenseNeuroVec object resulting from the arithmetic operation.

A DenseNeuroVec object resulting from the arithmetic operation.

## Details

**ROIVol contract:** element-wise arithmetic requires the same spatial
support. Both operands must share the same `NeuroSpace` and the same set
of voxel coordinates (order may differ). Missing voxels are not treated
as zero, and the result does not grow to the union of the two ROIs.
Values are aligned by linear voxel index before the operation.
