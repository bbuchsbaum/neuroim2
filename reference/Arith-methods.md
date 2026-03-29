# Arithmetic Operations

Methods for performing arithmetic operations on neuroimaging objects

This method performs arithmetic operations between two ROIVol objects
(`e1` and `e2`) using a generic arithmetic function. The dimensions of
both objects are checked for compatibility before performing the
operation.

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

- e1:

  A NeuroVol object.

- e2:

  A NeuroVec object.

## Value

A SparseNeuroVol object representing the result of the arithmetic
operation.

An ROIVol object resulting from the arithmetic operation.

An ROIVol object containing the result of the arithmetic operation
between `e1` and `e2`.

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
