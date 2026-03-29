# Comparison Operations

Methods for comparing neuroimaging objects. All volume comparisons
return
[`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
objects that preserve spatial metadata.

## Usage

``` r
# S4 method for class 'DenseNeuroVol,DenseNeuroVol'
Compare(e1, e2)

# S4 method for class 'DenseNeuroVol,numeric'
Compare(e1, e2)

# S4 method for class 'numeric,DenseNeuroVol'
Compare(e1, e2)

# S4 method for class 'SparseNeuroVol,numeric'
Compare(e1, e2)

# S4 method for class 'numeric,SparseNeuroVol'
Compare(e1, e2)

# S4 method for class 'NeuroVec,NeuroVec'
Compare(e1, e2)
```

## Arguments

- e1, e2:

  Neuroimaging objects or numeric values.

## Value

A
[`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
for volume comparisons.
