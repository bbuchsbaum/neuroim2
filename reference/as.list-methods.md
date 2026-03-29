# Convert FileBackedNeuroVec to List

Converts a FileBackedNeuroVec object to a list of DenseNeuroVol objects.

convert SparseNeuroVec to list of
[`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md)

## Usage

``` r
# S4 method for class 'FileBackedNeuroVec'
as.list(x)

# S4 method for class 'NeuroVec'
as.list(x)

# S4 method for class 'SparseNeuroVec'
as.list(x)
```

## Arguments

- x:

  the object

## Value

A list of DenseNeuroVol objects

## Details

This method creates a deferred list, where each element is a
DenseNeuroVol object representing a single volume from the
FileBackedNeuroVec.
