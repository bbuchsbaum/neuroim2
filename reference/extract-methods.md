# \[\[

This function extracts a single volume from a NeuroVec object.

Extracts a subset of data from a sparse four-dimensional brain image
based on provided indices.

## Usage

``` r
# S4 method for class 'DenseNeuroVol,numeric,missing,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'DenseNeuroVol,integer,missing,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'NeuroVec,numeric'
x[[i]]

# S4 method for class 'NeuroVec,character'
x[[i]]

# S4 method for class 'NeuroVol,ROICoords,missing,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'NeuroVol,ROIVol,missing,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'DenseNeuroVol,ROIVol,missing,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'SparseNeuroVol,numeric,numeric,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'ROIVol,numeric,missing,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'ROIVol,logical,missing,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'ROIVol,ROICoords,missing,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'ROIVol,matrix,missing,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'ROIVol,missing,missing,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'ROIVol,missing,numeric,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'ROIVol,numeric,numeric,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'ROIVol,logical,numeric,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'ROIVol,ROICoords,numeric,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'ROIVol,matrix,numeric,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'ROICoords,numeric,missing,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'ROIVol,numeric,missing,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'NeuroVol,ROICoords,missing,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'AbstractSparseNeuroVec,numeric,numeric,ANY'
x[i, j, k, m, ..., drop = TRUE]
```

## Arguments

- x:

  An object of class `AbstractSparseNeuroVec`

- i:

  Numeric vector specifying the indices for the first dimension

- j:

  Numeric vector specifying the indices for the second dimension

- k:

  Numeric vector specifying the indices for the third dimension
  (optional)

- ...:

  Additional arguments passed to methods

- drop:

  Logical indicating whether to drop dimensions of length one (default:
  TRUE)

- m:

  Numeric vector specifying the indices for the fourth dimension
  (optional)

## Value

a DenseNeuroVol object

A subset of the input object, with dimensions depending on the indexing
and the \`drop\` parameter.

An array containing the extracted subset
