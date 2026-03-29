# Extract full sparse rows across time.

This function extracts one or more rows from the sparse time-by-voxel
backing representation used by sparse neuroimaging vectors. It
complements
[`matricized_access()`](https://bbuchsbaum.github.io/neuroim2/reference/matricized_access-methods.md),
which is the column-oriented accessor.

## Usage

``` r
temporal_access(x, i, ...)

# S4 method for class 'SparseNeuroVec,integer'
temporal_access(x, i)

# S4 method for class 'SparseNeuroVec,numeric'
temporal_access(x, i)

# S4 method for class 'BigNeuroVec,integer'
temporal_access(x, i)

# S4 method for class 'BigNeuroVec,numeric'
temporal_access(x, i)
```

## Arguments

- x:

  a data source, typically a `SparseNeuroVec` object containing 4D
  neuroimaging data

- i:

  a numeric vector of temporal indices to extract

- ...:

  additional arguments to be passed to methods.

## Value

A matrix with one row per requested time index and one column per sparse
voxel in the backing representation.
