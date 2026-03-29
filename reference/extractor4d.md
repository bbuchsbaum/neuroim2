# Array-like access for 4-dimensional data structures

This generic function provides array-like access for 4-dimensional data
structures. It allows for flexible indexing and subsetting of 4D arrays
or array-like objects.

Provides array-like access to ClusteredNeuroVec objects, supporting
extraction patterns like x\[,,,t\] to get 3D volumes at specific time
points.

## Usage

``` r
# S4 method for class 'ArrayLike4D,matrix,missing,ANY'
x[i, j, k, m, ..., drop = TRUE]

# S4 method for class 'ArrayLike4D,numeric,numeric,ANY'
x[i, j, k, m, ..., drop = TRUE]

# S4 method for class 'ArrayLike4D,numeric,missing,ANY'
x[i, j, k, m, ..., drop = TRUE]

# S4 method for class 'ArrayLike4D,integer,missing,ANY'
x[i, j, k, m, ..., drop = TRUE]

# S4 method for class 'ArrayLike4D,missing,missing,ANY'
x[i, j, k, m, ..., drop = TRUE]

# S4 method for class 'ArrayLike4D,missing,numeric,ANY'
x[i, j, k, m, ..., drop = TRUE]

# S4 method for class 'ClusteredNeuroVec,missing,missing,ANY'
x[i, j, k, m, ..., drop = TRUE]

# S4 method for class 'ClusteredNeuroVec,missing,missing,ANY'
x[i, j, k, m, ..., drop = TRUE]

# S4 method for class 'ClusteredNeuroVec,numeric,numeric,ANY'
x[i, j, k, m, ..., drop = TRUE]
```

## Arguments

- x:

  The 4-dimensional object to be accessed.

- i:

  First index or dimension.

- j:

  Second index or dimension.

- k:

  Third index or dimension.

- m:

  Fourth index or dimension.

- ...:

  Additional arguments passed to methods.

- drop:

  Logical. If TRUE, the result is coerced to the lowest possible
  dimension.

## Value

A subset of the input object, with dimensions depending on the indexing
and the \`drop\` parameter.
