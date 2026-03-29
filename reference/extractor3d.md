# Array-like access for 3-dimensional data structures

This generic function provides array-like access for 3-dimensional data
structures. It allows for flexible indexing and subsetting of 3D arrays
or array-like objects.

## Usage

``` r
# S4 method for class 'ArrayLike3D,numeric,missing,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'ArrayLike3D,matrix,missing,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'ArrayLike3D,missing,missing,ANY'
x[i, j, k, ..., drop = TRUE]

# S4 method for class 'ArrayLike3D,missing,numeric,ANY'
x[i, j, k, ..., drop = TRUE]
```

## Arguments

- x:

  The 3-dimensional object to be accessed.

- i:

  First index or dimension.

- j:

  Second index or dimension.

- k:

  Third index or dimension.

- ...:

  Additional arguments passed to methods.

- drop:

  Logical. If TRUE, the result is coerced to the lowest possible
  dimension.

## Value

A subset of the input object, with dimensions depending on the indexing
and the \`drop\` parameter.
