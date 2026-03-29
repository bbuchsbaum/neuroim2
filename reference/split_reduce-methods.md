# Summarize Subsets of an Object by Splitting by Row and Applying a Summary Function

This function summarizes subsets of a numeric matrix or matrix-like
object by first splitting the object by row and then applying a summary
function.

## Usage

``` r
split_reduce(x, fac, FUN)

# S4 method for class 'matrix,integer,function'
split_reduce(x, fac, FUN)

# S4 method for class 'matrix,factor,missing'
split_reduce(x, fac)

# S4 method for class 'matrix,factor,function'
split_reduce(x, fac, FUN)

# S4 method for class 'NeuroVec,factor,function'
split_reduce(x, fac, FUN)

# S4 method for class 'NeuroVec,factor,missing'
split_reduce(x, fac, FUN)
```

## Arguments

- x:

  A numeric matrix or matrix-like object.

- fac:

  A factor to define subsets of the object.

- FUN:

  The summary function to apply to each subset. If not provided, the
  mean of each sub-matrix column is computed.

## Value

A `matrix` (or matrix-like object) containing the summarized values
after applying `FUN`.

## Details

If 'FUN' is supplied, it must take a vector and return a single scalar
value. If it returns more than one value, an error will occur.

If 'x' is a NeuroVec instance, voxels (dimensions 1:3) are treated as
columns and time-series (dimension 4) as rows. The summary function is
then applied to groups of voxels. However, if the goal is to apply a
function to groups of time-points.

## Examples

``` r
mat = matrix(rnorm(100*100), 100, 100)
fac = factor(sample(1:3, nrow(mat), replace=TRUE))
## Compute column means of each sub-matrix
ms <- split_reduce(mat, fac)
all.equal(row.names(ms), levels(fac))
#> [1] TRUE

## Compute column medians of each sub-matrix
ms <- split_reduce(mat, fac, median)

## Compute time-series means grouped over voxels.
## Here, 'length(fac)' must equal the number of voxels: 'prod(dim(bvec)[1:3])'
bvec <- NeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)), NeuroSpace(c(24,24,24,24), c(1,1,1)))
fac <- factor(sample(1:3, prod(dim(bvec)[1:3]), replace=TRUE))
ms <- split_reduce(bvec, fac)
ms2 <- split_reduce(bvec, fac, mean)
all.equal(row.names(ms), levels(fac))
#> [1] TRUE
all.equal(ms, ms2)
#> [1] TRUE
```
