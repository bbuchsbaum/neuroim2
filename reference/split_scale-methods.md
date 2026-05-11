# Center and/or Scale Row-subsets of a Matrix or Matrix-like Object

This function centers and/or scales the row-subsets of a numeric matrix
or matrix-like object.

## Usage

``` r
split_scale(x, f, center, scale)

# S4 method for class 'matrix,factor,logical,logical'
split_scale(x, f, center = TRUE, scale = TRUE)

# S4 method for class 'matrix,factor,missing,missing'
split_scale(x, f)

# S4 method for class 'DenseNeuroVec,factor,missing,missing'
split_scale(x, f)

# S4 method for class 'DenseNeuroVec,factor,logical,missing'
split_scale(x, f, center)

# S4 method for class 'DenseNeuroVec,factor,logical,logical'
split_scale(x, f, center, scale)
```

## Arguments

- x:

  A numeric matrix or matrix-like object.

- f:

  The splitting object, typically a factor or a set of integer indices.
  Must be equal to the number of rows in the matrix.

- center:

  Should values within each submatrix be centered? If TRUE, the mean is
  removed from each column of the submatrix.

- scale:

  Should values be scaled? If TRUE, the vector is divided by the
  standard deviation for each column of the submatrix.

## Value

An object of the same class as `x`, with row-subsets centered and/or
scaled according to `f`.

## Examples

``` r

M <- matrix(rnorm(1000), 10, 100)
fac <- factor(rep(1:2, each=5))
Ms <- split_scale(M, fac)

## Correctly centered
all(abs(apply(Ms[fac == 1,], 2, mean)) < .000001)
#> [1] TRUE
all(abs(apply(Ms[fac == 2,], 2, mean)) < .000001)
#> [1] TRUE

## Correctly scaled
all.equal(apply(Ms[fac == 1,], 2, sd), rep(1, ncol(Ms)))
#> [1] TRUE
all.equal(apply(Ms[fac == 2,], 2, sd), rep(1, ncol(Ms)))
#> [1] TRUE
```
