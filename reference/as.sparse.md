# Convert to from dense to sparse representation

Convert to from dense to sparse representation

## Usage

``` r
as.sparse(x, mask, ...)
```

## Arguments

- x:

  the object to make sparse, e.g. `DenseNeuroVol` or `DenseNeuroVec`

- mask:

  the elements to retain

- ...:

  additional arguments

## Value

A sparse representation of the input object, containing only the
elements specified by `mask`.

## Details

`mask` can be an integer vector of 1D indices or a mask volume of class
`LogicalNeuroVol`

## Examples

``` r
bvol <- NeuroVol(array(runif(24*24*24), c(24,24,24)), NeuroSpace(c(24,24,24), c(1,1,1)))
indmask <- sort(sample(1:(24*24*24), 100))
svol <- as.sparse(bvol, indmask)


mask <- LogicalNeuroVol(runif(length(indmask)), space=space(bvol), indices=indmask)
sum(mask) == 100
#> [1] TRUE
```
