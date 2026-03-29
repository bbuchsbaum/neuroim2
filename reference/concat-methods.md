# Concatenate two objects in the time dimension

Concatenate two objects in the time dimension

## Usage

``` r
concat(x, y, ...)

# S4 method for class 'NeuroVec,NeuroVol'
concat(x, y, ...)

# S4 method for class 'NeuroVol,NeuroVec'
concat(x, y, ...)

# S4 method for class 'NeuroVec,NeuroVec'
concat(x, y, ...)

# S4 method for class 'ROIVec,ROIVec'
concat(x, y, ...)

# S4 method for class 'DenseNeuroVol,missing'
concat(x, y, ...)

# S4 method for class 'DenseNeuroVol,DenseNeuroVol'
concat(x, y, ...)

# S4 method for class 'AbstractSparseNeuroVec,missing'
concat(x, y, ...)

# S4 method for class 'SparseNeuroVec,SparseNeuroVec'
concat(x, y, ...)
```

## Arguments

- x:

  the first object, typically `NeuroVol` or `NeuroVec`

- y:

  the second object, typically `NeuroVol` or `NeuroVec`

- ...:

  additional objects

## Value

A temporally concatenated object.

## Details

The `x` and `y` images must have compatible dimensions. A `NeuroVol` can
be concatenated to `NeuroVec`, and vice versa. See examples.

## Note

dimensions of x and y must be equal

## Examples

``` r
bv1 <- NeuroVol(rep(1,1000), NeuroSpace(c(10,10,10), c(1,1,1)))
bv2 <- NeuroVol(rep(2,1000), NeuroSpace(c(10,10,10), c(1,1,1)))
bv3 <- concat(bv1,bv2)
inherits(bv3, "NeuroVec")
#> [1] TRUE

bv4 <- concat(bv3, bv1)
dim(bv4)[4] == 3
#> [1] TRUE
bv5 <- concat(bv1, bv3)
dim(bv4)[4] == 3
#> [1] TRUE

bv6 <- concat(bv4,bv5)
dim(bv6)[4] == 6
#> [1] TRUE
```
