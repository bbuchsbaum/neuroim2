# convert 1d indices to n-dimensional real world coordinates

convert 1d indices to n-dimensional real world coordinates

## Usage

``` r
index_to_coord(x, idx)

# S4 method for class 'NeuroSpace,numeric'
index_to_coord(x, idx)

# S4 method for class 'NeuroSpace,integer'
index_to_coord(x, idx)

# S4 method for class 'NeuroVol,integer'
index_to_coord(x, idx)

# S4 method for class 'NeuroVec,integer'
index_to_coord(x, idx)
```

## Arguments

- x:

  the object

- idx:

  the 1D indices

## Value

A numeric `matrix` of real-world coordinates.

## Examples

``` r
bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
idx <- 1:10
g <- index_to_coord(bvol, idx)
idx2 <- coord_to_index(bvol, g)
all.equal(idx, idx2)
#> [1] TRUE
```
