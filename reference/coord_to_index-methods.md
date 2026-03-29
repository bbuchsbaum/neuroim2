# convert n-dimensional real world coordinates to 1D indices

convert n-dimensional real world coordinates to 1D indices

## Usage

``` r
coord_to_index(x, coords)

# S4 method for class 'NeuroSpace,matrix'
coord_to_index(x, coords)

# S4 method for class 'NeuroSpace,numeric'
coord_to_index(x, coords)

# S4 method for class 'NeuroVol,matrix'
coord_to_index(x, coords)
```

## Arguments

- x:

  the object

- coords:

  a matrix of real world coordinates

## Value

An integer `vector` of 1D indices corresponding to `coords`.

## Examples

``` r
bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
coords <- matrix(c(.5,.5,.5, 1.5,1.5,1.5), ncol=3, byrow=TRUE)
idx <- coord_to_index(bvol, coords)
coords2 <- index_to_coord(bvol, idx)
all.equal(coords, coords2)
#> [1] TRUE
```
