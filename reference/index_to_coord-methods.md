# convert 1d indices to n-dimensional real world coordinates

convert 1d indices to n-dimensional real world coordinates

## Usage

``` r
index_to_coord(x, idx)

# S4 method for class 'NeuroSpace,numeric'
index_to_coord(x, idx)

# S4 method for class 'NeuroObj,numeric'
index_to_coord(x, idx)
```

## Arguments

- x:

  the object: a `NeuroSpace`, or anything that has one, such as a
  `NeuroVol` or `NeuroVec`

- idx:

  the 1D indices, integer or double

## Value

A numeric `matrix` of real-world coordinates with one row per index,
including when there is only one.

## Examples

``` r
bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
idx <- 1:10
g <- index_to_coord(bvol, idx)
idx2 <- coord_to_index(bvol, g)
all.equal(idx, idx2)
#> [1] TRUE

# Where, in millimetres, is the largest value of a map?
vol <- read_vol(system.file("extdata", "mni_downsampled.nii.gz",
                            package = "neuroim2"))
index_to_coord(vol, which.max(vol))
#>           [,1]     [,2]      [,3]
#> [1,] -43.72921 -3.43872 -29.75004
```
