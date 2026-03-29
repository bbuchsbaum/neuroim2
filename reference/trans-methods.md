# Extract image coordinate transformation

Extract image coordinate transformation

Get transformation matrix

## Usage

``` r
trans(x)

# S4 method for class 'MetaInfo'
trans(x)

# S4 method for class 'NeuroObj'
trans(x)

# S4 method for class 'NeuroHyperVec'
trans(x)

# S4 method for class 'NeuroSpace'
trans(x)
```

## Arguments

- x:

  an object with a transformation

## Value

A numeric 4x4 matrix that maps from grid coordinates to real-world
coordinates.

## Details

This function returns a transformation that can be used to go from "grid
coordinates" to "real world coordinates" in millimeters. see
[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)

## Examples

``` r
bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
trans(bspace)
#>      [,1] [,2] [,3] [,4]
#> [1,]    2    0    0    0
#> [2,]    0    2    0    0
#> [3,]    0    0    2    0
#> [4,]    0    0    0    1
all.equal(dim(trans(bspace)), c(4,4))
#> [1] TRUE
```
