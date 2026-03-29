# Extract Voxel Dimensions of an Image

This function extracts the voxel dimensions of an image represented by
the input object.

## Usage

``` r
spacing(x)

# S4 method for class 'ROICoords'
spacing(x)

# S4 method for class 'NeuroObj'
spacing(x)

# S4 method for class 'NeuroHyperVec'
spacing(x)

# S4 method for class 'NeuroSpace'
spacing(x)
```

## Arguments

- x:

  The object representing the image.

## Value

A numeric vector specifying the voxel dimensions of `x`.

## Examples

``` r
bspace <- NeuroSpace(c(10, 10, 10), c(2, 2, 2))
all.equal(spacing(bspace), c(2, 2, 2))
#> [1] TRUE
```
