# Extract Geometric Properties of an Image

This function retrieves the geometric properties of a given image, such
as dimensions and voxel size.

Retrieves the NeuroSpace object associated with an IndexLookupVol
object.

## Usage

``` r
space(x, ...)

# S4 method for class 'ClusteredNeuroVec'
space(x)

# S4 method for class 'IndexLookupVol'
space(x)

# S4 method for class 'ROICoords'
space(x)

# S4 method for class 'NeuroObj'
space(x)

# S4 method for class 'NeuroHyperVec'
space(x)

# S4 method for class 'NeuroSpace'
space(x)
```

## Arguments

- x:

  An
  [`IndexLookupVol`](https://bbuchsbaum.github.io/neuroim2/reference/IndexLookupVol-class.md)
  object

- ...:

  Additional arguments, if needed.

## Value

A
[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
object representing the geometric space of `x`.

## Examples

``` r
# Create a NeuroSpace object with dimensions (10, 10, 10) and voxel size (1, 1, 1)
x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))

# Create a NeuroVol object with random data and the specified NeuroSpace
vol <- NeuroVol(rnorm(10 * 10 * 10), x)

# Retrieve the geometric properties of the NeuroVol object
identical(x, space(vol))
#> [1] TRUE

space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
ilv <- IndexLookupVol(space, c(1:100))
space(ilv)  # Get the associated NeuroSpace object
#> <NeuroSpace> [3D] 
#> ── Geometry ──────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 64
#>   Spacing       : 1 x 1 x 1 mm
#>   Origin        : 0, 0, 0
#>   Orientation   : RAS
#>   Voxels        : 262,144

```
