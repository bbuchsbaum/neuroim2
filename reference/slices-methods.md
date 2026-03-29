# Extract an ordered series of 2D slices from a 3D or 4D object

This function extracts an ordered series of 2D slices from a 3D or 4D
object. The returned slices are in the order they appear in the original
object.

## Usage

``` r
slices(x, ...)

# S4 method for class 'NeuroVol'
slices(x)
```

## Arguments

- x:

  A NeuroVol object

- ...:

  Additional arguments to be passed to the underlying methods

## Value

A `list` of 2D `matrices`, each containing a slice from the input `x`.

A deflist object containing functions that return 2D slices of the
volume along the z-axis. The length of the deflist equals the number of
slices in the z dimension.

## Examples

``` r
# Create a simple 3D volume
space <- NeuroSpace(c(10,10,10), c(1,1,1))
vol <- NeuroVol(array(rnorm(10*10*10), c(10,10,10)), space)

# Get all slices along the z-axis
slc <- slices(vol)

# Number of slices equals the z dimension
length(slc) == dim(vol)[3]
#> [1] TRUE

# Each slice is a 2D matrix
dim(slc[[1]]) == c(10,10)
#> [1] TRUE TRUE
```
