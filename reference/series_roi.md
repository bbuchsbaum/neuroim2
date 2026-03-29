# Extract time series from specific voxel coordinates and return as ROI object

Extracts time series data from a `NeuroVec` object at specified voxel
coordinates and returns it as an ROI object.

## Usage

``` r
series_roi(x, i, ...)
```

## Arguments

- x:

  The `NeuroVec` object

- i:

  Numeric index for the first dimension

- ...:

  Additional arguments

## Value

A `ROIVec` object containing the time series data for the specified
coordinates.

## See also

[`series`](https://bbuchsbaum.github.io/neuroim2/reference/series-methods.md)

## Examples

``` r
# Create a simple 4D neuroimaging vector
space <- NeuroSpace(c(10,10,10,20), c(1,1,1))
vec <- NeuroVec(array(rnorm(10*10*10*20), c(10,10,10,20)), space)

# Extract time series for first 100 voxels as ROI
roi1 <- series_roi(vec, 1:100)

# Extract time series using 3D coordinates
coords <- matrix(c(1,1,1, 2,2,2, 3,3,3), ncol=3, byrow=TRUE)
roi2 <- series_roi(vec, coords)
```
