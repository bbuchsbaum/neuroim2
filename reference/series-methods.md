# Extract one or more series from object

This function extracts time series data from specific voxel coordinates
in a 4D neuroimaging object. It supports multiple ways of specifying the
coordinates:

- Linear indices (1D)

- Grid coordinates (3D matrix)

- Individual x,y,z coordinates

## Usage

``` r
series(x, i, ...)

# S4 method for class 'ClusteredNeuroVec,numeric'
series(x, i, j, k, ...)

# S4 method for class 'NeuroHyperVec,ANY'
series(x, i, j, k, ...)

# S4 method for class 'NeuroVec,matrix'
series(x, i)

# S4 method for class 'NeuroVec,matrix'
series_roi(x, i)

# S4 method for class 'NeuroVec,ROICoords'
series(x, i)

# S4 method for class 'NeuroVec,ROICoords'
series_roi(x, i)

# S4 method for class 'NeuroVec,LogicalNeuroVol'
series(x, i)

# S4 method for class 'NeuroVec,NeuroVol'
series(x, i)

# S4 method for class 'NeuroVec,LogicalNeuroVol'
series_roi(x, i)

# S4 method for class 'NeuroVec,integer'
series(x, i, j, k, drop = TRUE)

# S4 method for class 'DenseNeuroVec,matrix'
series(x, i)

# S4 method for class 'DenseNeuroVec,integer'
series(x, i, j, k, drop = TRUE)

# S4 method for class 'NeuroVec,numeric'
series(x, i, j, k, drop = TRUE)

# S4 method for class 'NeuroVec,numeric'
series_roi(x, i, j, k, drop = TRUE)

# S4 method for class 'NeuroVecSeq,integer'
series(x, i, j, k, drop = TRUE)

# S4 method for class 'NeuroVecSeq,numeric'
series(x, i, j, k, drop = TRUE)

# S4 method for class 'NeuroVecSeq,matrix'
series(x, i)

# S4 method for class 'NeuroVecSeq,matrix'
series_roi(x, i)

# S4 method for class 'AbstractSparseNeuroVec,ROICoords'
series(x, i)

# S4 method for class 'AbstractSparseNeuroVec,matrix'
series(x, i)

# S4 method for class 'AbstractSparseNeuroVec,numeric'
series(x, i, j, k)

# S4 method for class 'AbstractSparseNeuroVec,integer'
series(x, i, j, k, drop = TRUE)
```

## Arguments

- x:

  A NeuroVecSeq object

- i:

  A matrix of ROI coordinates (n x 3)

- ...:

  Additional arguments (not used)

- j:

  second dimension index

- k:

  third dimension index

- drop:

  whether to drop dimension of length 1

## Value

A `list` or `array` containing the extracted series.

A 2D array with dimensions \[features x trials\]

A matrix where each column represents a voxel's time series

A ROIVec object containing the time series for the specified ROI

## Details

when x is a NeuroHyperVec object, the series method returns a 2D array
with dimensions \[features x trials\]

## See also

[`series_roi`](https://bbuchsbaum.github.io/neuroim2/reference/series_roi.md)

## Examples

``` r
# Create a simple 4D neuroimaging vector (10x10x10 volume with 20 timepoints)
space <- NeuroSpace(c(10,10,10,20), c(1,1,1))
vec <- NeuroVec(array(rnorm(10*10*10*20), c(10,10,10,20)), space)

# Extract time series using linear indices
ts1 <- series(vec, 1:10)  # Get time series for first 10 voxels

# Extract time series using 3D coordinates
coords <- matrix(c(1,1,1, 2,2,2, 3,3,3), ncol=3, byrow=TRUE)
ts2 <- series(vec, coords)  # Get time series for 3 specific voxel locations

# Extract single time series using x,y,z coordinates
ts3 <- series(vec, 5, 5, 5)  # Get time series from middle voxel
```
