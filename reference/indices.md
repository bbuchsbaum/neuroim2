# Extract indices

Extract indices

## Usage

``` r
indices(x)
```

## Arguments

- x:

  the object to extract indices

## Value

A `vector` of indices from `x`.

## Examples

``` r
# Create a NeuroSpace object with 3mm voxels
space <- NeuroSpace(c(10,10,10), spacing=c(3,3,3))

# Create ROI coordinates in voxel space
coords <- matrix(c(1,1,1, 2,2,2), ncol=3, byrow=TRUE)

# Create ROI volume
roi_vol <- ROIVol(space, coords, data=c(1,2))

# Get linear indices of ROI voxels
idx <- indices(roi_vol)
# These indices can be used to index into a 3D array of size 10x10x10
```
