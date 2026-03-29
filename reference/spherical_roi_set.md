# Create Multiple Spherical Regions of Interest

This function generates multiple spherical ROIs simultaneously, centered
at the provided voxel coordinates. It is more efficient than calling
`spherical_roi` multiple times when you need to create many ROIs.

## Usage

``` r
spherical_roi_set(bvol, centroids, radius, fill = NULL, nonzero = FALSE)
```

## Arguments

- bvol:

  A `NeuroVol` or `NeuroSpace` instance

- centroids:

  A matrix of voxel coordinates where each row represents a centroid
  (i,j,k)

- radius:

  The radius in real units (e.g. millimeters) of the spherical ROIs

- fill:

  Optional value(s) to store as data. If provided, must be either a
  single value or a vector with length equal to the number of ROIs

- nonzero:

  If `TRUE`, keep only nonzero elements from `bvol`

## Value

A list of `ROIVolWindow` objects, one for each centroid

## Examples

``` r
# Create a NeuroSpace object
sp1 <- NeuroSpace(c(10,10,10), c(1,2,3))

# Create multiple ROIs centered at different voxel coordinates
centroids <- matrix(c(5,5,5, 3,3,3, 7,7,7), ncol=3, byrow=TRUE)
rois <- spherical_roi_set(sp1, centroids, 3.5)

# Create ROIs with specific fill values
rois <- spherical_roi_set(sp1, centroids, 3.5, fill=c(1,2,3))
```
