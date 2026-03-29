# Create a Spherical Region of Interest

Creates a Spherical ROI based on a centroid.

## Usage

``` r
spherical_roi(
  bvol,
  centroid,
  radius,
  fill = NULL,
  nonzero = FALSE,
  use_cpp = TRUE
)
```

## Arguments

- bvol:

  an `NeuroVol` or `NeuroSpace` instance

- centroid:

  the center of the sphere in positive-coordinate (i,j,k) voxel space.

- radius:

  the radius in real units (e.g. millimeters) of the spherical ROI

- fill:

  optional value(s) to store as data

- nonzero:

  if `TRUE`, keep only nonzero elements from `bvol`

- use_cpp:

  whether to use compiled c++ code

## Value

an instance of class `ROIVol`

## See also

\[spherical_roi_set()\] for efficiently creating many spherical ROIs,
\[series_roi()\] and \[coords()\] for extracting time series and
coordinates from ROIs, and the vignette:
[`vignette("regionOfInterest", package = "neuroim2")`](https://bbuchsbaum.github.io/neuroim2/articles/regionOfInterest.md).

## Examples

``` r
 sp1 <- NeuroSpace(c(10,10,10), c(1,2,3))
 # create an ROI centered around the integer-valued positive voxel coordinate: i=5, j=5, k=5
 cube <- spherical_roi(sp1, c(5,5,5), 3.5)
 vox <- coords(cube)
 cds <- coords(cube, real=TRUE)
 ## fill in ROI with value of 6
 cube1 <- spherical_roi(sp1, c(5,5,5), 3.5, fill=6)
 all(cube1 == 6)
#> [1] TRUE

 ## Create multiple spherical ROIs at once (preferred):
 centers <- rbind(c(5,5,5), c(3,3,3), c(7,7,7))
 vols <- spherical_roi_set(bvol = sp1,
                          centroids = centers, radius = 3.5, fill = 1)
 length(vols)  # 3
#> [1] 3

 ## Equivalent, less efficient lapply variant:
 vols2 <- lapply(seq_len(nrow(centers)), function(i) {
   spherical_roi(sp1, centers[i,], radius = 3.5, fill = 1)
 })

 # create an ROI centered around the real-valued coordinates: x=5, y=5, z=5
 vox <- coord_to_grid(sp1, c(5, 5, 5))
 cube <- spherical_roi(sp1, vox, 3.5)
```
