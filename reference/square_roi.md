# Create a square region of interest

This function creates a square region of interest (ROI) in a 3D volume,
where the z-dimension is fixed at one voxel coordinate. The ROI is
defined within a given NeuroVol or NeuroSpace instance.

## Usage

``` r
square_roi(bvol, centroid, surround, fill = NULL, nonzero = FALSE, fixdim = 3)
```

## Arguments

- bvol:

  A `NeuroVol` or `NeuroSpace` instance representing the 3D volume or
  space.

- centroid:

  A numeric vector of length 3, representing the center of the square
  ROI in voxel coordinates.

- surround:

  A non-negative integer specifying the number of voxels on either side
  of the central voxel.

- fill:

  An optional value or values to assign to the data slot of the
  resulting ROI. If not provided, no data will be assigned.

- nonzero:

  A logical value indicating whether to keep only nonzero elements from
  `bvol`. If `bvol` is a `NeuroSpace` instance, this argument is
  ignored.

- fixdim:

  A logical value indicating whether the fixed dimension is the third,
  or z, dimension. Default is TRUE.

## Value

An instance of class `ROIVol` representing the square ROI.

## Examples

``` r
sp1 <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
square <- square_roi(sp1, c(5, 5, 5), 1)
vox <- coords(square)
## a 3 X 3 X 1 grid
nrow(vox) == 9
#> [1] TRUE
```
