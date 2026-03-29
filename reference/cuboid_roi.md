# Create A Cuboid Region of Interest

Create A Cuboid Region of Interest

## Usage

``` r
cuboid_roi(bvol, centroid, surround, fill = NULL, nonzero = FALSE)
```

## Arguments

- bvol:

  an `NeuroVol` or `NeuroSpace` instance

- centroid:

  the center of the cube in *voxel* coordinates

- surround:

  the number of voxels on either side of the central voxel. A `vector`
  of length 3.

- fill:

  optional value(s) to assign to data slot.

- nonzero:

  keep only nonzero elements from `bvol`. If `bvol` is A `NeuroSpace`
  then this argument is ignored.

## Value

An instance of class `ROIVol` representing the cuboid region of
interest, containing the coordinates and values of voxels within the
specified region.

## Examples

``` r
 sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
 cube <- cuboid_roi(sp1, c(5,5,5), 3)
 vox <- coords(cube)
 cube2 <- cuboid_roi(sp1, c(5,5,5), 3, fill=5)

```
