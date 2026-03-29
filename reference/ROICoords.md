# Create ROI Coordinates Object

Creates an
[`ROICoords`](https://bbuchsbaum.github.io/neuroim2/reference/ROICoords-class.md)
object from a matrix of coordinates representing points in 3D space.

## Usage

``` r
ROICoords(coords)
```

## Arguments

- coords:

  A matrix with 3 columns representing (x, y, z) coordinates

## Value

An
[`ROICoords`](https://bbuchsbaum.github.io/neuroim2/reference/ROICoords-class.md)
object

## Details

ROI Coordinates

## Examples

``` r
coords <- matrix(c(1,2,3, 4,5,6), ncol=3, byrow=TRUE)
roi_coords <- ROICoords(coords)
```
