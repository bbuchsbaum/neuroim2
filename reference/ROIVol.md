# Create ROI Volume Object

Creates an
[`ROIVol`](https://bbuchsbaum.github.io/neuroim2/reference/ROIVol-class.md)
object representing a set of values at specific 3D coordinates within a
spatial reference system.

## Usage

``` r
ROIVol(space, coords, data)
```

## Arguments

- space:

  A
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
  object defining the spatial reference

- coords:

  A matrix with 3 columns representing (x,y,z) coordinates

- data:

  A numeric vector of values corresponding to each coordinate

## Value

An
[`ROIVol`](https://bbuchsbaum.github.io/neuroim2/reference/ROIVol-class.md)
object

## Details

ROI Volume

## Examples

``` r
space <- NeuroSpace(c(64,64,64))
coords <- matrix(c(1,2,3, 4,5,6), ncol=3, byrow=TRUE)
data <- c(1.5, 2.5)
roi_vol <- ROIVol(space, coords, data)
```
