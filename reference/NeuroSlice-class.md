# NeuroSlice Class

Represents a two-dimensional brain image slice. This class extends both
the `array` class for data storage and the
[`NeuroObj`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroObj-class.md)
class for spatial information.

## Details

NeuroSlice objects are typically used to represent individual slices of
3D brain volumes or 2D projections of 3D data. They inherit the spatial
properties from NeuroObj and the data storage capabilities from array.

## See also

[`NeuroObj-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroObj-class.md),
[`NeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)

## Examples

``` r
# Create a simple 2D brain slice
slice_data <- matrix(rnorm(64*64), 64, 64)
slice_space <- NeuroSpace(dim=c(64L, 64L), origin=c(0, 0), spacing=c(1, 1))
brain_slice <- new("NeuroSlice", .Data=slice_data, space=slice_space)
```
