# AxisSet2D

A two-dimensional axis set representing an ordered pair of named axes.

## Slots

- `i`:

  The first axis, inherited from AxisSet1D

- `j`:

  The second axis, of class "NamedAxis"

## See also

[`AxisSet1D-class`](https://bbuchsbaum.github.io/neuroim2/reference/AxisSet1D-class.md),
[`AxisSet3D-class`](https://bbuchsbaum.github.io/neuroim2/reference/AxisSet3D-class.md)

## Examples

``` r
# Create an AxisSet2D object
axis1 <- new("NamedAxis", axis = "x", direction = 1)
axis2 <- new("NamedAxis", axis = "y", direction = 1)
axisSet2D <- new("AxisSet2D", i = axis1, j = axis2, ndim = 2L)
```
