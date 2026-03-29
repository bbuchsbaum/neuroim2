# AxisSet3D Class

A class representing a three-dimensional axis set, extending the
AxisSet2D class with an additional third axis.

## Slots

- `k`:

  A `NamedAxis` object representing the third axis.

## See also

[`AxisSet2D-class`](https://bbuchsbaum.github.io/neuroim2/reference/AxisSet2D-class.md),
[`NamedAxis-class`](https://bbuchsbaum.github.io/neuroim2/reference/NamedAxis-class.md)

## Examples

``` r
# Create NamedAxis objects for each dimension
x_axis <- new("NamedAxis", axis = "x", direction = 1)
y_axis <- new("NamedAxis", axis = "y", direction = 1)
z_axis <- new("NamedAxis", axis = "z", direction = 1)

# Create an AxisSet3D object
axis_set_3d <- new("AxisSet3D", i = x_axis, j = y_axis, k = z_axis, ndim = 3L)
```
