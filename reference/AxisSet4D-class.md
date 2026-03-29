# AxisSet4D Class

A class representing a four-dimensional axis set, extending the
AxisSet3D class with an additional fourth axis.

## Slots

- `l`:

  A `NamedAxis` object representing the fourth axis.

## See also

[`AxisSet3D-class`](https://bbuchsbaum.github.io/neuroim2/reference/AxisSet3D-class.md),
[`NamedAxis-class`](https://bbuchsbaum.github.io/neuroim2/reference/NamedAxis-class.md)

## Examples

``` r
# Create NamedAxis objects for each dimension
x_axis <- new("NamedAxis", axis = "x", direction = 1)
y_axis <- new("NamedAxis", axis = "y", direction = 1)
z_axis <- new("NamedAxis", axis = "z", direction = 1)
t_axis <- new("NamedAxis", axis = "t", direction = 1)

# Create an AxisSet4D object
axis_set_4d <- new("AxisSet4D", i = x_axis, j = y_axis, k = z_axis,
                   l = t_axis, ndim = 4L)
```
