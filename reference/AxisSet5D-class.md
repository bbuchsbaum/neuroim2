# AxisSet5D Class

A class representing a five-dimensional axis set, extending the
AxisSet4D class with an additional fifth axis.

## Slots

- `m`:

  A `NamedAxis` object representing the fifth axis.

## See also

[`AxisSet4D-class`](https://bbuchsbaum.github.io/neuroim2/reference/AxisSet4D-class.md),
[`NamedAxis-class`](https://bbuchsbaum.github.io/neuroim2/reference/NamedAxis-class.md)

## Examples

``` r
# Create NamedAxis objects for each dimension
x_axis <- new("NamedAxis", axis = "x", direction = 1)
y_axis <- new("NamedAxis", axis = "y", direction = 1)
z_axis <- new("NamedAxis", axis = "z", direction = 1)
t_axis <- new("NamedAxis", axis = "t", direction = 1)
v_axis <- new("NamedAxis", axis = "v", direction = 1)

# Create an AxisSet5D object
axis_set_5d <- new("AxisSet5D", i = x_axis, j = y_axis, k = z_axis,
                   l = t_axis, m = v_axis, ndim = 5L)
```
