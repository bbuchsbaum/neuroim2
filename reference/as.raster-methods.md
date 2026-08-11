# Convert neuroimaging objects to raster images

These methods convert `NeuroSlice` and `NeuroVol` objects into `raster`
images that can be displayed with
[`grid::grid.raster`](https://rdrr.io/r/grid/grid.raster.html) or base
`plot`.

## Usage

``` r
# S4 method for class 'NeuroSlice'
as.raster(
  x,
  cmap = gray(seq(0, 1, length.out = 255)),
  irange = range(x, na.rm = TRUE),
  ...
)

# S4 method for class 'NeuroVol'
as.raster(
  x,
  zlevel = ceiling(dim(x)[3]/2),
  cmap = gray(seq(0, 1, length.out = 255)),
  irange = range(x, na.rm = TRUE),
  ...
)
```

## Arguments

- x:

  object to convert.

- cmap:

  palette name or hex-color vector for the background (default
  `"grays"`). See
  [`resolve_cmap`](https://bbuchsbaum.github.io/neuroim2/reference/resolve_cmap.md).

- irange:

  numeric length-2 intensity range for the color scale.

- ...:

  additional arguments passed to methods

- zlevel:

  slice index for `NeuroVol` objects. Defaults to the middle slice.

## Value

A `raster` object

## See also

[`plot,NeuroSlice-method`](https://bbuchsbaum.github.io/neuroim2/reference/plot-methods.md),
[`plot,NeuroVol-method`](https://bbuchsbaum.github.io/neuroim2/reference/plot-methods.md)
