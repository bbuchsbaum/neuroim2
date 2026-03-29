# Plot a NeuroSlice

Display axial slices of a
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
as a faceted montage.

## Usage

``` r
# S4 method for class 'NeuroSlice,ANY'
plot(
  x,
  cmap = gray(seq(0, 1, length.out = 255)),
  irange = range(x, na.rm = TRUE),
  legend = TRUE
)

# S4 method for class 'NeuroVol,missing'
plot(
  x,
  y,
  cmap = "grays",
  zlevels = unique(round(seq(1, dim(x)[3], length.out = 9))),
  irange = range(x, na.rm = TRUE),
  thresh = c(0, 0),
  alpha = 1,
  legend = TRUE
)

# S4 method for class 'NeuroVol,NeuroVol'
plot(
  x,
  y,
  cmap = "grays",
  zlevels = unique(round(seq(1, dim(x)[3], length.out = 9))),
  ov_cmap = "inferno",
  ov_alpha = 0.5,
  ov_thresh = 0
)
```

## Arguments

- x:

  the background volume to display.

- cmap:

  palette name or hex-color vector for the background (default
  `"grays"`). See
  [`resolve_cmap`](https://bbuchsbaum.github.io/neuroim2/reference/resolve_cmap.md).

- irange:

  numeric length-2 intensity range for the color scale.

- legend:

  logical; show the colour bar?

- y:

  optional overlay volume (same dimensions as `x`). When supplied, the
  plot is rendered as a background + overlay composite.

- zlevels:

  integer slice indices to display. Default: 9 evenly-spaced slices (3
  \\\times\\ 3 grid).

- thresh:

  a 2-element vector indicating the lower and upper transparency
  thresholds.

- alpha:

  opacity for the background layer (0–1).

- ov_cmap:

  overlay palette name (default `"inferno"`).

- ov_alpha:

  overlay opacity (default 0.5).

- ov_thresh:

  overlay threshold; values with \\\|v\| \< \\ `ov_thresh` become
  transparent (default 0).

## Value

a ggplot2 object

## Details

The plot method uses `ggplot2` to create a raster visualization of the
slice data. The intensity values are mapped to colors using the
specified colormap and range.

when \`x\` is a NeuroSlice object, the plot method returns a `ggplot2`
object containing the raster visualization of the slice data. The plot
can be further customized using standard ggplot2 functions.

When a second volume `y` is supplied it is treated as an overlay (e.g.\\
a statistical map) composited on top of `x` with adjustable
transparency. This delegates to
[`plot_overlay`](https://bbuchsbaum.github.io/neuroim2/reference/plot_overlay.md).

## Examples

``` r
# Create example slice
slice_space <- NeuroSpace(c(100, 100))
slice_data <- matrix(rnorm(100*100), 100, 100)
slice <- NeuroSlice(slice_data, slice_space)
# \donttest{
# Basic plot
plot(slice)

# }


dat <- matrix(rnorm(100*100), 100, 100)
slice <- NeuroSlice(dat, NeuroSpace(c(100,100)))
# \donttest{
plot(slice)

# }
```
