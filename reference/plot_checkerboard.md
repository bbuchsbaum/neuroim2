# Checkerboard comparison of two registered volumes

Alternates tiles from a background volume and a comparison volume on
matched slices. This is useful for visual registration QC.

## Usage

``` r
plot_checkerboard(
  bgvol,
  overlay,
  zlevels = NULL,
  along = 3L,
  tile = 16L,
  cmap = "grays",
  bg_range = c("robust", "data"),
  ov_range = c("robust", "data"),
  probs = c(0.02, 0.98),
  ncol = 3L,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  draw = TRUE,
  style = c("light", "dark")
)
```

## Arguments

- bgvol:

  Background/reference 3D volume.

- overlay:

  Comparison 3D volume on the same NeuroSpace grid as \`bgvol\`.

- zlevels:

  Slices to plot. Defaults to nine evenly spaced slices.

- along:

  Native voxel-grid axis for slicing. Display orientation is inferred
  from the image affine.

- tile:

  Tile width in slice pixels.

- cmap:

  Palette used to render the normalized checkerboard image.

- bg_range, ov_range:

  "robust" or "data" intensity scaling.

- probs:

  Quantiles for robust scaling.

- ncol:

  Number of columns in the facet layout.

- title, subtitle, caption:

  Optional layout-level labels used when drawing.

- draw:

  Logical; if \`TRUE\`, draw the panels on the active graphics device.
  If \`FALSE\`, only return the ggplot objects invisibly.

- style:

  Visual style, either `"light"` or `"dark"`.
