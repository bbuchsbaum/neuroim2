# Overlay fixed and moving edge maps on a background volume

Displays a structural/reference background with two edge channels
rendered in distinct colors. This is intended for registration QC where
fixed/template edges and moving/result edges should coincide.

## Usage

``` r
plot_edge_overlay(
  bgvol,
  fixed_edges,
  moving_edges,
  zlevels = NULL,
  along = 3L,
  bg_cmap = "grays",
  fixed_color = "#00d5ff",
  moving_color = "#ff3b30",
  bg_range = c("robust", "data"),
  edge_range = c("robust", "data"),
  probs = c(0.02, 0.98),
  edge_thresh = 0,
  edge_alpha = 0.85,
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

  Background 3D volume.

- fixed_edges:

  Edge map for the fixed/reference image on the same NeuroSpace grid as
  \`bgvol\`.

- moving_edges:

  Edge map for the moving/result image on the same NeuroSpace grid as
  \`bgvol\`.

- zlevels:

  Slices to plot. Defaults to nine evenly spaced slices.

- along:

  Native voxel-grid axis for slicing. Display orientation is inferred
  from the image affine.

- bg_cmap:

  Background palette.

- fixed_color, moving_color:

  Overlay colors for the two edge maps.

- bg_range, edge_range:

  "robust" or "data" intensity scaling.

- probs:

  Quantiles for robust scaling.

- edge_thresh:

  Values below this edge magnitude are transparent.

- edge_alpha:

  Global alpha for edge overlays.

- ncol:

  Number of columns in the facet layout.

- title, subtitle, caption:

  Optional layout-level labels used when drawing.

- draw:

  Logical; if \`TRUE\`, draw the panels on the active graphics device.
  If \`FALSE\`, only return the ggplot objects invisibly.

- style:

  Visual style, either `"light"` or `"dark"`.
