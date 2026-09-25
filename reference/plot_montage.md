# Montage of slices through a volume

Draws a grid of slices through one volume (typically a structural image)
in the same style as
[`plot_overlay`](https://bbuchsbaum.github.io/neuroim2/reference/plot_overlay.md):
black tiles cropped to the head, world-coordinate labels, L/R markers on
the first tile, and a grid that fills the canvas. Also accepts a list of
`NeuroSlice` objects or plain matrices.

## Usage

``` r
plot_montage(
  x,
  zlevels = NULL,
  along = 3L,
  cmap = "grays",
  range = c("robust", "data"),
  probs = c(0.02, 0.98),
  ncol = NULL,
  downsample = 1L,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  style = c("light", "dark", "report"),
  crop = TRUE,
  interpolate = TRUE,
  cbar_title = "value",
  colorbar = NULL,
  unit = c("index", "mm"),
  annotate = TRUE,
  n_slices = 12L,
  draw = FALSE,
  canvas = NULL
)
```

## Arguments

- x:

  A 3D volume, or a list of `NeuroSlice` objects / matrices.

- zlevels:

  Slices to plot when `x` is a volume: indices along `along`
  (`unit = "index"`, the default) or world coordinates (`unit = "mm"`).
  `NULL` (default) picks `n_slices` slices spread over the brain.

- along:

  Native voxel-grid axis along which to slice. For canonically ordered
  images, 1 = sagittal, 2 = coronal, and 3 = axial. Display orientation
  is inferred from the image affine.

- cmap:

  Palette name or vector (see \[resolve_cmap()\]).

- range:

  "robust" (quantile-based), "data" (min/max), or an explicit numeric
  `c(lo, hi)`.

- probs:

  Quantiles for \`range="robust"\`.

- ncol:

  Number of columns in the facet layout.

- downsample:

  Integer decimation for speed.

- title, subtitle, caption:

  Optional ggplot labels.

- style:

  Visual style: `"light"`, `"dark"`, or `"report"` (light card, dark
  cropped tiles, typography, and a colorbar – matching
  [`plot_overlay`](https://bbuchsbaum.github.io/neuroim2/reference/plot_overlay.md)'s
  report look).

- crop, interpolate:

  Logical; crop to the head bounding box (volume input) / smooth the
  raster. Both default to `TRUE`.

- cbar_title:

  Character; the quantity label drawn above the colorbar. Supplying it
  explicitly also turns the colorbar on.

- colorbar:

  Logical or `NULL`. `NULL` (default) shows a slim colorbar only for
  non-grayscale palettes or when `cbar_title` is supplied; arbitrary
  structural intensity units carry no information.

- unit:

  `"index"` (default) or `"mm"`: how `zlevels` is interpreted for volume
  input. Panels are always labelled in world coordinates (mm).

- annotate:

  Logical; draw L/R (or A/P) orientation letters on the first panel.

- n_slices:

  Number of slices chosen automatically when `zlevels` is `NULL`.
  Automatic slices are spread over the extent of the brain (bright
  tissue), skipping neck, scalp-only and empty planes.

- draw:

  Logical; if `TRUE`, also print the figure immediately (and return it
  invisibly). By default it is returned visibly, like any ggplot.

- canvas:

  Optional `c(width, height)` in inches to fit the layout to (default:
  the open device, else 10 x 7.5 in).

## Value

A figure (class `neuro_fig`, a patchwork wrapping one faceted ggplot
with one facet per slice), returned visibly; invisibly when
`draw = TRUE`. Its layout is re-fitted to the device it is drawn on, so
`ggsave()` at any size gives a filled, centred grid.

## See also

Other plot_neuro:
[`plot_checkerboard()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_checkerboard.md),
[`plot_edge_overlay()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_edge_overlay.md),
[`plot_ortho()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_ortho.md),
[`plot_overlay()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_overlay.md)

## Examples

``` r
# \donttest{
bg <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
p <- plot_montage(bg, title = "MNI152 (downsampled)")
ggplot2::ggsave(tempfile(fileext = ".png"), p, width = 8, height = 6)
# }
```
