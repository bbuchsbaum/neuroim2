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
  tile = NULL,
  cmap = "grays",
  bg_range = c("robust", "data"),
  ov_range = c("robust", "data"),
  probs = c(0.02, 0.98),
  ncol = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  draw = FALSE,
  style = c("light", "dark", "report"),
  labels = c("fixed", "moving"),
  legend = TRUE,
  crop = TRUE,
  mask_background = TRUE,
  unit = c("index", "mm"),
  annotate = TRUE,
  assemble = TRUE,
  n_slices = 6L,
  match_intensity = TRUE,
  canvas = NULL,
  interpolate = TRUE,
  focus_brain = TRUE
)
```

## Arguments

- bgvol:

  Background/reference 3D volume.

- overlay:

  Comparison 3D volume on the same NeuroSpace grid as \`bgvol\`.

- zlevels:

  Slices to plot: indices along \`along\` (`unit = "index"`) or world
  coordinates (`unit = "mm"`). `NULL` (default) picks `n_slices` slices
  spread over the brain.

- along:

  Native voxel-grid axis for slicing. Display orientation is inferred
  from the image affine.

- tile:

  Tile width in voxels (the key states the equivalent in mm). `NULL`
  (default) picks a clean size in mm giving about ten tiles across the
  head.

- cmap:

  Palette used to render the normalized checkerboard image.

- bg_range, ov_range:

  Intensity scaling of each image: `"robust"` (computed over head
  voxels), `"data"`, or numeric `c(lo, hi)`. Each image is windowed
  independently so tissue contrast, not a brightness step, is what
  differs between tiles.

- probs:

  Quantiles for robust scaling.

- ncol:

  Number of columns (`NULL` = chosen to fill the canvas).

- title, subtitle, caption:

  Optional layout-level labels used when drawing.

- draw:

  Logical; if `TRUE`, also print the figure immediately (and return it
  invisibly). By default the figure is returned visibly.

- style:

  Visual style: `"light"`, `"dark"` or `"report"`.

- labels:

  Length-2 character vector naming the two images in the key.

- legend:

  Logical; draw the key under the panels.

- crop:

  Logical; crop panels to the head bounding box.

- mask_background:

  Logical; show the checker pattern only inside the head (union of both
  images' foreground), leaving air black.

- unit:

  `"index"` or `"mm"`: how `zlevels` is interpreted. Panels are labelled
  in world coordinates.

- annotate:

  Logical; draw orientation letters on the first panel.

- assemble:

  Logical; return one patchwork figure (default) or the list of panel
  ggplots.

- n_slices:

  Number of automatically chosen slices.

- match_intensity:

  Logical; map the comparison image's intensities onto the reference
  image's distribution (quantile matching over head voxels) before
  interleaving, so a brightness or contrast difference between the
  images does not masquerade as misalignment. Default `TRUE`.

- canvas:

  Optional `c(width, height)` in inches to fit the layout to (default:
  the open device).

- interpolate:

  Logical; smooth the rendered checkerboard (default `TRUE`), which also
  softens the head outline.

- focus_brain:

  Logical; interleave tiles only inside the (dilated) bright-tissue mask
  of `bgvol` and show the fixed image alone on the scalp, so the seams
  where registration is judged dominate the figure instead of the
  stair-stepped head outline. Default `TRUE`.

## Value

A figure (class `neuro_fig`, a patchwork whose layout is re-fitted to
the device it is drawn on) when `assemble = TRUE` or a named list of
panel ggplots (`assemble = FALSE`); invisibly when `draw = TRUE`.

## See also

Other plot_neuro:
[`plot_edge_overlay()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_edge_overlay.md),
[`plot_montage()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_montage.md),
[`plot_ortho()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_ortho.md),
[`plot_overlay()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_overlay.md)

## Examples

``` r
# \donttest{
fixed <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
moving <- fixed
moving[] <- sqrt(fixed[c(2:dim(fixed)[1], 1), , ])   # shifted, different contrast
p <- plot_checkerboard(fixed, moving)
ggplot2::ggsave(tempfile(fileext = ".png"), p, width = 8, height = 6)
# }
```
