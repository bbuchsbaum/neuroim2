# Orthogonal three-plane view with optional crosshairs and overlay

Draws sagittal, coronal and axial sections through one point, optionally
with a thresholded statistical map on top (the classic "stat map" view).
All three views share one physical scale and one head-bounding-box crop,
so the crosshair lines up across views.

## Usage

``` r
plot_ortho(
  vol,
  coord = NULL,
  unit = c("index", "mm"),
  cmap = "grays",
  range = c("robust", "data"),
  probs = c(0.02, 0.98),
  crosshair = TRUE,
  annotate = TRUE,
  downsample = 1L,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  draw = FALSE,
  style = c("light", "dark", "report"),
  enhance = FALSE,
  crop = TRUE,
  interpolate = TRUE,
  cbar_title = "value",
  colorbar = NULL,
  assemble = TRUE,
  overlay = NULL,
  ov_thresh = 0,
  ov_cmap = NULL,
  ov_range = c("robust", "data"),
  ov_alpha = 1,
  ov_alpha_mode = c("binary", "proportional", "ramp", "soft"),
  ov_symmetric = NULL,
  ov_cap = NULL,
  canvas = NULL
)
```

## Arguments

- vol:

  A 3D background volume.

- coord:

  Length-3 coordinate of the target point: voxel indices
  (`unit = "index"`, default) or world coordinates in mm
  (`unit = "mm"`). `NULL` (default) uses the peak absolute value of
  `overlay` when one is given, otherwise the volume centre.

- unit:

  `"index"` or `"mm"`: how `coord` is interpreted.

- cmap:

  Palette for the background.

- range:

  Background intensity limits shared by all panels: `"robust"` (computed
  over head voxels), `"data"`, or numeric `c(lo, hi)`.

- probs:

  Quantiles for robust scaling.

- crosshair:

  Logical; draw the (gapped) crosshair.

- annotate:

  Logical; draw orientation letters on every view.

- downsample:

  Integer decimation for speed.

- title, subtitle, caption:

  Optional figure labels.

- draw:

  Logical; if `TRUE`, also print the figure immediately (and return it
  invisibly). By default the figure is returned visibly.

- style:

  Visual style: `"light"`, `"dark"`, or `"report"`.

- enhance:

  Display-only enhancement of an unsmoothed statistical `vol`; see
  [`plot_overlay`](https://bbuchsbaum.github.io/neuroim2/reference/plot_overlay.md).

- crop:

  Logical; crop views to the head bounding box.

- interpolate:

  Logical; smooth the background raster (default `TRUE`).

- cbar_title:

  Character; the quantity label drawn above the colorbar. Supplying it
  explicitly also turns the colorbar on.

- colorbar:

  Logical or `NULL`. `NULL` (default) shows a colorbar when it carries
  information: an `overlay` is given, the background uses a
  non-grayscale palette, or `cbar_title` is supplied.

- assemble:

  Logical; if `TRUE` (default) return one assembled patchwork figure; if
  `FALSE` return the named list of the `axial`, `coronal` and `sagittal`
  ggplots.

- overlay:

  Optional 3D statistical volume on the same grid as `vol`, drawn over
  all three views.

- ov_thresh, ov_cmap, ov_range, ov_alpha, ov_alpha_mode, ov_symmetric,
  ov_cap:

  Overlay threshold, palette, scaling, opacity and opacity mode;
  identical in meaning to the same arguments of
  [`plot_overlay`](https://bbuchsbaum.github.io/neuroim2/reference/plot_overlay.md).

- canvas:

  Optional `c(width, height)` in inches to fit the layout to (default:
  the open device).

## Value

A figure (class `neuro_fig`, a patchwork whose layout is re-fitted to
the device it is drawn on) when `assemble = TRUE` or a named list of
ggplots (`assemble = FALSE`); invisibly when `draw = TRUE`.

## Details

The affine determines which native voxel axis is nearest each anatomical
plane and how that plane must be permuted or flipped for display.
Oblique images are shown on their regular native voxel planes; values
are not silently resampled. Use
[`deoblique()`](https://bbuchsbaum.github.io/neuroim2/reference/deoblique.md)
or
[`resample_to()`](https://bbuchsbaum.github.io/neuroim2/reference/resample_to.md)
first when true cardinal-plane sections are required.

Each view is labelled with the world coordinate of its plane (mm).

## See also

Other plot_neuro:
[`plot_checkerboard()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_checkerboard.md),
[`plot_edge_overlay()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_edge_overlay.md),
[`plot_montage()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_montage.md),
[`plot_overlay()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_overlay.md)

## Examples

``` r
# \donttest{
bg <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
p <- plot_ortho(bg, coord = c(24, 26, 26))
#> ℹ `coord` are read as voxel indices; panels are labelled in mm.
#>   Pass `unit = "mm"` to give positions in world coordinates.
#> This message is displayed once every 8 hours.
ggplot2::ggsave(tempfile(fileext = ".png"), p, width = 9, height = 3.5)
# }
```
