# Composite a statistical map on a structural background

Draws a grid of slices through a structural background (e.g. a T1) with
a thresholded statistical map on top, in the style of a journal figure:
black tiles cropped to the head, world-coordinate slice labels, L/R
markers, and a compact colorbar that marks the threshold.

## Usage

``` r
plot_overlay(
  bgvol,
  overlay,
  zlevels = NULL,
  along = 3L,
  bg_cmap = "grays",
  ov_cmap = NULL,
  bg_range = c("robust", "data"),
  ov_range = c("robust", "data"),
  probs = c(0.02, 0.98),
  ov_thresh = 0,
  ov_alpha = 1,
  ov_alpha_mode = c("binary", "proportional", "ramp", "soft"),
  ov_symmetric = NULL,
  alpha_gamma = NULL,
  ov_cap = NULL,
  ncol = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  draw = FALSE,
  style = c("light", "dark", "report"),
  enhance = FALSE,
  assemble = TRUE,
  colorbar = TRUE,
  legend = NULL,
  crop = TRUE,
  interpolate = TRUE,
  cbar_title = "value",
  unit = c("index", "mm"),
  annotate = TRUE,
  n_slices = 12L,
  canvas = NULL,
  alpha_knee = NULL,
  alpha_cap = NULL,
  alpha_floor = NULL,
  alpha_mid = 0.2,
  gamma_min = 1.5,
  gamma_max = 5
)
```

## Arguments

- bgvol:

  Background 3D volume.

- overlay:

  Overlay 3D volume on the same NeuroSpace grid as \`bgvol\`.

- zlevels:

  Slices to plot, as indices along \`along\` (`unit = "index"`, the
  default) or world coordinates in mm (`unit = "mm"`). `NULL` (default)
  chooses `n_slices` slices spread over the extent of the
  supra-threshold overlay.

- along:

  Native voxel-grid axis for slicing. Display orientation and anatomical
  plane labels are inferred from the image affine.

- bg_cmap:

  Background palette (e.g., "grays").

- ov_cmap:

  Overlay palette. `NULL` (default) chooses automatically: the two-sided
  `"cold_hot"` map for signed data, `"hot"` otherwise. Any name accepted
  by \[resolve_cmap()\] or a vector of colours.

- bg_range, ov_range:

  Background/overlay scaling. Either a mode string, `"robust"` or
  `"data"`, or an explicit numeric `c(lo, hi)` to pin the scale (e.g.
  `ov_range = c(-6, 6)`) for consistent colouring across figures and
  subjects. The robust background window is computed over head voxels
  only; the robust overlay scale is computed once from the whole overlay
  volume (supra-threshold values when a threshold is set) and rounded to
  two significant digits, so a given map always gets the same scale
  whichever slices are shown.

- probs:

  Quantiles for robust scaling.

- ov_thresh:

  Numeric threshold; values with \|v\| \< thresh are not drawn.

- ov_alpha:

  Global opacity of the overlay (0..1).

- ov_alpha_mode:

  One of `"binary"` (default: every supra-threshold voxel fully opaque),
  `"proportional"`, `"ramp"`, or `"soft"` (opacity rises nonlinearly
  with magnitude; the curve `alpha = floor + (1 - floor) * t^gamma`
  self-tunes `gamma` from the data). In the graded modes every voxel
  that passes a threshold keeps at least 60% opacity, and the colorbar
  is faded with the same curve so it matches the picture.

- ov_symmetric:

  Logical or `NULL`. `NULL` (default) uses symmetric limits around zero
  when the overlay has both signs.

- alpha_gamma:

  Optional exponent for `ov_alpha_mode = "soft"`. `NULL` (default)
  auto-tunes it from the data.

- ov_cap:

  Optional numeric; the magnitude at the upper end of the colour/opacity
  scale. Defaults to the data-driven limit.

- ncol:

  Number of columns. `NULL` (default) picks the layout that best fills
  the canvas.

- title, subtitle, caption:

  Optional figure labels, left-aligned with the tiles.

- draw:

  Logical; if `TRUE`, also print the figure immediately (and return it
  invisibly). By default the figure is returned visibly, like a ggplot,
  so it prints at the console or in a knitr chunk.

- style:

  Visual style: `"light"` (white card), `"report"` (warm off-white card
  with the key strip on), or `"dark"` (black card). Tiles are black in
  every style.

- enhance:

  Display-only enhancement of the (unsmoothed) statistical `overlay`.
  `FALSE` (default) leaves it untouched; `TRUE` applies
  [`enhance_stat_map`](https://bbuchsbaum.github.io/neuroim2/reference/enhance_stat_map.md)
  with defaults; a named `list` is forwarded as arguments to
  [`enhance_stat_map()`](https://bbuchsbaum.github.io/neuroim2/reference/enhance_stat_map.md).

- assemble:

  Logical; if `TRUE` (default), return one assembled patchwork figure.
  If `FALSE`, return the list of per-slice ggplots.

- colorbar:

  Logical; draw the colorbar (default `TRUE`).

- legend:

  Logical or `NULL`; add a one-line key under the tiles (plane,
  neurological convention, and what the threshold shows). `NULL`
  (default) shows it for `style = "report"` only.

- crop:

  Logical; crop every panel to the head bounding box (shared across
  slices, and always containing every supra-threshold voxel).

- interpolate:

  Logical; smooth the background raster (default `TRUE`). The overlay
  itself is always drawn with crisp voxels.

- cbar_title:

  Character; the quantity label drawn above the colorbar. Defaults to
  `"value"`; set it to the statistic actually shown (e.g. `"t"`,
  `"Semipartial r"`).

- unit:

  `"index"` or `"mm"`: how `zlevels` is interpreted. Panels are always
  labelled in world coordinates (mm).

- annotate:

  Logical; draw L/R orientation letters on the first panel.

- n_slices:

  Number of slices chosen when `zlevels` is `NULL`.

- canvas:

  Optional `c(width, height)` in inches to freeze the layout for one
  size. By default (`NULL`) the layout is re-fitted to whatever device
  the figure is drawn on, including `ggsave()`.

- alpha_knee, alpha_cap:

  Optional lower (non-negative) and upper (positive) magnitude anchors
  of the soft opacity curve, independent of the colour limits. `NULL`
  uses the threshold (or median magnitude) and the colour-scale cap. Fix
  both and `alpha_gamma` to reproduce one opacity mapping across
  figures; see
  [`soft_alpha_params`](https://bbuchsbaum.github.io/neuroim2/reference/soft_alpha_params.md).

- alpha_floor:

  Minimum soft opacity (0–1) above the knee, before `ov_alpha`; values
  below `ov_thresh` stay transparent. `NULL` (default) uses 0.6 when a
  threshold is set, otherwise 0.

- alpha_mid, gamma_min, gamma_max:

  Auto-gamma policy passed to
  [`soft_alpha_params`](https://bbuchsbaum.github.io/neuroim2/reference/soft_alpha_params.md)
  (0.2 at the median supra-knee magnitude, gamma clamped to \[1.5, 5\]).

## Value

A figure (class `neuro_fig`, a patchwork whose layout is re-fitted to
the device it is drawn on) when `assemble = TRUE` or a list of ggplots
(`assemble = FALSE`); invisibly when `draw = TRUE`.

## Details

**Signed maps.** For overlays with both signs (t/z/contrast maps), the
default palette is two-sided and the limits symmetric, so negative
values are as visible as positive ones. With a threshold, the colour
ramp starts at a saturated colour at \\\pm\\threshold and brightens
toward the cap; the colorbar shows the sub-threshold band in a neutral
tone and ticks the threshold and cap. When the data exceed the cap, the
end tick reads \\\ge\\ cap.

**Soft opacity.** In `ov_alpha_mode = "soft"` the resolved curve is
recorded in `attr(result, "soft_alpha")` for either return form; pass it
back through `alpha_knee`, `alpha_cap`, `alpha_gamma` and `alpha_floor`
to reuse it.

**Saving.** The figure's layout is fitted to the device it is drawn on:
`p <- plot_overlay(...); ggsave("fig.png", p, width = 6, height = 9)`
re-arranges the tiles for a 6 x 9 in page, and additions such as
`p + patchwork::plot_annotation(title = "...")` are kept. Pass `canvas`
only to freeze the layout for one size.

## See also

Other plot_neuro:
[`plot_checkerboard()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_checkerboard.md),
[`plot_edge_overlay()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_edge_overlay.md),
[`plot_montage()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_montage.md),
[`plot_ortho()`](https://bbuchsbaum.github.io/neuroim2/reference/plot_ortho.md)

## Examples

``` r
# \donttest{
bg <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
stat <- bg
stat[] <- rnorm(length(stat)) * (bg[] > stats::quantile(bg[], 0.6))
p <- plot_overlay(bg, stat, ov_thresh = 1.5, cbar_title = "z")
tf <- tempfile(fileext = ".png")
ggplot2::ggsave(tf, p, width = 7, height = 5)
# }
```
