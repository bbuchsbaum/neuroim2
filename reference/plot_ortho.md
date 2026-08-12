# Orthogonal three-plane view with optional crosshairs

Creates axial, coronal, and sagittal panels at a given coordinate with
harmonized aesthetics. Returns the three ggplot objects invisibly after
drawing, or without drawing when `draw = FALSE`.

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
  draw = TRUE,
  style = c("light", "dark", "report"),
  enhance = FALSE,
  crop = NULL,
  interpolate = NULL
)
```

## Arguments

- vol:

  A 3D volume handled by \`slice()\`.

- coord:

  Length-3 coordinate of the target point. Interpreted as voxel indices
  by default; set \`unit = "mm"\` to convert using \`coord_to_grid()\`
  if available in your environment.

- unit:

  "index" or "mm".

- cmap:

  Palette for the slices.

- range:

  Intensity limits shared by all panels: `"robust"`, `"data"`, or an
  explicit numeric `c(lo, hi)`.

- probs:

  Quantiles for robust range.

- crosshair:

  Logical; draw crosshair lines.

- annotate:

  Logical; add orientation glyphs.

- downsample:

  Integer decimation for speed.

- title, subtitle, caption:

  Optional layout-level labels used when drawing.

- draw:

  Logical; if \`TRUE\`, draw the panels on the active graphics device.
  If \`FALSE\`, only return the ggplot objects invisibly.

- style:

  Visual style: `"light"`, `"dark"`, or `"report"` (light card, dark
  cropped tiles, typography, and a colorbar – matching
  [`plot_overlay`](https://bbuchsbaum.github.io/neuroim2/reference/plot_overlay.md)'s
  report look).

- enhance:

  Display-only enhancement of an unsmoothed statistical `vol`. `FALSE`
  (default) leaves it untouched; `TRUE` applies
  [`enhance_stat_map`](https://bbuchsbaum.github.io/neuroim2/reference/enhance_stat_map.md)
  with defaults; a named `list` is forwarded as arguments to
  [`enhance_stat_map()`](https://bbuchsbaum.github.io/neuroim2/reference/enhance_stat_map.md).

- crop, interpolate:

  Logical or `NULL`; crop panels to the brain bounding box / smooth the
  raster. `NULL` (default) enables both for `style = "report"` only.

## Details

The affine determines which native voxel axis is nearest each anatomical
plane and how that plane must be permuted or flipped for display.
Oblique images are shown on their regular native voxel planes; values
are not silently resampled. Use
[`deoblique()`](https://bbuchsbaum.github.io/neuroim2/reference/deoblique.md)
or
[`resample_to()`](https://bbuchsbaum.github.io/neuroim2/reference/resample_to.md)
first when true cardinal-plane sections are required.
