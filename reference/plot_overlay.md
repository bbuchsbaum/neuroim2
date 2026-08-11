# Composite an overlay map on a structural background

Works without extra packages by colorizing both layers to rasters and
stacking them as grobs. Great for statistical maps over T1/T2
backgrounds.

## Usage

``` r
plot_overlay(
  bgvol,
  overlay,
  zlevels = NULL,
  along = 3L,
  bg_cmap = "grays",
  ov_cmap = "inferno",
  bg_range = c("robust", "data"),
  ov_range = c("robust", "data"),
  probs = c(0.02, 0.98),
  ov_thresh = 0,
  ov_alpha = 0.7,
  ov_alpha_mode = c("binary", "proportional", "ramp", "soft"),
  ov_symmetric = NULL,
  alpha_gamma = NULL,
  ov_cap = NULL,
  ncol = 3L,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  draw = TRUE,
  style = c("light", "dark", "report"),
  enhance = FALSE,
  assemble = TRUE,
  colorbar = TRUE,
  legend = NULL,
  crop = NULL,
  interpolate = NULL
)
```

## Arguments

- bgvol:

  Background 3D volume.

- overlay:

  Overlay 3D volume on the same NeuroSpace grid as \`bgvol\`.

- zlevels:

  Slices to plot (indices along the z/3rd axis by default).

- along:

  Axis for slicing (1 sagittal, 2 coronal, 3 axial).

- bg_cmap:

  Background palette (e.g., "grays").

- ov_cmap:

  Overlay palette (e.g., "inferno").

- bg_range, ov_range:

  Background/overlay scaling. Either a mode string, `"robust"` (quantile
  clip) or `"data"` (min/max), or an explicit numeric `c(lo, hi)` to pin
  the scale (e.g. `ov_range = c(-6, 6)`) for consistent coloring across
  panels and subjects.

- probs:

  Quantiles for robust scaling.

- ov_thresh:

  Numeric threshold; values with \|v\| \< thresh become transparent.

- ov_alpha:

  Global alpha for overlay (0..1).

- ov_alpha_mode:

  One of `"binary"` (default: pixels above threshold get full
  `ov_alpha`, others transparent), `"proportional"` (per-pixel alpha =
  \|v\| / cap), `"ramp"` (alpha ramps linearly from 0 at `ov_thresh` to
  1 at the cap), or `"soft"` (a nonlinear, self-tuning curve,
  `alpha = clamp((|v|-lo)/(hi-lo),0,1)^gamma`, where the knee `lo`
  defaults to `ov_thresh` or, if unset, the median in-mask magnitude,
  and `gamma` adapts to the value distribution so that opacity rises
  rapidly away from zero and the noisy bulk stays faint). The cap is
  shared across all panels so identical values get identical opacity
  everywhere.

- ov_symmetric:

  Logical or `NULL`. `NULL` (default) auto-selects symmetric limits
  around zero when the overlay has both positive and negative values;
  `TRUE`/`FALSE` forces the choice. Symmetric limits keep negative and
  positive values equally visible with a diverging palette.

- alpha_gamma:

  Optional exponent for `ov_alpha_mode = "soft"`. `NULL` (default)
  auto-tunes it from the data; larger values push more of the low-value
  range toward transparency.

- ov_cap:

  Optional numeric; the magnitude used as the upper end of the
  (symmetric) color/alpha scale. Defaults to the data-driven limit.

- ncol:

  Number of columns in the panel layout.

- title, subtitle, caption:

  Optional layout-level labels used when drawing.

- draw:

  Logical; if \`TRUE\`, draw on the active graphics device. If
  \`FALSE\`, return without drawing.

- style:

  Visual style, either `"light"` or `"dark"`.

- enhance:

  Display-only enhancement of the (unsmoothed) statistical `overlay`.
  `FALSE` (default) leaves it untouched; `TRUE` applies
  [`enhance_stat_map`](https://bbuchsbaum.github.io/neuroim2/reference/enhance_stat_map.md)
  with defaults; a named `list` is forwarded as arguments to
  [`enhance_stat_map()`](https://bbuchsbaum.github.io/neuroim2/reference/enhance_stat_map.md)
  (e.g. `enhance = list(detail_gain = 2, method = "bilateral")`).

- assemble:

  Logical; if `TRUE` (default), return a single assembled patchwork
  object (honoring `ncol` and the layout labels) suitable for
  `ggsave()`. If `FALSE`, draw a panel grid and return the per-slice
  ggplot list invisibly (useful for programmatic access to individual
  panels).

- colorbar:

  Logical; when `assemble = TRUE`, append a colorbar for the overlay
  statistic (with the threshold marked). Default `TRUE`.

- legend:

  Logical or `NULL`; when `assemble = TRUE`, add a bottom legend strip
  (positive/negative swatches, threshold, plane). `NULL` (default) shows
  it for `style = "report"` only.

- crop:

  Logical or `NULL`; crop every panel to the brain bounding box (shared
  across slices, so framing is consistent). `NULL` (default) crops for
  `style = "report"` only.

- interpolate:

  Logical or `NULL`; smooth the background raster. `NULL` (default)
  interpolates for `style = "report"` only.

## Details

**Return value.** By default (`assemble = TRUE`) the return value is a
single patchwork object that can be passed directly to `ggsave()`. With
`assemble = FALSE` the montage is drawn as a side effect and the return
value is the *list of per-slice ggplots* (invisibly); passing that list
to `ggsave()` saves only one panel.

**Signed maps.** For overlays with both signs (t/z/contrast maps), the
default palette switches to a diverging one and limits become symmetric
so negatives are as visible as positives; pass `ov_cmap`/`ov_symmetric`
to override.

**Report style.** `style = "report"` renders dark brain tiles on a light
card with bold/italic typography, a titled colorbar, a bottom legend
strip, brain-bbox cropping, and a smoothed background – a
publication-ready look. The individual features (`legend`, `crop`,
`interpolate`) can also be toggled on any style.
