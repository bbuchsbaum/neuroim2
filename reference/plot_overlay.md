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
  ov_alpha_mode = c("binary", "proportional"),
  ncol = 3L,
  title = NULL,
  subtitle = NULL,
  caption = NULL
)
```

## Arguments

- bgvol:

  Background 3D volume.

- overlay:

  Overlay 3D volume (same dims as bgvol).

- zlevels:

  Slices to plot (indices along the z/3rd axis by default).

- along:

  Axis for slicing (1 sagittal, 2 coronal, 3 axial).

- bg_cmap:

  Background palette (e.g., "grays").

- ov_cmap:

  Overlay palette (e.g., "inferno").

- bg_range, ov_range:

  "robust" or "data" for background/overlay scaling.

- probs:

  Quantiles for robust scaling.

- ov_thresh:

  Numeric threshold; values with \|v\| \< thresh become transparent.

- ov_alpha:

  Global alpha for overlay (0..1).

- ov_alpha_mode:

  Either `"binary"` (default, current behavior: pixels above threshold
  get full `ov_alpha`, others are transparent) or `"proportional"`
  (per-pixel alpha scaled by absolute value of the overlay, multiplied
  by `ov_alpha`).

- ncol:

  Number of columns in the facet layout.

- title, subtitle, caption:

  Optional labels.
