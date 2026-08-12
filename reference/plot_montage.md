# Plot a montage of axial (or any-plane) slices using facetting

This avoids extra dependencies by using a single ggplot with facets and
a shared colorbar. Supply a list of slice objects or a volume + indices.

## Usage

``` r
plot_montage(
  x,
  zlevels = NULL,
  along = 3L,
  cmap = "grays",
  range = c("robust", "data"),
  probs = c(0.02, 0.98),
  ncol = 6L,
  downsample = 1L,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  style = c("light", "dark", "report"),
  crop = NULL,
  interpolate = NULL
)
```

## Arguments

- x:

  Either a 3D volume object accepted by \`slice()\` or a list of slices.

- zlevels:

  Integer indices of slices to plot (if \`x\` is a volume).

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

  Logical or `NULL`; crop to the brain bounding box / smooth the raster.
  `NULL` (default) enables both for `style = "report"` only. (Cropping
  applies to the volume path.)
