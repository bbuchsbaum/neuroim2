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
  caption = NULL
)
```

## Arguments

- x:

  Either a 3D volume object accepted by \`slice()\` or a list of slices.

- zlevels:

  Integer indices of slices to plot (if \`x\` is a volume).

- along:

  Axis along which to slice (1 = sagittal, 2 = coronal, 3 = axial).

- cmap:

  Palette name or vector (see \[resolve_cmap()\]).

- range:

  "robust" (quantile-based) or "data" (min/max).

- probs:

  Quantiles for \`range="robust"\`.

- ncol:

  Number of columns in the facet layout.

- downsample:

  Integer decimation for speed.

- title, subtitle, caption:

  Optional ggplot labels.
