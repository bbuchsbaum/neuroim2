# Orthogonal three-plane view with optional crosshairs

Creates axial, coronal, and sagittal panels at a given coordinate with
harmonized aesthetics. Returns (invisibly) the three ggplot objects
after printing them in a single row using base grid (no extra deps).

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
  downsample = 1L
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

  "robust" or "data" for intensity limits shared by all panels.

- probs:

  Quantiles for robust range.

- crosshair:

  Logical; draw crosshair lines.

- annotate:

  Logical; add orientation glyphs.

- downsample:

  Integer decimation for speed.
