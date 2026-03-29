# Add L/R and A/P/S/I annotations (optional)

Add L/R and A/P/S/I annotations (optional)

## Usage

``` r
annotate_orientation(
  plane = c("axial", "coronal", "sagittal"),
  dims,
  gp = grid::gpar(col = "white", cex = 0.9, fontface = "bold")
)
```

## Arguments

- plane:

  "axial", "coronal", or "sagittal"

- dims:

  c(nrow, ncol) of the slice matrix

- gp:

  grid::gpar style

## Value

A ggplot2 layer with annotation_custom grobs
