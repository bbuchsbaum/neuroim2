# A ggplot2 fill scale with neuroimaging-friendly defaults

A ggplot2 fill scale with neuroimaging-friendly defaults

## Usage

``` r
scale_fill_neuro(
  cmap = "grays",
  range = c("robust", "data"),
  probs = c(0.02, 0.98),
  limits = NULL,
  na.value = "transparent",
  guide = "colorbar"
)
```

## Arguments

- cmap:

  Palette name or vector of colors. See \[resolve_cmap()\].

- range:

  Either "robust" (quantiles) or "data" (min/max) to determine the
  default scale limits when \`limits\` is not provided.

- probs:

  Two-length numeric vector of quantiles for \`range="robust"\`.

- limits:

  Optional numeric limits (min, max). Overrides \`range\`.

- na.value:

  Color for NA.

- guide:

  Legend guide (default "colorbar").

## Value

A ggplot2 scale object.
