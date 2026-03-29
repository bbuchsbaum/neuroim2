# Neuroimaging color palettes and helpers

Lightweight, perceptually-uniform color tools with safe fallbacks.

## Usage

``` r
resolve_cmap(name = "grays", n = 256)
```

## Arguments

- name:

  Palette name (e.g., "grays", "viridis", "inferno", "magma", "plasma",
  "turbo", "cividis"). Case-insensitive. If you pass a vector of colors,
  it's returned unchanged.

- n:

  Number of colors to generate.

## Value

A character vector of hex colors.
