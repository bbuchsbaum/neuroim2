# Neuroimaging color palettes and helpers

Lightweight, perceptually-uniform color tools with safe fallbacks.

## Usage

``` r
resolve_cmap(name = "grays", n = 256)
```

## Arguments

- name:

  Palette name (e.g., "grays", "viridis", "inferno", "magma", "plasma",
  "turbo", "cividis", "coldhot", "blue-red", "coolwarm"). Any palette in
  [`grDevices::hcl.pals()`](https://rdrr.io/r/grDevices/palettes.html)
  (e.g. "RdBu", "Spectral", "Reds") is also accepted. Case- and
  punctuation-insensitive. Unknown names emit a warning and fall back to
  a viridis-like ramp. If you pass a vector of colors, it's returned
  unchanged.

- n:

  Number of colors to generate.

## Value

A character vector of hex colors.
