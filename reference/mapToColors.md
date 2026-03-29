# Map intensity values to colors

Convert intensity values (e.g., a 2D slice) into a color representation
for plotting and overlays.

## Usage

``` r
mapToColors(
  imslice,
  col = heat.colors(128, alpha = 1),
  zero_col = "#00000000",
  alpha = 1,
  irange = range(imslice),
  threshold = c(0, 0)
)
```

## Arguments

- imslice:

  A numeric vector or array of intensities.

- col:

  A vector of colors used as a lookup table.

- zero_col:

  Color used for exactly-zero intensities (defaults to transparent).

- alpha:

  Global alpha multiplier applied to all colors when `alpha < 1`.

- irange:

  Intensity range used to normalize values before mapping to `col`.

- threshold:

  Optional length-2 numeric vector. If `diff(threshold) > 0`, values
  within `[threshold[1], threshold[2]]` are set to transparent.

## Value

If `alpha == 1`, returns a character vector/array of colors. If
`alpha < 1`, returns an array with an added RGBA channel (last dimension
length 4).
