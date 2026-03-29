# Laplacian Enhancement Filter for Volumetric Images

This function applies a multi-layer Laplacian enhancement filter to a
volumetric image (3D brain MRI data). The filter enhances details while
preserving edges using a non-local means approach with multiple scales.

## Usage

``` r
laplace_enhance(
  vol,
  mask,
  k = 2,
  patch_size = 3,
  search_radius = 2,
  h = 0.7,
  mapping_params = NULL,
  use_normalization_free = TRUE
)
```

## Arguments

- vol:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  object representing the image volume to be enhanced.

- mask:

  A
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
  object specifying the region to process. If not provided, the entire
  volume will be processed.

- k:

  An integer specifying the number of layers in the decomposition
  (default is 2).

- patch_size:

  An integer specifying the size of patches for non-local means. Must be
  odd (default is 3).

- search_radius:

  An integer specifying the radius of the search window (default is 2).

- h:

  A numeric value controlling the filtering strength. Higher values mean
  more smoothing (default is 0.7).

- mapping_params:

  An optional list of parameters for the enhancement mappings.

- use_normalization_free:

  Logical indicating whether to use normalization-free weights (default
  is TRUE).

## Value

A
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
object representing the enhanced image.
