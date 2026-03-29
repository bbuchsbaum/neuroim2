# Create a spherical random searchlight iterator

This function generates a spherical random searchlight iterator for
analyzing local neighborhoods of voxels within a given radius in a brain
mask.

## Usage

``` r
random_searchlight(mask, radius, nonzero = TRUE)
```

## Arguments

- mask:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  object representing the brain mask.

- radius:

  A numeric value specifying the radius of the searchlight sphere in
  voxel units.

- nonzero:

  Logical; if `TRUE` (default) discard zero-valued voxels in the mask
  when forming each searchlight.

## Value

A list of
[`ROIVolWindow`](https://bbuchsbaum.github.io/neuroim2/reference/ROIVolWindow-class.md)
objects, each representing a spherical searchlight region.

## Examples

``` r
# Create a simple brain mask
mask_data <- array(TRUE, c(10, 10, 10))
mask_data[1, 1, 1] <- FALSE
mask <- LogicalNeuroVol(mask_data, NeuroSpace(c(10,10,10)))

# Generate random searchlight iterator with a radius of 2 voxels
# \donttest{
searchlights <- random_searchlight(mask, radius = 6)
# }
```
