# Create an exhaustive searchlight iterator

This function generates an exhaustive searchlight iterator that returns
either voxel coordinates or ROIVolWindow objects for each searchlight
sphere within the provided mask. The iterator visits every non-zero
voxel in the mask as a potential center voxel.

## Usage

``` r
searchlight(mask, radius, eager = FALSE, nonzero = FALSE, cores = 0)
```

## Arguments

- mask:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  object representing the brain mask.

- radius:

  A numeric value specifying the radius (in mm) of the spherical
  searchlight.

- eager:

  A logical value specifying whether to eagerly compute the searchlight
  ROIs. Default is FALSE, which uses lazy evaluation.

- nonzero:

  A logical value indicating whether to include only coordinates with
  nonzero values in the supplied mask. Default is FALSE.

- cores:

  An integer specifying the number of cores to use for parallel
  computation. Default is 0, which uses a single core.

## Value

A `deferred_list` object containing either matrices of integer-valued
voxel coordinates or
[`ROIVolWindow`](https://bbuchsbaum.github.io/neuroim2/reference/ROIVolWindow-class.md)
objects, each representing a searchlight region.

## Examples

``` r
# Load an example brain mask
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Generate an exhaustive searchlight iterator with a radius of 6 mm

searchlights <- searchlight(mask, radius = 6, eager = FALSE)

```
