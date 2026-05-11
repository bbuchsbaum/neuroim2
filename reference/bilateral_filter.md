# Apply a bilateral filter to a volumetric image

This function smooths a volumetric image (3D brain MRI data) using a
bilateral filter. The bilateral filter considers both spatial closeness
and intensity similarity for smoothing. Only in-mask, in-bounds
neighbors contribute to each local weighted average.

## Usage

``` r
bilateral_filter(
  vol,
  mask,
  spatial_sigma = 2,
  intensity_sigma = 1,
  window = 1,
  range_scale = NULL
)
```

## Arguments

- vol:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  object representing the image volume to be smoothed.

- mask:

  An optional
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
  object representing the image mask that defines the region where the
  filtering is applied. If not provided, the entire volume is
  considered.

- spatial_sigma:

  A numeric value specifying the standard deviation of the spatial
  Gaussian kernel (default is 2).

- intensity_sigma:

  A numeric value specifying the standard deviation of the intensity
  Gaussian kernel (default is 1).

- window:

  An integer specifying the number of voxels around the center voxel to
  include on each side. For example, window=1 for a 3x3x3 kernel
  (default is 1).

- range_scale:

  Optional positive numeric range scale used by the intensity kernel. If
  `NULL`, the scale is estimated as the standard deviation of the
  current input values inside `mask`. Supply a fixed value to apply the
  same range bandwidth across observed and null maps.

## Value

A smoothed image of class
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md).

## Examples

``` r
brain_mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Apply bilateral filtering to the brain volume
filtered_vol <- bilateral_filter(brain_mask, brain_mask, spatial_sigma = 2,
intensity_sigma = 25, window = 1)
```
