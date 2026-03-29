# Apply a 4D bilateral filter to a NeuroVec

This function applies a full 4D bilateral filter to a `NeuroVec`,
smoothing jointly across space (x, y, z) and time (t). The filter uses
spatial, temporal, and intensity kernels to preserve edges while
reducing noise, leveraging a parallel C++ backend for performance.

## Usage

``` r
bilateral_filter_4d(
  vec,
  mask,
  spatial_sigma = 2,
  intensity_sigma = 1,
  temporal_sigma = 1,
  spatial_window = 1,
  temporal_window = 1,
  temporal_spacing = 1
)
```

## Arguments

- vec:

  A
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
  object (4D image).

- mask:

  An optional
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
  or
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  specifying the spatial region to process. If omitted, the entire
  spatial extent is processed.

- spatial_sigma:

  Numeric; standard deviation of the spatial Gaussian (default 2).

- intensity_sigma:

  Numeric; standard deviation of the intensity Gaussian (default 1).

- temporal_sigma:

  Numeric; standard deviation of the temporal Gaussian (default 1).

- spatial_window:

  Integer; half-width of the spatial window in voxels (default 1), e.g.,
  1 =\> 3x3x3 spatial neighborhood.

- temporal_window:

  Integer; half-width of the temporal window in frames (default 1),
  e.g., 1 =\> 3 timepoints (t-1, t, t+1).

- temporal_spacing:

  Numeric; spacing of the temporal dimension (e.g., TR in seconds).
  Default is 1. This sets the temporal scale used for the temporal
  kernel.

## Value

A
[`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
with filtered data.

## Details

Parameter guidance and units: - spatial_sigma: Measured in physical
units (millimeters). Distances are computed using `spacing(vec)[1:3]`,
so choose `spatial_sigma` relative to voxel size. As a rule of thumb,
set it to about 1-2 voxel sizes (e.g., 2-4 mm for 2 mm isotropic data)
for moderate smoothing. - intensity_sigma: Dimensionless multiplier of
the global intensity standard deviation. Internally, the filter uses
exp(-(dI)^2 / (2 \* (intensity_sigma \* sigma_I)^2)), where sigma_I is
the standard deviation of all finite voxel intensities within the mask
across time. Start with 1.0 for moderate smoothing; use 0.5-0.8 to
preserve more edges, or 1.5-2.0 for stronger smoothing. -
temporal_sigma: Measured in `temporal_spacing` units (e.g., seconds).
Typical values are 0.5-2 x TR. Larger values blend more across time.

Choosing the neighborhood window sizes: - spatial_window controls the
discrete spatial support. A common choice is
`ceiling(2 * spatial_sigma / min(spacing(vec)[1:3]))`, which covers
~95 - temporal_window similarly can be set to
`ceiling(2 * temporal_sigma / temporal_spacing)`.

Quick presets (typical fMRI with 2-3 mm voxels and TR~2s): - Light:
spatial_sigma = 1 x min(spacing), intensity_sigma = 0.8, temporal_sigma
= 0.5 x TR, windows = 1 - Moderate (default-ish): spatial_sigma = 1.5 x
min(spacing), intensity_sigma = 1.0, temporal_sigma = 1 x TR, windows =
1-2 - Strong: spatial_sigma = 2 x min(spacing), intensity_sigma = 1.5,
temporal_sigma = 1.5 x TR, windows = 2

Tip: If your time axis has known TR, pass it via `temporal_spacing`. For
NIfTI inputs, you can get TR via:

      hdr <- read_header(nifti_path)
      tr  <- hdr@header$pixdim[5]
      out <- bilateral_filter_4d(vec, mask, temporal_spacing = tr)

## See also

[`bilateral_filter`](https://bbuchsbaum.github.io/neuroim2/reference/bilateral_filter.md),
[`NeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md),
[`NeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)

## Examples

``` r
# \donttest{
vec <- read_vec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
out  <- bilateral_filter_4d(vec, mask,
                            spatial_sigma = 2, intensity_sigma = 1,
                            temporal_sigma = 1, spatial_window = 1,
                            temporal_window = 1, temporal_spacing = 1)
# }
```
