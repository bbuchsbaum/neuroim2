# Load a single 3D image volume from a file

Reads exactly one 3D volume from a neuroimaging file and returns it as a
[`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md).
Accepts both 3D files (where only `index = 1` is valid) and 4D files
(where `index` selects a single sub-volume along the 4th dimension).

## Usage

``` r
read_vol(file_name, index = 1)
```

## Arguments

- file_name:

  Path to a single image file (NIfTI `.nii` or `.nii.gz`). A character
  vector of length \> 1 is not supported — use
  [`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md)
  if you need to read multiple files, or call `read_vol` in a loop.

- index:

  Integer giving the index of the sub-volume to load. Must be `1` for a
  3D file. For a 4D file, must satisfy `1 <= index <= dim(file)[4]`.

## Value

A
[`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md)
(always 3D, always dense). The associated
[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
has three spatial dimensions even when the source file is 4D.

## See also

[`read_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md)
for loading 4D data as a
[`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md),
[`read_image`](https://bbuchsbaum.github.io/neuroim2/reference/read_image.md)
for automatic dimensionality-based dispatch, and
[`read_hyper_vec`](https://bbuchsbaum.github.io/neuroim2/reference/read_hyper_vec.md)
for 5D data.

## Examples

``` r
# Read the first volume from a 4D file
fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
x <- read_vol(fname)
print(dim(x))    # 3D
#> [1] 64 64 25
space(x)
#> <NeuroSpace> [3D] 
#> ── Geometry ──────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Spacing       : 3.5 x 3.5 x 3.7 mm
#>   Origin        : 112, -108, -46.2
#>   Orientation   : LAS
#>   Voxels        : 102,400

# Read the 3rd sub-volume from the same 4D file
x3 <- read_vol(fname, index = 3)
```
