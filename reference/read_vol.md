# Load an image volume from a file

Load an image volume from a file

## Usage

``` r
read_vol(file_name, index = 1)
```

## Arguments

- file_name:

  the name of the file to load

- index:

  the index of the volume (e.g. if the file is 4-dimensional)

## Value

an instance of the class
[`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md)

## Examples

``` r
fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
x <- read_vol(fname)
print(dim(x))
#> [1] 64 64 25
space(x)
#> <NeuroSpace> [3D] 
#> ── Geometry ──────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Spacing       : 3.5 x 3.5 x 3.7 mm
#>   Origin        : 112, -108, -46.2
#>   Orientation   : LAS
#>   Voxels        : 102,400
```
