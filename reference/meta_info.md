# Lightweight metadata for neuroimaging files

\`meta_info()\` provides a simple, CRAN-friendly way to retrieve
essential image metadata without teaching S4 details up front. It
accepts a file path or a \`FileMetaInfo\` object and returns a
normalized list containing common fields like dimensions, spacing,
origin, and transform.

The function does not read image data; it only parses header
information.

## Usage

``` r
meta_info(x)

# S4 method for class 'FileMetaInfo'
meta_info(x)

# S4 method for class 'character'
meta_info(x)
```

## Arguments

- x:

  A character file path (e.g., \`"image.nii.gz"\`) or an object of class
  [`FileMetaInfo`](https://bbuchsbaum.github.io/neuroim2/reference/FileMetaInfo-class.md).

## Value

A named list with the following elements:

- \`dim\` Integer vector of image dimensions.

- \`spacing\` Numeric voxel spacing (mm).

- \`origin\` Numeric coordinate origin.

- \`trans\` 4x4 transformation matrix mapping grid to world (mm).

- \`path\` Data file path.

- \`filename\` Basename of \`path\`.

- \`format\` File format label (e.g., "NIFTI", "AFNI").

- \`dtype\` Storage data type label.

- \`bytes_per_element\` Bytes per element.

- \`nvox\` Number of voxels in the spatial volume (prod of first 3
  dims).

- \`nvol\` Number of volumes (4th dim if present, else 1).

- \`size_bytes\` Approximate uncompressed size in bytes (\`nvox \* nvol
  \* bytes_per_element\`).

- \`time_step\` Time step (TR in seconds) if available for NIfTI, else
  \`NA_real\_\`.

## Details

Summarize Image Metadata

## See also

[`read_header`](https://bbuchsbaum.github.io/neuroim2/reference/read_header.md),
[`trans`](https://bbuchsbaum.github.io/neuroim2/reference/trans-methods.md),
[`FileMetaInfo`](https://bbuchsbaum.github.io/neuroim2/reference/FileMetaInfo-class.md),
[`NIFTIMetaInfo`](https://bbuchsbaum.github.io/neuroim2/reference/FileMetaInfo-class.md)

## Examples

``` r
# \donttest{
f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
mi <- meta_info(f)
mi$dim
#> [1] 64 64 25  4
mi$spacing
#> [1] 3.5 3.5 3.7
mi$origin
#> [1]  112.0 -108.0  -46.2
mi$filename
#> [1] "global_mask_v4.nii"
# 4x4 transform
mi$trans
#>      [,1] [,2] [,3]   [,4]
#> [1,] -3.5  0.0  0.0  112.0
#> [2,]  0.0  3.5  0.0 -108.0
#> [3,]  0.0  0.0  3.7  -46.2
#> [4,]  0.0  0.0  0.0    1.0
# }
```
