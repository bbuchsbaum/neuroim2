# Generic function to read image meta info given a file

Reads meta information from image files based on their format (NIFTI or
AFNI).

## Usage

``` r
read_meta_info(x, file_name)

# S4 method for class 'NIFTIFormat'
read_meta_info(x, file_name)

# S4 method for class 'AFNIFormat'
read_meta_info(x, file_name)
```

## Arguments

- x:

  A
  [FileFormat](https://bbuchsbaum.github.io/neuroim2/reference/FileFormat-class.md)
  object (either NIFTIFormat or AFNIFormat)

- file_name:

  A character string specifying the file name to read meta information
  from

## Value

A `list` containing the meta information read from the file.

An object of class
[NIFTIMetaInfo](https://bbuchsbaum.github.io/neuroim2/reference/FileMetaInfo-class.md)
or
[AFNIMetaInfo](https://bbuchsbaum.github.io/neuroim2/reference/FileMetaInfo-class.md),
depending on the input format

## Details

These methods use format-specific functions to read the header
information and create the appropriate meta information object. The
\`.read_meta_info\` helper function is used internally to streamline the
process for both formats.

## Examples

``` r
# Create a NIFTI format descriptor
fmt <- new("NIFTIFormat",
           file_format = "NIFTI",
           header_encoding = "raw",
           header_extension = "nii",
           data_encoding = "raw",
           data_extension = "nii")

# Read metadata from a NIFTI file
# \donttest{
fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
meta <- read_meta_info(fmt, fname)

# Access metadata properties
dim(meta)          # Image dimensions
#> [1] 64 64 25  4
trans(meta)        # Transformation matrix
#>      [,1] [,2] [,3]   [,4]
#> [1,] -3.5  0.0  0.0  112.0
#> [2,]  0.0  3.5  0.0 -108.0
#> [3,]  0.0  0.0  3.7  -46.2
#> [4,]  0.0  0.0  0.0    1.0
# }
```
