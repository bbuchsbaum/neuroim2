# read header information of an image file

read header information of an image file

## Usage

``` r
read_header(file_name)
```

## Arguments

- file_name:

  the name of the file to read

## Value

an instance of class
[`FileMetaInfo`](https://bbuchsbaum.github.io/neuroim2/reference/FileMetaInfo-class.md)

## Examples

``` r
hdr <- read_header(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
dim(hdr)                  # image dimensions
#> [1] 64 64 25  4
hdr@header$pixdim[5]      # TR in seconds
#> [1] 0
```
