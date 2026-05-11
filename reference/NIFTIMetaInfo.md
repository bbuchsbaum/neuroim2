# Create NIFTI Format Metadata Object

Creates a NIFTIMetaInfo object containing format-specific metadata for
NIFTI format neuroimaging files.

## Usage

``` r
NIFTIMetaInfo(descriptor, nifti_header)
```

## Arguments

- descriptor:

  NIFTIFormat object specifying file format details

- nifti_header:

  List containing NIFTI header information

## Value

A NIFTIMetaInfo object

## Details

Create NIFTIMetaInfo Object

The NIFTIMetaInfo object extends MetaInfo with NIFTI-specific features:

- NIFTI header fields (qform, sform matrices)

- Data scaling (slope, intercept)

- File organization (separate vs. single file)

- Orientation information

Validation ensures:

- Valid NIFTI format

- Consistent dimensions

- Valid transformation matrices

- Proper data scaling

## See also

[`MetaInfo`](https://bbuchsbaum.github.io/neuroim2/reference/MetaInfo.md)

## Examples

``` r

# Read NIFTI header
header <- read_header(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Create format descriptor
fmt <- new("NIFTIFormat",
           file_format = "NIFTI",
           header_encoding = "raw",
           header_extension = "nii",
           data_encoding = "raw",
           data_extension = "nii")

# Create metadata
meta <- NIFTIMetaInfo(fmt, header@header)

# Check dimensions
dim(meta)
#> [1] 64 64 25  4

```
