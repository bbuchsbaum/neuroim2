# Create AFNI Format Metadata Object

Creates an AFNIMetaInfo object containing format-specific metadata for
AFNI format neuroimaging files.

## Usage

``` r
AFNIMetaInfo(descriptor, afni_header)
```

## Arguments

- descriptor:

  AFNIFormat object specifying file format details

- afni_header:

  List containing AFNI header information

## Value

An AFNIMetaInfo object

## Details

Create AFNIMetaInfo Object

The AFNIMetaInfo object extends MetaInfo with AFNI-specific features:

- AFNI brick structure

- Sub-brick labels and scaling

- Space transformation

- Statistical parameters

The function handles:

- Dimension extraction and validation

- Label generation for sub-bricks

- Transformation from AFNI to NIFTI space

- Data type and scaling setup
