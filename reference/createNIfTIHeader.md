# Create an Empty NIfTI-1 Header List

Initializes a list of fields following the NIfTI-1 specification with
default or placeholder values. Users typically call this internally via
[`as_nifti_header`](https://bbuchsbaum.github.io/neuroim2/reference/as_nifti_header.md)
rather than using directly.

## Usage

``` r
createNIfTIHeader(oneFile = TRUE, file_name = NULL)
```

## Arguments

- oneFile:

  Logical; if `TRUE`, `magic` is set to `"n+1"` indicating a single-file
  (.nii) approach. Otherwise set to `"ni1"`.

- file_name:

  Optional character string to store in the header, usually referencing
  the intended output file name.

## Value

A named `list` containing approximately 30 fields that comprise the
NIfTI-1 header structure. Many of these are placeholders until filled by
downstream usage.

## Details

This function sets up the skeleton of a NIfTI-1 header, including fields
for `diminfo`, `pixdim`, `qform_code`, `magic`, etc. Most fields are
initialized to zero, empty characters, or standard placeholders. The
`oneFile` argument controls whether `"n+1"` or `"ni1"` is used for the
`magic` field.

## See also

[`as_nifti_header`](https://bbuchsbaum.github.io/neuroim2/reference/as_nifti_header.md)
for populating the returned header with actual data from a `NeuroVol`.
