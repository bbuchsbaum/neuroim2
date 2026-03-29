# Construct a Minimal NIfTI-1 Header from a NeuroVol

Given a
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol.md)
object (or similar), this function builds a basic NIfTI-1 header
structure, populating essential fields such as `dim`, `pixdim`,
`datatype`, the affine transform, and the quaternion parameters.

## Usage

``` r
as_nifti_header(
  vol,
  file_name,
  oneFile = TRUE,
  data_type = "FLOAT",
  extensions = NULL
)
```

## Arguments

- vol:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol.md)
  (or 3D array-like) specifying dimensions, spacing, and affine
  transform.

- file_name:

  A character string for the file name (used within the header but not
  necessarily to write data).

- oneFile:

  Logical; if `TRUE`, sets the NIfTI magic to `"n+1"`, implying a
  single-file format (`.nii`). If `FALSE`, uses `"ni1"` (header+image).

- data_type:

  Character specifying the data representation, e.g. `"FLOAT"`,
  `"DOUBLE"`. The internal code picks an integer NIfTI code.

- extensions:

  Optional
  [`NiftiExtensionList-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtensionList-class.md)
  object or list of
  [`NiftiExtension-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtension-class.md)
  objects to include in the header.

## Value

A `list` representing the NIfTI-1 header fields, containing elements
like `dimensions`, `pixdim`, `datatype`, `qform`, `quaternion`, `qfac`,
`extensions`, etc. This can be passed to other functions that write or
manipulate the header.

## Details

This is a convenience function that calls
[`createNIfTIHeader`](https://bbuchsbaum.github.io/neuroim2/reference/createNIfTIHeader.md)
first, then updates the fields (dimensions, `pixdim`, orientation, etc.)
based on the `vol` argument. The voxel offset is set to 352 bytes (or
larger if extensions are provided), and the quaternion is derived from
the transform matrix via
[`matrixToQuatern`](https://bbuchsbaum.github.io/neuroim2/reference/matrixToQuatern.md).

Note: This function primarily sets up a minimal header suitable for
writing standard single-file NIfTI-1. If you need a more comprehensive
or advanced usage, consider manually editing the returned list.

## See also

[`createNIfTIHeader`](https://bbuchsbaum.github.io/neuroim2/reference/createNIfTIHeader.md)
for the base constructor of an empty NIfTI header.
[`NiftiExtension`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtension.md)
for creating extensions.
