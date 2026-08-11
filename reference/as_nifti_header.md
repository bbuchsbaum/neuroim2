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
  data_type = NULL,
  extensions = NULL,
  values = NULL,
  version = NULL
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
  `"DOUBLE"`, `"SHORT"`. `NULL` (the default) keeps the source file's
  datatype when `vol` came from one and its values are still exactly
  representable in it, and uses `"FLOAT"` otherwise.

- extensions:

  Optional
  [`NiftiExtensionList-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtensionList-class.md)
  object or list of
  [`NiftiExtension-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtension-class.md)
  objects to include in the header.

- values:

  Optional numeric vector of the voxel values that will be written. Used
  to pick a datatype and to derive `scl_slope` / `scl_inter` for integer
  output; taken from `vol` when omitted.

- version:

  Either `1` or `2`, selecting the NIfTI version to describe. `NULL`
  (the default) uses NIfTI-1 unless a dimension exceeds what its 16-bit
  `dim` field can hold, in which case NIfTI-2 is used because NIfTI-1
  physically cannot represent the image.

## Value

A `list` representing the NIfTI header fields, containing elements like
`dimensions`, `pixdim`, `datatype`, `qform`, `quaternion`, `qfac`,
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

When `vol` was read from a NIfTI file, the fields that describe the
acquisition rather than the array – repetition time (`pixdim[5]`),
`xyzt_units`, `qform_code` and `sform_code`, `descrip`, `aux_file`, the
intent fields, `cal_min`/`cal_max`, the slice timing fields and
`toffset` – are carried over from that file. Only the fields the object
itself determines (geometry, dimensions, datatype, offsets) are
recomputed. Before this behaviour existed, a read-write round trip reset
all of them, which silently dropped the TR and relabelled an MNI-space
image as scanner-space.

## See also

[`createNIfTIHeader`](https://bbuchsbaum.github.io/neuroim2/reference/createNIfTIHeader.md)
for the base constructor of an empty NIfTI header.
[`NiftiExtension`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtension.md)
for creating extensions.
