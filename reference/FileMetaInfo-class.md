# FileMetaInfo Class

This class extends MetaInfo to include file-specific metadata for
neuroimaging data files.

This class extends FileMetaInfo with NIfTI-specific metadata.

This class extends FileMetaInfo with AFNI-specific metadata.

## Slots

- `header_file`:

  A `character` string specifying the name of the file containing meta
  information.

- `data_file`:

  A `character` string specifying the name of the file containing image
  data.

- `descriptor`:

  A
  [`FileFormat`](https://bbuchsbaum.github.io/neuroim2/reference/FileFormat-class.md)
  object describing the image file format.

- `endian`:

  A `character` string specifying the byte order of data ('little' or
  'big').

- `data_offset`:

  A `numeric` value indicating the number of bytes preceding the start
  of image data in the data file.

- `bytes_per_element`:

  An `integer` specifying the number of bytes per data element.

- `intercept`:

  A `numeric` vector of constant values added to image data (one per
  sub-image).

- `slope`:

  A `numeric` vector of multipliers for image data (one per sub-image).

- `header`:

  A `list` of format-specific attributes.

- `nifti_header`:

  A `list` of attributes specific to the NIfTI file format.

- `afni_header`:

  A `list` of attributes specific to the AFNI file format.

## See also

[`MetaInfo-class`](https://bbuchsbaum.github.io/neuroim2/reference/MetaInfo-class.md),
`NIFTIMetaInfo-class`, `AFNIMetaInfo-class`

`FileMetaInfo-class`

`FileMetaInfo-class`
