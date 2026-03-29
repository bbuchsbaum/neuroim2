# FileFormat Class

This class represents a neuroimaging file format descriptor, containing
information about the file format, encoding, and extensions for both
header and data components.

## Slots

- `file_format`:

  A `character` string specifying the name of the file format (e.g.,
  "NIfTI").

- `header_encoding`:

  A `character` string specifying the file encoding of the header file
  (e.g., "raw" for binary, "gzip" for gz compressed).

- `header_extension`:

  A `character` string specifying the file extension for the header file
  (e.g., "nii" for NIfTI single files).

- `data_encoding`:

  A `character` string specifying the file encoding for the data file.

- `data_extension`:

  A `character` string specifying the file extension for the data file
  (e.g., "nii" for NIfTI single files).

## Examples

``` r
# Create a FileFormat object for NIfTI format
nifti_format <- new("FileFormat",
                    file_format = "NIfTI",
                    header_encoding = "raw",
                    header_extension = "nii",
                    data_encoding = "raw",
                    data_extension = "nii")

```
