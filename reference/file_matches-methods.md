# Generic function to test whether a file name conforms to the given [`FileFormat`](https://bbuchsbaum.github.io/neuroim2/reference/FileFormat-class.md) instance. Will test for match to either header file or data file

Validates whether a file name conforms to the specified FileFormat and
verifies the existence of both header and data files.

## Usage

``` r
file_matches(x, file_name)

# S4 method for class 'FileFormat,character'
file_matches(x, file_name)
```

## Arguments

- x:

  A
  [FileFormat](https://bbuchsbaum.github.io/neuroim2/reference/FileFormat-class.md)
  object specifying the format requirements

- file_name:

  A character string specifying the file name to validate

## Value

`TRUE` for match, `FALSE` otherwise.

A logical value: `TRUE` if the file matches the format and both header
and data files exist, `FALSE` otherwise

## Details

The function performs the following validation steps:

1.  Checks if the file name matches either the header or data format

2.  Verifies the existence of the corresponding paired file

3.  Returns `FALSE` if either check fails

File names are validated using case-sensitive extension matching.

## See also

[`header_file_matches`](https://bbuchsbaum.github.io/neuroim2/reference/header_file_matches-methods.md),
[`data_file_matches`](https://bbuchsbaum.github.io/neuroim2/reference/data_file_matches-methods.md)
for individual file type checking

## Examples

``` r
# Create a FileFormat for NIFTI format
# \donttest{
fmt <- new("FileFormat",
  file_format = "NIFTI",
  header_encoding = "raw",
  header_extension = "nii",
  data_encoding = "raw",
  data_extension = "nii")

# Create temporary file
tmp <- tempfile("brainscan", fileext = ".nii")
file.create(tmp)
#> [1] TRUE

# Check if files exist and match format
file_matches(fmt, tmp)
#> [1] TRUE

# Clean up
unlink(tmp)
# }
```
