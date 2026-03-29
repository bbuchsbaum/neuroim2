# Generic function to strip extension from file name, given a [`FileFormat`](https://bbuchsbaum.github.io/neuroim2/reference/FileFormat-class.md) instance.

Removes the file extension from a given file name based on the
FileFormat specifications.

## Usage

``` r
strip_extension(x, file_name)

# S4 method for class 'FileFormat,character'
strip_extension(x, file_name)
```

## Arguments

- x:

  A
  [FileFormat](https://bbuchsbaum.github.io/neuroim2/reference/FileFormat-class.md)
  object specifying the format requirements

- file_name:

  A character string specifying the file name to strip the extension
  from

## Value

A `character` string `file_name` without its extension.

A character string representing the file name without the extension

## Details

The function performs the following steps:

1.  If the file_name matches the header file format, it removes the
    header extension.

2.  If the file_name matches the data file format, it removes the data
    extension.

3.  If the file_name doesn't match either format, it throws an error.

## See also

[`header_file`](https://bbuchsbaum.github.io/neuroim2/reference/header_file-methods.md),
[`data_file`](https://bbuchsbaum.github.io/neuroim2/reference/data_file-methods.md)
for related file name manipulation

## Examples

``` r
# Create a FileFormat for NIFTI files
fmt <- new("FileFormat",
           header_extension = "nii",
           data_extension = "nii")

# Strip extension from a NIFTI file
strip_extension(fmt, "brain_scan.nii")  # Returns "brain_scan"
#> [1] "brain_scan"
```
