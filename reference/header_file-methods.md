# Generic function to get the name of the header file, given a file name and a [`FileFormat`](https://bbuchsbaum.github.io/neuroim2/reference/FileFormat-class.md) instance.

Derives the header file name from a given file name based on the
FileFormat specifications.

## Usage

``` r
header_file(x, file_name)

# S4 method for class 'FileFormat,character'
header_file(x, file_name)
```

## Arguments

- x:

  A
  [FileFormat](https://bbuchsbaum.github.io/neuroim2/reference/FileFormat-class.md)
  object specifying the format requirements

- file_name:

  A character string specifying the file name to derive the header file
  name from

## Value

The correct header file name as a `character` string.

A character string representing the header file name

## Details

The function performs the following steps:

1.  If the input file_name already matches the header file format, it
    returns the file_name as is.

2.  If the file_name matches the data file format, it constructs and
    returns the corresponding header file name.

3.  If the file_name doesn't match either format, it throws an error.

## See also

[`data_file`](https://bbuchsbaum.github.io/neuroim2/reference/data_file-methods.md),
[`strip_extension`](https://bbuchsbaum.github.io/neuroim2/reference/strip_extension-methods.md)
for related file name manipulation

## Examples

``` r

fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
header_file(fmt, "brain_scan.hdr")  # Returns "brain_scan.hdr"
#> [1] "brain_scan.hdr"
header_file(fmt, "brain_scan.img")  # Returns "brain_scan.hdr"
#> [1] "brain_scan.hdr"

```
