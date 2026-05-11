# Generic function to get the name of the data file, given a file name and a [`FileFormat`](https://bbuchsbaum.github.io/neuroim2/reference/FileFormat-class.md) instance.

Derives the data file name from a given file name based on the
FileFormat specifications.

## Usage

``` r
data_file(x, file_name)

# S4 method for class 'FileFormat,character'
data_file(x, file_name)
```

## Arguments

- x:

  A
  [FileFormat](https://bbuchsbaum.github.io/neuroim2/reference/FileFormat-class.md)
  object specifying the format requirements

- file_name:

  A character string specifying the file name to derive the data file
  name from

## Value

The correct data file name as a `character` string.

A character string representing the data file name

## Details

The function performs the following steps:

1.  If the input file_name already matches the data file format, it
    returns the file_name as is.

2.  If the file_name matches the header file format, it constructs and
    returns the corresponding data file name.

3.  If the file_name doesn't match either format, it throws an error.

## See also

[`header_file`](https://bbuchsbaum.github.io/neuroim2/reference/header_file-methods.md),
[`strip_extension`](https://bbuchsbaum.github.io/neuroim2/reference/strip_extension-methods.md)
for related file name manipulation

## Examples

``` r

fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
data_file(fmt, "brain_scan.img")  # Returns "brain_scan.img"
#> [1] "brain_scan.img"
data_file(fmt, "brain_scan.hdr")  # Also Returns "brain_scan.img"
#> [1] "brain_scan.img"

```
