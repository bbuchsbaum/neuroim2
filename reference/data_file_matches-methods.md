# Generic function to test whether a file name conforms to the given a [`FileFormat`](https://bbuchsbaum.github.io/neuroim2/reference/FileFormat-class.md) instance. Will test for match to data file only

Validates whether a file name conforms to the data file format
specification.

## Usage

``` r
data_file_matches(x, file_name)

# S4 method for class 'FileFormat,character'
data_file_matches(x, file_name)
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

## Details

The function performs case-sensitive pattern matching to verify that the
file name ends with the specified data extension. The match is performed
using a regular expression that ensures the extension appears at the end
of the file name.

## See also

[`file_matches`](https://bbuchsbaum.github.io/neuroim2/reference/file_matches-methods.md),
[`header_file_matches`](https://bbuchsbaum.github.io/neuroim2/reference/header_file_matches-methods.md)
for related file format validation

## Examples

``` r

fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
data_file_matches(fmt, "brain_scan.img")  # TRUE
#> [1] TRUE
data_file_matches(fmt, "brain_scan.hdr")  # FALSE
#> [1] FALSE
data_file_matches(fmt, "brain.img.gz")    # FALSE
#> [1] FALSE

```
