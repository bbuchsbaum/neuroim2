# Generic function to test whether a file name conforms to the given [`FileFormat`](https://bbuchsbaum.github.io/neuroim2/reference/FileFormat-class.md) instance. Will test for match to header file only

Validates whether a file name conforms to the header file format
specification.

## Usage

``` r
header_file_matches(x, file_name)

# S4 method for class 'FileFormat,character'
header_file_matches(x, file_name)
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

A logical value: `TRUE` if the file name matches the header format,
`FALSE` otherwise

## Details

The function performs case-sensitive pattern matching to verify that the
file name ends with the specified header extension. The match is
performed using a regular expression that ensures the extension appears
at the end of the file name.

## See also

[`file_matches`](https://bbuchsbaum.github.io/neuroim2/reference/file_matches-methods.md),
[`data_file_matches`](https://bbuchsbaum.github.io/neuroim2/reference/data_file_matches-methods.md)
for related file format validation

## Examples

``` r

fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
header_file_matches(fmt, "brain_scan.hdr")  # TRUE
#> [1] TRUE
header_file_matches(fmt, "brain_scan.img")  # FALSE
#> [1] FALSE
header_file_matches(fmt, "brain.hdr.gz")    # FALSE
#> [1] FALSE

```
