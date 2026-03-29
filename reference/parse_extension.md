# Parse NIfTI Extension Data

Parses the raw data in a NIfTI extension based on its extension code.
Provides specialized parsing for known extension types.

## Usage

``` r
parse_extension(ext, ...)
```

## Arguments

- ext:

  A
  [`NiftiExtension-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtension-class.md)
  object.

- ...:

  Additional arguments passed to type-specific parsers.

## Value

Parsed data in an appropriate format:

- ecode 4 (AFNI): An XML document (if xml2 available) or character
  string

- ecode 6 (comment): Character string

- Other codes: Raw vector (unchanged)

## See also

[`parse_afni_extension`](https://bbuchsbaum.github.io/neuroim2/reference/parse_afni_extension.md)
for AFNI-specific parsing.

## Examples

``` r
# Parse a comment extension
ext <- NiftiExtension(ecode = 6L, data = "Test comment")
parse_extension(ext)  # Returns "Test comment"
#> [1] "Test comment"
```
