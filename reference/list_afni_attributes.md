# List AFNI Attributes in Extension

Returns a character vector of all attribute names in an AFNI extension.

## Usage

``` r
list_afni_attributes(ext)
```

## Arguments

- ext:

  A
  [`NiftiExtension-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtension-class.md)
  object with ecode = 4, or an xml_document from
  [`parse_afni_extension`](https://bbuchsbaum.github.io/neuroim2/reference/parse_afni_extension.md).

## Value

Character vector of attribute names.

## Examples

``` r
if (FALSE) { # \dontrun{
# List all attributes in an AFNI extension
attrs <- list_afni_attributes(afni_ext)
print(attrs)
} # }
```
