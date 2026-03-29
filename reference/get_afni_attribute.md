# Get AFNI Attribute from Extension

Extracts a specific attribute value from a parsed AFNI extension.

## Usage

``` r
get_afni_attribute(ext, attr_name)
```

## Arguments

- ext:

  A
  [`NiftiExtension-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtension-class.md)
  object with ecode = 4, or an xml_document from
  [`parse_afni_extension`](https://bbuchsbaum.github.io/neuroim2/reference/parse_afni_extension.md).

- attr_name:

  Character string specifying the attribute name to retrieve (e.g.,
  "HISTORY_NOTE", "BRICK_LABS").

## Value

The attribute value, or NULL if not found. The type depends on the
attribute's ni_type in the XML.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get the history note from an AFNI extension
history <- get_afni_attribute(afni_ext, "HISTORY_NOTE")
} # }
```
