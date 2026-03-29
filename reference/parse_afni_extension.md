# Parse AFNI Extension

Parses an AFNI extension (ecode = 4) containing XML-formatted
attributes.

## Usage

``` r
parse_afni_extension(ext, as_xml = TRUE)
```

## Arguments

- ext:

  A
  [`NiftiExtension-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtension-class.md)
  object with ecode = 4.

- as_xml:

  Logical; if TRUE (default) and xml2 is available, returns an
  xml_document object. Otherwise returns the raw XML string.

## Value

If `as_xml = TRUE` and xml2 is available, returns an xml_document.
Otherwise returns a character string containing the XML.

## Details

AFNI stores dataset attributes in an XML format within the NIfTI
extension. The XML contains elements like HISTORY_NOTE, volume labels,
tagged points, and other AFNI-specific metadata.

## See also

[`get_afni_attribute`](https://bbuchsbaum.github.io/neuroim2/reference/get_afni_attribute.md)
for extracting specific AFNI attributes.

## Examples

``` r
if (FALSE) { # \dontrun{
# Read a NIfTI file with AFNI extension
hdr <- read_nifti_header("afni_file.nii")
afni_ext <- hdr$extensions[[1]]
parsed <- parse_afni_extension(afni_ext)
} # }
```
