# Create a NIfTI Extension

Constructor function for creating a
[`NiftiExtension-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtension-class.md)
object with proper padding to ensure the size is a multiple of 16 bytes.

## Usage

``` r
NiftiExtension(ecode, data)
```

## Arguments

- ecode:

  Integer extension code. See
  [`NiftiExtensionCodes`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtensionCodes.md)
  for known codes. Common values: 4 (AFNI), 6 (comment), 32 (CIFTI).

- data:

  The extension data. Can be:

  - A character string (will be converted to raw with null terminator)

  - A raw vector (used as-is)

## Value

A
[`NiftiExtension-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtension-class.md)
object.

## Details

The function automatically handles padding to ensure the total extension
size (esize) is a multiple of 16 bytes, as required by the NIfTI
specification. The esize includes the 8-byte header (esize + ecode
fields).

## See also

[`NiftiExtension-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtension-class.md),
[`NiftiExtensionCodes`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtensionCodes.md)

## Examples

``` r
# Create a comment extension
ext <- NiftiExtension(ecode = 6L, data = "This is a comment")
ext@ecode
#> [1] 6
ext@esize
#> [1] 32

# Create an AFNI extension with XML data
afni_xml <- '<?xml version="1.0"?><AFNI_attributes></AFNI_attributes>'
afni_ext <- NiftiExtension(ecode = 4L, data = afni_xml)
```
