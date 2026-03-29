# NiftiExtension Class

Represents a single NIfTI header extension block. NIfTI extensions allow
additional metadata to be stored with the image file.

## Usage

``` r
# S4 method for class 'NiftiExtension'
show(object)
```

## Arguments

- object:

  A `NiftiExtension` object.

## Details

NIfTI-1.1 extensions follow this structure:

- Bytes 0-3: esize (int32) - total extension size, must be multiple of
  16

- Bytes 4-7: ecode (int32) - extension code identifying format

- Bytes 8-(esize-1): edata - the actual extension data

Extensions are chained sequentially after the NIfTI header (byte 352)
until the vox_offset is reached.

## Slots

- `ecode`:

  An `integer` extension code identifying the type of extension. See
  [`NiftiExtensionCodes`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtensionCodes.md)
  for known codes.

- `esize`:

  An `integer` giving the total size of the extension in bytes,
  including the 8-byte header (esize + ecode). Must be a multiple of 16.

- `edata`:

  A `raw` vector containing the extension data (length = esize - 8).

## See also

[`NiftiExtensionCodes`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtensionCodes.md)
for registered extension codes.
[`NiftiExtensionList-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtensionList-class.md)
for a collection of extensions.
[`parse_extension`](https://bbuchsbaum.github.io/neuroim2/reference/parse_extension.md)
for parsing extension data.

## Examples

``` r
# Create a simple comment extension
comment_text <- "This is a test comment"
ext <- NiftiExtension(ecode = 6L, data = comment_text)

# Access the extension code
ext@ecode
#> [1] 6
```
