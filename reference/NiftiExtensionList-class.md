# NiftiExtensionList Class

A validated list containing zero or more
[`NiftiExtension-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtension-class.md)
objects. This class ensures type safety when working with collections of
NIfTI extensions.

## Usage

``` r
# S4 method for class 'NiftiExtensionList'
show(object)
```

## Arguments

- object:

  A `NiftiExtensionList` object.

## Details

The class extends `list` and enforces that all elements must be
`NiftiExtension` objects. This provides a clean container for managing
multiple extensions attached to a NIfTI file.

## See also

[`NiftiExtension-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtension-class.md)
for individual extension objects.
[`extensions`](https://bbuchsbaum.github.io/neuroim2/reference/extensions.md)
for accessing extensions from image objects.

## Examples

``` r
# Create an empty extension list
ext_list <- new("NiftiExtensionList")

# Create a list with extensions
ext1 <- NiftiExtension(ecode = 6L, data = "Comment 1")
ext2 <- NiftiExtension(ecode = 6L, data = "Comment 2")
ext_list <- new("NiftiExtensionList", list(ext1, ext2))
```
