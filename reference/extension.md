# Get Extension by Code

Retrieve extensions with a specific extension code from a list.

## Usage

``` r
extension(x, ecode)

# S4 method for class 'NiftiExtensionList,numeric'
extension(x, ecode)
```

## Arguments

- x:

  A
  [`NiftiExtensionList-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtensionList-class.md)
  object.

- ecode:

  Integer extension code to filter by.

## Value

A
[`NiftiExtensionList-class`](https://bbuchsbaum.github.io/neuroim2/reference/NiftiExtensionList-class.md)
containing only extensions with the specified code.
