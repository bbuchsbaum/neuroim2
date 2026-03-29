# Get Indices from an IndexLookupVol Object

Retrieves the vector of indices that are included in the lookup volume.

## Usage

``` r
# S4 method for class 'IndexLookupVol'
indices(x)

# S4 method for class 'ROIVol'
indices(x)

# S4 method for class 'ROIVol'
indices(x)

# S4 method for class 'ROIVec'
indices(x)

# S4 method for class 'AbstractSparseNeuroVec'
indices(x)
```

## Arguments

- x:

  An
  [`IndexLookupVol`](https://bbuchsbaum.github.io/neuroim2/reference/IndexLookupVol-class.md)
  object

## Value

the indices of the lookup volume

## Examples

``` r
space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
ilv <- IndexLookupVol(space, c(1:100))
idx <- indices(ilv)  # Get included indices

```
