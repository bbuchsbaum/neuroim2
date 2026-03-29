# Lookup Values in an IndexLookupVol Object

Performs a lookup operation on an IndexLookupVol object.

## Usage

``` r
# S4 method for class 'IndexLookupVol,numeric'
lookup(x, i)

# S4 method for class 'AbstractSparseNeuroVec,numeric'
lookup(x, i)
```

## Arguments

- x:

  An
  [`IndexLookupVol`](https://bbuchsbaum.github.io/neuroim2/reference/IndexLookupVol-class.md)
  object

- i:

  A numeric vector of indices to look up

## Value

the values of the lookup volume

## Examples

``` r
space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
ilv <- IndexLookupVol(space, c(1:100))
lookup(ilv, c(1, 2, 3))  # Look up values for indices 1, 2, and 3
#> [1] 1 2 3

```
