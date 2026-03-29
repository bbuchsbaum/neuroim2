# Extract Coordinates from an IndexLookupVol Object

Extracts the coordinates from an IndexLookupVol object based on a given
index.

## Usage

``` r
# S4 method for class 'IndexLookupVol'
coords(x, i)

# S4 method for class 'ROIVol'
coords(x, real = FALSE)

# S4 method for class 'ROICoords'
coords(x, real = FALSE)

# S4 method for class 'AbstractSparseNeuroVec'
coords(x, i)
```

## Arguments

- x:

  An
  [`IndexLookupVol`](https://bbuchsbaum.github.io/neuroim2/reference/IndexLookupVol-class.md)
  object to extract coordinates from

- i:

  The index into the lookup volume

- real:

  if `TRUE`, return coordinates in real world units

## Value

A matrix of coordinates

## Examples

``` r
space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
ilv <- IndexLookupVol(space, c(1:100))
coords(ilv, 1)  # Extract coordinates for index 1
#>      [,1] [,2] [,3]
#> [1,]    1    1    1

```
