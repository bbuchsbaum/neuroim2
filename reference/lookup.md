# Index Lookup operation

Index Lookup operation

## Usage

``` r
lookup(x, i, ...)
```

## Arguments

- x:

  the object to query

- i:

  the index to lookup

- ...:

  additional arguments

## Value

The value(s) at the specified index/indices of `x`.

## Examples

``` r
# Create a 64x64x64 space
space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))

# Create a lookup volume with first 100 indices
ilv <- IndexLookupVol(space, 1:100)

# Look up values for indices 1, 2, and 3
# Returns their positions in the sparse representation
lookup(ilv, c(1, 2, 3))
#> [1] 1 2 3

# Look up values outside the included indices
# Returns 0 for indices not in the lookup volume
lookup(ilv, c(101, 102))
#> [1] 0 0
```
