# Extract inverse image coordinate transformation

Extract inverse image coordinate transformation

## Usage

``` r
inverse_trans(x)

# S4 method for class 'NeuroSpace'
inverse_trans(x)
```

## Arguments

- x:

  an object

## Value

A numeric 4x4 `matrix` that maps from real-world coordinates back to
grid coordinates.

## Examples

``` r
bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
itrans <- inverse_trans(bspace)
identical(trans(bspace) %*% inverse_trans(bspace), diag(4))
#> [1] TRUE
```
