# return the centroid of an object

return the centroid of an object

## Usage

``` r
centroid(x, ...)

# S4 method for class 'NeuroSpace'
centroid(x)

# S4 method for class 'ROICoords'
centroid(x)
```

## Arguments

- x:

  an object with a centroid

- ...:

  extra args

## Value

A numeric vector giving the centroid of `x`.

## Examples

``` r
bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
centroid(bspace)
#> [1] 10 10 10
```
