# Extract Spatial Bounds of an Image

This function extracts the spatial bounds (origin + dim \* spacing) of
an image represented by the input object.

## Usage

``` r
bounds(x)

# S4 method for class 'NeuroSpace'
bounds(x)
```

## Arguments

- x:

  The object with the \`bounds\` property, typically an image.

## Value

A numeric `matrix` with two columns specifying the min (column 1) and
max (column 2) bounds of each dimension of `x`.

## Examples

``` r
bspace <- NeuroSpace(c(10, 10, 10), c(2, 2, 2))
b <- bounds(bspace)
nrow(b) == ndim(bspace)
#> [1] TRUE
ncol(b) == 2
#> [1] TRUE
```
