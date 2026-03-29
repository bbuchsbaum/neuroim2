# Extract Image Axes

Extract Image Axes

## Usage

``` r
axes(x)

# S4 method for class 'NeuroSpace'
axes(x)
```

## Arguments

- x:

  an object with a set of axes

## Value

An object representing the `axes` of `x`.

## Examples

``` r
x <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
class(axes(x)) == "AxisSet3D"
#> [1] TRUE
```
