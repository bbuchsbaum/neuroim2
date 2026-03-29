# Find Dimensions of a Given Axis

This function returns the dimension of the specified axis for a given
object, such as a matrix or an array.

## Usage

``` r
which_dim(x, axis)

# S4 method for class 'NeuroSpace,NamedAxis'
which_dim(x, axis)
```

## Arguments

- x:

  The NeuroSpace object

- axis:

  The NamedAxis to find

## Value

An integer representing the dimension index of the specified `axis` for
the object `x`.

## Examples

``` r
x <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
which_dim(x, x@axes@j) == 2
#> [1] TRUE
```
