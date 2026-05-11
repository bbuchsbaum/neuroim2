# Get the length of a given dimension of an object

This function returns the length of a given axis (dimension) of an
object. The axis can be specified using its position or name.

## Usage

``` r
dim_of(x, axis)

# S4 method for class 'NeuroSpace,NamedAxis'
dim_of(x, axis)
```

## Arguments

- x:

  The NeuroSpace object

- axis:

  The NamedAxis to query

## Value

An integer representing the length of the specified `axis` of `x`.

## Examples

``` r

x <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
stopifnot(dim_of(x, x@axes@i) == 10)
```
