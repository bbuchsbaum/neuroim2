# Extract the number of dimensions of an object

Extract the number of dimensions of an object

Get number of dimensions in axis set

## Usage

``` r
ndim(x, ...)

# S4 method for class 'AxisSet'
ndim(x, ...)

# S4 method for class 'ClusteredNeuroVec'
ndim(x)

# S4 method for class 'NeuroObj'
ndim(x)

# S4 method for class 'NeuroHyperVec'
ndim(x)

# S4 method for class 'NeuroSpace'
ndim(x)
```

## Arguments

- x:

  An AxisSet object

- ...:

  Additional arguments (not used)

## Value

An integer representing the number of dimensions in `x`.

An integer representing the number of dimensions in `x`.

## Examples

``` r

x = NeuroSpace(c(10,10,10), spacing=c(1,1,1))
ndim(x) == 3
#> [1] TRUE
x = NeuroSpace(c(10,10,10,3), spacing=c(1,1,1))
ndim(x) == 4
#> [1] TRUE
```
