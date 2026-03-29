# Add a Dimension to an Object

This function adds a new dimension to a given object, such as a matrix
or an array.

## Usage

``` r
add_dim(x, n)

# S4 method for class 'NeuroSpace,numeric'
add_dim(x, n)
```

## Arguments

- x:

  The NeuroSpace object

- n:

  Numeric value specifying the size of the new dimension

## Value

An object of the same class as `x` with the new dimension added.

## Examples

``` r
# Create a NeuroSpace object
x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))

# Add a new dimension with size 10
x1 <- add_dim(x, 10)

# Check the new dimension
ndim(x1) == 4
#> [1] TRUE
dim(x1)[4] == 10
#> [1] TRUE
```
