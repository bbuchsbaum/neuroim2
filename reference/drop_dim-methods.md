# Drop a Dimension from an Object

This function removes a specified dimension from a given object, such as
a matrix or an array.

## Usage

``` r
drop_dim(x, dimnum)

# S4 method for class 'AxisSet2D,numeric'
drop_dim(x, dimnum)

# S4 method for class 'AxisSet2D,missing'
drop_dim(x, dimnum)

# S4 method for class 'AxisSet3D,numeric'
drop_dim(x, dimnum)

# S4 method for class 'AxisSet3D,missing'
drop_dim(x, dimnum)

# S4 method for class 'NeuroSpace,numeric'
drop_dim(x, dimnum)

# S4 method for class 'NeuroSpace,missing'
drop_dim(x)
```

## Arguments

- x:

  An AxisSet3D object

- dimnum:

  Numeric index of dimension to drop (optional)

## Value

An object of the same class as `x` with the specified dimension removed.

## Examples

``` r
# Create a NeuroSpace object with dimensions (10, 10, 10)
x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))

# Drop the first dimension
x1 <- drop_dim(x, 1)

# Check the new dimensions
ndim(x1) == 2
#> [1] TRUE
dim(x1)[1] == 10
#> [1] TRUE
```
