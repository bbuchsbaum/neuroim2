# Extract Image Origin

Extract Image Origin

## Usage

``` r
origin(x)

# S4 method for class 'NeuroHyperVec'
origin(x)

# S4 method for class 'NeuroSpace'
origin(x)

# S4 method for class 'NeuroVol'
origin(x)

# S4 method for class 'NeuroVec'
origin(x)
```

## Arguments

- x:

  an object with an origin

## Value

A numeric vector giving the origin of `x`.

## Examples

``` r
bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
stopifnot(origin(bspace) == c(0,0,0))
```
