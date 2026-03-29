# Summary Methods for Neuroimaging Objects

Methods for the `Summary` group generic (e.g., `sum`, `min`, `max`,
`range`, `prod`, `any`, `all`) applied to neuroimaging data objects.

## Usage

``` r
# S4 method for class 'SparseNeuroVec'
Summary(x, ..., na.rm = FALSE)

# S4 method for class 'SparseNeuroVol'
Summary(x, ..., na.rm = FALSE)

# S4 method for class 'DenseNeuroVol'
Summary(x, ..., na.rm = FALSE)

# S4 method for class 'DenseNeuroVol'
Summary(x, ..., na.rm = FALSE)
```

## Arguments

- x:

  A neuroimaging object (SparseNeuroVec, SparseNeuroVol, or
  DenseNeuroVol)

- ...:

  Additional arguments passed to methods

- na.rm:

  Logical indicating whether to remove NA values before computation

## Value

The result of the summary operation

## Examples

``` r
# Create a simple volume
vol <- DenseNeuroVol(array(1:27, c(3,3,3)),
                     NeuroSpace(c(3L,3L,3L), c(1,1,1)))
sum(vol)
#> [1] 378
range(vol)
#> [1]  1 27
```
