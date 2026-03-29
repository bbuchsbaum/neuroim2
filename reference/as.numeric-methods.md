# Convert SparseNeuroVol to numeric

Convert SparseNeuroVol to numeric

## Usage

``` r
# S4 method for class 'SparseNeuroVol'
as.numeric(x)

# S4 method for class 'ROIVol'
as.numeric(x)
```

## Arguments

- x:

  the object to convert

## Value

A numeric vector of length `nrow(x@coords)`
