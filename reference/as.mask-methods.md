# Convert NeuroVol to a mask

This method converts a NeuroVol object to a mask by setting all positive
values to TRUE and all non-positive values to FALSE.

This method converts a NeuroVol object to a mask by setting the
specified indices to TRUE and the remaining elements to FALSE.

## Usage

``` r
# S4 method for class 'NeuroVol,missing'
as.mask(x)

# S4 method for class 'NeuroVol,numeric'
as.mask(x, indices)
```

## Arguments

- x:

  A NeuroVol object to convert to a mask.

- indices:

  A numeric vector containing the indices of the input NeuroVol that
  should be set to TRUE in the resulting mask.

## Value

A LogicalNeuroVol object representing the mask created from the input
NeuroVol.

A LogicalNeuroVol object representing the mask created from the input
NeuroVol with specified indices.
