# Create a patch set from a NeuroVol object

This function creates a patch set from a NeuroVol object given specified
dimensions

This function creates a patch set from a NeuroVol object given specified
dimensions and a mask.

## Usage

``` r
# S4 method for class 'NeuroVol,numeric,missing'
patch_set(x, dims, mask, ...)

# S4 method for class 'NeuroVol,numeric,LogicalNeuroVol'
patch_set(x, dims, mask, ...)
```

## Arguments

- x:

  a NeuroVol object

- dims:

  the dimensions of the patch

- mask:

  the mask defining the valid patch centers

- ...:

  additional args

## Value

A deferred list of patches.

A deferred list of patches.
