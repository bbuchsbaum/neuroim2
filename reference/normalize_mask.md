# Normalize a mask to a logical array

Coerces various mask representations into a logical 3D array of the
requested dimensions. Accepts: logical arrays, integer index vectors,
`LogicalNeuroVol`, numeric `NeuroVol` (thresholded \> 0), logical
vectors, or `NULL` (all `TRUE`).

## Usage

``` r
normalize_mask(mask, target_dim)
```

## Arguments

- mask:

  The mask input (see description).

- target_dim:

  Integer vector of length 3 giving the target dimensions.

## Value

A logical array of dimension `target_dim`.
