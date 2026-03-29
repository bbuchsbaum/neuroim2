# Extract an ordered series of 3D volumes.

This function extracts an ordered series of 3D volumes from an object
that supplies volume data. The `indices` argument specifies the subset
of volumes to extract, and can be a vector of indices or a logical
vector. The return value is a list containing the extracted volumes in
the same order as the specified indices.

## Usage

``` r
vols(x, indices, ...)

# S4 method for class 'NeuroVec,numeric'
vols(x, indices)

# S4 method for class 'NeuroVec,missing'
vols(x)
```

## Arguments

- x:

  the object that supplies the volume data.

- indices:

  the subset of volumes to extract.

- ...:

  additional arguments to be passed to methods.

## Value

A `list` containing the extracted 3D volumes from `x` in the same order
as `indices`.

## Examples

``` r
vec <- read_vec(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
vs <- vols(vec)
length(vs) == dim(vec)[4]
#> [1] TRUE

vs <- vols(vec, indices=1:3)
length(vs) == 3
#> [1] TRUE
```
