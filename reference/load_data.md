# Read data from a data source.

This function loads data from a data source and returns it in a format
that is compatible with other functions in the neuroim2 package. The
format of the returned data depends on the type of data source used.

## Usage

``` r
load_data(x, ...)
```

## Arguments

- x:

  a data source.

- ...:

  additional arguments to be passed to methods.

## Value

An R object containing loaded data, in a format compatible with the
neuroim2 package.

## Examples

``` r
# Create a NeuroVolSource from a NIFTI file and load it
fname <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
src <- NeuroVolSource(fname)
vol <- load_data(src)
# The loaded volume is a DenseNeuroVol object
class(vol)
#> [1] "DenseNeuroVol"
#> attr(,"package")
#> [1] "neuroim2"
dim(vol)
#> [1] 64 64 25
```
