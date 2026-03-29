# Summary of Neuroimaging Objects

Provides a concise summary of neuroimaging volume and vector objects,
including spatial metadata and data statistics.

## Usage

``` r
# S4 method for class 'NeuroVol'
summary(object, ...)

# S4 method for class 'DenseNeuroVec'
summary(object, ...)

# S4 method for class 'SparseNeuroVec'
summary(object, ...)
```

## Arguments

- object:

  A neuroimaging object.

- ...:

  Additional arguments (currently ignored).

## Value

An object of class `"summary.NeuroVol"` or `"summary.NeuroVec"`
(invisibly), with a print method.

## Examples

``` r
vol <- DenseNeuroVol(array(rnorm(27), c(3,3,3)),
                     NeuroSpace(c(3L,3L,3L), c(1,1,1)))
summary(vol)
#> <DenseNeuroVol> [3 x 3 x 3]
#>   Spacing    : 1 x 1 x 1 mm
#>   Origin     : 0, 0, 0
#>   Orientation: RAS
#>   Range      : [-1.538, 1.985]
#>   Mean (SD)  : 0.27 (0.9236)
#>   Non-zero   : 27 / 27 voxels
```
